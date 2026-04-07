[CmdletBinding()]
param(
    [switch]$WithNotebook,
    [switch]$SmokeTest
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$WingetInstallCommand = "winget install --id Python.Python.3.11 --source winget"
$RepoRoot = Split-Path -Parent $PSScriptRoot
$BioAnalyzeRoot = Join-Path $RepoRoot "CURATED_SET\\BioAnalyze"
$RequirementsPath = Join-Path $BioAnalyzeRoot "requirements.txt"
$DevRequirementsPath = Join-Path $BioAnalyzeRoot "requirements-dev.txt"
$VenvDir = Join-Path $RepoRoot ".venv"
$VenvPython = Join-Path $VenvDir "Scripts\\python.exe"
$DefaultExternalStorageRoot = Join-Path (Split-Path -Parent $RepoRoot) "histonedb_external_storage"

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message" -ForegroundColor Cyan
}

function Invoke-Checked {
    param(
        [string]$Message,
        [scriptblock]$Action
    )

    Write-Step $Message
    & $Action
    if ($LASTEXITCODE -ne 0) {
        throw "Step failed: $Message"
    }
}

function Get-UsablePython {
    $candidates = New-Object System.Collections.Generic.List[string]

    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if ($pythonCommand -and $pythonCommand.Source) {
        $candidates.Add($pythonCommand.Source)
    }

    $candidateRoots = @(
        (Join-Path $env:LocalAppData "Programs\Python\Python311\python.exe"),
        (Join-Path $env:ProgramFiles "Python311\python.exe"),
        (Join-Path ${env:ProgramFiles(x86)} "Python311\python.exe")
    )
    foreach ($candidate in $candidateRoots) {
        if ($candidate -and (Test-Path $candidate)) {
            $candidates.Add($candidate)
        }
    }

    foreach ($candidate in $candidates | Select-Object -Unique) {
        try {
            $probe = & $candidate -c "import struct, sys; print(sys.executable); print(f'{sys.version_info[0]}.{sys.version_info[1]}.{sys.version_info[2]}'); print(struct.calcsize('P') * 8)" 2>$null
        } catch {
            continue
        }

        if ($LASTEXITCODE -ne 0 -or -not $probe) {
            continue
        }

        $lines = @($probe | Where-Object { $_ -ne "" })
        if ($lines.Count -lt 3) {
            continue
        }

        return [pscustomobject]@{
            Executable = $lines[0]
            Version = $lines[1]
            Bits = [int]$lines[2]
        }
    }

    return $null
}

function Invoke-PythonInline {
    param(
        [string]$Message,
        [string]$Code
    )

    Write-Step $Message
    $Code | & $VenvPython -
    if ($LASTEXITCODE -ne 0) {
        throw "Python inline check failed: $Message"
    }
}

$PythonInfo = Get-UsablePython
if (-not $PythonInfo) {
    Write-Error (
        "Python 3.11 x64 was not found on PATH. Install it with:`n" +
        "$WingetInstallCommand"
    )
}

$VersionParts = $PythonInfo.Version.Split(".")
if ($VersionParts.Count -lt 2 -or $VersionParts[0] -ne "3" -or $VersionParts[1] -ne "11") {
    Write-Error (
        "Detected Python $($PythonInfo.Version) at $($PythonInfo.Executable). " +
        "BioAnalyze uses Python 3.11 x64. Install it with:`n$WingetInstallCommand"
    )
}

if ($PythonInfo.Bits -ne 64) {
    Write-Error (
        "Detected a $($PythonInfo.Bits)-bit Python at $($PythonInfo.Executable). " +
        "BioAnalyze uses Python 3.11 x64. Install it with:`n$WingetInstallCommand"
    )
}

if (-not (Test-Path $RequirementsPath)) {
    throw "Missing requirements file: $RequirementsPath"
}

if (-not (Test-Path $VenvPython)) {
    Invoke-Checked "Creating repo-local virtual environment" {
        & $PythonInfo.Executable -m venv $VenvDir
    }
} else {
    Write-Step "Reusing existing virtual environment at $VenvDir"
}

Invoke-Checked "Upgrading pip, setuptools, and wheel" {
    & $VenvPython -m pip install --upgrade pip setuptools wheel
}

Invoke-Checked "Installing BioAnalyze runtime requirements" {
    & $VenvPython -m pip install -r $RequirementsPath
}

if ($WithNotebook) {
    Invoke-Checked "Installing notebook extras" {
        & $VenvPython -m pip install -r $DevRequirementsPath
    }
}

if ($SmokeTest) {
    $SmokeRoot = Join-Path $RepoRoot ("tmp\\bioanalyze-smoke\\" + (Get-Date -Format "yyyyMMdd-HHmmss"))
    $NtpmDataDir = Join-Path $SmokeRoot "expression_nTPM\\data"
    $NtpmFigDir = Join-Path $SmokeRoot "expression_nTPM\\figures"
    $CodonsOutDir = Join-Path $SmokeRoot "codons"
    $NtpmInput = Join-Path $DefaultExternalStorageRoot "BioAnalyze\\raw\\expression_nTPM\\human\\test.tsv"

    Invoke-PythonInline "Import smoke for runtime packages" @"
import requests
import numpy
import pandas
import matplotlib
import seaborn
import notebook
import ipykernel
from Bio import SeqIO

print("runtime imports ok")
"@

    $PreviousExternalStorage = [Environment]::GetEnvironmentVariable("HISTONEDB_EXTERNAL_STORAGE", "Process")
    Remove-Item Env:HISTONEDB_EXTERNAL_STORAGE -ErrorAction SilentlyContinue

    Invoke-PythonInline "Fallback path smoke without HISTONEDB_EXTERNAL_STORAGE" @"
from pathlib import Path
import sys

repo_root = Path(r"$RepoRoot")
scripts_root = repo_root / "CURATED_SET" / "BioAnalyze" / "scripts"
sys.path.insert(0, str(scripts_root))

from bioanalyze_paths import get_external_storage_root

expected = (repo_root.parent / "histonedb_external_storage").resolve()
actual = get_external_storage_root()
if actual != expected:
    raise SystemExit(f"Expected {expected}, got {actual}")

print(actual)
"@

    $env:HISTONEDB_EXTERNAL_STORAGE = $DefaultExternalStorageRoot
    Invoke-PythonInline "Environment override path smoke" @"
from pathlib import Path
import sys

repo_root = Path(r"$RepoRoot")
scripts_root = repo_root / "CURATED_SET" / "BioAnalyze" / "scripts"
sys.path.insert(0, str(scripts_root))

from bioanalyze_paths import get_external_storage_root

expected = Path(r"$DefaultExternalStorageRoot").resolve()
actual = get_external_storage_root()
if actual != expected:
    raise SystemExit(f"Expected {expected}, got {actual}")

print(actual)
"@

    if ($null -eq $PreviousExternalStorage) {
        Remove-Item Env:HISTONEDB_EXTERNAL_STORAGE -ErrorAction SilentlyContinue
    } else {
        $env:HISTONEDB_EXTERNAL_STORAGE = $PreviousExternalStorage
    }

    Invoke-Checked "CLI smoke: build_human_h2a_ntpm_gene_expression.py --help" {
        & $VenvPython (Join-Path $RepoRoot "CURATED_SET\\BioAnalyze\\scripts\\expression\\build_human_h2a_ntpm_gene_expression.py") --help
    }

    Invoke-Checked "CLI smoke: build_codon_heatmaps_all61.py --help" {
        & $VenvPython (Join-Path $RepoRoot "CURATED_SET\\BioAnalyze\\scripts\\codons\\build_codon_heatmaps_all61.py") --help
    }

    Invoke-Checked "CLI smoke: build_h2aj_synteny_plot.py --help" {
        & $VenvPython (Join-Path $RepoRoot "CURATED_SET\\BioAnalyze\\scripts\\h2aj_synteny\\build_h2aj_synteny_plot.py") --help
    }

    Invoke-Checked "Real smoke: human HPA nTPM pipeline into tmp/" {
        & $VenvPython `
            (Join-Path $RepoRoot "CURATED_SET\\BioAnalyze\\scripts\\expression\\build_human_h2a_ntpm_gene_expression.py") `
            --input-tsv $NtpmInput `
            --out-data-dir $NtpmDataDir `
            --out-fig-dir $NtpmFigDir
    }

    Invoke-Checked "Real smoke: codon heatmap pipeline into tmp/" {
        & $VenvPython `
            (Join-Path $RepoRoot "CURATED_SET\\BioAnalyze\\scripts\\codons\\build_codon_heatmaps_all61.py") `
            --dataset without-short `
            --output-dir $CodonsOutDir
    }

    Write-Step "Smoke artifacts written to $SmokeRoot"
}

Write-Host ""
Write-Host "Virtual environment ready:" -ForegroundColor Green
Write-Host "  $VenvPython"
