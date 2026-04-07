# BioAnalyze

BioAnalyze is the local Python pipeline for `CURATED_SET/BioAnalyze`.

## Windows setup

Recommended interpreter: `Python 3.11 x64`.

Install Python if it is missing:

```powershell
winget install --id Python.Python.3.11 --source winget
```

Create and populate the repo-local virtual environment:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\setup_bioanalyze_env.ps1
```

Add notebook tooling on top of the base environment:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\setup_bioanalyze_env.ps1 -WithNotebook
```

Run the bootstrap smoke checks:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\setup_bioanalyze_env.ps1 -SmokeTest
```

## Manual setup

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip setuptools wheel
.\.venv\Scripts\python.exe -m pip install -r CURATED_SET\BioAnalyze\requirements.txt
```

Notebook mode:

```powershell
.\.venv\Scripts\python.exe -m pip install -r CURATED_SET\BioAnalyze\requirements-dev.txt
```

## External storage

BioAnalyze resolves external data in this order:

1. `HISTONEDB_EXTERNAL_STORAGE`
2. Sibling folder next to the repo: `<repo-parent>\histonedb_external_storage`

Example override:

```powershell
$env:HISTONEDB_EXTERNAL_STORAGE = 'D:\data\histonedb_external_storage'
```

The raw-data subtree expected by the scripts is then:

```text
<external-storage>\BioAnalyze\raw
```

## Standard run style

Use the repo-local interpreter instead of a global Python:

```powershell
.\.venv\Scripts\python.exe CURATED_SET\BioAnalyze\scripts\expression\build_human_h2a_ntpm_gene_expression.py --help
```
