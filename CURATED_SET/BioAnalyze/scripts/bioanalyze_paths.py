#!/usr/bin/env python
"""Shared path helpers for BioAnalyze CLI scripts."""

from __future__ import annotations

import os
from pathlib import Path


EXTERNAL_STORAGE_ENV_VAR = "HISTONEDB_EXTERNAL_STORAGE"

SCRIPTS_ROOT = Path(__file__).resolve().parent
BIOANALYZE_ROOT = SCRIPTS_ROOT.parent
CURATED_SET_ROOT = BIOANALYZE_ROOT.parent
REPO_ROOT = CURATED_SET_ROOT.parent


def _normalize(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def get_repo_root() -> Path:
    return REPO_ROOT


def get_bioanalyze_root() -> Path:
    return BIOANALYZE_ROOT


def get_bioanalyze_path(*parts: str) -> Path:
    return BIOANALYZE_ROOT.joinpath(*parts)


def get_bioanalyze_data_root() -> Path:
    return get_bioanalyze_path("data")


def get_bioanalyze_figures_root() -> Path:
    return get_bioanalyze_path("figures")


def get_bioanalyze_stats_root() -> Path:
    return get_bioanalyze_path("stats")


def get_bioanalyze_audits_root() -> Path:
    return get_bioanalyze_path("audits")


def get_gene_tissue_stats_root() -> Path:
    return get_bioanalyze_stats_root() / "gene_tissue"


def get_gene_tissue_ranking_root() -> Path:
    return get_gene_tissue_stats_root() / "ranking"


def get_gene_tissue_barplot_root() -> Path:
    return get_gene_tissue_stats_root() / "barplot"


def get_gene_tissue_boxplot_root() -> Path:
    return get_gene_tissue_barplot_root()


def get_external_storage_root() -> Path:
    override = os.environ.get(EXTERNAL_STORAGE_ENV_VAR)
    if override:
        return _normalize(override)
    return _normalize(REPO_ROOT.parent / "histonedb_external_storage")


def get_bioanalyze_external_root() -> Path:
    return get_external_storage_root() / "BioAnalyze"


def get_bioanalyze_raw_root() -> Path:
    return get_bioanalyze_external_root() / "raw"


def get_bioanalyze_synteny_root() -> Path:
    return get_bioanalyze_external_root() / "synteny"


def get_documents_root() -> Path:
    return Path.home().expanduser() / "Documents"


def get_downloads_root() -> Path:
    return Path.home().expanduser() / "Downloads"
