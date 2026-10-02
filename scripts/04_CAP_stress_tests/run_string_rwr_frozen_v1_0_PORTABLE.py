# -*- coding: utf-8 -*-
"""
run_string_rwr_frozen_v1_0.py

Frozen-upstream STRING random walk with restart (RWR) for the syndrome target
projection workflow.

This script is deliberately disease-blind:
- it reads NO CAP transcriptomic data;
- it does NOT tune restart probability using disease alignment;
- it does NOT select a syndrome representation using downstream outcomes.

Primary upstream input
----------------------
LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv

Primary network definition
--------------------------
STRING v12.0 high-confidence human functional association network,
combined_score >= 700.

Primary RWR definition
----------------------
Let A be the symmetric weighted adjacency matrix and

    P = D^{-1} A

the ROW-stochastic transition matrix.

For a COLUMN probability vector p, the correct update is

    p_{k+1} = r * p0 + (1-r) * P.T @ p_k

where r is the restart probability.

IMPORTANT:
The transpose is intentional and required.  Using P @ p with a row-normalized
matrix and a column probability vector does not conserve probability mass.

Primary restart probability:
    r = 0.5

Predeclared sensitivity:
    r = 0.3 and r = 0.7

The script performs:
1. STRING edge loading / cleaning / thresholding.
2. Frozen LOPO seed mapping and mapped-weight renormalization.
3. Correct RWR with mass-conservation checks.
4. r=0.3/0.5/0.7 sensitivity analysis.
5. Toy-network unit test with an analytic solution.
6. Determinism test.
7. Node-order invariance test.
8. Connected-component and seed-mapping QA.
9. Pre/post-diffusion syndrome-similarity diagnostics.
10. Export of a target-database/network reference table for LATER matched-null
    analysis. Null testing itself is intentionally deferred until the disease
    alignment statistic is frozen.

No CAP data are read anywhere in this file.
"""

from __future__ import annotations

from pathlib import Path
import os
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import re
import shutil
import urllib.request
import warnings

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, diags, eye
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE
# =============================================================================

PROJECT_ROOT = Path(os.environ.get("DCSMI_PROJECT_ROOT", Path(__file__).resolve().parents[2]))

# Frozen upstream target seeds.
FROZEN_SEED_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
)

# Optional but strongly recommended.
# Used only to export the later matched-null reference universe.
TARGET_UBIQUITY_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "target_database_target_ubiquity_v1_1.csv"
)

# -------------------------------------------------------------------------
# STRING input mode
# -------------------------------------------------------------------------
# "gene_symbol_edges":
#     use an already prepared table such as:
#         geneA  geneB  score
#
# "string_raw":
#     use official STRING v12.0 files:
#         9606.protein.links.v12.0.txt.gz
#         9606.protein.info.v12.0.txt.gz
#     and map STRING protein IDs to preferred gene/protein display names.
#
# For continuity with the existing project, gene_symbol_edges is the default.
PPI_INPUT_MODE = "gene_symbol_edges"

PPI_GENE_EDGE_FILE = Path(os.environ.get("DCSMI_STRING_EDGE_FILE", PROJECT_ROOT / "external_data" / "string_ppi_edges_700.tsv"))
PPI_GENE_A_COLUMN = "geneA"
PPI_GENE_B_COLUMN = "geneB"
PPI_SCORE_COLUMN = "score"

# Raw STRING v12.0 mode.
STRING_RAW_DIR = PROJECT_ROOT / "STRING_v12_0"
STRING_RAW_LINKS_FILE = STRING_RAW_DIR / "9606.protein.links.v12.0.txt.gz"
STRING_RAW_INFO_FILE = STRING_RAW_DIR / "9606.protein.info.v12.0.txt.gz"

# Optional convenience download. Keep False if files are already available.
DOWNLOAD_STRING_V12_IF_MISSING = False

STRING_V12_LINKS_URL = (
    "https://stringdb-downloads.org/download/protein.links.v12.0/"
    "9606.protein.links.v12.0.txt.gz"
)
STRING_V12_INFO_URL = (
    "https://stringdb-downloads.org/download/protein.info.v12.0/"
    "9606.protein.info.v12.0.txt.gz"
)

# Output.
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "secondary_STRING_RWR"

# -------------------------------------------------------------------------
# Frozen analysis settings
# -------------------------------------------------------------------------

SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]

# STRING combined scores are distributed on the 0-1000 scale.
STRING_COMBINED_SCORE_THRESHOLD = 700

# Use STRING combined score as edge weight.
# False would make the thresholded graph unweighted; do not change this for
# the primary analysis after downstream outcomes are inspected.
USE_STRING_EDGE_WEIGHTS = True

# Restart probability r.
PRIMARY_RESTART = 0.5
SENSITIVITY_RESTARTS = [0.3, 0.7]

# Convergence.
RWR_TOL_L1 = 1e-12
RWR_MAX_ITER = 10000
MASS_TOL = 1e-10

# QA.
RANDOM_SEED = 20260921
NODE_ORDER_INVARIANCE_TOL = 1e-11
DETERMINISM_TOL = 1e-15
TOY_TEST_TOL = 1e-12

# Similarity diagnostics.
TOP_K_VALUES = [50, 100, 200]

# The upstream frozen seed weights must sum to 1 per syndrome before STRING
# mapping. This catches accidental use of a wrong file.
FROZEN_SEED_SUM_TOL = 1e-8

# Write a compact Excel QA workbook in addition to CSVs.
WRITE_EXCEL_SUMMARY = True

# Declared provenance. The SHA256 of actual input files is also written.
STRING_VERSION_DECLARED = "12.0"
STRING_SPECIES = "Homo sapiens"
STRING_TAXONOMY_ID = 9606
STRING_NETWORK_SEMANTICS = "functional association network"


# =============================================================================
# 2. BASIC HELPERS
# =============================================================================

def norm_gene(x) -> str:
    if pd.isna(x):
        return ""
    s = re.sub(r"\s+", "", str(x).strip())
    return s.upper()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def safe_sha256(path: Path) -> str | None:
    return sha256_file(path) if path.exists() else None


def require_columns(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def restart_label(r: float) -> str:
    return f"r{r:.2f}".replace(".", "p")


def read_delimited_auto(path: Path) -> pd.DataFrame:
    """
    Read common edge-table formats.
    """
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith(".csv"):
        return pd.read_csv(path)

    if suffixes.endswith(".tsv") or suffixes.endswith(".txt"):
        return pd.read_csv(path, sep="\t")

    if suffixes.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")

    if suffixes.endswith(".txt.gz"):
        # STRING raw links use whitespace; preprocessed tables are often tabs.
        # sep=r"\\s+" safely handles both.
        return pd.read_csv(path, sep=r"\s+", compression="gzip")

    raise ValueError(f"Unsupported delimited file type: {path}")


def maybe_download(url: str, destination: Path) -> None:
    if destination.exists():
        return
    if not DOWNLOAD_STRING_V12_IF_MISSING:
        raise FileNotFoundError(
            f"Missing STRING file:\n{destination}\n\n"
            "Download it from STRING v12.0 or set "
            "DOWNLOAD_STRING_V12_IF_MISSING=True."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading:\n  {url}\n-> {destination}")
    tmp = destination.with_suffix(destination.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(destination)


# =============================================================================
# 3. LOAD FROZEN SEEDS
# =============================================================================

def load_frozen_seeds(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Frozen seed file not found:\n{path}")

    x = pd.read_csv(path)
    require_columns(
        x,
        [
            "syndrome_code",
            "target_id",
            "ubiquity_corrected_sum1",
        ],
        path.name,
    )

    x = x.copy()
    x["syndrome_code"] = x["syndrome_code"].astype(str).str.strip()
    x["target_id"] = x["target_id"].map(norm_gene)
    x["ubiquity_corrected_sum1"] = pd.to_numeric(
        x["ubiquity_corrected_sum1"], errors="coerce"
    )

    if x["ubiquity_corrected_sum1"].isna().any():
        raise ValueError("Frozen seed file contains non-numeric seed weights.")

    x = x[
        x["syndrome_code"].isin(SYNDROME_CODES)
        & x["target_id"].ne("")
        & x["ubiquity_corrected_sum1"].gt(0)
    ].copy()

    # Frozen file should not contain duplicate syndrome-target rows.
    dup = x.duplicated(["syndrome_code", "target_id"], keep=False)
    if dup.any():
        bad = x.loc[dup, ["syndrome_code", "target_id"]].head(20)
        raise ValueError(
            "Frozen seed file contains duplicate syndrome-target rows.\n"
            f"{bad.to_string(index=False)}"
        )

    present = set(x["syndrome_code"])
    missing_syndromes = [s for s in SYNDROME_CODES if s not in present]
    if missing_syndromes:
        raise ValueError(
            f"Frozen seed file missing syndromes: {missing_syndromes}"
        )

    sums = x.groupby("syndrome_code")["ubiquity_corrected_sum1"].sum()
    bad_sums = sums[(sums - 1.0).abs() > FROZEN_SEED_SUM_TOL]
    if len(bad_sums):
        raise RuntimeError(
            "Frozen upstream seed weights do not sum to 1.\n"
            f"{bad_sums.to_string()}\n"
            "Do not silently renormalize the upstream frozen file."
        )

    return x


# =============================================================================
# 4. LOAD / PREPARE STRING NETWORK
# =============================================================================

def load_gene_symbol_edges() -> pd.DataFrame:
    path = PPI_GENE_EDGE_FILE
    if not path.exists():
        raise FileNotFoundError(
            f"PPI_GENE_EDGE_FILE not found:\n{path}\n\n"
            "Either point PPI_GENE_EDGE_FILE to the existing "
            "string_ppi_edges_700.tsv or switch PPI_INPUT_MODE='string_raw'."
        )

    raw = read_delimited_auto(path)
    require_columns(
        raw,
        [PPI_GENE_A_COLUMN, PPI_GENE_B_COLUMN, PPI_SCORE_COLUMN],
        path.name,
    )

    x = raw[
        [PPI_GENE_A_COLUMN, PPI_GENE_B_COLUMN, PPI_SCORE_COLUMN]
    ].copy()
    x.columns = ["geneA", "geneB", "score"]
    return x


def load_string_raw_edges() -> pd.DataFrame:
    maybe_download(STRING_V12_LINKS_URL, STRING_RAW_LINKS_FILE)
    maybe_download(STRING_V12_INFO_URL, STRING_RAW_INFO_FILE)

    print("Reading STRING protein info...")
    info = pd.read_csv(
        STRING_RAW_INFO_FILE,
        sep="\t",
        compression="gzip",
        dtype=str,
    )

    # STRING info file usually uses '#string_protein_id'.
    protein_col = None
    for c in info.columns:
        if str(c).lstrip("#").lower() == "string_protein_id":
            protein_col = c
            break

    if protein_col is None or "preferred_name" not in info.columns:
        raise ValueError(
            "Could not identify STRING protein-info columns. "
            f"Available: {list(info.columns)}"
        )

    id_to_gene = dict(
        zip(
            info[protein_col].astype(str),
            info["preferred_name"].map(norm_gene),
        )
    )

    print("Reading STRING protein links...")
    links = pd.read_csv(
        STRING_RAW_LINKS_FILE,
        sep=r"\s+",
        compression="gzip",
    )
    require_columns(
        links,
        ["protein1", "protein2", "combined_score"],
        STRING_RAW_LINKS_FILE.name,
    )

    x = links[["protein1", "protein2", "combined_score"]].copy()
    x["geneA"] = x["protein1"].map(id_to_gene).fillna("")
    x["geneB"] = x["protein2"].map(id_to_gene).fillna("")
    x["score"] = pd.to_numeric(x["combined_score"], errors="coerce")

    return x[["geneA", "geneB", "score"]]


def clean_string_edges(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    x = raw.copy()
    require_columns(x, ["geneA", "geneB", "score"], "STRING edges")

    n_raw = len(x)

    x["geneA"] = x["geneA"].map(norm_gene)
    x["geneB"] = x["geneB"].map(norm_gene)
    x["score"] = pd.to_numeric(x["score"], errors="coerce")

    x = x[
        x["geneA"].ne("")
        & x["geneB"].ne("")
        & x["score"].notna()
    ].copy()

    # Detect score scale.
    max_score = float(x["score"].max()) if len(x) else np.nan
    if not np.isfinite(max_score):
        raise ValueError("No valid STRING scores.")

    if max_score <= 1.000001:
        score_scale = "0_to_1"
        threshold = STRING_COMBINED_SCORE_THRESHOLD / 1000.0
        x["score_0_1"] = x["score"].astype(float)
    else:
        score_scale = "0_to_1000"
        threshold = float(STRING_COMBINED_SCORE_THRESHOLD)
        x["score_0_1"] = x["score"].astype(float) / 1000.0

    x = x[x["score"] >= threshold].copy()

    # Remove self loops.
    x = x[x["geneA"] != x["geneB"]].copy()

    # Canonical undirected pair. This prevents AB/BA duplicates from doubling
    # edge weight.
    a = x["geneA"].to_numpy(dtype=str)
    b = x["geneB"].to_numpy(dtype=str)
    x["u"] = np.where(a <= b, a, b)
    x["v"] = np.where(a <= b, b, a)

    # If duplicate evidence rows exist, keep the maximum combined score for
    # that undirected pair.
    x = (
        x.groupby(["u", "v"], as_index=False)
        .agg(
            string_score_raw=("score", "max"),
            string_score_0_1=("score_0_1", "max"),
        )
    )

    if USE_STRING_EDGE_WEIGHTS:
        x["edge_weight"] = x["string_score_0_1"]
    else:
        x["edge_weight"] = 1.0

    stats = {
        "raw_rows": int(n_raw),
        "clean_high_confidence_undirected_edges": int(len(x)),
        "score_scale_detected": score_scale,
        "threshold_on_detected_scale": threshold,
        "min_retained_raw_score": (
            float(x["string_score_raw"].min()) if len(x) else None
        ),
        "max_retained_raw_score": (
            float(x["string_score_raw"].max()) if len(x) else None
        ),
    }

    return x, stats


def load_string_edges() -> tuple[pd.DataFrame, dict, dict]:
    mode = PPI_INPUT_MODE.strip().lower()

    if mode == "gene_symbol_edges":
        raw = load_gene_symbol_edges()
        source_files = {
            "gene_symbol_edge_file": str(PPI_GENE_EDGE_FILE),
            "gene_symbol_edge_sha256": sha256_file(PPI_GENE_EDGE_FILE),
        }
    elif mode == "string_raw":
        raw = load_string_raw_edges()
        source_files = {
            "raw_links_file": str(STRING_RAW_LINKS_FILE),
            "raw_links_sha256": sha256_file(STRING_RAW_LINKS_FILE),
            "raw_info_file": str(STRING_RAW_INFO_FILE),
            "raw_info_sha256": sha256_file(STRING_RAW_INFO_FILE),
            "raw_links_url": STRING_V12_LINKS_URL,
            "raw_info_url": STRING_V12_INFO_URL,
        }
    else:
        raise ValueError(
            "PPI_INPUT_MODE must be 'gene_symbol_edges' or 'string_raw'."
        )

    clean, edge_stats = clean_string_edges(raw)

    if clean.empty:
        raise RuntimeError("No STRING edges remain after cleaning/thresholding.")

    return clean, edge_stats, source_files


# =============================================================================
# 5. BUILD SPARSE ROW-STOCHASTIC TRANSITION MATRIX
# =============================================================================

def build_network(
    edges: pd.DataFrame,
    node_order: list[str] | None = None,
) -> dict:
    if node_order is None:
        nodes = sorted(set(edges["u"]) | set(edges["v"]))
    else:
        nodes = list(node_order)

    node_set = set(nodes)
    edge_nodes = set(edges["u"]) | set(edges["v"])
    if node_set != edge_nodes:
        raise ValueError("Provided node_order does not match edge node universe.")

    index = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)

    i = edges["u"].map(index).to_numpy(dtype=int)
    j = edges["v"].map(index).to_numpy(dtype=int)
    w = edges["edge_weight"].to_numpy(dtype=float)

    # Add each undirected edge in both directions exactly once.
    A = coo_matrix(
        (
            np.concatenate([w, w]),
            (
                np.concatenate([i, j]),
                np.concatenate([j, i]),
            ),
        ),
        shape=(n, n),
        dtype=float,
    ).tocsr()

    weighted_degree = np.asarray(A.sum(axis=1)).ravel()
    unweighted_degree = np.diff(A.indptr).astype(int)

    if np.any(weighted_degree <= 0):
        raise RuntimeError(
            "Network contains zero-degree nodes after construction. "
            "This should not happen when nodes are derived from retained edges."
        )

    # Row-stochastic P = D^{-1} A.
    P = diags(1.0 / weighted_degree, offsets=0, format="csr") @ A
    row_sums = np.asarray(P.sum(axis=1)).ravel()
    max_row_sum_error = float(np.max(np.abs(row_sums - 1.0)))

    if max_row_sum_error > 1e-12:
        raise RuntimeError(
            "Transition matrix is not row-stochastic. "
            f"Max row-sum error={max_row_sum_error:.3e}"
        )

    n_components, component_labels = connected_components(
        A, directed=False, return_labels=True
    )
    component_sizes = np.bincount(component_labels, minlength=n_components)

    return {
        "nodes": nodes,
        "index": index,
        "A": A,
        "P": P.tocsr(),
        "weighted_degree": weighted_degree,
        "unweighted_degree": unweighted_degree,
        "component_labels": component_labels,
        "component_sizes": component_sizes,
        "n_components": int(n_components),
        "max_transition_row_sum_error": max_row_sum_error,
    }


# =============================================================================
# 6. SEED MAPPING
# =============================================================================

def build_seed_vectors(
    seeds: pd.DataFrame,
    network: dict,
) -> tuple[dict[str, np.ndarray], pd.DataFrame, pd.DataFrame]:
    nodes = network["nodes"]
    index = network["index"]
    node_set = set(nodes)

    vectors = {}
    qa_rows = []
    unmapped_parts = []

    for syndrome in SYNDROME_CODES:
        s = seeds[seeds["syndrome_code"].eq(syndrome)].copy()
        s["mapped_to_STRING"] = s["target_id"].isin(node_set)

        mapped = s[s["mapped_to_STRING"]].copy()
        unmapped = s[~s["mapped_to_STRING"]].copy()

        if mapped.empty:
            raise RuntimeError(
                f"No frozen {syndrome} targets map to the STRING network."
            )

        original_weight_sum = float(s["ubiquity_corrected_sum1"].sum())
        mapped_weight_sum = float(
            mapped["ubiquity_corrected_sum1"].sum()
        )

        p0 = np.zeros(len(nodes), dtype=float)
        mapped_norm = (
            mapped["ubiquity_corrected_sum1"].to_numpy(dtype=float)
            / mapped_weight_sum
        )
        idx = np.array(
            [index[g] for g in mapped["target_id"]],
            dtype=int,
        )
        p0[idx] = mapped_norm

        if abs(float(p0.sum()) - 1.0) > MASS_TOL:
            raise RuntimeError(
                f"{syndrome} mapped seed vector does not sum to 1."
            )

        vectors[syndrome] = p0

        mapped_components = network["component_labels"][idx]
        comp_weight = {}
        for c, wt in zip(mapped_components, mapped_norm):
            comp_weight[int(c)] = comp_weight.get(int(c), 0.0) + float(wt)

        largest_seed_component_weight = max(comp_weight.values())
        n_seed_components = len(comp_weight)

        qa_rows.append({
            "syndrome_code": syndrome,
            "frozen_seed_targets": int(len(s)),
            "mapped_seed_targets": int(len(mapped)),
            "unmapped_seed_targets": int(len(unmapped)),
            "target_mapping_fraction": len(mapped) / len(s),
            "frozen_seed_weight_sum": original_weight_sum,
            "mapped_frozen_weight_before_renormalization": mapped_weight_sum,
            "mapped_weight_fraction": (
                mapped_weight_sum / original_weight_sum
                if original_weight_sum > 0 else np.nan
            ),
            "mapped_seed_vector_sum_after_renormalization": float(p0.sum()),
            "seeded_network_components": n_seed_components,
            "largest_seed_component_weight_fraction": (
                largest_seed_component_weight
            ),
        })

        if not unmapped.empty:
            unmapped = unmapped.copy()
            unmapped["reason"] = "not_present_in_thresholded_STRING_network"
            unmapped_parts.append(unmapped)

    qa = pd.DataFrame(qa_rows)
    unmapped_df = (
        pd.concat(unmapped_parts, ignore_index=True)
        if unmapped_parts
        else pd.DataFrame()
    )

    return vectors, qa, unmapped_df


# =============================================================================
# 7. CORRECT RWR
# =============================================================================

def run_rwr(
    P: csr_matrix,
    p0: np.ndarray,
    restart: float,
    tol: float = RWR_TOL_L1,
    max_iter: int = RWR_MAX_ITER,
) -> tuple[np.ndarray, dict]:
    if not (0.0 < restart <= 1.0):
        raise ValueError("restart must be in (0, 1].")

    p0 = np.asarray(p0, dtype=float).reshape(-1)
    if len(p0) != P.shape[0]:
        raise ValueError("p0 length does not match transition matrix.")

    p0_sum = float(p0.sum())
    if abs(p0_sum - 1.0) > MASS_TOL:
        raise ValueError(
            f"p0 must sum to 1; observed {p0_sum:.16g}"
        )

    p = p0.copy()
    converged = False
    last_delta = np.nan
    max_mass_error = abs(float(p.sum()) - 1.0)

    # P is ROW-stochastic, while p is a COLUMN vector.
    # Therefore probability propagation is P.T @ p.
    for iteration in range(1, max_iter + 1):
        propagated = P.T @ p
        p_new = restart * p0 + (1.0 - restart) * propagated

        mass = float(p_new.sum())
        mass_error = abs(mass - 1.0)
        max_mass_error = max(max_mass_error, mass_error)

        if mass_error > MASS_TOL:
            raise RuntimeError(
                "RWR probability mass conservation failed at "
                f"iteration {iteration}: sum={mass:.16g}, "
                f"error={mass_error:.3e}. "
                "Check transition orientation/normalization."
            )

        last_delta = float(np.abs(p_new - p).sum())
        p = p_new

        if last_delta < tol:
            converged = True
            break

    if not converged:
        raise RuntimeError(
            f"RWR failed to converge within {max_iter} iterations. "
            f"Last L1 delta={last_delta:.3e}"
        )

    qa = {
        "restart_probability": float(restart),
        "iterations": int(iteration),
        "converged": True,
        "final_l1_delta": last_delta,
        "final_probability_sum": float(p.sum()),
        "final_mass_error": abs(float(p.sum()) - 1.0),
        "max_mass_error_during_iterations": float(max_mass_error),
    }

    return p, qa


# =============================================================================
# 8. UNIT / INVARIANCE TESTS
# =============================================================================

def toy_network_test() -> pd.DataFrame:
    """
    Analytic 3-node chain:
        A -- B -- C
    unweighted, p0 = [1,0,0], restart r=0.5.

    With row-stochastic P and column-vector update P.T @ p:
        exact solution = [7/12, 1/3, 1/12].
    """
    A = csr_matrix(
        np.array([
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ])
    )
    d = np.asarray(A.sum(axis=1)).ravel()
    P = diags(1.0 / d) @ A

    p0 = np.array([1.0, 0.0, 0.0])
    observed, qa = run_rwr(P.tocsr(), p0, restart=0.5)

    expected = np.array([7/12, 1/3, 1/12], dtype=float)
    abs_error = np.abs(observed - expected)
    max_error = float(abs_error.max())

    if max_error > TOY_TEST_TOL:
        raise RuntimeError(
            "Toy-network analytic RWR test failed. "
            f"Max error={max_error:.3e}"
        )

    return pd.DataFrame({
        "node": ["A", "B", "C"],
        "observed_rwr": observed,
        "expected_rwr": expected,
        "absolute_error": abs_error,
        "test_pass": max_error <= TOY_TEST_TOL,
        "iterations": qa["iterations"],
    })


def determinism_test(
    P: csr_matrix,
    seed_vectors: dict[str, np.ndarray],
) -> pd.DataFrame:
    rows = []
    for syndrome, p0 in seed_vectors.items():
        p1, _ = run_rwr(P, p0, PRIMARY_RESTART)
        p2, _ = run_rwr(P, p0, PRIMARY_RESTART)
        err = float(np.max(np.abs(p1 - p2)))
        rows.append({
            "syndrome_code": syndrome,
            "max_absolute_difference_two_identical_runs": err,
            "tolerance": DETERMINISM_TOL,
            "test_pass": err <= DETERMINISM_TOL,
        })

    out = pd.DataFrame(rows)
    if not out["test_pass"].all():
        raise RuntimeError("Determinism QA failed.")
    return out


def node_order_invariance_test(
    edges: pd.DataFrame,
    network: dict,
    seed_vectors: dict[str, np.ndarray],
) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)

    old_nodes = np.array(network["nodes"], dtype=object)
    perm = rng.permutation(len(old_nodes))
    shuffled_nodes = old_nodes[perm].tolist()

    shuffled_network = build_network(edges, node_order=shuffled_nodes)

    rows = []
    for syndrome, p0_old in seed_vectors.items():
        reference, _ = run_rwr(
            network["P"], p0_old, PRIMARY_RESTART
        )

        # New coordinate k corresponds to old coordinate perm[k].
        p0_shuffled = p0_old[perm]
        shuffled_result, _ = run_rwr(
            shuffled_network["P"], p0_shuffled, PRIMARY_RESTART
        )

        # Map shuffled result back to old node order.
        back = np.empty_like(shuffled_result)
        back[perm] = shuffled_result

        max_err = float(np.max(np.abs(reference - back)))
        l1_err = float(np.abs(reference - back).sum())

        rows.append({
            "syndrome_code": syndrome,
            "max_absolute_difference": max_err,
            "l1_difference": l1_err,
            "tolerance": NODE_ORDER_INVARIANCE_TOL,
            "test_pass": max_err <= NODE_ORDER_INVARIANCE_TOL,
        })

    out = pd.DataFrame(rows)
    if not out["test_pass"].all():
        raise RuntimeError("Node-order invariance QA failed.")
    return out


# =============================================================================
# 9. RUN ALL RESTARTS AND BUILD OUTPUT TABLES
# =============================================================================

def run_all_rwr(
    network: dict,
    seed_vectors: dict[str, np.ndarray],
    seeds: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    restarts = [PRIMARY_RESTART] + [
        r for r in SENSITIVITY_RESTARTS if r != PRIMARY_RESTART
    ]

    nodes = network["nodes"]
    score_parts = []
    qa_rows = []

    # Fast lookup for original frozen weights.
    original_weight = {
        (r["syndrome_code"], r["target_id"]):
            float(r["ubiquity_corrected_sum1"])
        for _, r in seeds.iterrows()
    }

    for restart in restarts:
        label = restart_label(restart)

        for syndrome in SYNDROME_CODES:
            p0 = seed_vectors[syndrome]
            p, qa = run_rwr(
                network["P"],
                p0,
                restart=restart,
            )

            qa_rows.append({
                "syndrome_code": syndrome,
                "restart_label": label,
                **qa,
            })

            part = pd.DataFrame({
                "gene": nodes,
                "syndrome_code": syndrome,
                "restart_probability": restart,
                "restart_label": label,
                "seed_weight_mapped_sum1": p0,
                "rwr_score": p,
                "network_degree_unweighted": network["unweighted_degree"],
                "network_degree_weighted": network["weighted_degree"],
                "component_id": network["component_labels"],
            })

            part["seed_weight_frozen_original"] = [
                original_weight.get((syndrome, g), 0.0)
                for g in nodes
            ]
            part["is_seed"] = (
                part["seed_weight_mapped_sum1"] > 0
            ).astype(int)
            part["rwr_rank"] = (
                part["rwr_score"]
                .rank(method="min", ascending=False)
                .astype(int)
            )
            score_parts.append(part)

    scores = pd.concat(score_parts, ignore_index=True)
    convergence = pd.DataFrame(qa_rows)

    return scores, convergence


# =============================================================================
# 10. SIMILARITY / RE-COLLAPSE DIAGNOSTICS
# =============================================================================

def pairwise_metrics(
    matrix: pd.DataFrame,
    stage: str,
    restart: float | None,
) -> pd.DataFrame:
    pairs = [("PDOL", "PHOL"), ("PDOL", "WHIL"), ("PHOL", "WHIL")]
    rows = []

    for a, b in pairs:
        va = matrix[a].to_numpy(dtype=float)
        vb = matrix[b].to_numpy(dtype=float)

        rho = float(spearmanr(va, vb).statistic)
        pearson = (
            float(np.corrcoef(va, vb)[0, 1])
            if np.std(va) > 0 and np.std(vb) > 0
            else np.nan
        )

        positive_a = set(matrix.index[matrix[a] > 0])
        positive_b = set(matrix.index[matrix[b] > 0])
        union = positive_a | positive_b
        inter = positive_a & positive_b

        rows.append({
            "stage": stage,
            "restart_probability": (
                np.nan if restart is None else float(restart)
            ),
            "syndrome_A": a,
            "syndrome_B": b,
            "spearman": rho,
            "pearson": pearson,
            "positive_overlap_n": len(inter),
            "positive_union_n": len(union),
            "positive_jaccard": (
                len(inter) / len(union) if union else np.nan
            ),
        })

    return pd.DataFrame(rows)


def build_pairwise_similarity(
    scores: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    # Mapped seed vectors before diffusion.
    primary = scores[
        scores["restart_probability"].eq(PRIMARY_RESTART)
    ].copy()
    seed_matrix = primary.pivot(
        index="gene",
        columns="syndrome_code",
        values="seed_weight_mapped_sum1",
    ).fillna(0)
    rows.append(
        pairwise_metrics(
            seed_matrix,
            stage="mapped_seed_before_RWR",
            restart=None,
        )
    )

    for r in [PRIMARY_RESTART] + SENSITIVITY_RESTARTS:
        x = scores[scores["restart_probability"].eq(r)]
        mat = x.pivot(
            index="gene",
            columns="syndrome_code",
            values="rwr_score",
        ).fillna(0)
        rows.append(
            pairwise_metrics(
                mat,
                stage="RWR",
                restart=float(r),
            )
        )

    return pd.concat(rows, ignore_index=True)


def build_topk_overlap(scores: pd.DataFrame) -> pd.DataFrame:
    pairs = [("PDOL", "PHOL"), ("PDOL", "WHIL"), ("PHOL", "WHIL")]
    rows = []

    # Seed-vector top-k before RWR.
    primary = scores[
        scores["restart_probability"].eq(PRIMARY_RESTART)
    ].copy()

    for stage, value_col, restart, frame in [
        (
            "mapped_seed_before_RWR",
            "seed_weight_mapped_sum1",
            None,
            primary,
        ),
    ]:
        for k in TOP_K_VALUES:
            sets = {}
            for syndrome in SYNDROME_CODES:
                z = frame[
                    frame["syndrome_code"].eq(syndrome)
                    & frame[value_col].gt(0)
                ].nlargest(k, value_col)
                sets[syndrome] = set(z["gene"])

            for a, b in pairs:
                inter = sets[a] & sets[b]
                union = sets[a] | sets[b]
                rows.append({
                    "stage": stage,
                    "restart_probability": (
                        np.nan if restart is None else float(restart)
                    ),
                    "top_k": k,
                    "syndrome_A": a,
                    "syndrome_B": b,
                    "overlap_n": len(inter),
                    "union_n": len(union),
                    "jaccard": len(inter) / len(union) if union else np.nan,
                })

    # RWR top-k.
    for r in [PRIMARY_RESTART] + SENSITIVITY_RESTARTS:
        frame = scores[scores["restart_probability"].eq(r)]
        for k in TOP_K_VALUES:
            sets = {}
            for syndrome in SYNDROME_CODES:
                z = frame[
                    frame["syndrome_code"].eq(syndrome)
                ].nlargest(k, "rwr_score")
                sets[syndrome] = set(z["gene"])

            for a, b in pairs:
                inter = sets[a] & sets[b]
                union = sets[a] | sets[b]
                rows.append({
                    "stage": "RWR",
                    "restart_probability": float(r),
                    "top_k": k,
                    "syndrome_A": a,
                    "syndrome_B": b,
                    "overlap_n": len(inter),
                    "union_n": len(union),
                    "jaccard": len(inter) / len(union) if union else np.nan,
                })

    return pd.DataFrame(rows)


def restart_sensitivity_stability(
    scores: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    primary = scores[
        scores["restart_probability"].eq(PRIMARY_RESTART)
    ][["gene", "syndrome_code", "rwr_score"]].rename(
        columns={"rwr_score": "primary_score"}
    )

    for r in SENSITIVITY_RESTARTS:
        sens = scores[
            scores["restart_probability"].eq(r)
        ][["gene", "syndrome_code", "rwr_score"]].rename(
            columns={"rwr_score": "sensitivity_score"}
        )

        merged = primary.merge(
            sens,
            on=["gene", "syndrome_code"],
            how="inner",
        )

        for syndrome in SYNDROME_CODES:
            x = merged[merged["syndrome_code"].eq(syndrome)]
            rho = float(
                spearmanr(
                    x["primary_score"],
                    x["sensitivity_score"],
                ).statistic
            )

            top_primary = set(
                x.nlargest(200, "primary_score")["gene"]
            )
            top_sens = set(
                x.nlargest(200, "sensitivity_score")["gene"]
            )
            inter = top_primary & top_sens
            union = top_primary | top_sens

            rows.append({
                "syndrome_code": syndrome,
                "primary_restart": PRIMARY_RESTART,
                "sensitivity_restart": float(r),
                "spearman_all_network_nodes": rho,
                "top200_overlap_n": len(inter),
                "top200_jaccard": len(inter) / len(union),
            })

    return pd.DataFrame(rows)


# =============================================================================
# 11. NETWORK / NULL-MATCHING REFERENCE
# =============================================================================

def load_target_ubiquity_optional() -> pd.DataFrame:
    if TARGET_UBIQUITY_FILE is None or not TARGET_UBIQUITY_FILE.exists():
        warnings.warn(
            "TARGET_UBIQUITY_FILE not found. "
            "RWR itself is unaffected, but later matched-null reference "
            "will lack target-database ubiquity variables."
        )
        return pd.DataFrame(
            columns=[
                "target_id",
                "target_herb_degree",
                "N_database_herbs",
                "target_idf",
            ]
        )

    x = pd.read_csv(TARGET_UBIQUITY_FILE)
    require_columns(
        x,
        ["target_id", "target_herb_degree", "target_idf"],
        TARGET_UBIQUITY_FILE.name,
    )
    x = x.copy()
    x["target_id"] = x["target_id"].map(norm_gene)
    return x.drop_duplicates("target_id")


def build_network_node_stats(
    network: dict,
    target_ubiquity: pd.DataFrame,
) -> pd.DataFrame:
    component_sizes = network["component_sizes"]
    labels = network["component_labels"]

    out = pd.DataFrame({
        "gene": network["nodes"],
        "network_degree_unweighted": network["unweighted_degree"],
        "network_degree_weighted": network["weighted_degree"],
        "component_id": labels,
        "component_size": component_sizes[labels],
    })

    if not target_ubiquity.empty:
        t = target_ubiquity.rename(columns={"target_id": "gene"})
        out = out.merge(t, on="gene", how="left")

    out["in_target_database_reference"] = (
        out.get("target_herb_degree", pd.Series(index=out.index, dtype=float))
        .notna()
        .astype(int)
    )

    return out


# =============================================================================
# 12. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("STRING RWR FROZEN-UPSTREAM v1.0")
    print("No CAP transcriptomic data are read or used.")
    print("=" * 78)

    print("\n[1/10] Loading frozen LOPO target seeds...")
    seeds = load_frozen_seeds(FROZEN_SEED_FILE)
    print(
        seeds.groupby("syndrome_code")
        .agg(
            frozen_targets=("target_id", "nunique"),
            weight_sum=("ubiquity_corrected_sum1", "sum"),
        )
        .to_string()
    )

    print("\n[2/10] Loading and cleaning STRING network...")
    edges, edge_stats, string_source_files = load_string_edges()
    print(
        f"    retained undirected edges: "
        f"{edge_stats['clean_high_confidence_undirected_edges']:,}"
    )

    print("\n[3/10] Building weighted row-stochastic transition matrix...")
    network = build_network(edges)
    print(f"    network nodes:      {len(network['nodes']):,}")
    print(f"    components:         {network['n_components']:,}")
    print(
        f"    largest component: "
        f"{int(network['component_sizes'].max()):,}"
    )
    print(
        f"    max row-sum error:  "
        f"{network['max_transition_row_sum_error']:.3e}"
    )

    print("\n[4/10] Mapping frozen syndrome seeds into STRING...")
    seed_vectors, mapping_qa, unmapped = build_seed_vectors(
        seeds, network
    )
    print(mapping_qa.to_string(index=False))

    print("\n[5/10] Running analytic toy-network QA...")
    toy_qa = toy_network_test()
    print(
        f"    max toy error: "
        f"{toy_qa['absolute_error'].max():.3e}"
    )

    print("\n[6/10] Running primary and sensitivity RWR...")
    scores, convergence_qa = run_all_rwr(
        network, seed_vectors, seeds
    )
    print(
        convergence_qa[
            [
                "syndrome_code",
                "restart_probability",
                "iterations",
                "final_l1_delta",
                "final_probability_sum",
                "max_mass_error_during_iterations",
            ]
        ].to_string(index=False)
    )

    print("\n[7/10] Running determinism and node-order invariance QA...")
    determinism_qa = determinism_test(
        network["P"], seed_vectors
    )
    order_qa = node_order_invariance_test(
        edges, network, seed_vectors
    )
    print("    determinism: PASS")
    print("    node-order invariance: PASS")

    print("\n[8/10] Quantifying pre/post-diffusion syndrome similarity...")
    similarity = build_pairwise_similarity(scores)
    topk = build_topk_overlap(scores)
    sensitivity = restart_sensitivity_stability(scores)

    primary_similarity = similarity[
        similarity["stage"].eq("RWR")
        & similarity["restart_probability"].eq(PRIMARY_RESTART)
    ]
    print(primary_similarity.to_string(index=False))

    print("\n[9/10] Building network and later-null reference tables...")
    target_ubiquity = load_target_ubiquity_optional()
    node_stats = build_network_node_stats(
        network, target_ubiquity
    )

    # Primary r=0.5 scores.
    primary_scores = scores[
        scores["restart_probability"].eq(PRIMARY_RESTART)
    ].copy()

    # Compact top-ranked output.
    top_primary = (
        primary_scores
        .sort_values(
            ["syndrome_code", "rwr_score"],
            ascending=[True, False],
        )
        .groupby("syndrome_code", group_keys=False)
        .head(500)
        .copy()
    )

    print("\n[10/10] Writing outputs...")

    # Core outputs.
    primary_scores.to_csv(
        OUTPUT_DIR / "RWR_PRIMARY_r0p50_scores_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    scores.to_csv(
        OUTPUT_DIR / "RWR_all_restart_scores_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    top_primary.to_csv(
        OUTPUT_DIR / "RWR_PRIMARY_top500_per_syndrome_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # QA.
    mapping_qa.to_csv(
        OUTPUT_DIR / "RWR_seed_mapping_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    unmapped.to_csv(
        OUTPUT_DIR / "RWR_unmapped_frozen_seeds_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    convergence_qa.to_csv(
        OUTPUT_DIR / "RWR_convergence_and_mass_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    toy_qa.to_csv(
        OUTPUT_DIR / "RWR_toy_network_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    determinism_qa.to_csv(
        OUTPUT_DIR / "RWR_determinism_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    order_qa.to_csv(
        OUTPUT_DIR / "RWR_node_order_invariance_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Diagnostics.
    similarity.to_csv(
        OUTPUT_DIR / "RWR_pairwise_similarity_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    topk.to_csv(
        OUTPUT_DIR / "RWR_topk_overlap_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    sensitivity.to_csv(
        OUTPUT_DIR / "RWR_restart_sensitivity_stability_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Network provenance and later null matching.
    edges.to_csv(
        OUTPUT_DIR / "STRING_clean_high_confidence_edges_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    node_stats.to_csv(
        OUTPUT_DIR / "STRING_network_node_stats_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    null_ref = node_stats[
        node_stats["in_target_database_reference"].eq(1)
    ].copy()
    null_ref.to_csv(
        OUTPUT_DIR / "LATER_null_matching_reference_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Human-readable QA workbook.
    if WRITE_EXCEL_SUMMARY:
        with pd.ExcelWriter(
            OUTPUT_DIR / "STRING_RWR_v1_0_QA_summary.xlsx",
            engine="openpyxl",
        ) as writer:
            mapping_qa.to_excel(
                writer, sheet_name="seed_mapping", index=False
            )
            convergence_qa.to_excel(
                writer, sheet_name="convergence_mass", index=False
            )
            similarity.to_excel(
                writer, sheet_name="pairwise_similarity", index=False
            )
            topk.to_excel(
                writer, sheet_name="topk_overlap", index=False
            )
            sensitivity.to_excel(
                writer, sheet_name="restart_sensitivity", index=False
            )
            toy_qa.to_excel(
                writer, sheet_name="toy_test", index=False
            )
            determinism_qa.to_excel(
                writer, sheet_name="determinism", index=False
            )
            order_qa.to_excel(
                writer, sheet_name="node_order", index=False
            )
            top_primary.to_excel(
                writer, sheet_name="primary_top500", index=False
            )

    # Metadata.
    metadata = {
        "analysis_name": "STRING RWR frozen-upstream v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"
        ),
        "status": "FREEZE_CANDIDATE_REQUIRES_OUTPUT_REVIEW",
        "CAP_transcriptomics_used": False,
        "GAT_used": False,
        "upstream_seed_file": str(FROZEN_SEED_FILE),
        "upstream_seed_sha256": sha256_file(FROZEN_SEED_FILE),
        "target_ubiquity_file": (
            str(TARGET_UBIQUITY_FILE)
            if TARGET_UBIQUITY_FILE is not None else None
        ),
        "target_ubiquity_sha256": (
            safe_sha256(TARGET_UBIQUITY_FILE)
            if TARGET_UBIQUITY_FILE is not None else None
        ),
        "STRING": {
            "declared_version": STRING_VERSION_DECLARED,
            "species": STRING_SPECIES,
            "taxonomy_id": STRING_TAXONOMY_ID,
            "network_semantics": STRING_NETWORK_SEMANTICS,
            "input_mode": PPI_INPUT_MODE,
            "combined_score_threshold": (
                STRING_COMBINED_SCORE_THRESHOLD
            ),
            "edge_weighting": (
                "combined_score"
                if USE_STRING_EDGE_WEIGHTS else "unweighted"
            ),
            "source_files": string_source_files,
            "edge_cleaning_stats": edge_stats,
        },
        "network": {
            "nodes": int(len(network["nodes"])),
            "undirected_edges": int(len(edges)),
            "connected_components": int(network["n_components"]),
            "largest_component_nodes": int(
                network["component_sizes"].max()
            ),
            "transition_definition": "P = D^{-1} A (row-stochastic)",
            "propagation_for_column_vector": "P.T @ p",
            "max_transition_row_sum_error": (
                network["max_transition_row_sum_error"]
            ),
        },
        "RWR": {
            "update": (
                "p_next = r*p0 + (1-r)*(P.T @ p)"
            ),
            "primary_restart_probability": PRIMARY_RESTART,
            "sensitivity_restart_probabilities": SENSITIVITY_RESTARTS,
            "convergence_tolerance_L1": RWR_TOL_L1,
            "max_iterations": RWR_MAX_ITER,
            "mass_tolerance": MASS_TOL,
            "mapped_seed_weights_renormalized_to_sum1": True,
        },
        "QA": {
            "toy_network_analytic_test": "PASS",
            "determinism_test": (
                "PASS" if determinism_qa["test_pass"].all()
                else "FAIL"
            ),
            "node_order_invariance_test": (
                "PASS" if order_qa["test_pass"].all()
                else "FAIL"
            ),
            "max_rwr_mass_error": float(
                convergence_qa[
                    "max_mass_error_during_iterations"
                ].max()
            ),
        },
        "matched_null_policy": {
            "generated_in_this_script": False,
            "reason": (
                "Matched-null sampling is deferred until the downstream "
                "disease-alignment statistic is frozen. This script exports "
                "network degree and target-database ubiquity variables for "
                "that later null."
            ),
            "later_matching_variables": [
                "network_degree_unweighted",
                "network_degree_weighted",
                "target_herb_degree",
                "target_idf",
                "seed-set size",
                "observed seed-weight distribution",
            ],
        },
        "freeze_policy": (
            "Do not freeze RWR merely because the code passes QA. "
            "First inspect mapping coverage, pre/post-RWR syndrome collapse, "
            "restart sensitivity, and network provenance. "
            "Do not use CAP alignment to choose r or alter the network."
        ),
    }

    with open(
        OUTPUT_DIR / "STRING_RWR_run_metadata_v1_0.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    statement = f"""STRING RWR v1.0 — FREEZE CANDIDATE, NOT YET FINAL

Upstream seed:
{FROZEN_SEED_FILE}

Primary network:
STRING {STRING_VERSION_DECLARED}, {STRING_SPECIES},
combined_score >= {STRING_COMBINED_SCORE_THRESHOLD},
{STRING_NETWORK_SEMANTICS}.

Primary RWR:
P = D^(-1) A
p_next = r*p0 + (1-r)*(P.T @ p)
restart r = {PRIMARY_RESTART}

Sensitivity:
r = {", ".join(map(str, SENSITIVITY_RESTARTS))}

Disease data used upstream:
NO

Required before final RWR freeze:
1. Confirm STRING source-file provenance / SHA256.
2. Confirm seed mapping and retained seed-weight fractions.
3. Confirm probability mass = 1 and convergence.
4. Confirm analytic toy test, determinism, and node-order invariance.
5. Inspect pre-RWR versus post-RWR syndrome similarity.
6. Inspect r=0.3/0.5/0.7 stability.
7. Do NOT use CAP alignment to change the above settings.

Matched-null testing is intentionally deferred to the disease-alignment stage.
"""
    (OUTPUT_DIR / "STRING_RWR_FREEZE_CANDIDATE_v1_0.txt").write_text(
        statement, encoding="utf-8"
    )

    print("\n" + "=" * 78)
    print("RWR computation and QA completed.")
    print("STATUS: FREEZE CANDIDATE — inspect collapse/provenance before final freeze.")
    print("=" * 78)
    print(f"\nOutputs:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
