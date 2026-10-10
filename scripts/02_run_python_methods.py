#!/usr/bin/env python3
"""Run Figure 4 Python-native analyses and retain every completed artifact.

This script never uses pathology labels as analysis input.  They remain in the
canonical input folder exclusively for later external validation.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None)


def append_ledger(path: Path, **row) -> None:
    frame = pd.DataFrame([row])
    exists = path.exists()
    frame.to_csv(path, sep="\t", mode="a", header=not exists, index=False)


def load_section(section_dir: Path):
    import anndata as ad

    spots = pd.read_csv(section_dir / "spots.tsv.gz", sep="\t", compression="gzip")
    if "spot_id" not in spots:
        raise ValueError(f"{section_dir}: spots.tsv.gz has no spot_id column")
    spots["spot_id"] = spots["spot_id"].astype(str)
    if (section_dir / "counts.tsv.gz").exists():
        counts = pd.read_csv(section_dir / "counts.tsv.gz", sep="\t", compression="gzip", index_col=0)
        counts.index = counts.index.astype(str)
        spots = spots.set_index("spot_id").loc[counts.index].reset_index()
        x = sparse.csr_matrix(counts.to_numpy(dtype=np.float32, copy=False))
        genes = counts.columns.astype(str).to_numpy()
    else:
        matrix = io.mmread(section_dir / "matrix.mtx.gz").tocsr().T
        barcodes = pd.read_csv(section_dir / "barcodes.tsv.gz", sep="\t", header=None).iloc[:, 0].astype(str).to_numpy()
        features = pd.read_csv(section_dir / "features.tsv.gz", sep="\t", header=None)
        genes = features.iloc[:, 1].astype(str).to_numpy() if features.shape[1] > 1 else features.iloc[:, 0].astype(str).to_numpy()
        spots = spots.set_index("spot_id").loc[barcodes].reset_index()
        x = matrix.astype(np.float32)
    if x.shape[0] != len(spots):
        raise ValueError(f"{section_dir}: matrix and spot rows disagree")
    adata = ad.AnnData(X=x, obs=spots.set_index("spot_id"), var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    x_col = next((c for c in ["x", "pixel_col", "array_col"] if c in adata.obs), None)
    y_col = next((c for c in ["y", "pixel_row", "array_row"] if c in adata.obs), None)
    if x_col is None or y_col is None:
        raise ValueError(f"{section_dir}: no usable coordinate columns")
    adata.obsm["spatial"] = adata.obs[[x_col, y_col]].astype(float).to_numpy()
    adata.uns["sparq_input"] = {"section_dir": str(section_dir), "coordinate_columns": [x_col, y_col]}
    return adata


def normalize_for_scanpy(adata, seed: int):
    import scanpy as sc

    out = adata.copy()
    out.layers["counts"] = out.X.copy()
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    n_top = min(2000, max(50, out.n_vars - 1))
    sc.pp.highly_variable_genes(out, n_top_genes=n_top, flavor="seurat")
    if int(out.var["highly_variable"].sum()) < 2:
        out.var["highly_variable"] = True
    sc.pp.pca(out, n_comps=min(30, out.n_obs - 1, int(out.var["highly_variable"].sum())), use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(out, n_neighbors=min(15, max(2, out.n_obs - 1)), n_pcs=min(30, out.obsm["X_pca"].shape[1]), random_state=seed)
    return out


def leiden(adata, seed: int):
    import scanpy as sc

    sc.tl.leiden(adata, resolution=1.0, key_added="leiden", random_state=seed, flavor="igraph", n_iterations=2)
    return pd.DataFrame({"spot_id": adata.obs_names.astype(str), "label": adata.obs["leiden"].astype(str).to_numpy()})


def save_partition(table: pd.DataFrame, folder: Path, name: str) -> None:
    write_tsv(table, folder / f"{name}_labels.tsv.gz")
    table.rename(columns={"label": "value"}).assign(result_id=name).to_csv(folder / f"{name}_standardized.tsv.gz", sep="\t", index=False, compression="gzip")


def run_squidpy(adata, folder: Path, tag: str, seed: int) -> None:
    import squidpy as sq

    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=min(6, max(2, adata.n_obs - 1)))
    result = sq.gr.nhood_enrichment(adata, cluster_key="leiden", n_perms=100, seed=seed, copy=True)
    if result is None:
        result = adata.uns.get("leiden_nhood_enrichment")
    z, counts = result
    labels = list(map(str, adata.obs["leiden"].cat.categories if hasattr(adata.obs["leiden"], "cat") else sorted(adata.obs["leiden"].unique())))
    rows = []
    for i, left in enumerate(labels):
        for j, right in enumerate(labels):
            rows.append({"source_label": left, "target_label": right, "z_score": float(z[i, j]), "observed_count": float(counts[i, j]), "edge_id": f"{left}--{right}", "selected": bool(z[i, j] > 0)})
    table = pd.DataFrame(rows)
    write_tsv(table, folder / f"{tag}_neighborhood_enrichment.tsv.gz")
    write_tsv(table.loc[table.selected, ["edge_id", "z_score", "observed_count", "selected"]], folder / f"{tag}_positive_edges.tsv.gz")
    table.to_pickle(folder / f"{tag}_neighborhood_enrichment.pkl")


def run_spatialdm(adata, folder: Path, tag: str, cpus: int) -> None:
    import spatialdm as sdm

    out = adata.copy()
    # Keep computational cost bounded but preserve the selected gene universe.
    if "highly_variable" in out.var and int(out.var.highly_variable.sum()) >= 10:
        out = out[:, out.var.highly_variable].copy()
    sdm.weight_matrix(out, n_neighbors=min(6, max(2, out.n_obs - 1)))
    sdm.extract_lr(out, species="human")
    sdm.spatialdm_global(out, method="z-score", n_perm=100, nproc=max(1, cpus))
    sdm.sig_pairs(out, method="z-score", fdr=True, threshold=0.10)
    global_result = out.uns.get("global_res")
    if global_result is None:
        raise RuntimeError("SpatialDM completed without adata.uns['global_res']")
    table = global_result.copy().reset_index().rename(columns={"index": "interaction_id"})
    score_col = next((c for c in ["z", "z_score", "stat", "score"] if c in table.columns), None)
    if score_col is None:
        numeric = table.select_dtypes(include=[np.number]).columns.tolist()
        if not numeric:
            raise RuntimeError("SpatialDM global result contains no numeric ranking column")
        score_col = numeric[0]
    fdr_col = next((c for c in ["fdr", "FDR", "padj", "p_adj"] if c in table.columns), None)
    table["ranking_score"] = pd.to_numeric(table[score_col], errors="coerce")
    table["selected"] = (pd.to_numeric(table[fdr_col], errors="coerce") <= 0.10) if fdr_col else np.isfinite(table["ranking_score"])
    table = table.sort_values("ranking_score", ascending=False, na_position="last").reset_index(drop=True)
    table["rank"] = np.arange(1, len(table) + 1)
    write_tsv(table, folder / f"{tag}_global_interactions.tsv.gz")
    write_tsv(table[["interaction_id", "ranking_score", "rank", "selected"]], folder / f"{tag}_ranked_interactions.tsv.gz")
    out.write_h5ad(folder / f"{tag}_spatialdm.h5ad", compression="gzip")


def run_spagcn_section(section_dir: Path, folder: Path, tag: str, seed: int) -> None:
    import SpaGCN as spg

    adata = load_section(section_dir)
    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]
    adjacency = spg.calculate_adj_matrix(x=x, y=y, histology=False)
    l_value = spg.search_l(0.5, adjacency, start=0.01, end=1000, tol=0.01, max_run=100)
    model = spg.SpaGCN()
    model.set_l(l_value)
    model.train(adata.X, adjacency, init_spa=True, init="louvain", res=0.7, tol=5e-3, lr=0.05, max_epochs=200)
    labels, probabilities = model.predict()
    table = pd.DataFrame({"spot_id": adata.obs_names.astype(str), "label": pd.Series(labels).astype(str), "max_probability": np.max(np.asarray(probabilities), axis=1)})
    save_partition(table[["spot_id", "label"]], folder, tag)
    write_tsv(table, folder / f"{tag}_probabilities.tsv.gz")
    np.save(folder / f"{tag}_adjacency.npy", adjacency)
    np.save(folder / f"{tag}_probabilities.npy", np.asarray(probabilities))


def method_folder(run_dir: Path, method: str, section: str) -> Path:
    path = run_dir / "raw" / "python" / method / section
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_core(run_dir: Path, iterations: int, retention: float, seed: int, cpus: int) -> None:
    import anndata

    canonical = sorted((run_dir / "inputs" / "canonical").glob("*/section.json"))
    ledger = run_dir / "manifests" / "python_method_ledger.tsv"
    for section_json in canonical:
        section_dir = section_json.parent
        section = section_dir.name
        try:
            raw = load_section(section_dir)
            write_tsv(raw.obs.reset_index(), run_dir / "standardized" / "inputs" / f"{section}_spots.tsv.gz")
            raw.write_h5ad(run_dir / "raw" / "python" / "input_objects" / f"{section}.h5ad", compression="gzip")
            rng = np.random.default_rng(seed + sum(map(ord, section)))
            plans = [("baseline", np.arange(raw.n_obs, dtype=int), seed)]
            perturb_dir = run_dir / "inputs" / "perturbations" / section
            perturb_dir.mkdir(parents=True, exist_ok=True)
            for iteration in range(1, iterations + 1):
                retained = np.sort(rng.choice(raw.n_obs, size=max(2, int(np.floor(raw.n_obs * retention))), replace=False))
                removed = np.setdiff1d(np.arange(raw.n_obs), retained)
                write_tsv(pd.DataFrame({"spot_id": raw.obs_names[retained].astype(str)}), perturb_dir / f"iteration_{iteration:03d}_retained_spot_ids.tsv.gz")
                write_tsv(pd.DataFrame({"spot_id": raw.obs_names[removed].astype(str)}), perturb_dir / f"iteration_{iteration:03d}_removed_spot_ids.tsv.gz")
                plans.append((f"iteration_{iteration:03d}", retained, seed + iteration))
            for tag, index, local_seed in plans:
                adata = raw[index].copy()
                started = now()
                try:
                    processed = normalize_for_scanpy(adata, local_seed)
                    labels = leiden(processed, local_seed)
                    p_folder = method_folder(run_dir, "pca_leiden", section)
                    save_partition(labels, p_folder, tag)
                    processed.write_h5ad(p_folder / f"{tag}_object.h5ad", compression="gzip")
                    q_folder = method_folder(run_dir, "squidpy_neighborhood", section)
                    run_squidpy(processed.copy(), q_folder, tag, local_seed)
                    d_folder = method_folder(run_dir, "spatialdm", section)
                    run_spatialdm(processed.copy(), d_folder, tag, cpus)
                    append_ledger(ledger, section=section, iteration=tag, method="pca_leiden+squidpy+spatialdm", status="ok", started_at=started, finished_at=now(), error="")
                except Exception as exc:
                    append_ledger(ledger, section=section, iteration=tag, method="pca_leiden+squidpy+spatialdm", status="failed", started_at=started, finished_at=now(), error=f"{type(exc).__name__}: {exc}")
                    (run_dir / "logs" / f"python_{section}_{tag}.traceback.txt").write_text(traceback.format_exc())
        except Exception as exc:
            append_ledger(ledger, section=section, iteration="section_setup", method="all_python", status="failed", started_at=now(), finished_at=now(), error=f"{type(exc).__name__}: {exc}")
            (run_dir / "logs" / f"python_{section}_setup.traceback.txt").write_text(traceback.format_exc())


def run_spagcn(run_dir: Path, iterations: int, retention: float, seed: int) -> None:
    ledger = run_dir / "manifests" / "spagcn_method_ledger.tsv"
    for section_json in sorted((run_dir / "inputs" / "canonical").glob("*/section.json")):
        section_dir = section_json.parent
        section = section_dir.name
        raw = load_section(section_dir)
        rng = np.random.default_rng(seed + sum(map(ord, section)))
        plans = [("baseline", np.arange(raw.n_obs, dtype=int), seed)]
        for iteration in range(1, iterations + 1):
            retained_path = run_dir / "inputs" / "perturbations" / section / f"iteration_{iteration:03d}_retained_spot_ids.tsv.gz"
            retained_ids = pd.read_csv(retained_path, sep="\t", compression="gzip")["spot_id"].astype(str)
            index = raw.obs_names.get_indexer(retained_ids)
            if (index < 0).any():
                raise RuntimeError(f"{section}: saved retained IDs do not match input")
            plans.append((f"iteration_{iteration:03d}", index, seed + iteration))
        for tag, index, local_seed in plans:
            started = now()
            try:
                temporary = run_dir / "raw" / "python" / "spagcn" / section / "_input"
                temporary.mkdir(parents=True, exist_ok=True)
                subset = raw[index].copy()
                subset.write_h5ad(temporary / f"{tag}.h5ad", compression="gzip")
                # SpaGCN loader reads the canonical source; an explicit index table keeps its exact input auditable.
                write_tsv(pd.DataFrame({"spot_id": subset.obs_names.astype(str)}), temporary / f"{tag}_spot_ids.tsv.gz")
                # Run directly from a transient canonical-like section by monkeypatching the loader input.
                original_loader = load_section
                def _loaded(_: Path): return subset.copy()
                globals()["load_section"] = _loaded
                try:
                    run_spagcn_section(section_dir, method_folder(run_dir, "spagcn", section), tag, local_seed)
                finally:
                    globals()["load_section"] = original_loader
                append_ledger(ledger, section=section, iteration=tag, method="spagcn", status="ok", started_at=started, finished_at=now(), error="")
            except Exception as exc:
                append_ledger(ledger, section=section, iteration=tag, method="spagcn", status="failed", started_at=started, finished_at=now(), error=f"{type(exc).__name__}: {exc}")
                (run_dir / "logs" / f"spagcn_{section}_{tag}.traceback.txt").write_text(traceback.format_exc())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--python-env", type=Path, required=True)
    parser.add_argument("--spagcn-env", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--retention", type=float, default=0.75)
    parser.add_argument("--seed", type=int, default=20261010)
    parser.add_argument("--only", choices=["core", "spagcn"], default="core")
    args = parser.parse_args()
    if not (1 <= args.iterations and 0 < args.retention < 1):
        raise ValueError("iterations must be positive and retention must lie in (0, 1)")
    args.run_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"started_at": now(), "python": sys.version, "platform": platform.platform(), "iterations": args.iterations, "retention": args.retention, "seed": args.seed, "mode": args.only}
    (args.run_dir / "manifests" / f"python_{args.only}_run.json").write_text(json.dumps(metadata, indent=2))
    if args.only == "core":
        run_core(args.run_dir, args.iterations, args.retention, args.seed, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
        command = ["conda", "run", "--no-capture-output", "-p", str(args.spagcn_env), "python", str(Path(__file__).resolve()), "--run-dir", str(args.run_dir), "--python-env", str(args.python_env), "--spagcn-env", str(args.spagcn_env), "--iterations", str(args.iterations), "--retention", str(args.retention), "--seed", str(args.seed), "--only", "spagcn"]
        result = subprocess.run(command, capture_output=True, text=True)
        (args.run_dir / "logs" / "spagcn_subprocess.out").write_text(result.stdout)
        (args.run_dir / "logs" / "spagcn_subprocess.err").write_text(result.stderr)
        if result.returncode:
            append_ledger(args.run_dir / "manifests" / "spagcn_method_ledger.tsv", section="all", iteration="subprocess", method="spagcn", status="failed", started_at=now(), finished_at=now(), error=f"Subprocess exit status {result.returncode}")
    else:
        run_spagcn(args.run_dir, args.iterations, args.retention, args.seed)
    (args.run_dir / "manifests" / f"python_{args.only}_completed.txt").write_text(now() + "\n")


if __name__ == "__main__":
    main()
