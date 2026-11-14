"""EFHC1 and CACNA1E in first-trimester thalamic and midbrain neurons.

Uses the first-trimester developing human brain atlas (Braun et al., Science 2023;
CELLxGENE dataset 13149914-ea03-4f01-8bf6-b793b667127b) through the CELLxGENE
Census (pinned release below). Neurons are taken from the dissections labeled
"dorsal plus ventral thalamus", "midbrain" and "cortex" (reference). Within the
thalamic dissection, GABAergic neurons (GAD2 > 0) include the prethalamic
lineage that forms the thalamic reticular nucleus; within the midbrain,
TH-positive neurons (TH > 0) approximate dopaminergic neurons.

For each group the script reports the number of cells, a pseudobulk value
(counts of the gene per million counts in the group, comparable in scale to the
nCPM of the adult single-nucleus data) and the percentage of cells in which the
gene was detected, overall and for two age bins (5-8 and 9-15 weeks after
fertilization).

Requires: pip install cellxgene-census (network access to the public Census).

Usage:
    python fetal_brain_single_cell.py <out_dir>
"""
import json
import sys
from pathlib import Path

import cellxgene_census
import numpy as np

CENSUS_VERSION = "2025-11-08"
DATASET = "13149914-ea03-4f01-8bf6-b793b667127b"
TISSUES = ["dorsal plus ventral thalamus", "midbrain", "cortex"]
GENES = ["EFHC1", "CACNA1E", "FOXJ1", "TH", "GAD2"]
REPORT = ["EFHC1", "CACNA1E", "FOXJ1"]


def age_bin(stage):
    return "5-8 weeks" if stage.startswith("Carnegie") else "9-15 weeks"


def summarize(adata, mask):
    sub = adata[mask]
    total = float(sub.obs["raw_sum"].sum())
    rec = {"cells": int(mask.sum())}
    for g in REPORT:
        x = sub[:, sub.var.feature_name == g].X
        counts = np.asarray(x.sum(axis=1)).ravel()
        rec[f"{g}_per_million"] = round(1e6 * counts.sum() / total, 1) if total else None
        rec[f"{g}_pct_cells"] = round(100 * (counts > 0).mean(), 1) if len(counts) else None
    return rec


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tissue_list = ", ".join(f"'{t}'" for t in TISSUES)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        adata = cellxgene_census.get_anndata(
            census, organism="Homo sapiens",
            obs_value_filter=(f"dataset_id == '{DATASET}' and cell_type == 'neuron' "
                              f"and tissue in [{tissue_list}]"),
            var_value_filter="feature_name in [" + ", ".join(f"'{g}'" for g in GENES) + "]",
            obs_column_names=["tissue", "development_stage", "raw_sum"])

    def gene(g):
        return np.asarray(adata[:, adata.var.feature_name == g].X.sum(axis=1)).ravel()

    tissue = adata.obs["tissue"].astype(str).to_numpy()
    ages = adata.obs["development_stage"].astype(str).map(age_bin).to_numpy()
    th, gad2 = gene("TH") > 0, gene("GAD2") > 0
    groups = {
        "thalamus, all neurons": tissue == "dorsal plus ventral thalamus",
        "thalamus, GAD2+ neurons": (tissue == "dorsal plus ventral thalamus") & gad2,
        "midbrain, TH+ neurons": (tissue == "midbrain") & th,
        "midbrain, other neurons": (tissue == "midbrain") & ~th,
        "cortex, all neurons": tissue == "cortex",
    }
    rows = []
    for name, mask in groups.items():
        for label, m in [("all ages", mask)] + [(b, mask & (ages == b)) for b in ("5-8 weeks", "9-15 weeks")]:
            rows.append({"group": name, "age": label} | summarize(adata, m))

    cols = list(rows[0])
    with open(out / "fetal_neurons_pseudobulk.tsv", "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    meta = {"source": "Braun et al. 2023, first-trimester developing human brain (CELLxGENE)",
            "dataset_id": DATASET, "census_version": CENSUS_VERSION,
            "cell_type": "neuron", "tissues": TISSUES,
            "markers": {"GABAergic": "GAD2 > 0", "dopaminergic (approximate)": "TH > 0"}}
    (out / "fetal_neurons_meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    for r in rows:
        print(r)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
