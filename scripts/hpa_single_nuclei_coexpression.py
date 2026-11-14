"""Cell-type co-expression of EFHC1 and CACNA1E in the adult human brain.

Question: is EFHC1 expressed in the neuronal cell types where Cav2.3
(CACNA1E) acts in the proposed mechanism (thalamic reticular-like GABAergic
neurons for absence seizures; midbrain dopaminergic neurons for dopamine
release), or mainly in ependymal / ciliated cells, as in mouse
(Suzuki et al., 2020)?

Input: Human Protein Atlas (v26) single-nuclei brain RNA data, which are
derived from the adult human brain atlas of Siletti et al. (2023). Values
are pseudobulk nCPM per cell type or per cluster. Co-detection in the same
cluster means the same cell population, not proven co-expression in single
cells; mRNA is not protein.

Usage:
    python hpa_single_nuclei_coexpression.py <hpa_dir> <out_dir>

<hpa_dir> holds the HPA v26 single-nuclei brain files downloaded by
download_hpa_single_nuclei.py. Outputs: cell_type_nCPM.tsv,
cluster_nCPM_selected.tsv, percentiles.json and summary.json.
"""
import csv
import io
import json
import statistics
import sys
import zipfile
from pathlib import Path

GENES = ["EFHC1", "CACNA1E", "TRPM2", "FOXJ1", "TH", "SLC6A3", "GAD2", "PVALB"]
DETECTED = 1.0  # nCPM
PERCENTILE_CELL_TYPES = ["dopaminergic", "thalamic projection neurons", "mge derived interneuron",
                         "midbrain derived inhibitory neurons", "forebrain excitatory",
                         "ependymal cells"]


def rows(path):
    z = zipfile.ZipFile(path)
    with z.open(z.namelist()[0]) as fh:
        yield from csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")


def main(hpa_dir, out_dir):
    hpa, out = Path(hpa_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = {"source": "Human Protein Atlas v26, single nuclei brain (Siletti et al., 2023)",
               "detection_threshold_nCPM": DETECTED}

    # 1. 48 brain cell types
    classes = {r["Cell type"]: r["Cell type class"]
               for r in rows(hpa / "rna_single_nuclei_cell_type_cell_types.tsv.zip")}
    ct = {}
    detected = {c: {} for c in PERCENTILE_CELL_TYPES}  # Ensembl gene ID -> nCPM
    for r in rows(hpa / "rna_single_nuclei_cell_type.tsv.zip"):
        if r["Gene name"] in GENES:
            ct.setdefault(r["Cell type"], {})[r["Gene name"]] = float(r["nCPM"])
        if r["Cell type"] in detected and float(r["nCPM"]) >= DETECTED:
            detected[r["Cell type"]][r["Gene"]] = float(r["nCPM"])

    # Rank of EFHC1 and CACNA1E among the genes detected in each cell type:
    # percentage of detected genes with lower expression.
    percentiles = {}
    for c, values in detected.items():
        n = len(values)
        percentiles[c] = {"genes_detected": n}
        for g in ("EFHC1", "CACNA1E"):
            x = ct[c][g]
            percentiles[c][g] = {"nCPM": x, "percentile_among_detected_genes":
                                 round(100 * sum(v < x for v in values.values()) / n, 1)}
    (out / "percentiles.json").write_text(json.dumps(percentiles, indent=2), encoding="utf-8")
    with open(out / "cell_type_nCPM.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cell_type", "class"] + GENES)
        for c in sorted(ct, key=lambda c: -ct[c].get("EFHC1", 0)):
            w.writerow([c, classes.get(c, "")] + [ct[c].get(g, "") for g in GENES])

    neuronal = {c: v for c, v in ct.items() if classes.get(c) == "neuronal cells"}
    top_neuron = max(neuronal, key=lambda c: neuronal[c]["EFHC1"])
    summary["EFHC1_ependymal"] = ct["ependymal cells"]["EFHC1"]
    summary["EFHC1_choroid_plexus"] = ct["choroid plexus epithelial cells"]["EFHC1"]
    summary["EFHC1_highest_neuronal_type"] = [top_neuron, neuronal[top_neuron]["EFHC1"]]
    summary["EFHC1_median_neuronal"] = round(statistics.median(v["EFHC1"] for v in neuronal.values()), 1)
    summary["EFHC1_ependymal_to_highest_neuron_ratio"] = round(
        ct["ependymal cells"]["EFHC1"] / max(neuronal[top_neuron]["EFHC1"], 1e-9), 1)
    for c in ["dopaminergic", "thalamic projection neurons", "mge derived interneuron",
              "midbrain derived inhibitory neurons", "forebrain excitatory", "ependymal cells"]:
        summary[f"{c}"] = {g: ct[c].get(g) for g in ["EFHC1", "CACNA1E", "TRPM2", "FOXJ1"]}
    summary["neuronal_types_EFHC1_detected"] = f"{sum(v['EFHC1'] >= DETECTED for v in neuronal.values())}/{len(neuronal)}"
    summary["neuronal_types_both_detected"] = f"{sum(v['EFHC1'] >= DETECTED and v['CACNA1E'] >= DETECTED for v in neuronal.values())}/{len(neuronal)}"

    # 2. Clusters: thalamic reticular-like and midbrain dopaminergic populations
    meta = {(r["Brain region"], r["Cluster"]): r
            for r in rows(hpa / "rna_single_nuclei_brain_clusters.tsv.zip")}
    cl = {}
    for r in rows(hpa / "rna_single_nuclei_brain_cluster.tsv.zip"):
        if r["Gene name"] in GENES:
            cl.setdefault((r["Brain region"], r["Cluster"]), {})[r["Gene name"]] = float(r["nCPM"])

    def describe(keys):
        return [{"region": k[0], "cluster": k[1], "cell_type": meta[k]["Cell type"],
                 "cells": int(meta[k]["Cell count"]), **{g: cl[k].get(g) for g in GENES}}
                for k in keys]

    thal = [k for k in cl if k[0] == "thalamus, neurons"]
    # Thalamic reticular nucleus neurons are GABAergic and parvalbumin-positive.
    trn_like = sorted([k for k in thal if cl[k]["GAD2"] >= 50 and cl[k]["PVALB"] >= 20],
                      key=lambda k: -cl[k]["PVALB"])
    da = sorted([k for k in cl if k[0] == "midbrain, neurons"
                 and cl[k]["TH"] >= 50 and cl[k]["SLC6A3"] >= 20],
                key=lambda k: -cl[k]["SLC6A3"])
    epend = sorted([k for k in cl if meta[k]["Cell type"] == "ependymal cells"],
                   key=lambda k: -cl[k]["EFHC1"])
    summary["thalamic_reticular_like_clusters"] = describe(trn_like)
    summary["midbrain_dopaminergic_clusters"] = describe(da)
    summary["ependymal_clusters"] = describe(epend)

    neuro_cl = [k for k in cl if meta[k]["Cell type class"] == "neuronal cells"]
    summary["neuronal_clusters_EFHC1_detected"] = f"{sum(cl[k]['EFHC1'] >= DETECTED for k in neuro_cl)}/{len(neuro_cl)}"
    summary["neuronal_clusters_both_detected"] = f"{sum(cl[k]['EFHC1'] >= DETECTED and cl[k]['CACNA1E'] >= DETECTED for k in neuro_cl)}/{len(neuro_cl)}"

    with open(out / "cluster_nCPM_selected.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "region", "cluster", "cell_type", "cells"] + GENES)
        for group, keys in [("thalamic_reticular_like", trn_like), ("midbrain_dopaminergic", da),
                            ("ependymal", epend)]:
            for d in describe(keys):
                w.writerow([group, d["region"], d["cluster"], d["cell_type"], d["cells"]]
                           + [d[g] for g in GENES])

    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                      encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
