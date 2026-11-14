"""EFHC1, CACNA1E and FOXJ1 in the developing human thalamus (BrainSpan).

BrainSpan's developmental transcriptome (RNA-seq, Gencode v10, summarized to
genes; RPKM) samples the mediodorsal nucleus of the thalamus (MD) from 13 weeks
after conception to 40 years. This script extracts the MD samples, groups them
by developmental period and reports, for each sample, the RPKM of EFHC1,
CACNA1E and the ciliated-cell marker FOXJ1, and the percentile of EFHC1 and
CACNA1E among the genes detected in that sample (RPKM >= 1), defined as for the
single-nucleus data (percentage of detected genes with lower expression).

Bulk tissue cannot separate cell types: FOXJ1 is reported to show how much
ependymal (ciliated) tissue may contribute to the EFHC1 signal.

Input: genes_matrix_csv.zip from
https://www.brainspan.org/api/v2/well_known_file_download/267666525

Usage:
    python brainspan_thalamus_by_age.py <genes_matrix_csv.zip> <out_dir>
"""
import csv
import io
import json
import statistics
import sys
import zipfile
from pathlib import Path

GENES = ["EFHC1", "CACNA1E", "FOXJ1"]
STRUCTURE = "MD"
DETECTED_RPKM = 1.0
PERIODS = [  # (label, lower bound in years after birth; prenatal ages are < 0)
    ("prenatal", None), ("infancy (0-1 y)", 0), ("childhood (2-11 y)", 2),
    ("adolescence (12-19 y)", 12), ("adulthood (20-40 y)", 20)]


def age_in_years(age):
    value, unit = age.split()
    value = float(value)
    if unit == "pcw":
        return value / 52 - 40 / 52          # before birth: negative
    if unit == "mos":
        return value / 12
    return value


def period(years):
    if years < 0:
        return PERIODS[0][0]
    label = PERIODS[1][0]
    for name, lower in PERIODS[1:]:
        if years >= lower:
            label = name
    return label


def percentile(value, detected):
    return round(100 * sum(v < value for v in detected) / len(detected), 1)


def main(zip_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    z = zipfile.ZipFile(zip_path)
    columns = list(csv.DictReader(io.TextIOWrapper(z.open("columns_metadata.csv"), encoding="utf-8")))
    rows = list(csv.DictReader(io.TextIOWrapper(z.open("rows_metadata.csv"), encoding="utf-8")))
    md = [i for i, c in enumerate(columns) if c["structure_acronym"] == STRUCTURE]
    gene_row = {r["gene_symbol"]: i for i, r in enumerate(rows) if r["gene_symbol"] in GENES}

    values = {i: [] for i in md}           # all genes, MD samples only
    selected = {g: {} for g in GENES}
    with io.TextIOWrapper(z.open("expression_matrix.csv"), encoding="utf-8") as fh:
        for r, line in enumerate(csv.reader(fh)):
            sample_values = [float(line[i + 1]) for i in md]   # first column is the row number
            for i, v in zip(md, sample_values):
                values[i].append(v)
            for g, gr in gene_row.items():
                if r == gr:
                    selected[g] = dict(zip(md, sample_values))

    table = []
    for i in md:
        c = columns[i]
        years = age_in_years(c["age"])
        detected = [v for v in values[i] if v >= DETECTED_RPKM]
        rec = {"donor": c["donor_name"], "age": c["age"], "period": period(years),
               "genes_detected": len(detected)}
        for g in GENES:
            rec[f"{g}_RPKM"] = round(selected[g][i], 2)
        for g in ("EFHC1", "CACNA1E"):
            rec[f"{g}_percentile"] = percentile(selected[g][i], detected)
        table.append((years, rec))
    table.sort(key=lambda t: t[0])

    with open(out / "thalamus_MD_by_age.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0][1]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for _, rec in table:
            w.writerow(rec)

    summary = {"source": "BrainSpan developmental transcriptome, RNA-seq Gencode v10 (RPKM)",
               "structure": "mediodorsal nucleus of thalamus (MD)",
               "detection_threshold_RPKM": DETECTED_RPKM, "periods": {}}
    for name, _ in PERIODS:
        recs = [rec for _, rec in table if rec["period"] == name]
        if not recs:
            continue
        summary["periods"][name] = {"samples": len(recs)} | {
            f"median_{k}": round(statistics.median(rec[k] for rec in recs), 2)
            for k in ("EFHC1_RPKM", "CACNA1E_RPKM", "FOXJ1_RPKM",
                      "EFHC1_percentile", "CACNA1E_percentile")}
    (out / "thalamus_MD_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
