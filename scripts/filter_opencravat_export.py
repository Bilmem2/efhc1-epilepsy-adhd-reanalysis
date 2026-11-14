"""Apply the variant-prioritization funnel to an OpenCRAVAT variant export.

OpenCRAVAT (v3.2.3) is used only for annotation; the export may be
unfiltered. Every rule is applied here so that the funnel is versioned and
reproducible.

Usage:
    python filter_opencravat_export.py <export_variant.tsv> <out_dir> [--track GENE,GENE,...]

The export is the OpenCRAVAT "variant" table exported as TSV (two header
rows: column groups and column names). --track lists genes whose variants
are traced through the steps in funnel_summary.json.

The outputs model_A_dominant.tsv and model_B_recessive.tsv contain the
individual's genotype data and must not be shared publicly.

Steps (each step's count is written to funnel_summary.json):
  0. Pre-filter: protein-altering or canonical splice-site consequence, and
     not classified benign/likely benign in ClinVar (absent = kept). The
     aggregate ClinVar classification is used: "Benign", "Likely benign" and
     "Benign/Likely benign" are removed; "Conflicting classifications" and
     VUS are kept.
  1. Technical QC: VCF FILTER == PASS, total reads >= 10,
     heterozygous variant allele fraction 0.20-0.80.
  2. Population frequency, applied with AND across sources; a missing
     value passes (it is not treated as rare, it is flagged):
       gnomAD v4 global AF   and   All of Us maximum ancestry-group AF.
     Model A (dominant / de novo): heterozygous, both <= 0.001.
     Model B (recessive): both <= 0.01, and homozygous or >= 2 qualifying
     variants in the same gene (phase unknown).
  3. Established disease gene: the gene has HPO disease annotations
     (OpenCRAVAT HPO module, gene level). The OpenCRAVAT OMIM column is
     variant level (OMIM allelic variants) and is not used for this.
  4. Neurodevelopmental relevance: the gene's HPO terms mention seizures,
     EEG, intellectual/developmental, language, learning or attention terms.
  5. Inheritance fit: for model A the gene's HPO terms include autosomal
     dominant or X-linked inheritance; for model B, recessive inheritance.
The last step (phenotype fit, variant-level evidence) is manual.
1000 Genomes and gnomAD v3 are not used: both are contained in gnomAD v4.
"""
import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

PROTEIN_ALTERING = {"missense_variant", "inframe_deletion", "inframe_insertion",
                    "stop_gained", "frameshift_truncation", "frameshift_elongation",
                    "splice_site_variant", "start_lost", "stop_lost",
                    "complex_substitution"}
QC_MIN_READS = 10
HET_VAF = (0.20, 0.80)
DOMINANT_MAX_AF = 0.001
RECESSIVE_MAX_AF = 0.01
NEURO_TERMS = re.compile(r"seizure|epilep|EEG|intellectual|developmental|language|"
                         r"speech|learning|attention|hyperactiv|autis", re.I)
# A heterozygous call fits a dominant, X-linked or de novo model; HPO gene terms
# carry the inheritance of the gene's known diseases.
DOMINANT_TERMS = re.compile(r"Autosomal dominant|X-linked|Somatic mosaicism", re.I)
RECESSIVE_TERMS = re.compile(r"Autosomal recessive|X-linked recessive", re.I)


def read_export(path):
    lines = [l for l in open(path, encoding="utf-8") if not l.startswith("#")]
    groups = lines[0].rstrip("\n").split("\t")
    names = lines[1].rstrip("\n").split("\t")
    cols, current = [], ""
    for i, name in enumerate(names):
        group = groups[i] if i < len(groups) else ""
        if group:
            current = group
        cols.append(f"{current}|{name}")
    rows = [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in lines[2:]]
    return cols, rows


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def at_most(value, threshold):
    return value is None or value <= threshold


def main(export_path, out_dir, tracked=()):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cols, rows = read_export(export_path)
    key = {
        "gene": "Variant_Annotation|Gene",
        "so": "Variant_Annotation|Sequence_Ontology",
        "protein": "Variant_Annotation|Protein_Change",
        "transcript": "Variant_Annotation|Transcript",
        "zyg": "VCF_Info|Zygosity",
        "filter": "VCF_Info|VCF_filter",
        "reads": "VCF_Info|Total_reads",
        "vaf": "VCF_Info|Variant_AF",
        "g4": "gnomAD4|Global_AF",
        "g4_hom": "gnomAD4|Global_Alt Homozygotes",
        "aou_max": "All_of Us|Max_AF",
        "aou_pop": "All_of Us|Max_Sub-Population",
        "clinvar": "ClinVar|Clinical_Significance",
        "hpo": "Human_Phenotype Ontology|HPO_Term",
        "revel": "REVEL|Score",
        "cadd": "CADD|Score",
    }
    missing = [k for k, c in key.items() if c not in cols]
    if missing:
        sys.exit(f"Missing expected columns: {missing}")

    summary = collections.OrderedDict()
    summary["input_rows"] = len(rows)

    def benign(sig):
        return sig.startswith("Benign") or sig.startswith("Likely benign")

    rows = [r for r in rows if r[key["so"]] in PROTEIN_ALTERING]
    summary["protein_altering_or_splice"] = len(rows)
    rows = [r for r in rows if not benign(r[key["clinvar"]])]
    summary["not_benign_in_ClinVar"] = len(rows)

    def qc(r):
        if r[key["filter"]] != "PASS":
            return False
        reads = num(r[key["reads"]])
        if reads is None or reads < QC_MIN_READS:
            return False
        if r[key["zyg"]] == "het":
            vaf = num(r[key["vaf"]])
            if vaf is None or not HET_VAF[0] <= vaf <= HET_VAF[1]:
                return False
        return True

    passed_qc = [r for r in rows if qc(r)]
    summary["after_QC"] = len(passed_qc)

    def freq_ok(r, t):
        return at_most(num(r[key["g4"]]), t) and at_most(num(r[key["aou_max"]]), t)

    model_a = [r for r in passed_qc
               if r[key["zyg"]] == "het" and freq_ok(r, DOMINANT_MAX_AF)]
    rec_pool = [r for r in passed_qc if freq_ok(r, RECESSIVE_MAX_AF)]
    per_gene = collections.Counter(r[key["gene"]] for r in rec_pool if r[key["zyg"]] == "het")
    model_b = [r for r in rec_pool
               if r[key["zyg"]] == "hom" or per_gene[r[key["gene"]]] >= 2]

    summary["model_A_dominant_het_AF_le_0.001"] = len(model_a)
    summary["model_A_with_gnomAD4_missing"] = sum(1 for r in model_a if num(r[key["g4"]]) is None)
    summary["model_B_recessive_candidates_AF_le_0.01"] = len(model_b)
    summary["model_B_genes"] = len({r[key["gene"]] for r in model_b})

    def disease_gene(r):
        return bool(r[key["hpo"]])

    def neuro(r):
        return bool(NEURO_TERMS.search(r[key["hpo"]]))

    def inheritance_fits(r, label):
        terms = r[key["hpo"]]
        if label == "A":
            return bool(DOMINANT_TERMS.search(terms))
        return bool(RECESSIVE_TERMS.search(terms))

    for label, data in [("A", model_a), ("B", model_b)]:
        summary[f"model_{label}_disease_gene"] = sum(1 for r in data if disease_gene(r))
        summary[f"model_{label}_disease_gene_neuro_terms"] = sum(1 for r in data if neuro(r))
        summary[f"model_{label}_neuro_and_inheritance_fits"] = sum(
            1 for r in data if neuro(r) and inheritance_fits(r, label))
        summary[f"model_{label}_no_disease_gene"] = sum(1 for r in data if not disease_gene(r))

    def tier(r, label):
        if neuro(r):
            fit = "inheritance fits" if inheritance_fits(r, label) else "inheritance does not fit"
            return f"neuro disease gene, {fit}"
        return "other disease gene" if disease_gene(r) else "no disease annotation"

    def where(r):
        if r in model_a:
            return f"model A, {tier(r, 'A')}"
        if r in model_b:
            return f"model B, {tier(r, 'B')}"
        if not qc(r):
            return "removed at QC"
        g4, mx = num(r[key["g4"]]), num(r[key["aou_max"]])
        return f"removed at frequency (gnomAD4 {g4}, AoU max {mx} {r[key['aou_pop']]})"

    if tracked:
        summary["tracked_genes"] = {}
        for r in rows:
            if r[key["gene"]] in tracked:
                summary["tracked_genes"].setdefault(r[key["gene"]], []).append(
                    f"{r[key['protein']] or r[key['so']]}: {where(r)}")
        for g in tracked:
            summary["tracked_genes"].setdefault(
                g, ["no call left after the consequence and ClinVar-benign steps"])

    report_cols = ["gene", "protein", "so", "transcript", "zyg", "g4", "g4_hom",
                   "aou_max", "aou_pop", "clinvar", "revel", "cadd"]
    for name, label, data in [("model_A_dominant.tsv", "A", model_a),
                              ("model_B_recessive.tsv", "B", model_b)]:
        with open(out / name, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(report_cols + ["gene_tier", "neuro_HPO_terms"])
            order = {"neuro disease gene, inheritance fits": 0,
                     "neuro disease gene, inheritance does not fit": 1,
                     "other disease gene": 2, "no disease annotation": 3}
            for r in sorted(data, key=lambda r: (order[tier(r, label)], r[key["gene"]])):
                terms = sorted({t for t in r[key["hpo"]].split(";") if NEURO_TERMS.search(t)})
                w.writerow([r[key[c]] for c in report_cols] + [tier(r, label), ";".join(terms)])
    (out / "funnel_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                             encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("export", help="OpenCRAVAT variant export (TSV)")
    parser.add_argument("out_dir", help="output directory")
    parser.add_argument("--track", default="",
                        help="comma-separated genes to trace through the steps")
    args = parser.parse_args()
    main(args.export, args.out_dir, [g for g in args.track.split(",") if g])
