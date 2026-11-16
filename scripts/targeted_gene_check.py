"""Search an OpenCRAVAT variant export directly for rare protein-altering
variants in selected genes.

Steps 3-5 of the prioritization (filter_opencravat_export.py) rely on gene
annotations, which are incomplete. This check therefore looks directly at two
gene groups named in the manuscript (Methods and Procedures):
  - the three ADHD genes identified by rare-variant analysis (Demontis et al.,
    2026): MAP1A, ANO8, ANK2;
  - genes whose pathogenic variants cause epilepsy with attention or learning
    problems: SLC6A1, SYNGAP1, CHD2, GRIN2A, SLC2A1, GABRA1, GABRG2, CACNA1A,
    STX1B.

Rules (the same definitions as filter_opencravat_export.py):
  - protein-altering or canonical splice-site consequence;
  - technical QC: FILTER == PASS, total reads >= 10, heterozygous variant
    allele fraction 0.20-0.80;
  - rare: gnomAD v4 global AF and All of Us maximum ancestry-group AF both
    <= 0.01 (the recessive-model threshold, so that both models are covered;
    a missing value passes).
Unlike the main funnel, variants classified as benign in ClinVar are NOT
removed, so that every rare protein-altering call is reported with its
ClinVar class.

Usage:
    python targeted_gene_check.py <export_variant.tsv> <out_dir>

Output: <out_dir>/targeted_gene_check.tsv, one row per gene. It lists, for
each gene, the number of calls of any consequence, the number of
protein-altering calls, the number of rare calls passing QC, and for those
calls the protein change, ClinVar class and population frequencies. Zygosity
and other genotype fields are not written.

A gene without any call cannot be judged from a VCF alone: the absence of
calls may mean that the person carries no variant there or that the region
was poorly covered. Such genes are named in the printed summary.
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_opencravat_export import (HET_VAF, PROTEIN_ALTERING, QC_MIN_READS,  # noqa: E402
                                      RECESSIVE_MAX_AF, at_most, num, read_export)

GENE_GROUPS = {
    "ADHD rare-variant gene": ["MAP1A", "ANO8", "ANK2"],
    "epilepsy with attention or learning problems": [
        "SLC6A1", "SYNGAP1", "CHD2", "GRIN2A", "SLC2A1", "GABRA1", "GABRG2",
        "CACNA1A", "STX1B"],
}
KEY = {
    "gene": "Variant_Annotation|Gene",
    "so": "Variant_Annotation|Sequence_Ontology",
    "protein": "Variant_Annotation|Protein_Change",
    "zyg": "VCF_Info|Zygosity",
    "filter": "VCF_Info|VCF_filter",
    "reads": "VCF_Info|Total_reads",
    "vaf": "VCF_Info|Variant_AF",
    "g4": "gnomAD4|Global_AF",
    "aou_max": "All_of Us|Max_AF",
    "clinvar": "ClinVar|Clinical_Significance",
}


def passes_qc(r):
    if r[KEY["filter"]] != "PASS":
        return False
    reads = num(r[KEY["reads"]])
    if reads is None or reads < QC_MIN_READS:
        return False
    if r[KEY["zyg"]] == "het":
        vaf = num(r[KEY["vaf"]])
        if vaf is None or not HET_VAF[0] <= vaf <= HET_VAF[1]:
            return False
    return True


def is_rare(r):
    return (at_most(num(r[KEY["g4"]]), RECESSIVE_MAX_AF)
            and at_most(num(r[KEY["aou_max"]]), RECESSIVE_MAX_AF))


def fmt_af(value):
    v = num(value)
    return "missing" if v is None else f"{v:.2g}"


def main(export_path, out_dir):
    cols, rows = read_export(export_path)
    missing = [k for k, c in KEY.items() if c not in cols]
    if missing:
        sys.exit(f"Missing expected columns: {missing}")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    table = []
    for group, genes in GENE_GROUPS.items():
        for gene in genes:
            all_calls = [r for r in rows if r[KEY["gene"]] == gene]
            calls = [r for r in all_calls if r[KEY["so"]] in PROTEIN_ALTERING]
            rare = [r for r in calls if passes_qc(r) and is_rare(r)]
            found = "; ".join(
                f"{r[KEY['protein']] or r[KEY['so']]} (ClinVar: {r[KEY['clinvar']] or 'no record'}; "
                f"gnomAD4 {fmt_af(r[KEY['g4']])}; AoU max {fmt_af(r[KEY['aou_max']])})"
                for r in rare)
            table.append([gene, group, len(all_calls), len(calls), len(rare), found or "none"])

    with open(out / "targeted_gene_check.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "group", "calls_any_consequence", "protein_altering_calls",
                    "rare_calls_passing_QC", "rare_calls"])
        w.writerows(table)

    print(f"Targeted check: {len(table)} genes; rare = AF <= {RECESSIVE_MAX_AF} "
          "in gnomAD v4 and All of Us; ClinVar-benign calls kept")
    for gene, group, n_all, n_calls, n_rare, found in table:
        print(f"  {gene:<8} {n_all:>2} calls, {n_calls} protein-altering, {n_rare} rare: {found}")
    n_genes = sum(1 for row in table if row[4])
    print(f"Genes with a rare protein-altering call: {n_genes} of {len(table)}")
    no_calls = [row[0] for row in table if row[2] == 0]
    if no_calls:
        print(f"No calls of any kind in {', '.join(no_calls)}: from a VCF alone, this cannot "
              "be told apart from low coverage of the region.")
    print(f"Written: {out / 'targeted_gene_check.tsv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("export", help="OpenCRAVAT variant export (TSV)")
    parser.add_argument("out_dir", help="output directory")
    args = parser.parse_args()
    main(args.export, args.out_dir)
