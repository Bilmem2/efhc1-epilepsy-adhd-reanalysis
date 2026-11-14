"""Join the manual review decisions (last prioritization step) with the
automated output of filter_opencravat_export.py.

The automated steps end with heterozygous variants in genes whose HPO
annotations include neurodevelopmental terms and a dominant or X-linked mode of
inheritance. Because HPO terms and inheritance are merged at the gene level, a
neurological term from one disease and a dominant mode from another disease of
the same gene can match automatically; the manual step therefore checks each
gene at the level of its individual disorders and records a reason for every
decision.

"Not prioritized" means not prioritized for this case's research question; it
is not a clinical classification of the variant.

Usage:
    python record_manual_review.py <model_A_dominant.tsv> <decisions.tsv> <out.tsv>
                                   [--vep <vep_review_set.tsv>]

decisions.tsv has the columns gene, decision, reason_category,
main_associated_disorders, reason (one row per gene). The optional VEP table
(columns gene and reference_protein_changes, protein changes separated by ";")
marks whether the same allele was also retained in an independent
re-annotation.
"""
import argparse
import csv
import sys

DECISION_ORDER = {"retained": 0, "not prioritized": 1}


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(model_a_path, decisions_path, out_path, vep_path=None):
    rows = [r for r in read_tsv(model_a_path)
            if r["gene_tier"] == "neuro disease gene, inheritance fits"]
    decisions = {d["gene"]: d for d in read_tsv(decisions_path)}

    genes = {r["gene"] for r in rows}
    if genes != set(decisions):
        sys.exit(f"Automated output and decisions differ: {sorted(genes ^ set(decisions))}")

    vep_alleles = None
    if vep_path:
        vep_alleles = {(r["gene"], p) for r in read_tsv(vep_path)
                       for p in r["reference_protein_changes"].split(";") if p}

    rows.sort(key=lambda r: (DECISION_ORDER[decisions[r["gene"]]["decision"]],
                             decisions[r["gene"]]["reason_category"], r["gene"]))
    cols = ["gene", "protein", "transcript", "gnomAD4_AF", "AoU_max_AF",
            "ClinVar_in_export", "decision", "reason_category",
            "main_associated_disorders", "reason"]
    if vep_alleles is not None:
        cols.append("same_allele_in_VEP_review_set")
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            d = decisions[r["gene"]]
            line = [r["gene"], r["protein"], r["transcript"], r["g4"] or "missing",
                    r["aou_max"] or "missing", r["clinvar"] or "not in ClinVar",
                    d["decision"], d["reason_category"], d["main_associated_disorders"],
                    d["reason"]]
            if vep_alleles is not None:
                line.append("yes" if (r["gene"], r["protein"]) in vep_alleles else "no")
            w.writerow(line)
    kept = sum(decisions[g]["decision"] == "retained" for g in genes)
    print(f"{len(rows)} reviewed, {kept} retained, {len(rows) - kept} not prioritized -> {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model_a", help="model_A_dominant.tsv from filter_opencravat_export.py")
    parser.add_argument("decisions", help="manual decisions, one row per gene")
    parser.add_argument("out", help="output table")
    parser.add_argument("--vep", help="optional independent VEP review set")
    args = parser.parse_args()
    main(args.model_a, args.decisions, args.out, args.vep)
