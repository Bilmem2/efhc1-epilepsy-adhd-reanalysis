"""EFHC1 variant landscape: all ClinVar records for EFHC1 and their position on
the MANE Select protein (ENST00000371068, Q5JVL4, 640 aa); optionally, the
EFHC1 calls of an individual exome.

Usage:
    python efhc1_variant_landscape.py <input_dir> <out_dir> [--vcf <exome.vcf>]

<input_dir> must contain clinvar_EFHC1_records.vcf (ClinVar GRCh38 VCF records
in chr6:52,355,000-52,535,000), mane_cds_exons.tsv (GRCh38 start and end of
each MANE Select CDS exon) and, optionally, gnomad41_exomes_EFHC1_region.tsv
(gnomAD v4.1 exomes; the first seven columns are chrom, pos, ref, alt, AC,
AN, AF).

Outputs clinvar_EFHC1_annotated.tsv and summary.json; with --vcf also
proband_EFHC1_calls.tsv. Outputs written with --vcf contain the individual's
genotype data and must not be shared publicly.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

CLINVAR_RELEASE = "ClinVar GRCh38 VCF, last modified 2 Nov 2026"
REGION = (52355000, 52535000)
DOMAINS = [("DM10-1", 93, 198), ("DM10-2", 239, 359), ("DM10-3", 416, 520),
           ("EF-hand", 574, 609)]
SIG_GROUPS = [("P/LP", ("Pathogenic", "Likely_pathogenic")),
              ("Conflicting", ("Conflicting",)),
              ("VUS", ("Uncertain_significance",)),
              ("B/LB", ("Benign", "Likely_benign"))]
PHENO_FLAGS = ("attention", "adhd", "hyperactiv", "behavio", "autis", "intellectual")


def info_dict(field):
    return dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in field.split(";"))


def sig_group(sig):
    for name, keys in SIG_GROUPS:
        if any(sig.startswith(k) for k in keys):
            return name
    return "Other/none"


def cds_codon(pos, exons):
    """1-based codon number on the + strand MANE CDS, or None if outside CDS."""
    offset = 0
    for start, end in exons:
        if start <= pos <= end:
            return (offset + pos - start) // 3 + 1
        offset += end - start + 1
    return None


def domain(codon):
    if codon is None:
        return ""
    for name, s, e in DOMAINS:
        if s <= codon <= e:
            return name
    return "outside domains"


def load_gnomad(path):
    freq = {}
    if path.exists():
        for line in open(path, encoding="utf-8"):
            c, p, r, a, ac, an, af = line.rstrip("\n").split("\t")[:7]
            freq[(int(p), r, a)] = (int(ac), int(an), float(af))
    return freq


def main(d, out, proband_vcf=None):
    d, out = Path(d), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    exons = [tuple(map(int, l.split("\t")[:2])) for l in open(d / "mane_cds_exons.tsv")]
    cds_len = sum(e - s + 1 for s, e in exons)
    gnomad = load_gnomad(d / "gnomad41_exomes_EFHC1_region.tsv")

    proband = {}
    if proband_vcf:
        for line in open(proband_vcf, encoding="utf-8"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] in ("chr6", "6") and REGION[0] <= int(f[1]) <= REGION[1]:
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
                for alt in f[4].split(","):
                    proband[(int(f[1]), f[3], alt)] = fmt

    records = []
    for line in open(d / "clinvar_EFHC1_records.vcf", encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        info = info_dict(f[7])
        pos, ref, alt = int(f[1]), f[3], f[4]
        sig = info.get("CLNSIG", "")
        mc = info.get("MC", "")
        consequence = mc.split("|")[1].split(",")[0] if "|" in mc else ""
        codon = cds_codon(pos, exons) if consequence in (
            "missense_variant", "nonsense", "frameshift_variant", "synonymous_variant",
            "inframe_deletion", "inframe_insertion", "inframe_indel") else None
        g = gnomad.get((pos, ref, alt))
        record = {
            "variation_id": f[2], "pos": pos, "ref": ref, "alt": alt,
            "significance": sig, "group": sig_group(sig),
            "review": info.get("CLNREVSTAT", ""), "consequence": consequence,
            "codon": codon or "", "domain": domain(codon),
            "condition": info.get("CLNDN", ""),
            "gnomAD_exome_AC_AN": f"{g[0]}/{g[1]}" if g else "",
            "gnomAD_exome_AF": g[2] if g else "",
        }
        if proband_vcf:
            record["in_proband"] = (pos, ref, alt) in proband
        records.append(record)

    cols = list(records[0])
    with open(out / "clinvar_EFHC1_annotated.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(records, key=lambda r: r["pos"]))

    clinvar_by_allele = {(r["pos"], r["ref"], r["alt"]): r for r in records}
    if proband_vcf:
        with open(out / "proband_EFHC1_calls.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pos", "ref", "alt", "GT", "DP", "GQ", "gnomAD_exome_AC_AN",
                        "gnomAD_exome_AF", "clinvar_id", "clinvar_significance"])
            for k in sorted(proband):
                g, cv, fmt = gnomad.get(k), clinvar_by_allele.get(k), proband[k]
                w.writerow([*k, fmt.get("GT"), fmt.get("DP"), fmt.get("GQ"),
                            f"{g[0]}/{g[1]}" if g else "", g[2] if g else "",
                            cv["variation_id"] if cv else "", cv["significance"] if cv else ""])

    summary = collections.OrderedDict()
    summary["clinvar_snapshot"] = CLINVAR_RELEASE
    summary["mane_cds_length_nt"] = cds_len
    summary["clinvar_records"] = len(records)
    summary["by_significance"] = collections.Counter(r["group"] for r in records)
    summary["by_significance_and_consequence"] = {
        f"{g} | {c or 'non-coding/other'}": n for (g, c), n in collections.Counter(
            (r["group"], r["consequence"]) for r in records).most_common()}
    summary["P_LP_records"] = [{k: r[k] for k in ("variation_id", "pos", "ref", "alt",
                                "significance", "review", "consequence", "codon",
                                "domain", "condition", "gnomAD_exome_AF")}
                               for r in records if r["group"] == "P/LP"]
    ef = [r for r in records if r["domain"] == "EF-hand"]
    summary["EF_hand_records"] = collections.Counter(f"{r['group']} | {r['consequence']}" for r in ef)
    summary["missense_by_domain_and_group"] = {
        f"{dom} | {g}": n for (dom, g), n in sorted(collections.Counter(
            (r["domain"], r["group"]) for r in records
            if r["consequence"] == "missense_variant").items())}
    summary["records_with_neurodevelopmental_condition_terms"] = [
        {k: r[k] for k in ("variation_id", "codon", "significance", "condition")}
        for r in records if any(t in r["condition"].lower() for t in PHENO_FLAGS)]
    summary["conditions"] = collections.Counter(
        c for r in records for c in r["condition"].split("|") if c)
    if proband_vcf:
        summary["proband_calls"] = len(proband)
        summary["proband_calls_in_clinvar"] = [
            {"pos": k[0], "ref": k[1], "alt": k[2],
             "variation_id": clinvar_by_allele[k]["variation_id"],
             "significance": clinvar_by_allele[k]["significance"]}
            for k in sorted(proband) if k in clinvar_by_allele]
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                      encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_dir", help="directory with the public ClinVar, MANE and gnomAD inputs")
    parser.add_argument("out_dir", help="output directory")
    parser.add_argument("--vcf", help="optional individual exome VCF (outputs are then private)")
    args = parser.parse_args()
    main(args.input_dir, args.out_dir, args.vcf)
