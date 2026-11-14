"""Look up the final candidates in the Turkish Variome (TR-Variome).

TR-Variome (Kars et al., PNAS 2021; 3,362 unrelated individuals, 2,589 exomes and
773 genomes) is distributed as one allele-frequency table on figshare:

    https://figshare.com/articles/dataset/15147642  (file TurkishVariome.tsv.gz)

The table has a GRCh38 position column, so variants are matched on chromosome,
GRCh38 position, reference and alternate allele. Variants that are not in the
table are reported as "not listed"; the table has no allele number for absent
sites, so no frequency is given for them.

Usage:
    python turkish_variome_lookup.py <TurkishVariome.tsv.gz> <variants.tsv> <out.tsv>

<variants.tsv> needs the columns gene, protein and variant_id_GRCh38
(chrom-pos-ref-alt), as in results/population/tgd_candidate_check.tsv.
"""
import csv
import gzip
import sys

SOURCE = "TR-Variome (Kars et al. 2021), figshare 15147642"
OUT_COLUMNS = ["gene", "protein", "variant_id_GRCh38", "source", "status",
               "AC", "AN", "AF", "homozygotes", "sequencing"]


def main(variome_path, variants_path, out_path):
    with open(variants_path, encoding="utf-8", newline="") as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))
    wanted = {}
    for v in variants:
        chrom, pos, ref, alt = v["variant_id_GRCh38"].split("-")
        wanted[(chrom, pos, ref, alt)] = v

    found = {}
    positions = {(c, p) for c, p, _, _ in wanted}
    with gzip.open(variome_path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            first = line.split("\t", 3)
            if (first[0], first[2]) not in positions:
                continue
            row = line.rstrip("\n").split("\t")
            key = (row[col["CHROM"]], row[col["GRCh38Pos"]], row[col["REF"]], row[col["ALT"]])
            if key in wanted:
                found[key] = row

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        out.writerow(OUT_COLUMNS)
        for key, v in wanted.items():
            row = found.get(key)
            if row is None:
                out.writerow([v["gene"], v["protein"], v["variant_id_GRCh38"], SOURCE,
                              "not listed", "", "", "", "", ""])
            else:
                out.writerow([v["gene"], v["protein"], v["variant_id_GRCh38"], SOURCE, "listed",
                              row[col["AC"]], row[col["AN"]], row[col["AF"]],
                              row[col["Hom"]], row[col["SequencingMethod"]]])
    print(f"{len(found)} of {len(wanted)} variants listed in TR-Variome; written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
