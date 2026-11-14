"""Download the Human Protein Atlas (v26) single-nuclei brain RNA files used by
hpa_single_nuclei_coexpression.py and check them against the SHA-256 sums of
the copies used for the manuscript (downloaded 13 Nov 2026).

Usage:
    python download_hpa_single_nuclei.py <out_dir>

The cluster-level file is about 140 MB.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

BASE_URL = "https://v26.proteinatlas.org/download/tsv/"
FILES = {
    "rna_single_nuclei_brain_cluster.tsv.zip":
        "6bf4803512ab8cd39d851b11e938ff2d9bc16ca2f1e9b9934968b6dba4bbc1d6",
    "rna_single_nuclei_brain_clusters.tsv.zip":
        "1df48fc7a65e50c47a948cf6e3f557e6b7081b719acd32692f7616cfbe7b71b0",
    "rna_single_nuclei_cell_type.tsv.zip":
        "e3b09b5d8feab271eb4e82489f6d20d91ef515295bdbdd6ebb33983a4784be1c",
    "rna_single_nuclei_cell_type_cell_types.tsv.zip":
        "cd481dc86dc055c6cc656ef643373369e6f1d5d0689c48acd1978200111c2796",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mismatches = []
    for name, expected in FILES.items():
        target = out / name
        if not target.exists():
            print(f"downloading {name}")
            urllib.request.urlretrieve(BASE_URL + name, target)
        ok = sha256(target) == expected
        print(f"{name}: {'checksum ok' if ok else 'CHECKSUM DIFFERS'}")
        if not ok:
            mismatches.append(name)
    if mismatches:
        sys.exit("Some files differ from the versions used for the manuscript; "
                 "results may not match exactly.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
