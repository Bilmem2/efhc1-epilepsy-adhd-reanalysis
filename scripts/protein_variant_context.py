"""Structural and annotation context of an amino acid substitution, from public
protein resources only.

For a UniProt accession and a substitution such as E604D, the script downloads
the UniProt record, the AlphaFold model and the precomputed AlphaMissense
scores for the whole protein, and reports for the residue:
  - UniProt domains, regions, motifs and sites that contain it;
  - AlphaFold per-residue confidence (pLDDT);
  - relative solvent accessibility in the AlphaFold model (Shrake-Rupley,
    normalized by the maximum values of Tien et al., 2013; requires Biopython);
  - the AlphaMissense score and class of the substitution, and the scores of
    all other substitutions at the same position.
Only the protein accession is sent to the services; the substitution is looked
up locally.

Usage:
    python protein_variant_context.py <uniprot_accession> <substitution> <out_dir>
    e.g. python protein_variant_context.py Q5JVL4 E604D results/protein_context
"""
import csv
import io
import json
import sys
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
FEATURE_TYPES = {"Domain", "Region", "Motif", "Site", "Binding site", "Calcium binding",
                 "Repeat", "Zinc finger", "Modified residue"}
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
                "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
                "TYR": "Y", "VAL": "V"}
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223, "G": 104,
           "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155,
           "T": 172, "W": 285, "Y": 263, "V": 174}  # Tien et al. 2013, theoretical


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "protein-variant-context"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def main(acc, substitution, out_dir):
    ref, pos, alt = substitution[0], int(substitution[1:-1]), substitution[-1]
    out = Path(out_dir)
    cache = out / "inputs"
    cache.mkdir(parents=True, exist_ok=True)

    uniprot_path = cache / f"uniprot_{acc}.json"
    if not uniprot_path.exists():
        uniprot_path.write_bytes(fetch(UNIPROT.format(acc=acc)))
    uniprot = json.loads(uniprot_path.read_text(encoding="utf-8"))
    seq = uniprot["sequence"]["value"]
    if seq[pos - 1] != ref:
        sys.exit(f"Residue {pos} of {acc} is {seq[pos - 1]}, not {ref}")
    features = [
        {"type": f["type"], "start": f["location"]["start"]["value"],
         "end": f["location"]["end"]["value"], "description": f.get("description", "")}
        for f in uniprot.get("features", [])
        if f["type"] in FEATURE_TYPES
        and f["location"]["start"]["value"] <= pos <= f["location"]["end"]["value"]]

    record = json.loads(fetch(ALPHAFOLD_API.format(acc=acc)))[0]
    pdb_path = cache / Path(record["pdbUrl"]).name
    am_path = cache / Path(record["amAnnotationsUrl"]).name
    for path, url in ((pdb_path, record["pdbUrl"]), (am_path, record["amAnnotationsUrl"])):
        if not path.exists():
            path.write_bytes(fetch(url))

    plddt = None
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == pos:
            plddt = float(line[60:66])
            if THREE_TO_ONE.get(line[17:20]) != ref:
                sys.exit(f"AlphaFold model residue {pos} is {line[17:20]}, not {ref}")

    rsa = None
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
        structure = PDBParser(QUIET=True).get_structure(acc, pdb_path)
        ShrakeRupley().compute(structure[0], level="R")
        residue = structure[0]["A"][pos]
        rsa = round(residue.sasa / MAX_ASA[ref], 2)
    except ImportError:
        print("Biopython not installed; solvent accessibility skipped")

    scores = {}
    for row in csv.DictReader(io.StringIO(am_path.read_text())):
        v = row["protein_variant"]
        if int(v[1:-1]) == pos:
            scores[v] = {"score": float(row["am_pathogenicity"]), "class": row["am_class"]}
    key = f"{ref}{pos}{alt}"

    result = {
        "uniprot_accession": acc, "protein_length": len(seq), "substitution": key,
        "alphafold_model": Path(record["pdbUrl"]).name,
        "uniprot_features_containing_residue": features,
        "sequence_context": seq[max(0, pos - 11):pos + 10],
        "plddt": plddt,
        "relative_solvent_accessibility": rsa,
        "alphamissense": scores.get(key),
        "alphamissense_all_substitutions_at_position": dict(sorted(scores.items())),
    }
    (out / f"{acc}_{key}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items()
                      if k != "alphamissense_all_substitutions_at_position"}, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
