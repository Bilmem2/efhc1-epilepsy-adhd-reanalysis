#!/usr/bin/env bash
# Reruns every step that uses only public data (steps 3-7 in README.md;
# step 7 only if the TR-Variome table has been downloaded, see README.md).
# Steps 1 and 2 need the proband's exome data, which is not shared.
set -euo pipefail
cd "$(dirname "$0")"

# Use $PYTHON if set, otherwise the first python3/python that is 3.10 or later
if [ -z "${PYTHON:-}" ]; then
  for candidate in python3 python; do
    if "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 10))' 2>/dev/null; then
      PYTHON="$candidate"
      break
    fi
  done
fi
if [ -z "${PYTHON:-}" ]; then
  echo "Python 3.10 or later not found; set PYTHON to its path." >&2
  exit 1
fi

echo "== 3. EFHC1 in ClinVar"
"$PYTHON" scripts/efhc1_variant_landscape.py data/public results/efhc1_clinvar

echo "== 4. Protein context"
"$PYTHON" scripts/protein_variant_context.py Q5JVL4 E604D results/protein_context
"$PYTHON" scripts/protein_variant_context.py P31321 R96W results/protein_context

echo "== 5. Single-nucleus expression"
"$PYTHON" scripts/download_hpa_single_nuclei.py data/hpa_v26
"$PYTHON" scripts/hpa_single_nuclei_coexpression.py data/hpa_v26 results/single_nucleus

echo "== 6. Figures"
"$PYTHON" scripts/plot_prioritization_funnel.py results/prioritization/funnel_summary.json results/figures/Fig1_prioritization
"$PYTHON" scripts/plot_efhc1_protein.py results/protein_context/inputs/uniprot_Q5JVL4.json results/figures/Fig2_efhc1_protein
"$PYTHON" scripts/plot_single_nuclei_coexpression.py results/single_nucleus results/figures/Fig3_single_nuclei_coexpression

echo "== 7. Turkish population frequencies (TR-Variome)"
VARIOME=data/turkish_variome/TurkishVariome.tsv.gz
if [ -f "$VARIOME" ]; then
  "$PYTHON" scripts/turkish_variome_lookup.py "$VARIOME" results/population/tgd_candidate_check.tsv results/population/trvariome_candidate_check.tsv
else
  echo "skipped: $VARIOME not found (download command in README.md)"
fi
