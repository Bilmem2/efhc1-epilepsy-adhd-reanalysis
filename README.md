# EFHC1 and epilepsy–ADHD comorbidity: single-case exome reanalysis

Analysis code and summary results for the manuscript

> *EFHC1 and Epilepsy–ADHD Comorbidity: A Gene-Level Pleiotropy Hypothesis from a Single-Case Exome Reanalysis* (under review, American Journal of Undergraduate Research).

The repository documents the analysis: the variant prioritization rules, the number of variants left after each step, the reason for each manual decision, and the scripts and public inputs behind the protein, ClinVar and single-nucleus results and the figures.

## What is and is not included

**Included:** the analysis scripts; the prioritization rules (in `scripts/filter_opencravat_export.py`); the counts after each step; the manual review decisions; the six retained candidates; the frequency checks in two Turkish population databases (TR-Variome and TGD); the public input subsets; the protein, ClinVar and single-nucleus results; the figures.

**Not included:** the proband's exome VCF and the full OpenCRAVAT annotation table. An exome is identifying genetic information even without names, and the written consent covers publication of anonymized findings, not distribution of the genome. The case-specific steps (1 and 2 below) therefore cannot be rerun from this repository; all steps based on public data can. One gene set aside in the manual review is a cancer-predisposition gene unrelated to the research question; it is not named in the public decision table.

## Requirements

Python 3.10 or later.

```bash
pip install -r requirements.txt
```

Biopython is needed only for the solvent accessibility in step 4. Without the manuscript's typeface (garamondx, converted for matplotlib by `scripts/convert_garamondx_ttf.py`) the figures are drawn in a default serif font; their content is the same.

## Workflow

| Step | Script | Input | Output | Manuscript |
|---|---|---|---|---|
| 1. Variant prioritization | `filter_opencravat_export.py` | OpenCRAVAT export (private) | `results/prioritization/funnel_summary.json` (counts only) | Figure 1, Results |
| 2. Manual review (20 to 6) | `record_manual_review.py` | step 1 output (private), decision table | `results/manual_review/` | Table 2, Results |
| 3. *EFHC1* in ClinVar | `efhc1_variant_landscape.py` | `data/public/` | `results/efhc1_clinvar/` | Results |
| 4. Protein context | `protein_variant_context.py` | UniProt, AlphaFold DB, AlphaMissense | `results/protein_context/` | Results, Discussion |
| 5. Single-nucleus expression | `download_hpa_single_nuclei.py`, `hpa_single_nuclei_coexpression.py` | Human Protein Atlas v26 | `results/single_nucleus/` | Figure 3, Results |
| 6. Figures | `plot_prioritization_funnel.py`, `plot_efhc1_protein.py`, `plot_single_nuclei_coexpression.py` | results of steps 1, 4, 5 | `results/figures/` | Figures 1–3 |
| 7. Turkish population frequencies | `turkish_variome_lookup.py` (TR-Variome); manual lookup in the TGD portal, which has no download interface | TR-Variome table (figshare); TGD portal | `results/population/` | Table 2, Results |

## Reproducing the public steps

```bash
bash run_public_steps.sh
```

The script runs steps 3–7 in order and overwrites the files in `results/`. It needs an internet connection: step 4 queries AlphaFold DB, and step 5 downloads about 140 MB from the Human Protein Atlas and checks the files against the SHA-256 sums of the copies used for the manuscript.

Step 7 runs only if the TR-Variome table (2 GB) is present, because the script does not download it:

```bash
mkdir -p data/turkish_variome
curl -L -o data/turkish_variome/TurkishVariome.tsv.gz https://figshare.com/ndownloader/files/30739108
```

The copy used for the manuscript has the SHA-256 sum `8067a71c9b589074a7fdbb82804190c746dd4768b11be3037734c21096736b54`.

## Case-specific steps (require the private data)

```bash
python scripts/filter_opencravat_export.py <export_variant.tsv> <private_dir>
python scripts/record_manual_review.py <private_dir>/model_A_dominant.tsv <decisions.tsv> <private_dir>/manual_review.tsv
```

The OpenCRAVAT (3.2.3) export is the unfiltered variant table with the gnomAD v4, All of Us, ClinVar, Human Phenotype Ontology, REVEL and CADD annotations. Step 1 writes the counts after each filter to `funnel_summary.json`; its tables `model_A_dominant.tsv` and `model_B_recessive.tsv` contain the individual's genotypes and must not be shared. Step 2 checks that the decision table covers exactly the genes in `model_A_dominant.tsv`, so it needs the private decision table in which all 20 genes are named; the public `results/manual_review/manual_review_decisions.tsv` leaves one gene unnamed and does not pass this check.

## Data sources and versions

| Resource | Version or date |
|---|---|
| OpenCRAVAT | 3.2.3 (annotation only) |
| gnomAD | v4 (OpenCRAVAT module); v4.1 exomes for `data/public/gnomad41_exomes_EFHC1_region.tsv` |
| All of Us | OpenCRAVAT module |
| ClinVar | GRCh38 VCF, 2 November 2026 (`data/public/clinvar_EFHC1_records.vcf`, chr6:52,355,000–52,535,000) |
| Reference transcripts | MANE Select (*EFHC1*: ENST00000371068) |
| Ensembl VEP | release 115, used only as an independent annotation check; its per-variant output is not shared |
| Human Protein Atlas | v26, single-nuclei brain data (Siletti *et al.*, 2023) |
| UniProt, AlphaFold DB, AlphaMissense | records retrieved in November 2026 (copies in `results/protein_context/inputs/`) |
| Turkish Variome (TR-Variome) | Kars *et al.*, 2021 (3,362 individuals); figshare dataset 15147642, file `TurkishVariome.tsv.gz` |
| Turkish Genome Project database (TGD) | release 3.1, checked 13 November 2026 |

## Licence

The code is released under the MIT licence (`LICENSE`). Public data subsets and records retrieved from ClinVar, gnomAD, UniProt, AlphaFold DB, the Human Protein Atlas and TR-Variome remain subject to the terms of their sources.
