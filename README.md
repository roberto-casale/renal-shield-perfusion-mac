# A consensus transcriptomic signature of renal ischemia–reperfusion injury

Code, results and figures of the article *A Consensus Transcriptomic Signature of Renal
Ischemia–Reperfusion Injury to Guide Ex-Vivo Graft Reconditioning*.

Nine public transcriptomic series of renal ischemia–reperfusion injury (human, mouse and
rat, from the Gene Expression Omnibus) are analysed one by one, mapped to human orthologs
and combined by robust rank aggregation into a consensus signature. The signature is
characterised by over-representation analysis and tested on an independent mouse series,
GSE267650, that took no part in its construction.

## From the article to the files

| Article | File | Notebook |
|---|---|---|
| Table 1 | `results/Table1_datasets.csv` | 1 |
| Figure 1 | `figures/Figure1_dataset_selection.png` (counts in `results/geo_identification.csv` and `results/geo_corpus_coverage.csv`) | 1 and 2 |
| Table 2 | `results/Table2_signature.csv` | 2 |
| Table 3 | `results/Table3_pathways.csv` | 2 |
| Table S1 | `results/TableS1_enrichment.csv` | 2 |
| Figure 2 | `figures/Figure2_interventions.png` | 2 |
| Figure 3 | `figures/Figure3_workflow.png` | 2 |
| Section 2.4: concordance, median oriented log2 fold change, permutation p-value | `results/validation_summary.csv` (gene by gene: `results/validation_genes_vs_controls.csv`, `results/validation_genes_vs_ischemia.csv`) | 3 |

In `results/Table2_signature.csv` the signature rule uses `FDR_posthoc` (Section 4.4); the
columns `FDR_up`, `FDR_down`, `meta_FDR` and `consensus` are intermediate values of the
rank aggregation.

Intermediate tables: `results/DE_<series>.csv` (one per corpus series, notebook 1);
`results/validation_permutations_*.csv` and `results/validation_design.csv` (notebook 3).
Sensitivity analysis on one-to-one orthologs: `results/orthologs_per_dataset.csv` and
`results/signature_one_to_one_sensitivity.csv` (notebook 2).

## Contents

```
notebook/    01_datasets.ipynb     GEO series -> per-dataset differential expression
             02_analysis.ipynb     signature, enrichment, figures
             03_validation.ipynb   external validation on GSE267650
R/           the code the notebooks run, one module per step
results/     tables (CSV)
figures/     figures (PNG)
data/cache/  cached answers of GEO and WebGestalt, so that a re-run gives the same result
renv.lock    versions of the R packages
```

## How to run

Requirements: R 4.4.3 with Bioconductor 3.20, and Python 3 with Jupyter.

```
Rscript R/00_install.R          # R packages into env/rlib (versions in renv.lock)
./run_all_macos.sh setup        # register the R kernel for Jupyter
./run_all_macos.sh all          # notebooks 1, 2 and 3
```

On Windows: `.\run_all.ps1` (setup, then the three notebooks). Notebook 1 downloads the
GEO series into `data/`; notebooks 2 and 3 can be run on their own once the files in
`results/` exist.

## Licence

Code: MIT. Results, figures and this README: CC BY 4.0. See `LICENSE`.
