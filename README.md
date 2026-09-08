# microbial-MR-and-host-transcriptomic-integration
Codes for paper "Species- and function-resolved gut microbial associations with colorectal cancer: Mendelian randomization and host transcriptomic integration"

# Reproducibility workflow

Run all commands from the repository root after loading `.env`. Intermediate
and result files are written only to repository-relative locations.

## 1. Input preparation and exposure inventory

The frozen exposure inventory is included in `config/`. To reconstruct it from
source metadata, place `swedish_gws_eligible_mr_exposures.tsv` under
`config/source_metadata/`, place the HUNT supplementary workbook at the path
used by `build_v2_exposure_inventory.py`, and run:

```bash
python3 scripts/build_v2_exposure_inventory.py
python3 scripts/download_missing_candidate_gwas.py --workers 4
python3 scripts/audit_gwas_coordinates.py
```

## 2. Descriptive LDSC and exposure selection

```bash
python3 scripts/run_candidate_ldsc.py --workers 2
python3 scripts/select_ldsc_exposures.py
```

LDSC estimates are descriptive. The candidate rule determines eligibility;
heritability estimates are used only to resolve biologically and
phenotype-model-equivalent duplicate GWAS.

## 3. Mendelian randomization

```bash
Rscript legacy_sources/run_strict_ref1_crc_mr.R
python3 scripts/finalize_v2_mr_results.py
Rscript scripts/build_v2_mr_artifacts.R
Rscript scripts/annotate_gmrgs_gprofiler.R
```

The MR script reads the relative `source_gwas_path` values in the frozen
manifest. PLINK and the European LD reference prefix are supplied through
`.env`.

## 4. Bulk expression and feature selection

```bash
Rscript scripts/run_v2_bulk_deg.R
Rscript scripts/run_v2_dynamic_v1_stage.R prepare87211
Rscript scripts/run_v2_patient_grouped_ml.R
Rscript scripts/run_v2_ml_stability_audit.R 200
```

## 5. External tissue-cohort model, GSEA and immune context

```bash
Rscript scripts/run_ridge_nomogram_external_validation.R
Rscript scripts/plot_ridge_nomogram_figure3.R
Rscript scripts/run_v2_dynamic_v1_stage.R gsea
Rscript scripts/finalize_v2_figure4_tables.R
CIBERSORT_PERMUTATIONS=100 CIBERSORT_CORES=8 \
  Rscript scripts/run_v2_dynamic_v1_stage.R immune
```

## 6. GSE178341 broad cell-type analysis

```bash
python3 scripts/build_gse178341_label_free_sketch.py
Rscript scripts/cluster_gse178341_label_free_sketch.R
python3 scripts/annotate_gse178341_label_free_sketch.py
python3 scripts/map_gse178341_label_free_cell_types.py
python3 scripts/prepare_gse178341_label_free_pseudobulk.py
python3 scripts/export_gse178341_label_free_edger.py
Rscript scripts/run_gse178341_label_free_edger.R
python3 scripts/build_gse178341_label_free_qc_sources.py
python3 scripts/build_gse178341_label_free_panels.py
python3 scripts/build_gse178341_figure6.py
python3 scripts/build_gse178341_wfdc2_supplementary.py
```

## 7. GSE178341 epithelial-state analysis

```bash
python3 scripts/build_gse178341_label_free_epithelial_sketch.py
Rscript scripts/cluster_gse178341_label_free_epithelial.R
python3 scripts/annotate_map_gse178341_label_free_epithelial.py
python3 scripts/prepare_gse178341_label_free_figure7.py
Rscript scripts/run_gse178341_label_free_cellchat.R
python3 scripts/postprocess_gse178341_label_free_figure7.py
Rscript scripts/fit_gse178341_biomarker_pseudotime_gamm.R
python3 scripts/build_gse178341_label_free_figure7.py
python3 scripts/build_gse178341_label_free_figure_s5.py
python3 scripts/export_label_free_single_cell_supplementary_tables.py
```

CellChat and graph-pseudotime results are descriptive tissue-context analyses.
Generated tables and figures are written to `results/` and `manuscript/`.
