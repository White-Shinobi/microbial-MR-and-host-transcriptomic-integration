# Public data sources and expected locations

Raw data are not redistributed in this repository. Download them from the
original public repositories and retain the accession identifiers in all local
records.

| Resource | Accession or identifier | Role |
|---|---|---|
| GWAS Catalog | GCST90013862 | Overall colorectal cancer outcome GWAS |
| GWAS Catalog | Accessions in `config/v2_candidate_exposure_inventory.tsv` | Swedish and HUNT microbial exposure GWAS |
| NCBI GEO | GSE44076 | Bulk-tissue discovery cohort |
| NCBI GEO | GSE39582 | Bulk-tissue processing dependency retained by the source workflow |
| NCBI GEO | GSE87211 | Independent expression-direction validation |
| NCBI GEO | GSE41258 | External tissue-cohort classifier evaluation |
| NCBI GEO | GSE37364 | External tissue-cohort classifier evaluation |
| NCBI GEO | GSE178341 | Paired colorectal single-cell analysis |
| 1000 Genomes Project | Phase 3 EUR | LD clumping reference |
| MSigDB | KEGG collection used in the completed run | Gene-set enrichment input |
| CIBERSORT/LM22 | LM22 signature matrix | Immune deconvolution input |

Exact relative input paths, file sizes, and checksums from the completed run
are in `config/input_manifest.tsv` and
`config/analysis_input_checksums.sha256`. Public-source licenses and access
conditions continue to apply; this repository does not grant redistribution
rights for third-party data.
