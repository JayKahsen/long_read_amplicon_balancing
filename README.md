# long_read_amplicon_balancing

Reproducible R code and supporting datasets for the manuscript:

**A rapid and flexible method for long-read amplicon library preparation and balancing**

This repository contains all scripts and data required to reproduce the processed matrices and Figures 2–8 from the manuscript. All analyses can be run directly from the included files after unzipping the required inputs.

---

## Contents

### R scripts
- Importing QIIME2 feature tables
- Building standardized working matrices
- Rarefaction and filtering
- Taxonomy summaries
- Alpha diversity statistics
- PCA ordinations
- PERMANOVA / PERMDISP testing
- Pipeline retention plots
- Differential-abundance–style visualizations (Figure 8)

### Data
- Metadata tables
- Matrix registry (`matrix_names.csv`)
- Zipped QIIME2 exports used to rebuild working matrices
- Zipped taxonomy and ASV reference files
- Pre-filtered, analysis-ready matrices
- Output folders for generated figures and tables

---

## Directory structure

long_read_amplicon_balancing/
├─ long_read_amplicon_balancing.Rproj
├─ _globalStuff.R
│
├─ 1_QIIME Tables to working files.R
├─ Figure_2_taxonomy.R
├─ Figure_3_6_alpha_diversity_stats_ANOVA.R
├─ Figure_4_7_PCA.R
├─ Figure_5_pipeline.R
├─ Figure_8_abundance.R
└─ helper.R
│
├─ data_tables/
│  ├─ meta.csv
│  ├─ meta_ALL.csv
│  ├─ matrix_names.csv
│  ├─ sample_reads.csv
│  ├─ R_package_citations.txt
│  ├─ ASV_sequence.zip
│  └─ ASV_taxa.zip
│  │
│  ├─ original/
│  │  ├─ Doubleton.zip
│  │  ├─ Feces.zip
│  │  ├─ Skin.zip
│  │  └─ Soil.zip
│  │
│  └─ filtered_matrix/
│     ├─ filtered_matrix_ASV_Set1.csv
│     ├─ filtered_matrix_Genus_Set1.csv
│     ├─ filtered_matrix_Phylum_Set1.csv
│     └─ filtered_matrix_Species_Set1.csv
│
├─ output_data/
└─ output_plot/

---

## Setup

### 0) Unzip required inputs (required)
Unzip the following before running any scripts:

- `data_tables/ASV_taxa.zip`

To rebuild working matrices from QIIME2 exports, also unzip:

- all files in `data_tables/original/`
- then run `1_QIIME Tables to working files.R`

The repository already includes pre-filtered matrices, so rebuilding is optional.

---

## Running the analyses

Each figure can be reproduced independently using the corresponding script:

1. `Figure_2_taxonomy.R`
2. `Figure_3_6_alpha_diversity_stats_ANOVA.R`
3. `Figure_4_7_PCA.R`
4. `Figure_5_pipeline.R`
5. `Figure_8_abundance.R`

Outputs are written automatically to:

- `output_plot/` (figures)
- `output_data/` (tables and intermediate results)

---

## Expected outputs

- Figures: `output_plot/`
- Tables and intermediate results: `output_data/`
- Intermediate matrices: `data_tables/raw_matrix/`
- Final working matrices: `data_tables/filtered_matrix/`

---

## Analysis notes

- Community analyses use CLR-transformed ASV data with Euclidean distance.
- CLR is applied after adding a small pseudo-count.
- PERMANOVA: `vegan::adonis2`
- PERMDISP: `vegan::betadisper` with permutation tests via `permutest`

---

## Data availability

This repository includes all processed tables required to reproduce the figures, along with the original QIIME2-derived tables (zipped). Metadata required for plotting and statistical testing is included.

Raw PacBio reads are available through NCBI SRA under the BioProject listed in the manuscript.

---

## Citation

Wu LYA#, Kahsen J#, Kunstman K, Naqib A, Green SJ.  
**A rapid and flexible method for long-read amplicon library preparation and balancing.**  
(Manuscript in preparation / under review)

---

## Contact

For questions or reuse requests:  
Jeremy Kahsen (Rush University Medical Center)  
Jeremy_Kahsen@Rush.edu
