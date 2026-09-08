# real_data — Reproduction Code for the mixFunMap Manuscript

This directory contains the reproduction scripts for the five analyses in the
mixFunMap manuscript. Each script has two parts:

1. **Analysis (`full` mode)**: reruns model fitting and the genome-wide scans
   from the analysis-ready inputs in `inputs/`, reproducing the results in
   `results/`.
2. **Figures (default mode)**: redraws every manuscript panel (Figures 1–5 and
   Supplementary Figures S1–S5, PDF + TIFF) directly from the frozen results
   in `results/`, without refitting any model.

## Directory layout

```
real_data/
├── 01_simulation.R            # Formal V15 logistic simulation (Figure 1, S1–S5)
├── 02_wheat_height.R          # FIP1 wheat plant-height GWAS (Figure 2)
├── 03_staph_mic0.R            # S. aureus MIC=0 growth-curve GWAS (Figure 3)
├── 04_rice_diversity.R              # Rice RIL two-environment analysis (Figure 4)
├── 05_wheat_canopy_cover.R    # FIP1 wheat canopy-cover LOP analysis (Figure 5)
├── inputs/                    # Analysis-ready inputs (rds / csv)
│   ├── wheat_height/          #   FPWW012 height phenotypes + genotypes
│   │                          #   (incl. 0/1/2 hard calls for Functional Mapping)
│   ├── staph_mic0/            #   growdata_mic0.csv / snpdata.csv / FASTSTRUCTURE.csv
│   ├── rice_diversity/              #   RiceCGM/RDP1 two-environment joint inputs
│   └── wheat_canopy_cover/    #   FPWW012 canopy-cover LOP inputs
├── results/                   # Frozen analysis results (read by the figure code;
│   │                          #   can be rebuilt by full mode)
│   ├── simulation/summary/    #   Simulation metrics summary (17,200 datasets)
│   ├── wheat_height/          #   Three-method scans, thresholds, top-4 rds
│   ├── staph_mic0/            #   Three-method scans, top-4 interpretability rds
│   ├── rice_diversity/              #   control / lowwater / gxe and benchmark subdirs
│   ├── wheat_cc/              #   LOP degree selection, scans, comparators,
│   │                          #   region annotation
│   └── annotation/            #   Table S significant-loci annotation (v2 layout,
│                              #   used for gene labels in the figures)
└── figures/                   # All panels written by the figure code (PDF + TIFF)
```

## Usage

By default every script only redraws the figures (seconds):

```bash
Rscript 01_simulation.R
Rscript 02_wheat_height.R
Rscript 03_staph_mic0.R
Rscript 04_rice_diversity.R
Rscript 05_wheat_canopy_cover.R
```

To rerun the full analyses (time-consuming; requires the mixFunMap package
and GMMAT):

```bash
Rscript 01_simulation.R full [cores] [n_rep_null] [n_rep_alt]
Rscript 02_wheat_height.R full [cores] [output_dir]
Rscript 03_staph_mic0.R  full [cores] [output_dir]
Rscript 04_rice_diversity.R    full [cores] [output_dir]
Rscript 05_wheat_canopy_cover.R full [cores] [output_dir]
```

Note: the formal simulation design comprises 28 null cells × 300 replicates
plus 44 alternative cells × 200 replicates = 17,200 datasets. Pass small
replicate counts for a pilot run (e.g. `full 4 2 1`).

## Requirements

- **Figure mode**: R ≥ 4.4 with `ggplot2`, `ggrepel`, `readxl`
  (01 additionally needs `gridExtra`).
- **Full mode**: additionally the mixFunMap package (repository root
  `./mixFunMap`, installed or on `.libPaths`) and the `GMMAT` package
  (minP benchmark).

## Model and threshold summary

| Analysis | Mean curve | mixFunMap setting | Thresholds |
|---|---|---|---|
| Simulation (V15) | Three-parameter logistic (1, 0.65, 7) | 3-df P3D Wald, fixed-τ² null | 0.05/m |
| Wheat height | Three-parameter logistic | Q = PC1–5 + VanRaden K, 3-df P3D Wald | Bonferroni 0.05/M_eff = 1.53e-5 (M_eff = 3,270); suggestive 1/M = 5.38e-5 (M = 18,583) |
| S. aureus MIC=0 | Standard logistic | Q = PC1–3 + K, 3-df P3D Wald | Bonferroni 0.05/M; suggestive 1/M |
| Rice RIL | Three-parameter logistic (single environment) / joint G×E | 3-df Wald per environment; joint `fit_mixfunmap_joint` 3-df Wald | Bonferroni 0.05/M_eff (M_eff = 6,904); suggestive 1/M (M = 33,697) |
| Wheat canopy cover | LOP (degree 4 by BIC) | Q = PC1–5 + K, 5-df P3D Wald | Same as wheat height |

The header of each script documents the exact model specification, QC steps
and data provenance for that analysis.

## Data sources

- **Wheat height / canopy cover**: FIP1 (GABI-WHEAT) FPWW012 population,
  90K SNP array.
- **S. aureus**: 99 strains, growth curves at 14 time points under MIC = 0
  (raw growdata / snpdata / FASTSTRUCTURE csv files).
- **Rice**: RiceCGM/RDP1 RIL population (349 lines, 21 days, 33,697 markers,
  control and low-water environments).
- **Simulation**: manuscript V15 configuration (P3D-Wald, relative-QTL
  direction (1, 1, −1), K-neutral causal markers, n = 100, target family
  size 5, Q/noise ratio 0.10, K/noise ratio 0.50, Henderson REML).
