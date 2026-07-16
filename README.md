# mETHYLotest

[![R package](https://img.shields.io/badge/R%20package-1.0.0-blue)](https://github.com/OMICShub-Imagine/mETHYLotest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Unified interactive pipeline for DNA methylation analysis.**

mETHYLotest is a wrapper package providing two analysis modules with
Shiny-based interfaces, automated QC, differential methylation detection,
signature validation, and HTML reporting.

It offers two main pipelines for your data and exposes individual
functions that you can integrate into your own workflows for automated
QC, reporting, interactive exploration, and more.

---

## Disclaimer

mETHYLotest has been developed and tested internally on small datasets
at Institut Imagine. As it is still in early development, updates will
be published as issues are identified and resolved.

Bug reports, feature requests, and feedback are welcome:
**[Open an issue](https://github.com/OMICShub-Imagine/mETHYLotest/issues)**

*mETHYLotest is a free R package and comes with ABSOLUTELY NO WARRANTY.*

## Modules

---

| Module | Status | Input | Wrapper |
|--------|--------|-------|---------|
| **EPIC** | Stable | Illumina 450K / EPICv1 / EPICv2 / Mouse IDAT files | [ChAMP](https://bioconductor.org/packages/ChAMP/) |
| **NGS** | Stable | WGBS, RRBS, Nanopore (modkit, f5c), PacBio (pb-CpG-tools) | [methylKit](https://bioconductor.org/packages/methylKit/) |
| **Episignatures** | Experimental | Existing EPIC project | |

---

## Quick Start

### Installation

#### 1. Prerequisites

mETHYLotest requires **R >= 4.2** and Bioconductor packages.

#### 2. methylKit

mETHYLotest has been created using [methylKit](https://bioconductor.org/packages/methylKit/) **1.36.0**, you can simply install it using :

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("methylKit"))
```

#### 3. ChAMP and ChAMPdata

ChAMP and ChAMPdata require specific patched versions to work with
mETHYLotest. We provide pre-built packages in the
[releases](https://github.com/OMICShub-Imagine/mETHYLotest/releases)
section.

**Download the following files** from the latest release:
- `ChAMPdata_2.31.1.tar.gz`
- `ChAMP_2.29.1.tar.gz`
- `kpmt_0.1.0.tar.gz`

Then install in this order:

```r
# 1. ChAMPdata (no dependencies)
install.packages("ChAMPdata_2.31.1.tar.gz",
                 repos = NULL, type = "source")

# 2. kpmt (ChAMP dependency)
if (!require("remotes")) install.packages("remotes")
remotes::install_local("kpmt_0.1.0.tar.gz",
                       dependencies = TRUE,
                       upgrade = "never")

# 3. ChAMP
remotes::install_local("ChAMP_2.29.1.tar.gz",
                       dependencies = TRUE,
                       upgrade = "never")
```

These packages include minor patches required for mETHYLotest compatibility.
If you rather install them from another technique here are the versions :
- ChAMP **2.29.1**
- ChAMPdata **2.31.1**

#### 4. Additional Bioconductor dependencies

```r
BiocManager::install(c(
  # Annotation (required for WGBS-LR annotation step)
  "annotatr",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",

  # Optional but recommended
  "Rtsne",
  "pROC"
))
```

#### 5. mETHYLotest

```r
if (!require("remotes")) install.packages("remotes")

remotes::install_github("OMICShub-Imagine/mETHYLotest")
```

## Quick Demo

### EPIC Arrays (450K / EPICv1 / EPICv2)

Test the full EPIC pipeline using example 450K IDATs from [ChAMPdata](https://bioconductor.org/packages/ChAMPdata/) (8 samples: 4 controls, 4 cases — lung dataset):

```r
library(mETHYLotest)
mETHYLotest.EPIC.pipeline("demo")
```

Opens the EPIC Project Setup UI with IDAT directory and sample sheet **pre-filled**.
Click **"Load & Validate Sample Sheet"**, configure parameters, and run.

> See the [EPIC Tutorial](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/mETHYLotest_EPIC-Tutorial.html) for full documentation.

### NGS (WGBS / Nanopore / PacBio)

Test the full NGS pipeline with built-in synthetic data (10 samples, no setup required):

```r
library(mETHYLotest)
mETHYLotest.NGS.pipeline("demo")
```

Opens the Project Setup UI with a pre-loaded dataset (5 controls + 5 cases).
Click **"Load & Validate Metadata"**, configure parameters, and run.

> See the [NGS Tutorial](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/mETHYLotest_NGS-Tutorial.html) for full documentation.

## Pipelines

### EPIC Array Analysis

```r
library(mETHYLotest)

results <- mETHYLotest.EPIC.pipeline()

# Or on an existing project
results <- mETHYLotest.EPIC.pipeline("/path/to/My_EPIC_Project")
```

#### Example Report

Reports generated from the EPIC demo dataset:

| Report | Description | Link |
|--------|-------------|------|
| **QC Report** | PCA, t-SNE, Beta Distribution, Sample Correlation, SVD, Detection P-val, Sex Prediction | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/EPIC_DEMO_QC_report.html) |
| **Validation Report** | SVM, PCA, heatmap, silhouette, top loci per signature | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/EPIC_DEMO_Validation_Report.html) |
| **Final Report** | Project summary, diff methylation, Signatures and Validation | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/EPIC_DEMO_Report.html) |


[View example EPIC Report](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/EPIC_DEMO_Report.html){target="_blank"}

### NGS Methylation Analysis

```r
library(mETHYLotest)

results <- mETHYLotest.NGS.pipeline()

# Or on an existing project
results <- mETHYLotest.NGS.pipeline("/path/to/My_NGS_Project")
```

#### Example Reports

Reports generated from the NGS demo dataset:

| Report | Description | Link |
|--------|-------------|------|
| **QC Report** | Coverage stats, control genomes, global methylation, QC flags | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/NGS_DEMO_QC_Report.html) |
| **Validation Report** | SVM, PCA, heatmap, silhouette, top loci per signature | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/NGS_DEMO_Validation_Report.html) |
| **Final Report** | Project summary, diff methylation, annotation, file inventory | [View](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/NGS_DEMO_Final_Report.html) |

### Episignature Scoring on EPIC Array (Experimental)

```r
# Run on an existing EPIC project
scores <- mETHYLotest.EPIC.Episignatures("/path/to/My_EPIC_Project")
```

---

## EPIC Module

```
champ.load()           Import IDAT files
      |
  Harmonize             Mixed-array support (450K + EPICv1 + EPICv2)
      |
  Pre-QC                Quality control
      |
champ.impute()         Imputation
      |
champ.norm()           Normalisation (BMIQ / PBC / SWAN / illumina)
      |
  Post-QC               Quality control
      |
champ.refbase()        Cell type correction (blood, optional)
      |
champ.runCombat()      Batch correction (optional)
      |
  Post-Combat QC        Quality control
      |
champ.DMP()            Differentially Methylated Positions
      |
      +-- Default signatures (top 500 CpGs by score)
      +-- DMP Explorer UI (interactive Shiny)
      +-- User signatures (created in UI)
      +-- Signature validation (SVM / PCA / Silhouette)
      +-- Validation UI (interactive Shiny)
      +-- HTML reports
      |
champ.DMR()            Differentially Methylated Regions
      |
champ.Block()          Block methylation
      |
champ.GSEA()           Gene set enrichment
      |
champ.CNA()            Copy number aberrations
```

### Input

**Metadata** (Excel `.xlsx`):

| Column | Required | Description |
|--------|----------|-------------|
| `Sample_Name` | Yes | Unique identifier |
| `Sample_Plate` | Yes | `450K` / `EPICv1` / `EPICv2` / `Mouse` |
| `Sample_Group` | Yes | Biological group (CTL, CASE...) |
| `Sentrix_ID` | Yes | BeadChip barcode |
| `Sentrix_Position` | Yes | Position on chip (R06C01) |

Additional columns (`Batch`, `Sex`, `Age`...) are available for batch
correction and covariate selection.

**+ IDAT directory** containing all `.idat` files.

---

## NGS Module

```
methRead()              Import (Bismark / modkit / f5c / pb-CpG-tools)
      |
  Offset correction      0-based BED to 1-based (optional)
      |
  Iterative QC           Coverage filtering, sample/chr exclusion
      |
  normalizeCoverage()    Coverage normalization (optional)
      |
  unite()                Intersect CpGs across samples
      |
  Clustering             Correlation + hierarchical clustering
      |
  calculateDiffMeth()    Multi-scenario differential methylation
      |                  (Basic / Batch / Covariate / Full)
      |
  tileMethylCounts()     Tiling windows (optional, region-level DMRs)
      |                  DMP-based confidence filtering
      |                  (High / Medium / Low / Unsupported)
      |
  methSeg()              Segmentation (optional, UMR/LMR/PMD/HMR)
      |
  Signatures             Default (top 500) + user-created
      |
  Validation             SVM / PCA / Silhouette
      |
  Annotation             Genomic annotation (annotatr)
      |
  Final Report           HTML
```

### Input

**Metadata** (Excel `.xlsx`):

| Column | Required | Description |
|--------|----------|-------------|
| `Sample_ID` | Yes | Unique identifier |
| `Path` | Yes | Absolute path to methylation call file |
| `Treatment` | Yes | `0` (Control) or `1` (Case) |

Additional columns available for batch/covariate adjustment.

### Supported Formats

| Technology | Tool | Preset |
|---|---|---|
| Short-read WGBS/RRBS | Bismark | Standard pipelines |
| Oxford Nanopore | modkit | `modkit` preset |
| Oxford Nanopore | f5c / Nanopolish | `f5c` preset |
| PacBio HiFi | pb-CpG-tools | `pbcpg` preset |

---

## Interactive UIs

Six Shiny interfaces with consistent styling:

| UI | Module | Purpose |
|----|--------|---------|
| Project Setup | EPIC / NGS | Configuration, sample selection |
| DMP Explorer | EPIC | Volcano, ChromoMap, heatmap, signatures |
| DiffMeth Explorer | NGS | Volcano, DMC table, signatures |
| QC UI | NGS | Iterative QC with live adjustments |
| Signature Validator | Both | SVM metrics, PCA, benchmarking |

---

## Project Structure

```
My_Project/
├── data/
│   └── selected_samples.xlsx
└── Results/
    ├── project_config.R
    ├── interim/                   R objects (.rds)
    ├── QC*/                       Quality control plots
    ├── DMP/ or Differential_Analysis/
    ├── DMR/ or Tiling_Windows/
    │   ├── DMR_tiles_{scenario}.xlsx      Supported DMRs (≥1 DMP)
    │   ├── DMR_tiles_{scenario}all.xlsx  All DMRs with confidence
    │   ├── DMR_tiles{scenario}.bed       Supported DMRs for genome browser
    │   └── Full_tiles_{scenario}.xlsx     All tested regions
    ├── Episignatures/ or Signatures/
    ├── Validation/
    ├── DMP_Reports/ or Annotation/
    ├── Block/ | GSEA/ | CNA/     (EPIC only)
    ├── Segmentation/              (NGS only)
    └── *.html                     Reports
```

---

## Key Features

- **Automatic array detection** from `Sample_Plate`
- **Mixed-array harmonization** (EPICv1 + EPICv2 + 450K)
- **Smart loading** — skips import if cached `.rds` exists
- **Confounding detection** before batch correction
- **Multi-scenario analysis** (Basic / Batch / Covariate / Full)
- **DMR confidence filtering** — cross-references DMRs with DMPs to assign confidence levels (High / Medium / Low / Unsupported), reducing false positives especially with low sample sizes
- **Signature validation** with SVM, PCA, Silhouette
- **Automated HTML reports** at every stage

---

## DMR Confidence Levels (NGS Module)

When tiling windows are enabled, each DMR is cross-referenced with
individual DMPs (same scenario, same thresholds) to assess reliability.

| Confidence | Criteria | Interpretation |
|---|---|---|
| 🟢 **High** | ≥ 3 DMPs in region AND ≥ 25% difference | Strong signal confirmed by individual CpGs |
| 🟡 **Medium** | ≥ 1 DMP AND ≥ 15% difference | Likely real, consider for validation |
| 🟠 **Low** | ≥ 1 DMP | Weak but DMP-supported signal |
| 🔴 **Unsupported** | No DMP overlap | Likely false positive — excluded from main export |

**Why is this needed?**

The tiling approach aggregates read counts across CpGs within each
window, which increases statistical power but can inflate significance
— especially with few biological replicates.
DMRs without any individually significant CpG (DMP) lack independent
confirmation and should be interpreted with caution.

The main export (`DMR_tiles_*.xlsx`) contains only supported DMRs.
The complete set including unsupported regions is available in
`DMR_tiles_*_all.xlsx` for exploratory analysis.

---

## Tutorials

[How to use mETHYLotest with EPIC (.idat) data ?](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/mETHYLotest_EPIC-Tutorial.html)

[How to use mETHYLotest with NGS (WGBS, nanopore, pacbio...) data ?](https://htmlpreview.github.io/?https://github.com/OMICShub-Imagine/mETHYLotest/blob/main/inst/tutorials/mETHYLotest_NGS-Tutorial.html)

```r
library(mETHYLotest)

# EPIC tutorial
mETHYLotest.tutorial("EPIC")

# NGS tutorial
mETHYLotest.tutorial("NGS")
```

---

## Citation

If you use mETHYLotest in your research, please cite:

> Doldi N, Nitschké P, de Dieuleveult M, Puig Lombardi ME (2026).
> *mETHYLotest: a unified toolkit for reproducible multi-platform DNA methylation analysis.*
> R package version 1.0.0. Institut Imagine, Paris, France.

```r
citation("mETHYLotest")
```

### Please also cite the underlying tools

mETHYLotest wraps several published methods. Please acknowledge them
depending on which module(s) you used:

#### NGS Module

| Tool | Citation |
|------|----------|
| **methylKit** | Akalin A et al. (2012). *methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles.* Genome Biology, 13:R87. [doi:10.1186/gb-2012-13-10-r87](https://doi.org/10.1186/gb-2012-13-10-r87) |
| **annotatr** | Cavalcante RG, Sartor MA (2017). *annotatr: genomic regions in context.* Bioinformatics, 33(15):2381–2383. [doi:10.1093/bioinformatics/btx183](https://doi.org/10.1093/bioinformatics/btx183) |

#### EPIC Module

| Tool | Citation |
|------|----------|
| **ChAMP** | Tian Y et al. (2017). *ChAMP: updated methylation analysis pipeline for Illumina BeadChips.* Bioinformatics, 33(24):3982–3984. [doi:10.1093/bioinformatics/btx513](https://doi.org/10.1093/bioinformatics/btx513) |
| **minfi** | Aryee MJ et al. (2014). *Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays.* Bioinformatics, 30(10):1363–1369. [doi:10.1093/bioinformatics/btu049](https://doi.org/10.1093/bioinformatics/btu049) |

#### Validation (both modules)

| Tool | Citation |
|------|----------|
| **caret** | Kuhn M (2008). *Building predictive models in R using the caret package.* Journal of Statistical Software, 28(5):1–26. [doi:10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05) |
| **pROC** | Robin X et al. (2011). *pROC: an open-source package for R and S+ to analyze and compare ROC curves.* BMC Bioinformatics, 12:77. [doi:10.1186/1471-2105-12-77](https://doi.org/10.1186/1471-2105-12-77) |

---

## Authors

- **Nicolas Doldi** — Author, maintainer | Bioinformatics Core, SFR Necker – Institut Imagine, INSERM U1163, Paris, France | [![ORCID](https://img.shields.io/badge/ORCID-0009--0000--2816--9543-green)](https://orcid.org/0009-0000-2816-9543)
- **Patrick Nitschké** — Contributor | Bioinformatics Core, SFR Necker – Institut Imagine, INSERM U1163, Paris, France
- **Maud de Dieuleveult** — Contributor | INSERM U1163, Université de Paris, Imagine Institute, Paris, France
- **Maria Emilia Puig Lombardi** — Supervisor | Bioinformatics Core, SFR Necker – Institut Imagine, INSERM U1163, Paris, France

Funded by **Institut Imagine**.

---

## License

GPL (>= 3)
