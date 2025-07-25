
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GencoDymo2: Comprehensive Analysis of GENCODE Annotations and Splice Site Motifs <img src="man/figures/GencoDymo2_logo.png" align="right" height="140"/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Licence: GPL
v3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/monahton)
[![GitHub
issues](https://img.shields.io/github/issues/monahton/GencoDymo2)](https://github.com/monahton/GencoDymo2/issues)
[![Platform](https://img.shields.io/badge/platform-all-green)](https://cran.r-project.org/)
[![Website](https://img.shields.io/badge/docs-website-blue)](https://monahton.github.io/GencoDymo2/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15302316.svg)](https://doi.org/10.5281/zenodo.15302316)
[![](https://www.r-pkg.org/badges/version/GencoDymo2?color=green)](https://cran.r-project.org/package=GencoDymo2)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GencoDymo2?color=green)](https://cran.r-project.org/package=GencoDymo2)  
[![name status
badge](https://monahton.r-universe.dev/badges/:name)](https://monahton.r-universe.dev/)
[![registry status
badge](https://monahton.r-universe.dev/badges/:registry)](https://monahton.r-universe.dev/)

[![Last
commit](https://img.shields.io/github/last-commit/monahton/GencoDymo2)](https://github.com/monahton/GencoDymo2/commits/main)

------------------------------------------------------------------------

## 📦 Overview

**GencoDymo2** is an R package tailored for dynamic extraction,
exploration, and comparison of gene annotations from the
[GENCODE](https://www.gencodegenes.org) database for human and mouse
genomes. This tool facilitates:

- Automated retrieval of the latest or specific GENCODE releases
- Visualization and quantification of annotation changes across versions
- Extraction of introns, exons, splice motifs
- Generation of splice site FASTA files for tools like MaxEntScan

It bridges bioinformatics workflows and annotation dynamics, enhancing
reproducibility and comparative studies in transcriptome and splicing
research.

------------------------------------------------------------------------

## 💻 Installation

You can install the **stable version** of `GencoDymo2` from
[CRAN](https://cran.r-project.org/package=GencoDymo2):

``` r
# Install the stable version from CRAN
install.packages("GencoDymo2")
```

Or you can install the **development version** from GitHub for the
latest features:

``` r
#Install the development version from GitHub
install.packages("pak")
pak::pkg_install("monahton/GencoDymo2")
```

``` r
# Load the package
library(GencoDymo2)
```

------------------------------------------------------------------------

## 🚀 Getting Started

To get started, view the vignette:

``` r
vignette("GencoDymo2")
```

Or visit the documentation website:  
👉 <https://monahton.github.io/GencoDymo2/>  
👉
<https://monahton.github.io/GencoDymo2/articles/GencoDymo2_vignette.html>

------------------------------------------------------------------------

## 📁 Functions Highlights

| Function | Description |
|----|----|
| `get_latest_release()` | Retrieves the latest available GENCODE release per species |
| `compare_release()` | compare annotation statistics between releases |
| `extract_introns()` | Extracts and processes introns from annotation |
| `assign_splice_sites()` | Assign the donor and acceptor splice sites |
| `extract_ss_motif()` | Extract splicing motifs for MaxEntScan tool |

------------------------------------------------------------------------

## 🛠️ Development & Contributing

**GencoDymo2** is actively developed. Contributions and suggestions are
welcome!

- 🔧 Open issues: <https://github.com/monahton/GencoDymo2/issues>
- 📬 Email: <aboualezz.monah@hsr.it>
- 🤝 Pull requests encouraged!

------------------------------------------------------------------------

## :writing_hand: Author

**Monah Abou Alezz, PhD** – <aboualezz.monah@hsr.it>.

San Raffaele Telethon Institute for Gene Therapy (SR-TIGET)  
IRCCS San Raffaele Scientific Institute, Milan, Italy

🌍 [Personal website](https://monahton.github.io)

[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/monahton)
[![](https://img.shields.io/badge/follow%20me%20on-LinkedIn-blue.svg)](https://linkedin.com/in/monah-abou-alezz-phd-06a948ba)

------------------------------------------------------------------------

## 🧪 Related Projects

- [GENCODE](https://www.gencodegenes.org)
- MaxEntScan

------------------------------------------------------------------------

## 📣 Acknowledgments

Developed as part of ongoing research on lncRNA splicing and gene
annotation evolution.  
Special thanks to colleagues at IGM-CNR and collaborators across
splicing research projects.
