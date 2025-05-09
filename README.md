
<!-- README.md is generated from README.Rmd. Please edit that file -->

# atacama.waxes

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/584162411.svg)](https://doi.org/10.5281/zenodo.15354457)
[![License: GPL (\>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![Dependencies](https://img.shields.io/badge/dependencies-2/94-green?style=flat)](#)
<!-- badges: end -->

Research Compendium for the paper **Hydroclimate variations over the
last 17,000 years as estimated by leaf waxes in rodent middens from the
south-central Atacama Desert, Chile**

### How to cite

Please cite this compendium as:

> Frugone-Álvarez, M., Contreras, S., Meseguer-Ruiz, O., Tejos, E.,
> Delgado-Huertas, A., Valero-Garcés, B., Díaz, F.P., Briceño, M.,
> Bustos-Morales, M., Latorre, C., 2023. Hydroclimate variations over
> the last 17,000 years as estimated by leaf waxes in rodent middens
> from the south-central Atacama Desert, Chile. Quat. Sci. Rev. 311,
> 108084. <https://doi.org/10.1016/j.quascirev.2023.108084>

## Abstract

Leaf cuticular waxes are one of the most important environment-plant
interaction structural systems that enable desert plants to withstand
extreme climatic conditions. We present a long chain n-alkyl lipids
study in fresh plant leaves and rodent palaeomiddens collected along an
elevational gradient in the south-central Atacama Desert of Chile,
covering six different vegetation belts: Steppe (4500-4000 m asl), Puna
(4000-3300 m asl), pre-Puna (3300-2400 m asl), Absolute Desert
(2400-1000 m asl) and Coastal Desert (1000-0 m asl). The 28 rodent
palaeomiddens analyzed from Quebrada Incahuasi (25.6 °S, 3600 m asl)
span the last 17,000 years. Modern-day distribution of long-chain
n-alkanes and n-alkanoic acids varies among the dominant plant
associations of the Atacama Desert. These plants show a species-specific
chemotaxonomy linked to the climatic conditions. Furthermore,
differences in average chain length (ACL) and carbon preference index
(CPI) suggest that these plant communities are highly adapted to extreme
environmental conditions. The sum of leaf wax n-alkanes was highest
under wet conditions, while n-alkanoic acids (between n-C24 and n-C28)
increased with hyperaridity. Similarly, analysis of n-alkane time series
from palaeomiddens showed that the greatest changes in leaf wax n-alkane
distributions (ACL and CPI) corresponded to the greatest increases in
moisture during the Central Andean Pluvial Event (CAPE; between 18 and 9
ka cal BP) and the Late Holocene. The shift in the palaeomidden n-alkane
distributions is corroborated by the relative abundance of
rainfall-dependent extra-local taxa. This is the first study to report
leaf wax content obtained from ancient rodent middens, and shows
promising results as a robust hydroclimate proxy for the Atacama Desert
region.

### Content

This repository is structured as follow:

- [`data/`](https://github.com/mat1506/atacama.waxes/tree/master/data):
  contains all raw data required to perform analyses

- [`analyses/`](https://github.com/mat1506/atacama.waxes/tree/master/analyses/):
  contains R scripts to run each step of the workflow

- [`R/`](https://github.com/mat1506/atacama.waxes/tree/master/R):
  contains R functions developed especially for this project

- [`man/`](https://github.com/mat1506/atacama.waxes/tree/master/man):
  contains help files of R functions

- [`DESCRIPTION`](https://github.com/mat1506/atacama.waxes/tree/master/DESCRIPTION):
  contains project metadata (author, date, dependencies, etc.)

- [`make.R`](https://github.com/mat1506/atacama.waxes/tree/master/make.R):
  main R script to run the entire project by calling each R script
  stored in the `analyses/` folder

### Usage

- Clone this repository
- Open a terminal
- Build the Docker image with:

``` sh
docker build -t "atacama.waxes" .
```

- Start a container based on this image:

``` sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true atacama.waxes
```

- On a web browser enter this URL: `127.0.0.1:8787`. A new RStudio
  Server instance will be available.
- To run the analysis:

``` r
source("make.R")
```

- Clone this repository
- Open a terminal
- Build the Docker image with:

``` sh
docker build -t "atacama.waxes" .
```

- Start a container based on this image:

``` sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true atacama.waxes
```

- On a web browser enter this URL: `127.0.0.1:8787`. A new RStudio
  Server instance will be available.
- To run the analysis:

### Notes

- All required packages, listed in the `DESCRIPTION` file, will be
  installed (if necessary)
- All required packages and R functions will be loaded
- Some analyses listed in the `make.R` might take time
