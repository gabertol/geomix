# GEOMIX

GeoMix is an R package for multi-proxy geological data unmixing, developed aiming to sedimentary provenance but with potential to work in any multiproxy geological dataset.

It allows the integration of compositional and distributional datasets (e.g., detrital zircon ages, heavy minerals, petrography, grain size) within a unified unmixing framework using Non-Negative Least Squares (NNLS).

The package is designed to work with multiple data blocks, allowing simultaneous inversion of different provenance proxies.

Traditional provenance analysis often treats each proxy independently.

GeoMix allows the integration of multiple proxies simultaneously:

- Detrital zircon age distributions (KDE)
- Heavy mineral proportions
- Petrography (QFL)
- Grain size distributions
- Other compositional datasets

The model estimates:

<b>A matrix </b> → mixture proportions per sample

<b>B matrices </b> → endmember signatures for each proxy

This approach allows a consistent interpretation of sediment sources across multiple datasets.

# Installation

To install geomix directly from GitHub:

install.packages("devtools")
library(devtools)

devtools::install_github("gabertol/geomix")

Then load the package:

library(geomix)
