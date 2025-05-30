<!-- README.md is generated from README.Rmd. Please edit that file -->

# mapa <img src="man/figures/mapa_logo.png" align="right" alt="" width="120" />
<!-- badges: start -->
[![R-CMD-check](https://github.com/jaspershen-lab/mapa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jaspershen-lab/mapa/actions/workflows/R-CMD-check.yaml)
[![](https://www.r-pkg.org/badges/version/mapa?color=green)](https://cran.r-project.org/package=mapa)
[![](https://img.shields.io/github/languages/code-size/jaspershen/mapa.svg)](https://github.com/jaspershen/mapa)
[![Dependencies](https://tinyverse.netlify.com/badge/mapa)](https://cran.r-project.org/package=mapa)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->
------

# About

`mapa` is a package used to merge redundant pathways and GO terms for transcriptomics and proteomics data.


# Installation

You can install `mapa` from [GitHub](https://github.com/jaspershen-lab/mapa)

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

remotes::install_github(
  "jaspershen-lab/mapa",
  dependencies = TRUE,
  repos        = BiocManager::repositories(),
  upgrade      = "ask" 
)
```

More information can be found [here](https://jaspershen.github.com/mapa).

# Get started

Please see the `Help documents`.

# Need help?

If you have any questions about `mapa`, please don’t hesitate to
email me (<shenxt1990@outlook.com>) or reach out me via the social medias below.

<i class="fa fa-weixin"></i>
[shenxt1990](https://www.shenxt.info/files/wechat_QR.jpg)

<i class="fa fa-envelope"></i> <shenxt@stanford.edu>

<i class="fa fa-twitter"></i>
[Twitter](https://twitter.com/JasperShen1990)

<i class="fa fa-map-marker-alt"></i> [M339, Alway Buidling, Cooper Lane,
Palo Alto, CA
94304](https://www.google.com/maps/place/Alway+Building/@37.4322345,-122.1770883,17z/data=!3m1!4b1!4m5!3m4!1s0x808fa4d335c3be37:0x9057931f3b312c29!8m2!3d37.4322345!4d-122.1748996)

# Citation

If you use `mapa` in your publications, please cite this paper:

Shen, X., Yan, H., Wang, C. et al. TidyMass an object-oriented reproducible analysis framework for LC–MS data. Nat Commun 13, 4365 (2022). 

[Weblink](https://www.nature.com/articles/s41467-022-32155-w)

Thanks very much!
