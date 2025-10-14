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
email me (<xiaotao.shen@outlook.com>) or reach out me via the social medias below.

<i class="fa fa-weixin"></i>
[jaspershen1990](https://www.jaspershenlab.com/image/wechat_QR.jpg)

📩 <shenxt@stanford.edu>

<i class="fa fa-twitter"></i>
[Twitter](https://twitter.com/xiaotaoshen1990)

<i class="fa fa-home"></i>
[Shen Lab website](https://www.shen-lab.org/)

<i class="fa fa-map-marker-alt"></i> [59 Nanyang Dr, Singapore, Singapore 636921](https://www.google.com/maps/place/Alway+Building/@37.4322345,-122.1770883,17z/data=!3m1!4b1!4m5!3m4!1s0x808fa4d335c3be37:0x9057931f3b312c29!8m2!3d37.4322345!4d-122.1748996)

# Citation

If you use `mapa` in your publications, please cite this paper:

Yifei Ge, Feifan Zhang, Yijiang Liu, Chao Jiang, Peng Gao, Nguan Soon Tan, Sai Zhang, Yuchen Shen, Qianyi Zhou, Xin Zhou, Chuchu Wang, Xiaotao Shen. Leveraging Large Language Models for Redundancy-Aware Pathway Analysis and Deep Biological Interpretation. 

[Weblink](https://www.biorxiv.org/content/10.1101/2025.08.23.671949v1.abstract)

Thanks very much!
