<h1 align="center">Joint modeling of response accuracy and time in process data: A multilevel latent class IRT model</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Francesco Bartolucci &middot; Fulvia Pennoni &middot; Giorgio Vittadini</strong></em></span> </p>
<br>

The `RTMLC` package contains functions to specify and estimate Response Time Multilevel Latent Class (RTMLC) model for the joint analysis of response accuracy and response time in item-response data. The model is designed for dichotomously scored items administered to individuals nested within higher-level units, such as students within school classes.

The proposed framework relies on discrete latent variables at both the individual and cluster levels, allowing unobserved heterogeneity to be represented through latent classes. Response accuracy is modeled through a normal-ogive IRT sub-model, whereas response times are modeled by a log-normal sub-model. Individual- and cluster-level covariates can be included in both sub-models.

The code implements maximum likelihood estimation through the Expectation-Maximization algorithm. It also allows for structurally missing item responses arising from the random assignment of item subsets from a larger item pool.

To install the `RTMLC` package directly from GitHub:
```r
# install.packages("devtools")
require(devtools)
devtools::install_github("LB1304/RTMLC")
```

To download the .tar.gz file (for manual installation), use [this link](https://github.com/LB1304/RTMLC/archive/main.tar.gz).
