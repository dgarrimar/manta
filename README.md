# MANTA <img src='man/figures/logo.png' align="right" height="139"/>

[![R-CMD-check](https://github.com/dgarrimar/manta/actions/workflows/check-full.yaml/badge.svg)](https://github.com/dgarrimar/manta/actions/workflows/check-full.yaml)
[![Codecov test coverage](https://codecov.io/gh/dgarrimar/manta/branch/master/graph/badge.svg)](https://app.codecov.io/gh/dgarrimar/manta)

The Multivariate Asymptotic Non-parametric Test of Association (MANTA) enables non-parametric, asymptotic P-value computation for multivariate linear models. 

## Installation 

```r
# install.packages("devtools")
devtools::install_github("dgarrimar/manta")
```

**R 3.3.2 or higher** is required.

## Usage

```r
library(manta)
manta(biomarkers ~ ., data = patients)
#>
#> Call:
#> manta(formula = biomarkers ~ ., data = patients)
#> 
#> Type II Sum of Squares
#> 
#>           Df Sum Sq Mean Sq F value      R2    Pr(>F)    
#> age        1  400.6  400.63  7.3566 0.04242  0.001685 ** 
#> gender     1   34.3   34.28  0.6295 0.00363    0.5144    
#> status     2 2152.7 1076.33 19.7643 0.22793 1.348e-12 ***
#> Residuals 91 4955.7   54.46                              
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 4 observations deleted due to missingness
```

## Cite MANTA

If you find MANTA useful in your research please cite the related publication:

Garrido-Martín, D., Calvo, M., Reverter, F., Guigó, R. A fast non-parametric test of association for multiple traits. *bioRxiv* (2022). [https://doi.org/10.1101/2022.06.06.493041](https://doi.org/10.1101/2022.06.06.493041)
