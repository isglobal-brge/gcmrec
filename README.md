# gcmrec

<!-- badges: start -->
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->

General class of models for recurrent event data.

`gcmrec` fits the general class of semiparametric models for recurrent event
data of Peña and Hollander (2004). Unlike standard survival models, it
accounts for what an event does to the subject through an *effective age*
function, for the effect of accumulating occurrences, and for unobserved
heterogeneity through gamma frailties. It applies to settings such as repeated
hospital readmissions, cancer relapses or successive failures of a machine,
and it also fits the extension for cancer relapses of González, Peña and
Slate (2005).

## Installation

```r
install.packages("gcmrec")
```

Development version:

```r
# install.packages("remotes")
remotes::install_github("isglobal-brge/gcmrec")
```

## Usage

```r
library(gcmrec)
data(readmission)

fit <- gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
              data = readmission, s = 3000)

summary(fit)      # hazard ratios with confidence intervals
plotForest(fit)   # and the corresponding figure
plot(fit)         # baseline survivor function
```

Beyond fitting, the package provides the mean cumulative function of the
events (`mcf()`), event charts (`graph.caltimes()`), predictions for covariate
profiles (`predict()`, `plotPredict()`) and likelihood ratio tests between
nested models (`anova()`).

The vignette walks through a complete analysis:

```r
vignette("gcmrec")
```

## References

Peña, E.A. and Hollander, M. (2004). Models for recurrent events in
reliability and survival analysis. In *Mathematical Reliability: An Expository
Perspective*, chapter 6, pp. 105–123. Kluwer Academic Publishers.

González, J.R., Peña, E.A. and Slate, E.H. (2005). Modelling intervention
effects after cancer relapses. *Statistics in Medicine*, 24(24), 3959–3975.
[doi:10.1002/sim.2410](https://doi.org/10.1002/sim.2410)

## License

GPL (>= 2)
