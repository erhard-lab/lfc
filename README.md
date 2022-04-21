# lfc
R package for modelling count ratio data

## How to get started

Install the R package using the following commands on the R console:

```
install.packages("lfc")
library(lfc)
```

Alternatively, you can install from github:
```
install.packages("devtools")
devtools::install_github("erhard-lab/lfc")
library(lfc)
```


## Brief introduction

Digital expression measurements (e.g. RNA-seq) are often used to determine the change of quantities upon some treatment or stimulus. The resulting value of interest is the *fold change* (often logarithmized).

This effect size of the change is often treated as a value that can be computed as lfc(A,B) = log2 A/B. However, due to the probabilistic nature of the experiments, the effect size rather is a random variable that must be estimated. This fact becomes obvious when considering that A or B can be 0, even if the true abundance is non-zero.

We have shown that this can be modelled in a [Bayesian framework](https://dx.doi.org/10.1093/nar/gkv696). The intuitively computed effect size is the maximum likelihood estimator of a binomial model, where the effect size is not represented as fold change, but as proportion (of note, the log fold change simply is the logit transformed proportion). The Bayesian prior corresponds to pseudocounts frequently used to prevent infinite fold changes by A or B being zero. Furthermore, the Bayesian framework offers more advanced estimators (e.g. interval estimators or the posterior mean, which is the optimal estimator in terms of squared errors).

This R package offers the implementation to harness the power of this framework.

