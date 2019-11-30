# rstanmsm
rstanmsm: Markov-Switching Models using Stan

## Introduction

This package implements Markov Switching Models in Stan. A Markov-Switching Model (MSM) considers the data generating process (DGP) as composed of a discrete model taking the form of a Markov process and a continuous model whose parameters and distribution are partially or completely conditional on the realization of the Markov process. The MSM applies to cases in which the discrete process is hidden and only its realizations in the form of the observed outcome are available. Given a sufficient correspondence to the true DGP, the model can recover the distributions of the state-specific parameters as well as the state sequence for K states.


## Emphasis

This package focuses primarily on longitudinal data and thereby also covers the univariate time-series case. This leads to a variety of options when specifying the model. Specifically, we can consider the state-process and/or the transition probabilities between the states as being shared across time-series.


## Current state

The current version supports balanced panels for data that is Gaussian. Ordering constraints on the parameters are arbitrary and not restricted to the intercept term to provide greater freedom. Transition probabilities, while still time-constant, can vary or be shared. The same applies to the discrete process.


## Future developments

* Extension to unbalanced panels
* Implementation for missing values within the time series
* User-defined priors
* Posterior Predictive Checks, e.g. integration with ```bayesplot```
* Summary statistics
