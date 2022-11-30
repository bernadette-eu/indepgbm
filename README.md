
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bayesian analysis of diffusion-driven multi-type epidemic models with application to COVID-19

Lampros Bouranis^(1,\*), Nikolaos Demiris^1, Konstantinos
Kalogeropoulos^2, Ioannis Ntzoufras^1

^1 Department of Statistics, Athens University of Economics and
Business, Athens, Greece

^2 The London School of Economics and Political Science, London, United
kingdom

^\* Corresponding author (<bouranis@aueb.gr>)

arXiv preprint: [link](https://arxiv.org/abs/2211.15229)

## Summary

We consider a flexible Bayesian evidence synthesis approach to model the
age-specific transmission dynamics of COVID-19 based on daily
age-stratified mortality counts. The temporal evolution of transmission
rates in populations containing multiple types of individual are
reconstructed via an appropriate dimension-reduction formulation driven
by independent diffusion processes assigned to the key epidemiological
parameters. A suitably tailored Susceptible-Exposed-Infected-Removed
(SEIR) compartmental model is used to capture the latent counts of
infections and to account for fluctuations in transmission influenced by
phenomena like public health interventions and changes in human
behaviour. We analyze the outbreak of COVID-19 in Greece and Austria and
validate the proposed model using the estimated counts of cumulative
infections from a large-scale seroprevalence survey in England.

## Installation

You can install the development version of the **Bernadette** R library
from GitHub following these [installation
instructions](https://github.com/bernadette-eu/Bernadette/).

## Libraries and source files

The workflow is initialised by loading the required libraries and
sourcing .R files via
[/R/1\_Libraries.R](https://github.com/bernadette-eu/indepgbm/blob/main/R/1_Libraries.R).

## Case Study: COVID-19 in Greece

-   Single Brownian motion (SBM) model: see
    [/Sampling/SingleBM\_GR.R](https://github.com/bernadette-eu/indepgbm/blob/main/Sampling/SingleBM_GR.R).

-   Multi Brownian motion (MBM) model: see
    [/Sampling/MultiBM\_GR.R](https://github.com/bernadette-eu/indepgbm/blob/main/Sampling/MultiBM_GR.R).

## Case Study: COVID-19 in Austria

-   Single Brownian motion (SBM) model: see
    [/Sampling/SingleBM\_AT.R](https://github.com/bernadette-eu/indepgbm/blob/main/Sampling/SingleBM_AT.R).

-   Multi Brownian motion (MBM) model: see
    [/Sampling/MultiBM\_AT.R](https://github.com/bernadette-eu/indepgbm/blob/main/Sampling/MultiBM_AT.R).

## Model validation: COVID-19 in England

See
[/Sampling/MultiBM\_ENG.R](https://github.com/bernadette-eu/indepgbm/blob/main/Sampling/MultiBM_ENG.R).

## Estimation of model information criteria

-   Estimation of the Deviance information criterion (DIC) and the
    respective effective number of model parameters via
    [/R/6\_MCMC\_diagnostics\_SingleBM.R](https://github.com/bernadette-eu/indepgbm/blob/main/R/6_MCMC_diagnostics_SingleBM.r)
    for the SBM model and via
    [/R/6\_MCMC\_diagnostics\_MultiBM.R](https://github.com/bernadette-eu/indepgbm/blob/main/R/6_MCMC_diagnostics_SingleBM.r)
    for the MBM model.

-   Estimation of the Pareto smoothed importance sampling Leave-One-Out
    information criterion and the respective effective number of model
    parameters via
    [/R/5\_LooIC.R](https://github.com/bernadette-eu/indepgbm/blob/main/R/5_LooIC.R).

## Graphical outputs

Information about the non-pharmaceutical interventions implemented by
European governments during the study period is available at
[/Data/](https://github.com/bernadette-eu/indepgbm/tree/main/Data).

The graphs available in the main text and the supplementary material can
be reproduced by executing the code in
[/R/7\_Model\_Fit\_graphs\_generation.R](https://github.com/bernadette-eu/indepgbm/blob/main/R/7_Model_Fit_graphs_generation.R).
