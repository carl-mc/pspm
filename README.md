# pspm: Probabilistic Spatial Partition Model

## Introduction

This packages contains a series of functions and classes to set up and
fit a Probabilistic Spatial Partition Model and sample partitionings
following the the methodology introduced by Müller-Crepon, Schvitz and
Cederman (2023). In short, the model is based on a simplification of
continuous geographic space as a (planar) graph. The observed
partitioning is encoded on the graphs’ vertices such that each vertex is
a member of one and only one partition. The graphs partitioning is
modeled as the result of *attractive* and *repulsive* force active on
its edges which affect the probability that the vertices they connect
belong to the same or to different partitions. The model allows to
estimate the effects of edge-elvel attributes on the
repulsion/attraction between vertices, thus estimating their impact on
the overall partitioning of the graph.

Please refer to the original publication for all technical details
beyond the summary provided below. It can be found [here]() and as
ungated version
[here](http://www.carlmueller-crepon.org/publication/state_shape/).

When using the pspm package, please cite:

Müller-Crepon, Carl, Guy Schvitz, Lars-Erik Cederman (2023). Shaping
States into Nations: The Effects of Ethnic Geography on State Borders.
*American Journal of Political Science*, conditionally accepted for
publication.

## The model

We model the distribution over all possible partitionings $P$ of lattice
$G$ as a Boltzmann distribution:

``` math
Pr(P = p_{i}) =  {e^{-\epsilon_{i}}\over\displaystyle\sum_{i = 1}^{|\mathbb{P}|}e^{-\epsilon_{i}} }
```

here, a partitioning $i$’s chance of realization decreases with its
energy $\epsilon_i$. The energy of a partitioning is defined as the sum
of energies on the graph’s edges that run between nodes $j$ and $k$
where $j$ and $k$ are located in the same partition ($s_{j,k} = 1$).
Edges energies are not realized where $j$ and $k$ do not belong to the
same partition ($s_{j,k} = 0$):

``` math
\epsilon_{i} = \displaystyle\sum_{j,k \in L}\, \epsilon_{j,k}\,s_{j,k}
```

Potential and realized edge-level energies are in turn affected by a
constant ($\beta_0$) and a set of edge-level attributes
$\textbf{x}_{j,k}$ that are weighted by a parameter vector $\beta$

``` math
\epsilon_{j,k} = \beta_0 + \beta\, \textbf{x}_{j,k}
```

The empirical goal of the PSPM is it to estimate parameters $\beta_0$
and $\beta$ from observed data. To that intent, the PSPM uses a maximum
composite likelihood approach. Uncertainty estimates can be derived via
a parametric bootstrap that samples partitionings with a given set of
$\beta_0$ and $\beta$ parameters and then fits the model on the sampled
partitionings.

## Installation

You can directly download and install the pspm package from GitHub.
Before doing so, please make sure that you have
[Python3](https://www.python.org/downloads/) installed. Upon
installation, the package automatically installs necessary python
dependencies via the
[reticulate](https://cran.r-project.org/web/packages/reticulate/index.html)
R-package.

``` r
library(devtools)
install_github(repo = "carl-mc/pspm")
```

## Getting started

``` r
# Packages
library(pspm)
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
# Set seeds in R and python
pspm_set_seed(1)
```

## Setting up and handling a toy lattice

``` r
# Make mock PSPM Object with a sampled partitioning
sl <- generate_grid_data(N_sqrd = 10, ## 10 x 10 lattice
                        beta0 = -2, ## Negative constant = baseline attraction between nodes
                        beta = c(2,1), ## Include two repulsive edge-level predictors
                        dep_structure = "von_neumann", ## Each node connects to 4 neighbors
                        burnin = 10 ## Sample with a 10 burn-in periods
                        )
sl$plot_partitioning(edge_predictor = 1, 
                     edge.width = 5, vertex.size = 10,
                     main = "A Toy Lattice")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## PSPM to igraph conversion

``` r
# Transform PSPM to igraph
graph <- PSPM2igraph(sl)
edge_attr_names(graph)
```

    ## [1] "x1" "x2"

``` r
vertex_attr_names(graph)
```

    ## [1] "X1" "Y"

``` r
# Transform igraph back to PSPM
sl.from.g <- igraph2PSPM(g = graph, outcome_name = "Y",
                         edge_pred_names = c("x1", "x2"))
```

## Fitting a PSPM Model

``` r
## The easy way

### Estimate
m.simple <- fit_pspm_model(formula = Y ~ x1 + x2, 
                           g_ls = list(graph),
                           return_pspm = TRUE)

# ### Bootstrap CIs
# bs.simple <- bootstrap_pspm(m.simple, 
#                             n_boot_iter = 10, ## Should be > 100
#                             burnin = 10, ## Could be higher, depending on complexite of graph and model
#                             cl = 10L, ## Number of CPUs for parallelization
#                             return_sims = TRUE, ## Return full distribution of estimates
#                             ci_level = .95)
# 
# ### Summary
# summary(m.simple)
# print(bs.simple$ci_mat)
# plot(density(bs.simple$beta_boot[,"x1"]),
#      main = "Distribution of bootstrapped estimates of x1")
```

``` r
## The complicated way (inside the wrapper)

### Initiate PSPMLearn Object
learn_obj <- PSPMLearn$new(list(sl.from.g))

# ### Fit
# m.compl <- learn_obj$fit_composite_log_likelihood(beta_init = c(0,0,0))
# 
# ### Bootstrap
# bs.compl <- learn_obj$par_bootstrap_composite_log_likelihood(n_boot_iter = 10, burnin = 10, 
#                                                  cl = 10L, return_sims = FALSE, ci_level = .95)
# 
# ### Summary
# summary(m.compl)
# print(bs.compl)
```

## Producing nice tables

``` r
# Integration with texreg
pspm2table(list(m.simple),
           # bootci = list(bs.simple$ci_mat), boottype = "percentile",
           type = "text", add.stats = c("Edges" = "N_edges", "Vertices" = "N"))
```

    ## 
    ## ==========================
    ##                 Model 1   
    ## --------------------------
    ## Constant         -2.51 ***
    ##                  (0.59)   
    ## x1                1.95    
    ##                  (1.84)   
    ## x2                0.45    
    ##                  (0.93)   
    ## --------------------------
    ## Edges           180       
    ## Vertices        100       
    ## Log-Likelihood  -13.45    
    ## Num. obs.       100       
    ## ==========================
    ## *** p < 0.01; ** p < 0.05; * p < 0.1
