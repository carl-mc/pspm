When using the pspm package, please cite: Müller-Crepon, Carl, Guy
Schvitz, Lars-Erik Cederman (2023). Shaping States into Nations: The
Effects of Ethnic Geography on State Borders. *American Journal of
Political Science*, conditionally accepted for publication.

## Installation

You can directly download and install the pspm package from GitHub.
Before doing so, please make sure that you have
[Python3](https://www.python.org/downloads/) installed. Upon
installation, the package automatically installs necessary python
dependencies via the
[reticulate](https://cran.r-project.org/web/packages/reticulate/index.html)
R-package.

    library(devtools)
    install_github(repo = "carl-mc/pspm")

## Getting started

    # Packages
    library(pspm)
    library(igraph)

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    # Set seeds in R and python
    pspm_set_seed(1)

## Handling a toy lattice

    # Make mock PSPM Object with a d partitioning
    sl <- generate_grid_data(N_sqrd = 10,
                            beta0 = -2,
                            beta = 2,
                            dep_structure = "von_neumann",
                            burnin = 10)
    sl$plot_partitioning(edge_predictor = 1, edge.width = 5, vertex.size = 10,
                         main = "A Toy Lattice")

![](README_files/figure-markdown_strict/unnamed-chunk-3-1.png)

    # Transform PSPM to igraph
    graph <- PSPM2igraph(sl)
    edge_attr_names(graph)

    ## [1] "x1"

    vertex_attr_names(graph)

    ## [1] "X1" "Y"

    # Transform igraph back to PSPM
    sl.from.g <- igraph2PSPM(g = graph, outcome_name = "Y",
                             edge_pred_names = "x1")

## Fitting a PSPM Model

    ## The easy way

    ### Estimate
    m.simple <- fit_pspm_model(formula = Y ~ x1, 
                               g_ls = list(graph),
                               return_pspm = TRUE)

    ### Bootstrap CIs
    bs.simple <- bootstrap_pspm(m.simple, n_boot_iter = 10, burnin = 10, 
                                cl = 10L, return_sims = TRUE, ci_level = .95)

    ## [1] "Load Learn Object on cluster"
    ## [1] "Run Bootstrap"

    ### Summary
    summary(m.simple)

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 21 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -26.18528 
    ## 2  free parameters
    ## Estimates:
    ##          Estimate Std. error t value  Pr(> t)    
    ## Constant  -2.1262     0.4044  -5.258 1.46e-07 ***
    ## x1         1.0165     1.1015   0.923    0.356    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

    print(bs.simple$ci_mat)

    ##            Constant         x1
    ## LB_Basic -2.4998667 0.05452731
    ## UB_Basic -0.6604442 1.89724697
    ## LB_Perc  -3.5920536 0.13574815
    ## UB_Perc  -1.7526311 1.97846781

    plot(density(bs.simple$beta_boot[,"x1"]),
         main = "Distribution of bootstrapped estimates of x1")

![](README_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    ## The complicated way (inside the wrapper)

    ### Initiate PSPMLearn Object
    learn_obj <- PSPMLearn$new(list(sl.from.g))

    ### Fit
    m.compl <- learn_obj$fit_composite_log_likelihood(beta_init = c(0,0))

    ### Bootstrap
    bs.compl <- learn_obj$par_bootstrap_composite_log_likelihood(n_boot_iter = 10, burnin = 10, 
                                                     cl = 10L, return_sims = FALSE, ci_level = .95)

    ## [1] "Load Learn Object on cluster"
    ## [1] "Run Bootstrap"

    ### Summary
    summary(m.compl)

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 21 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -26.18528 
    ## 2  free parameters
    ## Estimates:
    ##      Estimate Std. error t value  Pr(> t)    
    ## [1,]  -2.1262     0.4044  -5.258 1.46e-07 ***
    ## [2,]   1.0165     1.1015   0.923    0.356    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## --------------------------------------------

    print(bs.compl)

    ##            Constant         x1
    ## LB_Basic -2.4998667 0.05452731
    ## UB_Basic -0.6604442 1.89724697
    ## LB_Perc  -3.5920536 0.13574815
    ## UB_Perc  -1.7526311 1.97846781

## Producing nice tables

    # Integration with texreg
    pspm2table(list(m.simple),
               bootci = list(bs.simple$ci_mat), boottype = "percentile",
               type = "text", add.stats = c("Edges" = "N_edges", "Vertices" = "N"))

    ## 
    ## ==============================
    ##                 Model 1       
    ## ------------------------------
    ## Constant         -2.13 *      
    ##                 [-3.59; -1.75]
    ## x1                1.02 *      
    ##                 [ 0.14;  1.98]
    ## ------------------------------
    ## Edges           180           
    ## Vertices        100           
    ## Log-Likelihood  -26.19        
    ## Num. obs.       100           
    ## ==============================
    ## * 0 outside the confidence interval.
