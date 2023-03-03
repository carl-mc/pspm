
########################################
# Bootstrap helper functions
########################################

#' @name get_basic_ci
#' @title Get basic CI from bootstrapped estimates
#' @param bmat Bootstrapped coefficient matrix
#' @param bpoint Point estimate
#' @param alpha Confidence level
#' @keywords internal
get_basic_ci <- function(bmat, bpoint, alpha) {
  ci_mat <- matrix(NA, 2, ncol(bmat))
  for (k in 1:ncol(bmat)) {
    p <- bpoint[k]
    b <- bmat[,k]
    ci_mat[,k] <- c(2*p - quantile(b, 1 - alpha/2),
                    2*p - quantile(b, alpha/2))
  }
  return(ci_mat)
}


#' @name get_percentile_ci
#' @title Get percentile-based CI from bootstrapped estimates
#' @param bmat Bootstrapped coefficient matrix
#' @param bpoint Point estimate
#' @param alpha Confidence level
#' @keywords internal
get_percentile_ci <- function(bmat, alpha) {
  ci_mat <- matrix(NA, 2, ncol(bmat))
  for (k in 1:ncol(bmat)) {
    b <- bmat[,k]
    ci_mat[,k] <- c(quantile(b, alpha/2),
                    quantile(b, 1-alpha/2))
  }
  return(ci_mat)
}


########################################
# Fast cluster Export via disk
########################################

#' @name clusterExport_fast
#' @title Fast Export of R-Objects to CPU Node Cluster
#' @description Runs significantly faster than \code{parallel:clusterExport}, 
#' as it uses a temporary file that is written to disc and then read in parallel. 
#' This comes at the risk that the temporary file may survive if the function breaks
#' and that some types of objects, in particular when they include pointers, 
#' are broken upon arrival on the cluster. 
#' 
#' @param cl cluster object, created by \code{parallel::makeCluster}
#' @param varlist Character vector of object names
#' @param envir Environment from which to export objects. Defaults to \code{.GlobalEnv}.
#' 
#' @returns \code{TRUE}
#' @export
clusterExport_fast <- function(cl, varlist, envir = .GlobalEnv){
  # Save var.list to tempfile
  tmp.file <- tempfile()
  save(list = varlist, file = tmp.file, envir = envir)
  clusterExport(cl, list("tmp.file"), envir = environment())
  
  # Load data
  load <- parLapply(cl, seq_along(cl), function(i){
    load(tmp.file, envir = .GlobalEnv)
    TRUE
  })
  
  # Cleanup
  clusterEvalQ(cl = cl, expr = {
    rm(tmp.file)
    TRUE
  })
  unlink(tmp.file)
  
  # Return
  return(TRUE)
}

############################
# Reset functions
############################

#' @name reset_pspm
#' @title Reset PSPM object
#' @description Reset PSPM object to restore pointers and internal Python objects. This function
#' has to be called after reloading a PSPM object from disc or exporting it to a CPU cluster. 
#' 
#' @param x A PSPM object to be reset. 
#' 
#' @returns The reset PSPM object.
#' @export
reset_pspm <- function(x){
  PSPM$new(x$Y, x$Z, 
              x$X, x$A, x$net_stats,
              x$beta, x$points, 
              force_contiguous = x$force_contiguous, 
              verbose = x$verbose)
}

#' @name reset_pspmlearn
#' @title Reset PSPMLearn object
#' @description Reset PSPMLearn object to restore pointers and internal Python objects. This function
#' has to be called after reloading a PSPMLearn object from disc or exporting it to a CPU cluster. 
#' 
#' @param x A PSPMLearn object to be reset. 
#' 
#' @returns The reset PSPMLearn object.
#' @export
reset_pspmlearn <- function(x){
  x$reset()
  return(x)
}

######################
# Spatial Lattice Class
# For storing training/testing instances and simulation
######################
#' @name PSPM
#' @title Probabilistic Spatial Partition Model Data Container
#' 
#' @description Data container for a spatial network to be analysed with a Probabilistic Spatial Partition Model. 
#' Objects of this class are the main input to create a PSPMLearn object, 
#' which allows to fit the Probabilistic Partition Model to derive parameter estimates. 
#' 
#' @details Write about details here.
#' 
#' 
#' @import R6
#' reticulate
#' @export
PSPM <- R6Class("PSPM",
                   
                   public = list(
                     
                     ## PUBLIC DATA
                     #' @field A 
                     #' NxN adjacency matrix of baseline lattice
                     A = NULL, 
                     
                     #' @field g 
                     #' A python object of class SpatialLattice
                     g = NULL, 
                     
                     #' @field pm
                     #' A python object of class PartitionModel
                     pm = NULL, 
                     
                     #' @field N 
                     #' Number of nodes in spatial lattice
                     N = NULL, 
                     
                     #' @field L 
                     #' Number of edge-wise covariates
                     L = NULL, 
                     
                     #' @field Y 
                     #' Vector of partition labels
                     Y = NULL, 
                     
                     #' @field Z 
                     #' List of length L of NxN edge-wise covariates
                     Z = NULL, # 
                     
                     #' @field  Z.labs
                     #' Vector of edge-wise covariate labels
                     Z.labs = NULL, 
                     
                     #' @field  P
                     #' Number of vertex predictors. Currently 1, not used. Only here for future extensions.
                     P = NULL, 
                     
                     #' @field X 
                     #' NxP Matrix of vertex predictors, currently contains only one column full of 1s. 
                     #' Currently not used. Only here for future extensions.
                     X = NULL, 
                     
                     #' @field points 
                     #' A matrix of xy coordinates or a SpatialPoints* object of nrow/length = N; 
                     #' corresponds to nodes.
                     points = NULL, 
                     
                     #' @field beta 
                     #' Vector of beta parameters. Defaults to 0s
                     beta = NULL,
                     
                     
                     #' @field net_stats 
                     #' List of network statistics. 
                     net_stats = c("LatticeNeighborsNS"), 
                     
                     #' @field verbose 
                     #' Level of verbosity, binary flag. 
                     verbose = FALSE, 
                     
                     #' @field force_contiguous 
                     #' Flag on whether to enforce contiguous partitions. Input by user. 
                     force_contiguous = TRUE, 
                     
                     ## PUBLIC METHODS
                     
                     ### Initialize object
                     #' @description Initialize a new PSPM object
                     #' 
                     #' @param Y Integer vector of outcome partitioning. 
                     #' @param Z List of weighted adjacency matrices, each encoding one edge-level predictor. 
                     #' Can be named, names are then recycled as variable names.
                     #' @param X Matrix of Vertex attributes. Not in use at the moment, 
                     #' for future extensions and backwards compatibility. 
                     #' @param A Adjacency matrix
                     #' @param net_stats List of network statistics. 
                     #' Defaults to \code{c("LatticeNeighborsNS")} and included for 
                     #' for future extensions and backwards compatibility. 
                     #' @param beta Preset beta parameters, if any. 
                     #' @param points SpatialPoints object that give geolocation of vertices.
                     #' @param force_contiguous Flag to enforce contiguity of the graph 
                     #' and throw an error if that is not given. 
                     #' @param verbose Flag for verbosity. 
                     #' 
                     #' @examples
                     #' 
                     initialize = function(Y, Z, X = NULL, A, 
                                           net_stats = c("LatticeNeighborsNS"),
                                           beta = NULL, points = NULL, force_contiguous=TRUE, verbose = FALSE) {
                       
                       # Replace X if NULL with empty placeholder
                       if(is.null(X)){
                         X <- as.matrix(rep(1, length(Y)), ncol = 1)
                         colnames(X) <- "X"
                       }
                       
                       # Constructor
                       self$A = A
                       self$N = nrow(A)
                       self$L = length(Z)
                       self$Z = Z
                       self$Z.labs = names(Z)
                       names(Z) <- NULL
                       self$P = ncol(X)
                       self$X = X
                       self$Y = Y
                       self$beta = beta
                       self$points = points
                       self$verbose = verbose
                       self$force_contiguous = force_contiguous
                       self$net_stats = net_stats

                       # Initialize SpatialLattice
                       self$g <-  SpatialLattice(A, c(list(A), Z), 
                                                 X)
                       
                       # Initialize partition model
                       self$pm <- PartitionModel(self$g, self$Y, 
                                                 force_contiguous = self$force_contiguous,
                                                 net_stats = self$net_stats)
                       
                       # Initialize beta (also sets D)
                       if(!is.null(self$beta)){
                         self$set_beta(self$beta)
                       }

                     },
            

                      #' @description Set beta parameters on PSPM object. 
                      #' 
                      #' @param beta Numeric vector with a constant and 
                      #' one entry per edge-level predictor, in this order
                      #'
                      #'
                      #' @examples
                     set_beta = function(beta) {
                       # Sets beta and updates D (the predictor matrix)
                       stopifnot(length(beta) == length(unlist(self$pm$parameters)))
                       self$beta <- beta
                       
                       # Set name
                       if(!is.null(self$Z.labs)){
                         names(self$beta) <- c("Constant", self$Z.labs)
                       }
                       
                       # Update pm object
                       self$pm$update_parameters_theta(beta)
                     },


                     #' @description
                     #' Get beta parameters on PSPM object. 
                     #' 
                     #' @returns Parameter vector, a numeric vector with a constant and 
                     #' one entry per edge-level predictor, in this order
                     #'
                     #' @examples
                     get_beta = function() {
                       return(self$beta)
                     },


                     #' @description Plot partitioning of PSPM object as a labelled graph. 
                     #'
                     #' @param edge_predictor ID or name of edge predictor to include
                     #'  to plot different edge types. Only binary edge predictors supported at the moment.
                     #'  
                     #' @examples
                     plot_partitioning = function(edge_predictor = NULL) {
                       # Plots the current lattice partitioning
                       
                       plot_adj_mat <- self$A
                       G <- (outer(self$Y, self$Y, function(a, b) {a==b})*1)*plot_adj_mat
                       plot_adj_mat[G==1] <- 2
                       plot_g <- igraph::graph_from_adjacency_matrix(plot_adj_mat, 
                                                                     mode = 'undirected', 
                                                                     weighted = TRUE)
                       
                       labels <- self$Y
                       n_partitions <- length(unique(labels))
                       cols <- rainbow(n_partitions)
                       
                       if (!is.null(edge_predictor)) {
                         if(is.character(edge_predictor)){
                           edge_predictor <- which(self$Z.labs == edge_predictor)
                         } 
                         pred_mat <- self$Z[[edge_predictor]]
                         elist <- as_edgelist(plot_g)
                         predictor <- rep(NA, nrow(elist))
                         for (i in 1:nrow(elist)) {
                           predictor[i] <- pred_mat[elist[i,1], elist[i,2]]
                         }
                         
                         if(any(!unique(predictor) %in% c(0,1))){
                           warning("Plotting works only with binary edge predictor")
                         } else {
                           edge_attr(plot_g, 'weight') <- predictor
                         }
                         
                       }
                       
                       # Encode edge types if any
                       edge_types <- E(plot_g)$weight
                       
                       # Plot
                       plot(plot_g,
                            vertex.color = cols[as.numeric(as.factor(self$Y))],
                            vertex.label = labels,
                            layout = layout_on_grid(plot_g),
                            edge.lty = ifelse(edge_types == 1, 2, 1),
                            edge.width = ifelse(edge_types == 1, 0.5, 2))
                     },
                     
                     
                      #' @description Sample partitioning 
                      #'
                      #' @param burnin Length of burn-in period
                      #' @param return_full Flag to return all partitionings, one for each burn-in iteration.
                      #'
                      #' @return If \code{return_full = TRUE}, a matrix of partition IDs. 
                      #' If \code{return_full = FALSE}, the partitioning of the PSPM object
                      #' is updated internally. 
                      #' 
                      #'
                      #' @examples
                       sample = function(burnin = 0, return_full = FALSE) {
                       # Samples a new lattice partitioning via Gibbs sampling.
                       
                       # Set up result set if requested
                       if(return_full){
                         sample.mat <- matrix(NA, self$N, burnin+1)
                       }
                       
                       # Sample
                       for (m in 1:(burnin+1)) {
                         self$pm$sample_all(M = 1L)
                         if(return_full){
                           sample.mat[ ,m] <- as.integer(unlist(self$pm$partitioning))
                         }
                       }
                       
                       # Update visible labels
                       self$Y = as.integer(unlist(self$pm$partitioning))
                       
                       # Return if requested
                       if(return_full){
                         return(sample.mat)
                       }
                     },

                     #' @description Compute the log composite likelihood
                     #'
                     #' @param beta Beta parameters for which to query log composite likelihood
                     #'
                     #' @return Returns log composite likelihood
                     #' 
                     #' @examples
                     get_composite_log_likelihood = function(beta) {
                       # Returns log pseudo likelihood (see Sutton and McCullum 2011: p. 345)
                       
                       # Set beta (also sets D)
                       self$pm$update_parameters_theta(beta)
                       
                       # Get the log likelihood
                       ll <- self$pm$get_composite_log_likelihood()
                       
                       return(ll)
                     }
                   )
)

######################
# PSPMLearn Class
# For parameter learning on collections of PSPM objects
######################
#' @name PSPMLearn
#' @title Probabilistic Spatial Partition Model Fitter
#' 
#' @description Class to fit the parameters of a Probabilistic Spatial Partition Model across
#' a set of PSPM objects. 
#' 
#' @details 
#' 
#' @return Returns a PSPMLearn object.
#' 
#' @details 
#' 
#' @import R6
#' reticulate
#' maxLik
#' foreach
#' doParallel
#' parallel
#' foreach
#' @export
PSPMLearn <- R6Class("PSPMLearn",
                        public = list(
                          
                          # Class variables
                          
                          #' @field pspm_ls 
                          #' List of PSPM objects. Input by user. 
                          pspm_ls = list(),
                          #' @field M 
                          #' Number of objects in \code{pspm_ls}. Set internally.
                          M = NULL,
                          #' @field  sigma2 
                          #' Penalization parameter. Input by user.
                          sigma2 = NULL,
                          # current_instance = NULL,
                          
                          #' @field samples_cache 
                          #' Cache of samples saved during parametric bootstrapping. 
                          samples_cache = list(),
                          
                          #' @description Initialize a new PSPMLearn object
                          #'
                          #' @param pspm_ls List of PSPM objects
                          #' @param sigma2 Penalty parameter for penalized maximum composite likelihood estimation
                          #'
                          #' @return A PSPMLearn object
                          #' @export
                          #'
                          #' @examples
                          initialize = function(pspm_ls, sigma2 = 10) {
                            
                            # Set the list of instances
                            self$pspm_ls = pspm_ls
                            self$M <- length(pspm_ls)
                            
                            # Set the penalty parameter
                            self$sigma2 = sigma2
                            
                            # # Set the 'current_instance' iterator
                            # self$current_instance = 1
                            
                          },
                          
                          
                          #' @description Reset PSPMLearn object
                          #'
                          #' 
                          #' @export
                          #'
                          #' @examples
                          reset = function(){
                            
                            # Re-Initialize partition graphs in pspm.ls
                            for(s in seq_along(self$pspm_ls)){
                              self$pspm_ls[[s]] <- reset_pspm(self$pspm_ls[[s]])
                            }
                            
                          }, 
                          
                          #' @description Reset beta parameters to 0
                          #'
                          #' @param beta Can be NULL -- then chooses (default) betas set in 
                          #' first PSPM object -- 
                          #' or a vector of beta parameters (constant and one for each edge-predictor)
                          #' 
                          #' @export
                          #'
                          #' @examples
                          reset_beta = function(beta = NULL) {
                            ## Resets the beta parameter vector of all instances
                            # beta: Optional argument; beta set to zero vector if missing
                            
                            for (m in 1:self$M) {
                              if (is.null(beta)) {
                                # L <- self$pspm_ls[[m]]$L
                                # P <- self$pspm_ls[[m]]$P
                                # NS <- length(self$pspm_ls[[m]]$net_stats)
                                beta <- rep(0, length(unlist(self$pspm_ls[[m]]$pm$parameters)))
                              }
                              self$pspm_ls[[m]]$set_beta(beta)
                            }
                          },
                          
                          # @description Reset beta parameters to 0
                          #
                          #
                          # @export
                          #
                          # @examples
                          # get_current_instance = function() {
                          #   ## Current instance iterator
                          # 
                          #   rv <- self$current_instance
                          #   if (self$current_instance == self$M) {
                          #     self$current_instance <- 1
                          #   } else {
                          #     self$current_instance <- rv+1
                          #   }
                          # 
                          #   return(rv)
                          # },


                          #' @description  Get composite log likelihood 
                          #' 
                          #' Returns sum of composite log likelihoods over
                          #' all PSPM objects in the PSPMLearn object for a set of beta parameters. 
                          #'
                          #' @param beta Set of Beta parameters to compute composite log likelihood for. 
                          #'
                          #' @return Numeric composite log likelihood penalized by \code{sum(beta^2)/(2*sigma2)}
                          #' @export
                          #'
                          #' @examples
                          composite_log_likelihood = function(beta) {
                            ## Returns composite likelihood
                            
                            # For each instance: Get log-likelihood
                            ll = 0
                            for (m in 1:self$M) {
                              ll <- ll + self$pspm_ls[[m]]$get_composite_log_likelihood(beta = beta)
                            }
                            
                            # Penalize
                            if(!is.null(self$sigma2)){
                              ll <- ll - sum(beta^2)/(2*self$sigma2)
                            }

                            # Return
                            return(ll)
                          },
                          
                          

                          #' @description Fit beta parameters to minimize composite log likelihood
                          #' 
                          #' Returns set of beta parameters the minimizes (penalized) composite log likelihood
                          #' across PSPM objects in PSPMLearn object. 
                          #'
                          #' @param beta_init Set of Beta parameters to initialize fitting with
                          #'
                          #' @param set_cache (Re)set cache to speed up fitting
                          #'
                          #' @param clear_cache Clear cache after fitting (reduces memory demands). 
                          #'
                          #' @return Numeric composite log likelihood penalized by \code{sum(beta^2)/(2*sigma2)}
                          #' 
                          #' @import maxLik
                          #' 
                          #' @export
                          #'
                          #' @examples
                          fit_composite_log_likelihood = function(beta_init = NULL, 
                                                                    set_cache = TRUE, 
                                                                    clear_cache = TRUE) {
                            ## Estimation via maximum composite likelihood
                            # beta_init: Optional argument; beta_init is copied from first instance if missing
                            
                            # Set up cached membership options
                            if(set_cache){
                              for(m in 1:self$M){
                                self$pspm_ls[[m]]$pm$generate_likelihood_cache()
                              }
                            }
                            
                            # Set initial beta
                            if (is.null(beta_init)) {
                              beta_init <- self$pspm_ls[[1]]$get_beta()
                            }
                            
                            # Fit model
                            mle.fit <- maxLik(self$composite_log_likelihood, start = beta_init, method = 'BFGS')
                            
                            # Set beta vectors of instances to estimated beta
                            beta_est <- coef(mle.fit)
                            self$reset_beta(beta_est)
                            
                            # Clear cached membership options
                            if(clear_cache){
                              for(m in 1:self$M){
                                self$pspm_ls[[m]]$pm$likelihood_cache = dict()
                              }
                            }
                            
                            # Return
                            return(mle.fit)
                          },
                          
                          
                          

                          #' @description Parametric bootstrap
                          #' 
                          #' Returns confidence intervals and full set of 
                          #' bootstrapped estimates from parametric bootstrap. 
                          #'
                          #' @param n_boot_iter Number of bootstrap iterations. 
                          #' @param burnin Number of burn-in period for each sampling.
                          #' @param ci_level Confidence interval level.
                          #' @param return_sims Return all sampled partitionings?
                          #' @param report_every Periodicity of reporting after bootstrap iterations.
                          #' @param cache_samples Cache the sampled partitionings in the PSPM objects inside the PSPMLearn object. 
                          #'
                          #' @return A list with estimated parameters, confidence intervals, and 
                          #' possibly sampled partitionings. 
                          #' 
                          #' @export
                          #'
                          #' @examples
                          bootstrap_composite_log_likelihood = function(n_boot_iter = 100, 
                                                                        burnin = 1e1,
                                                                     ci_level = 0.95,
                                                                     return_sims = FALSE, 
                                                                     report_every = 50,
                                                                     cache_samples = T) {
                            ## Estimates CIs for beta using a parametric bootstrap
                            
                            # Clear previous cache
                            if(cache_samples){
                              self$samples_cache <- list()
                            }
                            
                            # Get the point estimates from the first pspm instance
                            beta_point <- self$pspm_ls[[1]]$get_beta()
                            
                            # Bootstrap iterations
                            beta_boot_ls <- vector('list', n_boot_iter)
                            for (j in 1:n_boot_iter) {
                              # Init caching of samples
                              if(cache_samples){
                                self$samples_cache[[length(self$samples_cache) + 1]] <- list()
                              }
                              
                              # Generate new instances
                              boot_pspm_ls <- vector('list', self$M)
                              for (m in 1:self$M) {
                                A <- self$pspm_ls[[m]]$A
                                Z <- self$pspm_ls[[m]]$Z
                                X <- self$pspm_ls[[m]]$X
                                Y <- 1:self$pspm_ls[[m]]$N
                                net_stats <- self$pspm_ls[[m]]$net_stats
                                fc <- 1:self$pspm_ls[[m]]$force_contiguous
                                beta <- self$pspm_ls[[m]]$get_beta()
                                
                                pspm_m <- PSPM$new(Y = Y, Z = Z, X =  X, 
                                                         A = A, 
                                                         net_stats = net_stats,
                                                         force_contiguous=fc, 
                                                         verbose = FALSE)
                                pspm_m$set_beta(beta)
                                pspm_m$sample(burnin = burnin)
                                boot_pspm_ls[[m]] <- pspm_m
                                
                                # Cache this sample
                                if(cache_samples){
                                  self$samples_cache[[length(self$samples_cache)]][[m]] <- pspm_m$Y
                                }
                              }
                              
                              # Generate new learning obj, fit
                              learn <- PSPMLearn$new(pspm_ls = boot_pspm_ls, sigma2 = self$sigma2)
                              learn$reset_beta()
                              llfit <- learn$fit_composite_log_likelihood()
                              beta_boot_j <- coef(llfit)
                              beta_boot_ls[[j]] <- beta_boot_j
                              
                              # Report progress
                              if (j %% report_every == 0) {
                                print(paste('Boostrap Iteration', j))
                                flush.console()
                              }
                            }
                            
                            # Compute CIs
                            beta_boot_mat <- do.call('rbind', beta_boot_ls)
                            alpha <- 1 - ci_level
                            ci_basic_mat <- get_basic_ci(beta_boot_mat, beta_point, alpha)
                            ci_perc_mat <- get_percentile_ci(beta_boot_mat, alpha)
                            
                            # Labels
                            if(!is.null(self$pspm_ls[[1]]$Z.labs)){
                              colnames(beta_boot_mat) <- c("Constant", self$pspm_ls[[1]]$Z.labs)
                              colnames(ci_basic_mat) <- c("Constant", self$pspm_ls[[1]]$Z.labs)
                              colnames(ci_perc_mat) <- c("Constant", self$pspm_ls[[1]]$Z.labs)
                              rownames(ci_basic_mat) <- c("LB_Basic", "UB_Basic")
                              rownames(ci_perc_mat) <- c("LB_Perc", "UB_Perc")
                            }
                     
                            # Prep return value
                            if (return_sims) {
                              out <- list(beta_boot = beta_boot_mat,
                                          ci_mat = rbind(ci_basic_mat,ci_perc_mat))
                            } else {
                              out <- list(ci_basic = ci_basic_mat,
                                          ci_mat = rbind(ci_basic_mat,ci_perc_mat))
                            }
                            
                            return(out)
                          },
                          

                          #' @description Parallelized parametric bootstrap
                          #' 
                          #' Returns confidence intervals and full set of 
                          #' bootstrapped estimates from parametric bootstrap. 
                          #'
                          #' @param cl CPU cluster object (created through \code{parallel::makeCluster}) 
                          #' or integer number of core. If the latter, a cluster is setup internally and
                          #' closed upon completion.                           #' 
                          #' @param n_boot_iter Number of bootstrap iterations. 
                          #' @param burnin Number of burn-in period for each sampling.
                          #' @param ci_level Confidence interval level.
                          #' @param return_sims Return all sampled partitionings?
                          #' @param report_every Periodicity of reporting after bootstrap iterations.
                          #' @param cache_samples Cache the sampled partitionings in the PSPM objects inside the PSPMLearn object. 
                          #'
                          #' @return A list with estimated parameters, confidence intervals, and 
                          #' possibly sampled partitionings. 
                          #' 
                          #' @export
                          #'
                          #' @examples
                          par_bootstrap_composite_log_likelihood = function(cl = NULL, n_boot_iter = 100, burnin = 1e2,
                                                                         ci_level = 0.95, 
                                                                         return_sims = FALSE, cache_samples = T){
                            
                            if(is.null(cl)){
                              out <- self$bootstrap_composite_log_likelihood(n_boot_iter = n_boot_iter, burnin = burnin,
                                                                          ci_level = ci_level,
                                                                          return_sims = return_sims, report_every = 100,
                                                                          cache_samples = cache_samples)
                            } else {
                              # Setup cluster
                              if(is.numeric(cl)){
                                cl <- parallel::makeCluster(getOption("cl.cores", cl))
                                close.cluster = T
                              } else {
                                close.cluster = F
                              }
                              
                              ## Init Cluster
                              init.ls <- parallel::clusterEvalQ(cl = cl, expr = {
                                # Packages 
                                library(pspm)
                                library(reticulate)
                              })
                              
                              # Register for doParallel 
                              doParallel::registerDoParallel(cl)
                              
                              # Load self to cluster
                              print("Load Learn Object on cluster")
                              learn_obj <- self
                              clusterExport_fast(cl = cl, varlist = c("learn_obj"),
                                                 envir = environment())
                              
                              # Reset learn object
                              load_learnobj <- parallel::parLapply(cl, seq_along(cl), 
                                                        function(i){
                                                         learn_obj <<- reset_pspmlearn(learn_obj)
                                                         TRUE
                                                       })
                              
                              # Run parallel bootstrap
                              print("Run Bootstrap")
                              bootse <- foreach::foreach(i = 1:n_boot_iter, .noexport = c("learn_obj")) %dopar% {
                                set.seed(i)
                                reticulate::py_set_seed(i, disable_hash_randomization = TRUE)
                                bse <- learn_obj$bootstrap_composite_log_likelihood(n_boot_iter = 1, burnin = burnin,
                                                                                 ci_level = ci_level,
                                                                                 return_sims = T, 
                                                                                 report_every = 5, 
                                                                                 cache_samples = T)
                                return(list(bse = bse, sample = learn_obj$samples_cache))
                              }
                              
                              # Close cluster
                              if(close.cluster){
                                stopCluster(cl)
                              }
                              
                              # Save samples to cache
                              if(cache_samples){
                                self$samples_cache <- lapply(bootse, function(bse){bse$sample[[1]]})
                              }
                              
                              # Compute CIs
                              beta_boot_mat <- do.call('rbind', lapply(bootse, function(bse){bse$bse$beta_boot}))
                              alpha <- 1 - ci_level
                              ci_mat <- matrix(NA, 4, ncol(beta_boot_mat))
                              
                              ## Basic CI
                              for (k in 1:ncol(beta_boot_mat)) {
                                if(inherits(learn_obj, "PSPMLearn")){
                                  p <- learn_obj$pspm_ls[[1]]$get_beta()[k]
                                } 
                                
                                b <- beta_boot_mat[,k]
                                ci_mat[1:2,k] <- c(2*p - quantile(b, 1 - alpha/2),
                                                   2*p - quantile(b, alpha/2))
                              }
                              
                              ## Percentile based
                              for (k in 1:ncol(beta_boot_mat)) {
                                b <- beta_boot_mat[,k]
                                ci_mat[3:4,k] <- c(quantile(b, alpha/2),
                                                   quantile(b, 1-alpha/2))
                              }
                              
                              # Labels
                              if(!is.null(self$pspm_ls[[1]]$Z.labs)){
                                colnames(beta_boot_mat) <- c("Constant", self$pspm_ls[[1]]$Z.labs)
                                colnames(ci_mat) <- c("Constant", self$pspm_ls[[1]]$Z.labs)
                                rownames(ci_mat) <- c("LB_Basic", "UB_Basic",
                                                      "LB_Perc", "UB_Perc")
                                
                              }
                              
                              # Prep return value
                              if (return_sims) {
                                out <- list(beta_boot = beta_boot_mat, ci_mat = ci_mat)
                              } else {
                                out <- ci_mat
                              }
                            }
                            
                            # Return
                            return(out)
                          }
                        )
)
