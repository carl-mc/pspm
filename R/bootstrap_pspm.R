
#' @description Parametric bootstrap
#' 
#' Returns confidence intervals and full set of 
#' bootstrapped estimates from a parametric bootstrap of a fitted PSPM model. 
#'
#' @param model A model as returned by fit_pspm_model(..., return_pspm = TRUE). 
#' Must include a PSPMLearn object called "learn_obj". Can be NULL if PSPMLearn object is provided directly.
#' @param PSPMLearn  A PSPMLearn object.  Can be NULL if model is provided instead.
#' @param cl CPU cluster object (created through \code{parallel::makeCluster}) 
#' or integer number of core. If the latter, a cluster is setup internally and
#' closed upon completion.    
#' @param n_boot_iter Number of bootstrap iterations. 
#' @param burnin Number of burn-in period for each sampling.
#' @param ci_level Confidence interval level.
#' @param return_sims Return all sampled partitionings?
#' @param cache_samples Cache the sampled partitionings in the PSPM objects inside the PSPMLearn object. 
#' @param beta Defaults to NULL. Can be a vector of beta parameters to reset 
#' beta parameters in the model and associated PSPMLearn object. If fitted, the latter
#' beta parameters in the PSPMLearn objects are the fitted parameters. 
#' The function resets beta parameters to its old values when finished. 
#' To avoid inconsistencies inside PSPMLearn objects, 
#' samples cannot be cached if beta parameter is specified. 
#'
#' @return A list with estimated parameters, confidence intervals, and 
#' possibly sampled partitionings. 
#' 
#' @export
#'
#' @examples
bootstrap_pspm <- function(model = NULL, 
                           PSPMLearn = NULL,
                           cl = NULL,
                           n_boot_iter = 100, 
                           burnin = 1e2,
                           ci_level = 0.95, 
                           return_sims = FALSE, 
                           cache_samples = T,
                           beta = NULL){
  
  # Checks
  if(!((!is.null(model) & is.null(PSPMLearn)) |
       (is.null(model) & !is.null(PSPMLearn)))){
    stop("Either model or PSPMLearn must be specified but not both.")
  }
  
  # Run bootstrap
  if(!is.null(model)){
    if(! "learn_obj" %in% names(model)){
      stop("Modelmust include a PSPMLearn object called 'learn_obj'. Set return_pspm = TRUE when running fit_pspm_model().")
    }
    
    ## Reset Beta to old value
    if(!is.null(beta)){
      old_beta <- model$PSPMLearn$pspm_ls[[1]]$get_beta()
      model$PSPMLearn$reset_beta(beta)
      cache_samples <- FALSE
    }
    
    ## Run bootstrap
    bs <- model$PSPMLearn$par_bootstrap_composite_log_likelihood(cl = cl, n_boot_iter = n_boot_iter, burnin = burnin,
                                                                 ci_level = ci_level, 
                                                                 return_sims = return_sims, cache_samples = cache_samples)
    
    ## Reset Beta to new value
    if(!is.null(beta)){
      model$PSPMLearn$reset_beta(old_beta)
    }
  } else {
    
    ## Reset Beta to new value
    if(!is.null(beta)){
      old_beta <- PSPMLearn$pspm_ls[[1]]$get_beta()
      PSPMLearn$reset_beta(beta)
      cache_samples <- FALSE
    }
    
    ## Run bootstrap
    bs <- PSPMLearn$par_bootstrap_composite_log_likelihood(cl = cl, n_boot_iter = n_boot_iter, burnin = burnin,
                                                           ci_level = ci_level, 
                                                           return_sims = return_sims, cache_samples = cache_samples)
    
    ## Reset Beta to old value
    if(!is.null(beta)){
      PSPMLearn$reset_beta(old_beta)
    }
  }
  
  # Return
  return(bs)
}