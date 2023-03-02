##################################
# SCCRF FIT OUTPUT TO LATEX / HTML
##################################

require(texreg)

# 
#' Transfer fitted PSPM Model output to tex table
#' 
#' @description Builds directly on Philip Leifeld's texreg package. 
#'
#' @param model.ls List of models returned from \code{fit_pspm_model} or \code{PSPMLearn$fit_composite_log_likelihood}
#' @param type One of latex, text or html
#' @param add.stats Statistics to add. 
#' @param custom.gof.rows Add custom rows.
#' @param stars Stars to add (only if native standard errors are used)
#' @param bootci List of bootstrapped CIs as returns from \code{bootstrap_pspm_model} or  or \code{PSPMLearn$bootstrap_composite_log_likelihood}
#' @param boottype Type of bootstrapped standard error (basic or percentile)
#' @param ... Additional inputs taken by texreg (latex) / screenreg (text) / htmlreg (html), 
#' depending on \code{type}
#'
#' @return
#' @export
#' @import texreg
#'
#' @examples
pspm2table <- function(model.ls, type = "latex", 
                        add.stats = c("N_groups", "N_edges", "N_instances","MutInfo"),
                        custom.gof.rows = NULL, 
                        stars = c(.01, .05, .1), 
                        bootci = F, boottype = "percentile", ...){
  
  
  # Extract info
  
  ## Add Mutual Information if asked for
  if("MutInfo" %in% add.stats){
    model.ls <- lapply(model.ls, function(m){
      m$MutInfo = get_fit_statistic(m$learn_obj, fit_mutualinfo, return_boot = FALSE)
      m
    })
  }

  
  ## Add to Table
  all.stats <- c("N_groups", "N_edges", "N_instances","MutInfo")
  for(i in rev(add.stats)){
    custom.gof.rows <- c(list(unlist(lapply(model.ls, function(m){sum(m[[i]])}))),
                         custom.gof.rows)
    name <- names(rev(add.stats)[rev(add.stats) == i])
    if(is.null(name)){
      if(type == "latex"){
        name <- rev(c("N_{groups}", "N_{edges}", "N_{instances}", "Mutual Info"))[
          rev(add.stats) == i]
      } else {
        name <- rev(add.stats)[rev(add.stats) == i]
      }
    }
    names(custom.gof.rows)[1] <- name
  }
  
  # Bootstrapped standard errors
  if(boottype == "percentile"){
    b.rows <- 3:4
  } else {
    b.rows <- 1:2
  }
  
  if(is.list(bootci)){
    override.ci.low = lapply(bootci, function(m){m[b.rows[1],]})
    override.ci.up = lapply(bootci, function(m){m[b.rows[2],]})
  } else if(bootci){
    override.ci.low = lapply(model.ls, function(m){m$bootci[b.rows[1],]})
    override.ci.up = lapply(model.ls, function(m){m$bootci[b.rows[2],]})
  } else {
    override.ci.low = 0
    override.ci.up = 0
  }
  
  # To Table
  if(type == "latex"){
    table <- texreg(model.ls, 
              stars = stars,custom.gof.rows = custom.gof.rows,
              caption.above = T, 
              override.ci.low = override.ci.low, override.ci.up = override.ci.up, 
              ...)
    
  } else if(type == "text"){
    table <- screenreg(model.ls, 
              stars = stars,custom.gof.rows = custom.gof.rows,
              caption.above = T, 
              override.ci.low = override.ci.low, override.ci.up = override.ci.up, 
              ...)
    
  } else if(type == "html"){
    table <- htmlreg(model.ls, 
                       stars = stars,custom.gof.rows = custom.gof.rows,
                       caption.above = T, 
                     override.ci.low = override.ci.low, override.ci.up = override.ci.up, 
                     ...)
  } else {
    stop("input valid type.")
  }
  return(table)
}


# extend texreg packages
extract.pspmfit <- function(model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]
  loglik.value <- s$loglik
  n <- sum(model$N)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value)
    gof.names <- c(gof.names, "Log-Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}
