
################################
# FIT PSPM MODEL: WRAPPER AROUND PSPM CLASS
###############################



#' Reduce graph to data necessary to estimate PSPM
#'
#' @param g igraph object
#' @param formula Regression formula\code{y ~ x1 + x2}. y must refer to a vertex-attribute that 
#' encodes partition memberships. x1, x2, etc. must refer to edge-level attributes that 
#' encode edge-wise predictors. 
#' @param network_stat (Empty) list of network statistics. Not used, for future extensions and backwards compatibility. 
#' Type of likelihood used. Only \code{"composite_log_likelihood"} 
#' @param force_contiguous Whether to force partitions encoded by y to be contiguous. 
#' @param na.rm Flag to remove vertices and edges with missing data from the graph. 
#'
#' @return A clean graph to estimate a PSPM
#'
#' @examples
model_graphframe <- function(g, formula, 
                             network_stat = list(),
                             force_contiguous = T, na.rm = F){
  
  # Check args
  stopifnot(inherits(g, 'igraph'))
  stopifnot(inherits(formula, 'formula'))
  
  # Get raw edge variables from igraph
  vars <- c(edge_attr_names(g), vertex_attr_names(g))
  vars <- vars[vars %in% all.vars(formula) & !vars %in% all.vars(formula)[1]]
  names(vars) <- vars
  
  # Remove missing parts from graph
  if(na.rm){
    # Delete verteces
    
    ## Get missings
    keep.vertex <- !is.na(vertex_attr(g, name = all.vars(formula)[1]))
    for(v in unique(unlist(network_stat))){
      keep.vertex <- keep.vertex & !is.na(vertex_attr(g, name = v))
    }
    
    ## Delete
    if(sum(!keep.vertex) > 0){
      warning(paste("Deleting ", sum(!keep.vertex), "vertices"))
      g <- induced_subgraph(g, vids = which(keep.vertex))
    }
    
    # Delete edges
    
    ## Get Edge data
    edge.df <- data.frame(do.call(cbind, lapply(vars, function(v){edge_attr(g, v)})))
    
    ## Delete
    has.predictor <- apply(edge.df, 1, function(x){all(!is.na(x))})
    if(sum(!has.predictor) > 0){
      warning(paste("Deleting ", sum(!has.predictor), "edges"))
      g <- delete.edges(g, which(!has.predictor))
      edge.df <- edge.df[has.predictor, , drop = F] 
    }
  } else {
    edge.df <- data.frame(do.call(cbind, lapply(vars, function(v){edge_attr(g, v)})))
  }
  if(length(vars) > 0){
    stopifnot(length(E(g)) == nrow(edge.df))
  } else {
    edge.df <- data.frame(matrix(vector(), nrow=length(E(g)), ncol=0))
  }
  
  if(length(E(g)) == 0 ){
    return(NULL)
  }
  
  # Create temporary edge 'outcome'
  edge.df[,all.vars(formula)[1]] <- 1
  
  # Check presence of all vars
  stopifnot(all.vars(formula) %in% colnames(edge.df))
  
  # Disentangle formula
  model.df <- model.frame(formula, edge.df)
  
  # Transfer back to graph edges
  for(v in colnames(model.df)[-1]){
    edge_attr(g, v) <- model.df[,v]
  }
  
  # Make the outcome variable
  
  ## Enforce contiguous partitioning
  if(force_contiguous){
    vertex_attr(g, all.vars(formula)[1]) <- 
      recode_contig_partitioning(g, all.vars(formula)[1])
  }
  
  ## Finalize
  y.df <- data.frame(vertex_attr(g, all.vars(formula)[1]))
  colnames(y.df) <- all.vars(formula)[1]
  y.df <-  model.frame(as.formula(paste("~", as.character(formula)[2])), y.df)
  vertex_attr(g, as.character(formula)[2]) <- y.df[, as.character(formula)[2]]
  
  # Save outcome and predictors on graph
  g$outcome <- as.character(formula)[2]
  g$predictor <- colnames(model.df)[-1]
  g$network_stat <- network_stat
  
  # Return
  return(g)
}



#' @title Recode non-contiguous partitioning to be contiguous
#'
#' @description Recodes partition memberships such that every partition is contiguous on the graph. 
#' @param graph igraph object
#' @param part_attr Vertex attribute that encodes partition membership.
#'
#' @return Numeric vector that encodes vertices partition IDs for contiguous partitions. 
#' @import igraph
#'
#' @examples
recode_contig_partitioning <- function(graph, part_attr){
  
  # Input partitionings
  part.input <- vertex_attr(graph, part_attr)
  part.output <- as.character(part.input)
  
  # Vertex temp.ids
  V(graph)$tempid <- seq_along(V(graph))
  
  # Loop to detect connected components by partition
  for(i in na.omit(unique(part.input))){
    # Subset graph
    temp.g <- induced_subgraph(graph, vid = which(part.input == i & !is.na(part.input)))
    
    # Connected components
    concomp <- clusters(temp.g)
    
    # Save
    part.output[V(temp.g)$tempid] <-
      paste0(part.output[V(temp.g)$tempid],
             ".", concomp$membership)
  }
  
  # Return
  return(as.numeric(as.factor(part.output)))
}


#' @title Fit a PSPM Model on a list of graphs
#' @name fit_pspm_model
#' 
#' @return Returns a MaxLik model object. 
#' 
#' @param formula Regression formula\code{y ~ x1 + x2}. y must refer to a vertex-attribut that 
#' encodes partition memberships. x1, x2, etc. must refer to edge-level attributes that 
#' encode edge-wise predictors. 
#'
#' @param g_ls List of igraph objects
#' @param network_stat (Empty) list of network statistics. Not used, for future extensions and backwards compatibility. 
#' @param model_type Type of likelihood used. Only \code{"composite_log_likelihood"} 
#' is implemented at the moment but future extensions are possible. 
#' @param sigma2 Penalization parameter. Defaults to 10. 
#' @param na.rm Remove missings from the data (both vertices and edges). 
#' @param force_contiguous Enforce partitions' contiguity by resetting partition value such that 
#' non-contiguous partitions are split. 
#' @param split_connect_comp Split graph into connected components. This maybe necessary for graphs with unconnected components. 
#' @param min_component_size Minimum number of vertices per connected components 
#' (can be used to delete disconnected "islands" with only one vertex which 
#' by definition form their own partition if partitions are enforces to be contiguous)
#' @param vertex_coords Input spatial vertex coordinates to keep alongside the model. 
#' As matrix or SpatialPoints* object. If \code{na.rm = T}, 
#' coordinates will be deleted where vertex information is missing. 
#' @param return_pspm Return PSPMLearn object alongside MaxLik model. Necessary to bootstrap standard errors.
#' @param return_g_ls Return (reduced) graph objects alongside MaxLik model. Useful for plotting purposes. 
#'

#' 
#' @details 
#'      
#' 
#' @import infotheo
#' @export
fit_pspm_model <- function(formula, g_ls, 
                            network_stat = list(),
                            model_type = "composite_log_likelihood", 
                            sigma2 = 10, na.rm = TRUE, 
                            force_contiguous = TRUE,
                            split_connect_comp = FALSE,
                            min_component_size = 2,
                            vertex_coords = NULL,
                            return_pspm = FALSE, return_g_ls = FALSE){

  # If single igraph, make it a list
  if(inherits(g_ls, "igraph")){
    g_ls <- list(g_ls)
  } else {
    stopifnot(inherits(g_ls, 'list'))
  }
  
  # Check for duplicate networks stats
  stopifnot(all(!duplicated(names(network_stat))))
  
  # Prepare graphs
  g_ls <- lapply(g_ls, function(g){
    model_graphframe(g = g, 
                     formula = formula, 
                     network_stat = network_stat,
                     na.rm = na.rm, 
                     force_contiguous = force_contiguous)
  })
  
  # Delete NULL graphs
  g_ls[sapply(g_ls, is.null)] <- NULL
  
  # Check
  if(length(g_ls) == 0){
    stop("All data is missing.")
  }
  
  
  # Drop small connected components
  g_ls <- lapply(g_ls, function(g){
    # Find clusters
    clusters <- clusters(g)
    
    # Drop small ones
    induced_subgraph(g, vids = which(clusters$membership %in% which(clusters$csize >= min_component_size)))
  })
  
  # Split into connected components
  if(split_connect_comp){
    g_ls <- unlist(lapply(g_ls, function(g){
      # Find clusters
      clusters <- clusters(g)
      
      # Split
      lapply(unique(clusters$membership), function(c){
        induced_subgraph(g, vids = which(clusters$membership == c))
      })
    }), recursive = F)
  }
  
  
  # Make PSPM Object
  pspm_ls <- lapply(g_ls, function(g){
    igraph2PSPM(g = g, 
                   outcome_name = g$outcome, 
                   edge_pred_names = g$predictor, 
                   network_stat = network_stat,
                   vertex_coords = vertex_coords, 
                   verbose = F, force_contiguous = force_contiguous)
  })
  
  # Make learning object
  learn_obj <- PSPMLearn$new(pspm_ls, sigma2 = sigma2)
  learn_obj$reset_beta()
  beta_init <- pspm_ls[[1]]$get_beta()
  
  # Estimate model
  if(model_type == "composite_log_likelihood"){
    llfit <- learn_obj$fit_composite_log_likelihood(beta = beta_init,
                                                      set_cache = TRUE, 
                                                      clear_cache = TRUE)
  } else {
    stop("Estimation method not available")
  }
  
  # Name parameters decently
  var.names <- c("Constant", 
                 g_ls[[1]]$predictor, 
                 paste0(rep(names(network_stat), sapply(network_stat, length)),
                        rep(".", length(unlist(network_stat))), 
                        unlist(network_stat)))
  names(llfit$estimate) <- var.names
  names(llfit$gradient) <- var.names
  colnames(llfit$hessian) <- var.names
  rownames(llfit$hessian) <- var.names
  names(llfit$estimate) <- var.names
  
  # Add valuable information
  llfit$N <- unlist(lapply(pspm_ls, function(s){s$N}))
  llfit$N_edges <- unlist(lapply(g_ls, function(s){length(E(s))}))
  llfit$N_groups <- unlist(lapply(g_ls, function(s){length(unique(vertex_attr(s, s$outcome)))}))
  llfit$N_instances <- length(g_ls)
  
  # Return
  if(return_pspm){
    llfit$pspm_ls <- pspm_ls
    llfit$learn_obj <- learn_obj
  } 
  if(return_g_ls){
    llfit$g_ls <- g_ls
  }
  return(llfit)
}


#' @name get_fit_statistic
#' @title Compute any fit statistic for predicted partitionings stored in a learn object
#' (produced during parametric bootstrap). 
#' 
#' @return If \code{return_boot = TRUE}, the full distribution across all bootstrap iterations, 
#' otherwise their mean.
#' 
#' @details Inputs are defined as
#'      \code{learn_obj}: As returned from \code{learn_obj$bootstrap_composite_log_likelihood} or 
#'      \code{learn_obj$par_bootstrap_composite_log_likelihood}
#'      
#'      \code{fun}: Function that computes fit statistic for partitionings, e.g. \code{fit_mutualinfo}
#'      
#' 
#' @import infotheo
#' @export
get_fit_statistic <- function(learn_obj, 
                              fun, 
                              return_boot = FALSE){
  ## Weights according to network size
  weights <- sapply(learn_obj$pspm_ls, function(s){length(s$Y)})
  weights <- weights / sum(weights)
  
  ## Fit statistics
  boot.stat <- c()
  for(i in seq_len(n_boot_iter)){
    stat_vec <- c()
    for (m in 1:learn_obj$M) {
      stat_vec <- fun(learn_obj$samples_cache[[i]][[m]],
                      learn_obj$pspm_ls[[m]]$Y)
    }
    boot.stat <- c(boot.stat, sum(stat_vec * weights))
  }
  
  if(return_boot){
    return(boot.stat)
  } else {
    return(mean(boot.stat))
  }
}

#' @name fit_mutualinfo
#' @title Computes the normalized mutual information
#' 
#' @return A number between 0 and 1.
#' 
#' @details Inputs are defined as
#'      \code{X}: Predicted partitioning
#'      
#'      \code{Y}: Observed partitioning
#'      
#' 
#' @import infotheo
#' @export
fit_mutualinfo <- function(X, Y){
  infotheo::mutinformation(X, Y) /  sqrt(infotheo::entropy(Y)*infotheo::entropy(X))
}
