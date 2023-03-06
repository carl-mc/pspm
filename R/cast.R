########################################
# Casting functions
########################################
#' @name PSPM2igraph
#' @title Cast PSPM as igraph object
#' 
#' @return Returns an igraph object that encodes all the information stored in the 
#' PSPM object as vertex and edge-level attributes. 
#' 
#' @param pspm A PSPM object to be cast as an igraph. 
#' 
#' @import sp
#' igraph
#' @export
PSPM2igraph <- function(pspm) {

  # Check args
  stopifnot(inherits(pspm, 'PSPM'))
  
  # Graph from adjacency
  A <- pspm$A
  g <- igraph::graph_from_adjacency_matrix(A, mode = 'undirected')
  
  # Add edge attributes
  edge_ind <- igraph::get.edgelist(g)
  Z <- pspm$Z
  L <- length(Z)
  for (l in 1:L) {
    X <- Z[[l]]
    x_edges <- X[edge_ind]
    if (is.null(pspm$Z.labs)) {
      name <- paste0('x', l)
    } else {
      name <- pspm$Z.labs[l]
    }
    g <- igraph::set.edge.attribute(g, name, value = x_edges)
  }
  
  # Add other vertex attributes
  for(c in seq_len(ncol(pspm$X))){
    if(is.null(colnames(pspm$X))){
      x.name <- paste0("X", c)
    } else {
      x.name <- colnames(pspm$X)[c]
    }
    g <- igraph::set.vertex.attribute(g, name = x.name, value = pspm$X[, c])
  }
  
  # Add outcome as vertex attribute
  g <- igraph::set.vertex.attribute(g, name = 'Y', value = pspm$Y)
  
  # Add spatial info
  if(!is.null(pspm$points)){
    g <- igraph::set.vertex.attribute(g, name = 'x', value = pspm$points@coords[,1])
    g <- igraph::set.vertex.attribute(g, name = 'y', value = pspm$points@coords[,2])
  }
  
  # Return
  return(g)
}


#' @name igraph2PSPM
#' @title Cast igraph as PSPM object
#' 
#' @return Returns a PSPM object.
#' 
#' @param g An igraph
#' @param outcome_name Vertex-level variable names as strings that encodes partition IDs, the main outcome. 
#' @param edge_pred_names Edge-level variable names as strings that encode the predictors
#' @param network_stat List of supra-edge network statistics. 
#'      Currently not used. For future backward compatibility. 
#' @param vertex_pred_names List of vertex-level predictors. 
#'      Currently not used. For future backward compatibility. 
#' @param vertex_coords Vector of x and y coordinates of vertices. 
#'      If provided, saves coordinates of vertices with PSPM object as \code{sp::SpatialPoints}.
#' @param verbose Verbosity
#' @param force_contiguous Force contiguous partitioning.
#' 
#' @import sp
#' igraph
#' @export
igraph2PSPM <- function(g, outcome_name, edge_pred_names, 
                           network_stat = list(),
                           vertex_pred_names = NULL, 
                           vertex_coords = NULL,
                           verbose = FALSE,
                           force_contiguous = TRUE){

  # Check args
  stopifnot(inherits(g, 'igraph'))
  
  # Adjacency matrix
  A <- igraph::as_adjacency_matrix(g, type = "both",
                                   attr = NULL, 
                                   names = FALSE, sparse = FALSE)
  
  # Outcome vector
  Y <- vertex_attr(g, outcome_name)
  
  # Predictor matrices
  Z <- lapply(edge_pred_names, function(nm){
    igraph::as_adjacency_matrix(g, type = "both", attr = nm,
                                names = FALSE, sparse = FALSE)
  })
  names(Z) <- edge_pred_names
  
  # Vertex attributes
  if(length(network_stat) > 0){
    vertex_pred_names <- unique(unlist(network_stat))
    X <- do.call(cbind, lapply(vertex_pred_names, function(c){
      vertex_attr(g, c)
    }))
    colnames(X) <- vertex_pred_names
  } else {
    X <- as.matrix(rep(1, length(Y)), ncol = 1)
    colnames(X) <- "X"
  }
  
  # Coordinates
  if(!is.null(vertex_coords)){
    points <- sp::SpatialPoints(cbind(vertex_attr(g, vertex_coords[1]),
                                  vertex_attr(g, vertex_coords[2])))
  } else {
    points <- NULL
  }
  
  # Make object and return
  pspm <- PSPM$new(Y = Y,
                         Z = Z,
                         X = X,
                         A = A, 
                         net_stats = c("LatticeNeighborsNS",names(network_stat)),
                         points = points,
                         verbose = verbose,
                         force_contiguous = force_contiguous)
  return(pspm)
}

