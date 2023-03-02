
#' @name is_von_neumann_neighbor
#' @title Checks whether nodes i and j are von Neumann neighbors on a square grid
#' with side length N_sqrd.
#' @keywords internal
is_von_neumann_neighbor <- function(i, j, N_sqrd) {
  # Checks whether nodes i and j are vN neighbors on a square grid
  # with side length N_sqrd.
  
  # no self-neighbors
  if (i == j) {
    return(FALSE)
  }
  
  irow <- (i-1)%%N_sqrd + 1
  jrow <- (j-1)%%N_sqrd + 1
  icol <- (i-1)%/%N_sqrd + 1
  jcol <- (j-1)%/%N_sqrd + 1
  
  # check vertical neighbors
  if (abs(irow - jrow) == 1 & icol == jcol) {
    return(TRUE)
  }
  
  # check horizontal neighbors
  if (irow == jrow & abs(icol-jcol) == 1) {
    return(TRUE)
  }
  
  # else not neighbors
  return(FALSE)
}

#' @name is_moore_neighbor
#' @title Checks whether nodes i and j are Moore neighbors on a square grid
#' with side length N_sqrd.
#' @keywords internal
is_moore_neighbor <- function(i, j, N_sqrd) {
  # Checks whether nodes i and j are Moore neighbors on a square grid
  # with side length N_sqrd.
  
  # no self-neighbors
  if (i == j) {
    return(FALSE)
  }
  
  irow <- (i-1)%%N_sqrd + 1
  jrow <- (j-1)%%N_sqrd + 1
  icol <- (i-1)%/%N_sqrd + 1
  jcol <- (j-1)%/%N_sqrd + 1
  
  # whether manhattan distance is at most 1
  if (abs(irow - jrow) <= 1 & abs(icol - jcol) <= 1) {
    return(TRUE)
  }
  
  # else not neighbors
  return(FALSE)
}

#' @name is_hexagonal_neighbor
#' @title Checks whether nodes i and j are neighbors on a square grid
#' with hexagonally sampled points
#' with side length N_sqrd.
#' @keywords internal
is_hexagonal_neighbor <- function(i, j, N_sqrd) {
  # Checks whether nodes i and j are neighbors on a square grid with hexagonally sampled points
  # with side length N_sqrd.
  
  # no self-neighbors
  if (i == j) {
    return(FALSE)
  }
  
  irow <- (i-1)%%N_sqrd + 1
  jrow <- (j-1)%%N_sqrd + 1
  icol <- (i-1)%/%N_sqrd + 1
  jcol <- (j-1)%/%N_sqrd + 1
  
   # check vertical neighbors
  if (abs(irow - jrow) == 1 & icol == jcol) {
    return(TRUE)
  }
  
  # check horizontal neighbors
  if (irow == jrow & abs(icol-jcol) == 1) {
    return(TRUE)
  }
  # check diagonal neighbors 
  if ( ((irow %% 2) == 0 & (icol %% 2) == 1 & irow-jrow == 1 & abs(icol-jcol) == 1) | 
      ((irow %% 2) == 1 & (icol %% 2) == 0 & irow-jrow == -1 & abs(icol-jcol) == 1) |
        ((irow %% 2) == 1 & (icol %% 2) == 1 & irow-jrow == 1 & abs(icol-jcol) == 1) | 
      ((irow %% 2) == 0 & (icol %% 2) == 0 & irow-jrow == -1 & abs(icol-jcol) == 1)) {
    return(TRUE)
  }
  
  # else not neighbors
  return(FALSE)
}

#' @name is_separated_by_col
#' @title Checks whether units i and k are separated by sep_col in a 
#' grid with side length N_sqrd
#' @keywords internal
is_separated_by_col <- function(i, j, N_sqrd, sep_col) {
  # Checks whether units i and k are separated by sep_col in a
  # grid with side length N_sqrd
  
  # no self-neighbors
  if (i == j) {
    return(FALSE)
  }
  
  icol <- (i-1)%/%N_sqrd + 1
  jcol <- (j-1)%/%N_sqrd + 1
  
  # check seperation
  if (icol <= sep_col & jcol > sep_col) {
    return(TRUE)
  }
  
  if (icol > sep_col & jcol <= sep_col) {
    return(TRUE)
  }
  
  # else not separated
  return(FALSE)
}

#' @name generate_square_grid_lattice
#' @title Returns adjacency matrix for square lattice of side-length \code{N_sqrd}
#' @export
generate_square_grid_lattice <- function(N_sqrd, dep_structure = c("von_neumann", "moore", "hexagonal")) {
  # Returns adjacency matrix for square grid lattice
  
  dep_structure <- match.arg(dep_structure)
  
  # Make NxN dependency matrix
  N <- N_sqrd^2
  pairs <- expand.grid(1:N, 1:N)
  p1 <- pairs[,1]
  p2 <- pairs[,2]
  
  if (dep_structure == "von_neumann") {
    dep_function <- is_von_neumann_neighbor
  } else if (dep_structure == "moore") {
    dep_function <- is_moore_neighbor
  } else if (dep_structure == "hexagonal") {
    dep_function <- is_hexagonal_neighbor
  }
  
  is_neighbor <- mapply(function(a, b) {
    dep_function(a, b, N_sqrd)
  }, p1, p2)
  A <- matrix(is_neighbor, N, N)*1
  
  return(A)
}

#' @name generate_grid_data
#' @title Generates square lattice with sampled partitioning
#' 
#' @return A PSPM object based on sqare lattice of side-length \code{N_sqrd}.
#' 
#' @details Inputs are defined as
#'      \code{N_sqrd}: Integer; side length of the lattice. 
#'      
#'      \code{beta0}: Baseline repulsion (positive) or attraction (negative) among nodes.  
#'      
#'      \code{beta}: Effect of predictors. 
#'      Each predictor separates the lattice into a left and right part. Parts do not align.
#'      
#'      \code{dep_structure}: Dependency structure. 
#'      One of \code{c("von_neumann", "moore", "hexagonal")}
#'      
#'      \code{burnin}: Length of burn-in period during sampling. 
#'      
#'      \code{temperatures}: Not used. For potential parallel tempering.
#'      
#'      \code{swap_iter}: Not used. For potential parallel tempering with swapping. 
#'      
#' 
#' @export
generate_grid_data <- function(N_sqrd,
                               beta0,
                               beta = NULL,
                               dep_structure = c("von_neumann", "moore", "hexagonal"),
                               burnin = 10,
                               temperatures = NULL, swap_iter = NULL) {
  
  dep_structure <- match.arg(dep_structure)
  
  # Make NxN dependency matrix
  N <- N_sqrd^2
  pairs <- expand.grid(1:N, 1:N)
  p1 <- pairs[,1]
  p2 <- pairs[,2]
  A <- generate_square_grid_lattice(N_sqrd, dep_structure = dep_structure)
  
  # Generate distance covariates (column separators, one for each beta)
  if (!is.null(beta)) {
    L <- length(beta)
    stopifnot(L < N_sqrd - 1)
    separating_cols <- floor(seq(1, N_sqrd, length.out = L+2))[c(-1, -(L+2))]
    X.ls <- vector('list', L)
    for (l in 1:L) {
      is_distant <- mapply(function(a, b) {
        is_separated_by_col(a, b, N_sqrd, separating_cols[l])
      }, p1, p2)
      D <- A*matrix(is_distant, N, N)*1
      X.ls[[l]] <- D
    }
  } else {
    X.ls <- list()
  }
  
  # Placeholder Y
  Y <- 1:N
  
  # Create the pspm object
  pspm <- PSPM$new(Y = Y, Z = X.ls, X =  matrix(Y, ncol = 1), 
                    A = A, force_contiguous=TRUE, verbose = FALSE)
  pspm$set_beta(c(beta0, beta))
  
  # Sample a new lattice partitioning
  if(!is.null(temperatures) | !is.null(swap_iter)){
    stop("Parallel tempering not (yet) implemented")
  } else {
    pspm$sample(burnin = burnin)
  }
  
  # Return
  return(pspm)
}
