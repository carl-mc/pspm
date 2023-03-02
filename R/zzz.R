
########################
# On Load
########################
SpatialLattice <- NULL
PartitionModel <- NULL
abstractmethod <- NULL
.onLoad <- function(libname, pkgname){
  
  # Check whether reticulate works on python 3
  if(basename(py_config()$python) != "python3"){
    stop("Must use Python 3.")
  }
  
  ## Source python code
  source_python(file.path(system.file("python", package = pkgname),
                          "sl3.py"))

  SpatialLattice <<- SpatialLattice
  PartitionModel <<- PartitionModel
  abstractmethod <<- abstractmethod

}

