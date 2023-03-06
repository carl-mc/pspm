
########################
# On Load
########################
SpatialLattice <- NULL
PartitionModel <- NULL
abstractmethod <- NULL
.onLoad <- function(libname, pkgname){

  ## Source python code
  source_python(file.path(system.file("python", package = pkgname),
                          "sl3.py"))
  SpatialLattice <<- SpatialLattice
  PartitionModel <<- PartitionModel
  abstractmethod <<- abstractmethod
  

}

