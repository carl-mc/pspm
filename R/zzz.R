
########################
# On Load
########################
SpatialLattice <- NULL
PartitionModel <- NULL
abstractmethod <- NULL
.onLoad <- function(libname, pkgname){
  
  # Check whether reticulate works on python 3
  if(basename(py_config()$python) != "python3"){
    stop("You must use Python 3.")
  }
  
  # Message to install python modules if not yet installed
  modules <- c("scipy",  "abc", "numpy", "networkx", "typing", "typing", "collections")
  for(m in modules){
    have_m <- py_module_available(m)
    if(!have_m){
      print(warning(paste("Please install python module", 
                          m)))
    }
  }
  
  
  ## Source python code
  source_python(file.path(system.file("python", package = pkgname),
                          "sl3.py"))
  SpatialLattice <<- SpatialLattice
  PartitionModel <<- PartitionModel
  abstractmethod <<- abstractmethod
  

}

