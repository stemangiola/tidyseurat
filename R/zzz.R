#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
	version = packageDescription(pkgname, fields = "Version")
	
	msg = paste0("========================================
", pkgname, " version ", version, "
If you use TIDYSEURAT in published research, please cite:

Mangiola et al. Interfacing Seurat with the R tidy universe. Bioinformatics 2021.

This message can be suppressed by:
  suppressPackageStartupMessages(library(tidyseurat))
  
To restore the Seurat default display use options(\"restore_Seurat_show\" = TRUE) 
========================================
")	
	
	packageStartupMessage(msg)
    # Attach tidyverse
    attached <- tidyverse_attach()
}

# rv = R.Version()

# if(getRversion() >= "4.0.0" && as.numeric(rv$`svn rev`) >= 77889) {
# 	unitType = get("unitType", envir = asNamespace("grid"))
# } else {
# 	unitType = function(x, recurse = TRUE) attr(x, "unit")
# }