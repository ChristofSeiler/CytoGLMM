#' Check status of jobs on sherlock cluster
#'
#' @export
#'
get_status = function(user = "cseiler",registry = "registry_2017_10_26_113306") {

  server = "login.sherlock.stanford.edu"

  cd = "cd /scratch/users/cseiler/Biogen"
  module = "module load R"
  library = "R -e 'library(batchtools)"
  load_registry = paste0("loadRegistry(as.character(substitute(",registry,")))")
  status = "getStatus()'"
  code = paste(cd,module,library,load_registry,status,sep = "; ")

  command = paste0("ssh -t ",user,"@",server," ",'"',code,'"')
  system(command)

}
