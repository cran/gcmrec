
############ First.lib ###############

.First.lib <- function(lib, pkg){
   require(survrec)
   library.dynam("gcmrec", pkg, lib)
   
}
############ End of .First.lib ###############


