# Get all the functions
source("function_lib.R")

# define all parameters
#param <- expand.grid(
#  nodes = c(3,5,10,15,20,25,30),
#  samples = c(100,1000,5000),
#  sp = seq(0.1,1,by=0.1),
#  conf = seq(0.1,1,by=0.1),
#  pl = seq(0.1,1,by=0.1))


# Get the chunk and chunk size from the command line
args <- commandArgs(T)
job_id <- as.numeric(args[1])
job_size <- as.numeric(args[2])
datadir <- args[3]

param <- expand.grid(
  nodes = c(10),
  samples = c(1000),
  sp = seq(0.1,1,by=0.1),
  conf = seq(0,1,by=0.1),
  pl = seq(0,1,by=0.1))

#save(param, file=paste0(datadir,"/parameters.rdata"))

# Static parameters
iter=100
edges=0
cycles=0
cycle_size=0
broke=FALSE

# Define the section of params to run for this job
start <- (job_id - 1) * job_size + 1
end <- min(job_id * job_size, nrow(param))
param <- param[start:end, ]

# Run only the section for this job
l <- list()
for(i in 1:nrow(param))
{
	out <- do_test(iter, param[i, ]$nodes, param[i, ]$samples, edges, cycles, cycle_size, edgeset=param[i, ]$nodes,broken=broke,sparsity=param[i, ]$sp,cf=param[i, ]$conf,pl=param[i, ]$pl)
	for(j in names(param))
	{
		out[[j]] <- param[[j]][i]
	}
	l[[i]] <- out
}
l <- dplyr::bind_rows(l)
save(l, file=paste0(datadir,"/sim_data/out", job_id, ".rdata"))
