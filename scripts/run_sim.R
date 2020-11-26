# Get all the functions
source("function_lib.R")

# define all parameters
#param <- expand.grid(
#  nodes = c(3,5,10,15,20,25,30),
#  samples = c(100,1000,5000),
#  sp = seq(0.1,1,by=0.1),
#  conf = seq(0.1,1,by=0.1),
#  pl = seq(0.1,1,by=0.1))


param <- expand.grid(
  nodes = c(10,20,30),
  samples = c(1000),
  sp = seq(0.5,0.6,by=0.1),
  conf = seq(0.5,0.6,by=0.1),
  pl = seq(0.5,0.6,by=0.1))

# Static parameters
iter=10
observations=2000
edges=0
cycles=0
cycle_size=0
broke=FALSE

# Get the chunk and chunk size from the command line
args <- commandArgs(T)
job_id <- as.numeric(args[1])
job_size <- as.numeric(args[2])
datadir <- args[3]

# Define the section of params to run for this job
start <- (job_id - 1) * job_size + 1
end <- min(job_id * job_size, nrow(param))
param <- param[start:end, ]

# Run only the section for this job
l <- list()
for(i in 1:nrow(param))
{
	out <- as.data.frame(do_test(iter, param[i, ]$nodes, observations, edges, cycles, cycle_size, edgeset=param[i, ]$nodes,broken=broke,sparsity=param[i, ]$sp,cf=param[i, ]$conf,pl=param[i, ]$pl))
	for(j in names(param))
	{
		out[[j]] <- param[[j]][i]
	}
	l[[i]] <- out
	#l[[i]] <- run_tests_subgr(base=param[i, ]$nodes,param[i, ]$nodes,iter,broke=FALSE,sparsity=param[i, ]$sp,cf=param[i, ]$conf,pl=param[i, ]$pl)
}
l <- dplyr::bind_rows(l)
save(l, file=paste0(datadir,"/sim_data/out", job_id, ".rdata"))
