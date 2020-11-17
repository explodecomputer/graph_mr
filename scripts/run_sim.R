# Get all the functions
source("Deconvolution_test_mod.R")

# define all parameters
param <- expand.grid(
	1:100,
	1:100
)

# Get the chunk and chunk size from the command line
args <- commandArgs(T)
job_id <- as.numeric(args[1])
job_size <- as.numeric(args[2])

# Define the section of params to run for this job
start <- (job_id - 1) * job_size + 1
end <- min(job_id * job_size, nrow(param))
param <- param[start:end, ]

# Run only the section for this job
l <- list()
for(i in 1:nrow(param))
{
	l[[i]] <- run_simulations(param[i, ])
}


save(l, file=paste("../results/output_", job_id, ".rdata"))

