source("function_lib.R")


run_tests_subgr(base=8,9,10,broke=FALSE,sparsity=0.5,cf=0.2,pl=0.5)
# The largest base size for the simulations to be finished in a reasonable time
run_tests_subgr(base=10,10,10,broke=FALSE,sparsity=0.5)
#run_tests_sparse(5,5, broke = TRUE)
#do_test(4,10,100,0,0,0,edgeset=100,broken=FALSE,sparsity=-1)

dat <- init_data(50000,4)
edg <- cbind(c(1,2),c(1,3),c(1,4),c(2,4),c(2,3))
#edg <- getRandomDag(5,0.5)
dat <- single_test(4, 50000, 0, 0, 0, prRes=TRUE, edgeset = edg, confset = conf, data = dat, broken=FALSE)

# slow!
dat <- init_data(50000,8)
edg <- cbind(c(1,2),c(1,3),c(1,4),c(2,4),c(2,3))
#edg <- getRandomDag(5,0.5)
dat <- single_test(8, 50000, 0, 0, 0, prRes=TRUE, edgeset = edg, confset = conf, data = dat, broken=FALSE)


dat <- init_data(50000,100)
edg <- cbind(c(1,2),c(1,3),c(1,4),c(2,4),c(2,3))
#edg <- getRandomDag(5,0.5)
dat <- single_test(100, 50000, 0, 0, 0, prRes=TRUE, edgeset = edg, confset = conf, data = dat, broken=FALSE)

