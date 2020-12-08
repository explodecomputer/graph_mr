suppressPackageStartupMessages(suppressWarnings({
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(gtools)
  library(gridExtra)
  library(network)
  library(sna)
  library(ggnetwork)
  library(dagitty)
  #library(MRPC)
  library(glasso)
  library(corpcor)
  library(pROC)
  library(ggpubr)
}))

#data initialisation functions
options(warn=-1)


# empty matrix dims p x p
# set to ident matrix
# n - number of samples
# p - dimension, number of traits
# for any particular phenotype d_i, it's going to be:
# d_ij = b_g*g_ij + SUM(h_k*d_ik) + e_ij,
# where i - samples, j - phenotype
# rbinom <- number of observations, size of trials, probability of success
# scale() forces column centering by subtracting mean of each column from their data points
# d <- split each column into the list and add randomly distributed numbers to each element
# do.call(cbind, d) concatenates the list elements to reform the previous matrix representation

init_data <- function(n, p, pl){
  r <- matrix(0, p, p)
  diag(r) <- 1
  g <- matrix(rbinom(n*p, 2, 0.5), n)
  g <- scale(g)
  d <- list()
  for(i in 1:p){
    d[[i]] <- g[,i] + rnorm(n)
  }
  d <- do.call(cbind, d)
  #which snps are pleiotropic
  plei_snps <- sample(1:p, p * pl, replace=FALSE)
  for(i in plei_snps) {
    rand_trait <- sample(1:p, 1)
    d[,rand_trait] <- d[,rand_trait] + g[,i]
    }
  return(list(r=r, d=d, g=g, p=p))
}


# This function allows for an association between 1 trait and another trait, i has an effect on j
# r is a p x p ident matrix, this sets location j, i to the effect
# d is the scaled data representation, this adds the value of d[,i] * the effect to d[,j]
# conf effect is applied to the node level weightings, not the true graph edges
# conf-r influences both i and j:
# d_ij = b_g*g_ij + SUM(h_k*d_ik) + e_ij + SUM(U)
make_edge <- function(i, j, effect, data)
{
  data$r[j, i] <- effect
  u <- rnorm(nrow(data$d))
  data$d[,i] <- data$d[,i]
  data$d[,j] <- data$d[,j] + data$d[,i] * effect
  return(data)
}


normalise_r <- function(dat)
{
  s <- apply(dat$d, 2, sd)
  o <- dat$r
  for(i in 1:dat$p)
  {
    for(j in 1:dat$p)
    {
      o[i,j] <- o[i,j] * s[j] / s[i]
    }
  }
  dat$r_norm <- o
  return(dat)
}


make_conf <- function(i, j, data, conf_effect=0)
{
  u <- rnorm(nrow(data$d))
  data$d[,i] <- data$d[,i]+ u*conf_effect
  data$d[,j] <- data$d[,j] + u*conf_effect
  return(data)
}

#Work on generating good data variants, cycles etc.
generate_cycle_set <- function(nodes, ncycles=-1, scycle=3){
  if (ncycles == -1) {
    ncycles <- sample(0:nodes$l, 1)
  }
  cycle_lst <- NULL 
  x <- 0
  while (x != ncycles){
    cycle_size <- scycle#sample(3:nodes$l, 1)
    samp <- c(1:cycle_size)#sample(nodes$n, cycle_size)
    cycle_lst <- cbind(cycle_lst, samp)
    x <- x + 1
  }
  return(unique(cycle_lst))
}

# choosing which edges relate to each other (just connections - adjacency matrix)
# randomly generates the graph based on how many nodes, edges etc
generate_edge_set <- function(nodes, nedges=-1){
  if (nedges == -1) {
    nedges <- sample(0:choose(nodes$l,2), 1)
  }
  edge_lst <- NULL
  x <- 0
  while (x != nedges){
    samp <- sample(nodes$n, 2)
    edge_lst <- cbind(edge_lst, c(samp))
    x <- x + 1
  }
  return(unique(edge_lst))
}


#Takes a n x 2 edgeset and augments the data to represent it
set_to_edges <- function(set, conf, data, edge=FALSE){
  if (is.null(set) | length(set) < 1){
    return(data)
  }
  l <- ncol(set)
  n <- ncol(conf)
  
  for (i in 1:l){
    cycle <- set[,i]
    if (edge){
      effect <- runif(1,min=-4,max=4)
      if(effect == 0){
        effect <- 1
      }
      data <- make_edge(cycle[[1]],cycle[[2]], effect, data)
    }else{
      for (j in 1:k){
        effect <- runif(1,min=-4,max=4)
        if(effect == 0){
          effect <- 1
        }
        nxt <- (j+1)%%(k+1)
        if (nxt == 0){
          nxt <- 1
        }
        data <- make_edge(cycle[j],cycle[nxt], effect, data)
        
      }
    }
  }
  
  if (is.null(conf) | length(conf) < 1){
    return(data)
  }

  for (i in 1:n){
    cycle <- conf[,i]
    if (edge){
      conf_ef <- runif(1,min=-2,max=2)
      data <- make_conf(cycle[[1]],cycle[[2]], data, conf_effect = conf_ef)
    }else{
      for (j in 1:k){
        nxt <- (j+1)%%(k+1)
        if (nxt == 0){
          nxt <- 1
        }
        conf_ef <- runif(1,min=-2,max=2)
        data <- make_conf(cycle[j],cycle[nxt], data, conf_effect = conf_ef)
        
      }
    }
  }
  
  return(data)
}


# A wrapper for the above functions
# Either generates edge set and makes edges, or makes edges for a given edge set
graph_gen <- function(ncycles, scycle, nedges, data, edgeset = 0, confset = 0){
  nodes <- list(n=c(1:data$p),l=data$p)
  if(length(edgeset)>0){
    if (edgeset != 0){
      data <- set_to_edges(edgeset, confset, data, TRUE)
    }else{
      # recommended to keep relatively sparse, as edges can accidentally create cycles
      cycles <- generate_cycle_set(nodes, ncycles, scycle=scycle)
      edges <- generate_edge_set(nodes, nedges)
      if (length(cycles) != 0){data <- set_to_edges(cycles, confset, data)}
      if (length(edges) != 0){data <- set_to_edges(edges, confset, data, TRUE)}
    }
  }
  
  data <- normalise_r(data)
  return(data)
}

# MSE

meanSquareError <- function(gr, tGr){
  MSE <- sum((tGr - gr)^2)/length(tGr)
  cor(as.numeric(gr), as.numeric(tGr))^2
}



# After the data is generated, process it


# fast linear regression
# is.finite() & is.finite() pairwise AND comparison of integer finiteness, Inf, Inf, NaN, Na all result in False
# sum tallies the number of True occurences
# y[index] strips y to only finite values
# rsq = correlation^2
# this is the square of the correlation coefficient
# tval = pearson product-moment correlation coefficient fval is the square of tval
# this is identical to the built-in R cor() pearson method
# pf is a cumulative distribution function(values, numerator dgf, denom dgf, tail or not)

fastAssoc <- function(y, x){
  index <- is.finite(y) & is.finite(x)
  n <- sum(index)
  y <- y[index]
  x <- x[index]
  
  vx <- var(x)
  bhat <- cov(y, x) / vx
  ahat <- mean(y) - bhat * mean(x)
  
  
  rsq <- (cov(y, x))^2 / (vx * var(y))
  fval <- rsq * (n-2) / (1-rsq)
  tval <- sqrt(fval)
  se <- abs(bhat / tval)
  
  p <- pf(fval, 1, n-2, lowe=F)
  return(
    list(ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p)
  )
}

# Another level of regression
# Simple linear regression applied to data based on least square estimates bhat and ahat
# bhat = conv(y, x)/var(x)
# ahat = Ex(y) - bhat * Ex(x)
# fitted values Y : Yi = ahat + bhat * xi

getFittedVals <- function(y, x){
  n <- length(x)
  bhat <- cov(y, x) / var(x)
  ahat <- mean(y) - bhat * mean(x)
  fitted <- ahat + x * bhat
  return(fitted)
}

#tsls gets the fitted values for some data, and then applies fastAssoc to the resulting data with the original data
# 2-stage regression:
# 1. x on y
# 2. association between prediction of x on y
tsls <- function(y, x, g){
  xhat <- getFittedVals(x, g)
  res <- fastAssoc(y, xhat)
  return(res)
}

# get pairwise causal effects
# Output matrices: causal effect matrix and standard error
# ncol(data) number of columns in data
# make matrix of size p x p, all ones
# same but all zeroes
# for all none diagonal matrices{
# b non-diagonals set to bhat value, cov(y, x)/var(x), diags = cov(x, x)/var(x) = 1
# se non-diagonals set to se value, abs(bhat/tval)
# }

graph_mr <- function(data){
  p <- ncol(data$d)
  b <- matrix(1, p, p)
  se <- matrix(0, p, p)
  pval <- matrix(0, p, p)
  for(i in 1:p){
    for(j in 1:p){
      if(i != j){
        a <- tsls(data$d[, i], data$d[, j], data$g[, j])
        b[i, j] <- a$bhat
        se[i, j] <- a$se
        pval[i, j] <- a$pval
      }
    }
  }
  return(list(b=b, se=se, pval=pval))
}

# make correlation matrix
# obfuscating function call for my own readability
make_cor_mat <- function(d){
  b <- cor(d)
  se <- sqrt((1-b^2)/(nrow(d)-2))
  pval <- pnorm(abs(b/se), lower.tail=FALSE)
  return(list(b=b, se=se, pval=pval))
}

# Methods:

# a = negative of the inversion of res$b, with diag set to 1, as trait self-relations are trivially related
inversion_method <- function(res){
  #print(dim(res$b))
  a <- -solve(res$b)
  diag(a) <- 1
  return(a)
}

# set the diag of mat to zero, then get the inverse of (ident matrix (size nrow(mat)) + original matrix)
# then matrix multiply with mat
deconvolution_method <- function(res){
  mat <- res$b
  diag(mat) <- 0
  out <- mat %*% solve(diag(nrow(mat)) + mat)
  diag(out) <- 1
  return(out)
}

#' Estimate the bigge
#'
#' <full description>
#'
#' @param res <what param does>
#' @param nboot=100 <what param does>
#' @param minp=1e-300 <what param does>
#' @param fn <what param does>
#' @param ... <what param does>
#'
#' @export
#' @return
bootstrap_deconv <- function(res, nboot=100, minp=1e-300, fn=deconv, ...)
{
  out <- list()
  b <- fn(res$b, ...)

  out$b <- b
  if(!is.null(nboot))
  {
    n <- ncol(res$b)
    l <- list()
    for(i in 1:nboot)
    {
      newb <- matrix(rnorm(n * n, mean=res$b, sd=res$se), n, n)
      temp <- fn(newb, ...)
      l[[i]] <- temp
    }

    addsq <- function(x, y) {
      return(x^2 + y^2)
    }

    m <- Reduce('+', l)
    m2 <- Reduce('+', lapply(l, function(x) x^2))
    
    out$se <- sqrt((m2 - m^2 / nboot) / nboot)
    out$pval <- pnorm(abs(out$b/out$se), lower=FALSE)
    out$pval[out$pval < minp] <- minp
    diag(out$pval) <- 1
  }
  return(out)
}


# function that matrix multiplies x by the inverse of itself plus a same size ident matrix 
deconv <- function(x, remove_diag=TRUE)
{
  if(remove_diag) diag(x) <- 0
  d <- determinant(x)$modulus
  if(is.finite(d))
  {
    o <- x %*% solve(x + diag(nrow(x)))
  } else {
    o <- x %*% corpcor::pseudoinverse(x + diag(nrow(x)))
  }
  if(remove_diag) diag(o) <- 1
  return(o)
}


# normalizes values by a value alpha, dependent on the relative scales of max_e and min_e, being the largest and smallest eigenvalues
normalize_beta <- function(x, beta=0.95){
  alpha <- tryCatch({
    xe <- eigen(x)
    max_e <- max(xe$values)
    min_e <- min(xe$values)
    alpha <- min(beta / ((1-beta) * max_e), -beta / ((1+beta) * min_e))
    alpha
  },
  error=function(cond) {
    return(1)
  })
  return(alpha)
}

#deconv correlation matrix, normalize determines whether to normalize, alpha act as an adjustment, the diagonal 1s are stripped out for deconv()
deconv_corr <- function(corr, normalize=TRUE, alpha=1){
  diag(corr) <- 0
  inp <- corr
  if(normalize && alpha == 1){
    alpha <- normalize_beta(inp)
    deconv_matrix <- deconv(inp*alpha)   
  }else if (alpha == 1)
    deconv_matrix <- deconv(inp)
  else
    deconv_matrix <- deconv(alpha * (inp))
  diag(corr) <- 1
  return(deconv_matrix)
}

# visualisation functions

plot_from_matrix <- function(mat, title="", MSE)
{
  diag(mat) <- 0
  net <- graph.adjacency(round(t(mat), 1), weighted=TRUE, mode="directed")
  layout=layout.circle(net)
  plot(net, edge.label = E(net)$weight, main=title, sub=MSE, layout=layout)
}

plot_from_matrix_clean <- function(mat, title="")
{
  diag(mat) <- 0
  n <- network(round(t(mat), 1))
  p <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    geom_nodes() +
    theme_blank() +
    labs(title=title)
  return(p)
}

# # dagitty facilitator


# Gets a random dag of given nodes and sparsity
getRandomDag <- function(nodes, sparsity){
  dag <- NULL
  dag <- randomDAG(nodes, sparsity)
  
  edg <- makeEdgeList(edges(dag))
  
  return(edg)
}

# takes a DAGitty edge list and reformats it for use
makeEdgeList <- function(edges){
  lst <- apply(edges,1, function(x)(return(c(strtoi(sub('.', '', x['v'])), strtoi(sub('.', '', x['w']))))))
  return(lst)
}

# Data checks

edge_matrix <- function(mat)
{
  tr <- diag(nrow(mat))
  tr[lower.tri(tr)] <- mat[lower.tri(mat)] != 0 | t(mat)[lower.tri(t(mat))] != 0
  tr[upper.tri(tr)] <- t(tr)[upper.tri(tr)]
  tr
}


#' Use p-value to predict edge detection
#'
#' Asymmetric matrix has two opportunities to detect an edge - upper and lower triangles
#' Only compare the off-diagonal elements
#' Get the union of edges from lower and upper triangles from the truth adjacency matrix
#'
#' @param mat Deconv list (includes pval)
#' @param truth True adjacency matrix
#'
#' @export
#' @return
auc_edge_detection <- function(res, truth)
{
  # just use lower triangle - it's symmatrical now
  mat <- -log10(res$pval)
  mat1 <- mat[lower.tri(mat)]
  mat2 <- t(mat)[lower.tri(t(mat))]
  mat <- pmax(mat1, mat2)
  truth <- edge_matrix(truth)
  truth <- truth[lower.tri(truth)]
  tryCatch(
      pROC::roc(truth, mat, quiet=TRUE) %>% pROC::auc(),
    error=function(cond) NA
  )  
}


#' Use -log10 p-values to predict oriented edges
#'
#' @param res Deconv list (includes pval)
#' @param truth True adjacency matrix
#'
#' @export
#' @return
auc_edge_orientation <- function(res, truth)
{
  mat <- -log10(res$pval)
  truth <- truth != 0

  # remove diagonal
  mat <- c(mat[lower.tri(mat)], mat[upper.tri(mat)])
  truth <- c(truth[lower.tri(truth)], truth[upper.tri(truth)])
  tryCatch(
    pROC::roc(truth, mat, quiet=TRUE) %>% pROC::auc(),
    error=function(cond) NA
  )
}


# cor

#' Estimate bias of estimates
#'
#' Similar to MSE but uses correlation - scale free so can compare MR and correlation matrices
#'
#' @param mat estimated graph
#' @param truth True adjacency matrix
#'
#' @export
#' @return
corComparison <- function(mat, truth){
  mat <- c(mat[lower.tri(mat)], mat[upper.tri(mat)])
  truth <- c(truth[lower.tri(truth)], truth[upper.tri(truth)])
  cor(mat, truth)^2
}


sig_graph <- function(res, fdrthresh=0.05)
{
  s <- res$pval %>% p.adjust(., "fdr") %>% {. < fdrthresh} %>% matrix(., nrow(res$pval))
  b <- res$b
  b[!s] <- 0
  return(b)
}


# does everything
single_test <- function(dat, prRes=FALSE, broken=FALSE, save_comparison=FALSE){
  # For missing node tests, strips a data column
  if (broken == TRUE && length(edgeset) > 0){
    rm <- sample(edgeset, 1)[[1]][1]
    
    dat$d <- dat$d[,-rm]
    dat$g <- dat$g[,-rm]
    dat$r <- dat$r[,-rm]
    dat$r <- dat$r[-rm,]
    dat$r_norm <- dat$r_norm[,-rm]
    dat$r_norm <- dat$r_norm[-rm,]
  }
 
  method_list <- c(
    "Total - MR", 
    "Total - Cor", 
    "ND - MR", 
    "ND - Cor", 
    "ND+norm - Cor"
  )
  res1 <- graph_mr(dat)
  res2 <- make_cor_mat(dat$d)
  res3 <- bootstrap_deconv(res1, fn=deconv)
  res4 <- bootstrap_deconv(res2, fn=deconv_corr, normalize=FALSE)
  res5 <- bootstrap_deconv(res2, fn=deconv_corr)

  resall <- list()
  resall$performance <- tibble(
    method = method_list,
    accuracy = c(
      corComparison(res1$b, dat$r),
      corComparison(res2$b, dat$r_norm),
      corComparison(res3$b, dat$r),
      corComparison(res4$b, dat$r_norm),
      corComparison(res5$b, dat$r_norm)
    ),
    edge_detection = c(
      auc_edge_detection(res1, dat$r),
      auc_edge_detection(res2, dat$r_norm),
      auc_edge_detection(res3, dat$r),
      auc_edge_detection(res4, dat$r_norm),
      auc_edge_detection(res5, dat$r_norm)
    ),
    edge_orientation = c(
      auc_edge_orientation(res1, dat$r),
      auc_edge_orientation(res2, dat$r_norm),
      auc_edge_orientation(res3, dat$r),
      auc_edge_orientation(res4, dat$r_norm),
      auc_edge_orientation(res5, dat$r_norm)
    )
  )

  m <- expand.grid(i=1:dat$p, j=1:dat$p)
  
  if(save_comparison)
  {
    resall$comparison <- list(
      tibble(
        method = method_list[1],
        i = m$i,
        j = m$j,
        estimate = c(res1$b),
        truth = c(dat$r)
      ),
      tibble(
        method = method_list[2],
        i = m$i,
        j = m$j,
        estimate = c(res2$b),
        truth = c(dat$r_norm)
      ),    
      tibble(
        method = method_list[3],
        i = m$i,
        j = m$j,
        estimate = c(res3$b),
        truth = c(dat$r)
      ),
      tibble(
        method = method_list[4],
        i = m$i,
        j = m$j,
        estimate = c(res4$b),
        truth = c(dat$r_norm)
      ),
      tibble(
        method = method_list[5],
        i = m$i,
        j = m$j,
        estimate = c(res5$b),
        truth = c(dat$r_norm)
      )
    ) %>% bind_rows()

  }
  
  if(prRes){
    par(mfrow=c(2,3))
    plot_from_matrix(dat$r, "True graph", "")
    plot_from_matrix(sig_graph(res1), method_list[1] , resall$performance$accuracy[1])
    plot_from_matrix(sig_graph(res2), method_list[2] , resall$performance$accuracy[2])
    plot_from_matrix(sig_graph(res3), method_list[3] , resall$performance$accuracy[3])
    plot_from_matrix(sig_graph(res4), method_list[4] , resall$performance$accuracy[4])
    plot_from_matrix(sig_graph(res5), method_list[5] , resall$performance$accuracy[4])
  }
  return(resall)
}

# does lots of different levels of sparsity
do_test <- function(iter, nodes, observations, edges, cycles, cycle_size, edgeset = 0, broken = FALSE, sparsity = 0, cf = 0, pl = 0){
  #print(pl)
  avgRes <- rbind(
    data.frame(mseTot=0,mseAvg=0,method="Feizi Method",aucTot=0,aucAvg=0,aucSd=0),
    data.frame(mseTot=0,mseAvg=0,method="ND Correlation Mat",aucTot=0,aucAvg=0,aucSd=0)
  )
  InvSd <- list()
  FeiSd <- list()
  NDSd <- list()
  avgCompRes <- list()
  for (x in 1:(edgeset+1)){
    avgCompRes <- c(avgCompRes,list(avgRes))
  }
  totalres <- list()
  j <- 1
  for(it in 1:iter){
    #print(it)
    InvSdTmp <- list()
    FeiSdTmp <- list()
    NDSdTmp <- list()
    pr = FALSE
    if (iter == 1){
      pr = FALSE
    }
    dat <- init_data(observations, nodes, pl)
    if(sparsity == -1){
      for(i in seq(0.01,1,length.out=100)){
        sp <- i
        edgeset2 <- getRandomDag(nodes, sp)
        confset2 <- getRandomDag(nodes, cf)
        
        dat <- graph_gen(cycles, cycle_size, edges, data=dat, edgeset = edgeset2, confset = confset2)
        edgeRes <- single_test(dat, prRes = pr, broken=broken, save_comparison=FALSE)$performance
        edgeRes$it <- it
        edgeRes$sp <- i
        totalres[[j]] <- edgeRes
        j <- j + 1
      }
      
    }else{
      
      for(edge_lim in 1:edgeset){
        if (edge_lim == 0) {
          edgeset2 = NULL
        }else if (sparsity > 0){
          edgeset2 <- getRandomDag(edge_lim, sparsity)
        }else{
          edgeset2 = list()
          for(ed in 1 : edge_lim){
            edgeset2 <- cbind(edgeset2, c(ed, ifelse(((ed+1)%%(edgeset+1) == 0), 1, (ed+1))))
          }
        }
        confset2 <- getRandomDag(nodes, cf)
        if(edge_lim < 3){
          broken=FALSE
        }
        dat <- graph_gen(cycles, cycle_size, edges, data=dat, edgeset = edgeset2, confset = confset2)
        edgeRes <- single_test(dat, prRes = pr, broken=broken, save_comparison=FALSE)$performance
        edgeRes$it <- it
        edgeRes$sp <- edge_lim
        totalres[[j]] <- edgeRes
        j <- j + 1
      }
    }
  }

  out <- totalres %>% bind_rows %>%
    tidyr::gather(key="measure", value="value", accuracy, edge_detection, edge_orientation) %>%
    group_by(method, sp, measure) %>%
    summarise(
      m = mean(value, na.rm=TRUE),
      s = sd(value, na.rm=TRUE),
      lci = quantile(value, 0.025, na.rm=TRUE),
      uci = quantile(value, 0.975, na.rm=TRUE)
    )
  return(out)
}


plot_Data <- function(out)
{
  ggplot(out, aes(x=sp, y=m)) +
  geom_point(aes(colour=method)) +
  geom_line(aes(colour=method)) +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour=method), width=0) +
  facet_grid(. ~ measure) +
  scale_colour_brewer(type="qual")
}


# for subgraph test of network sizes base -> limit, very slow
run_tests_subgr <- function(base=3,limit, iter,broke=FALSE,sparsity,cf,pl)
{
  for(i  in base : (limit)){
    #print(i)
    dis <- do_test(iter, i, 2000, 0, 0, 0, edgeset=i,broken=broke,sparsity=sparsity,cf=cf,pl=pl)
    #print(head(dis))
    plot_Data(dis, i,FALSE)
  }
  print("Done")
  
}

# for sparsity tests, for given network size 
run_tests_sparse <- function(nwork, iter, broke=FALSE){
  print("sparsity test")
  dis <- do_test(iter,nwork,2000,0,0,0,edgeset = 99,broken=broke,sparsity=-1,cf=0.5)
  plot_Data(dis, 99,TRUE)
}

