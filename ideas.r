---
title: Deconvoluting causal networks
author: Gibran Hemani
date: "`r Sys.Date()`"
output: pdf_document
---



```{r, echo=FALSE}
# http://www.mathpages.com/home/kmath198/2-1/2-1.htm
# http://videolectures.net/gpip06_mackay_gpb/

library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, echo=FALSE)
```

```{r, echo=FALSE}
set.seed(1234)
read_chunk("~/repo/cit_measurement_error/scripts/mr_directionality_analysis.R")
opts_chunk$set(warning=FALSE, message=FALSE, cache=TRUE, echo=FALSE)
```

source("functions.r")

## Introduction

Summary data for GWASs on phenotypes is becoming increasingly available, which is leading to the situation where we can test 'everything against everything'. With this capability does arise a number of problems, such as 

- Multiple testing
- Decomposing many indirect effects into a terse set of direct effects

The latter problem can be summarised as follows. Suppose there are five variables of interest, 1-5, and the causal relationships are

```
1 -> 2
2 -> 3
3 -> 4
4 -> 5
```

This can be depicted in graph form as:

```{r }

dat1 <- init_dat(500000, 5)
dat1 <- make_edge(1, 2, -1, dat1)
dat1 <- make_edge(2, 3, -2, dat1)
dat1 <- make_edge(3, 4, 3, dat1)
dat1 <- make_edge(4, 5, 4, dat1)

plot_from_matrix(dat1$r)

```

If, however, we performed MR of 1 against 3, 1 against 4, 1 against 5, 2 against 4, etc, we would identify associations because they exist indirectly. Hence, our graph would look like:

```{r }

res1 <- graph_mr(dat1)
plot_from_matrix(res1$b)

```

The task is to decompose 

```{r eval=FALSE}

dat1 <- init_dat(500000, 5)
dat1 <- make_edge(1, 2, -1, dat1)
dat1 <- make_edge(1, 3, -2, dat1)
dat1 <- make_edge(2, 3, 3, dat1)
dat1 <- make_edge(3, 4, 4, dat1)
dat1 <- make_edge(3, 5, 5, dat1)
res1 <- graph_mr(dat1)

res1

plot_from_matrix(res1$b)

res1o <- get_orthogonal_graph(res1)

res1o
plot_from_matrix(round(res1o,1))
plot_from_matrix(dat1$r)
dat1
names(dat1)
dat1$r

p <- 15
dat2 <- init_dat(500000, p)
for(i in 1:(p-1))
{
	dat2 <- make_edge(i, i + 1, rnorm(1), dat2)
}
res2 <- graph_mr(dat2)
res2o <- get_orthogonal_graph(res2)

round(res2o,1)
round(dat2$r,1)
plot(res2o ~ dat2$r)



dat3 <- init_dat(500000, 4)
dat3 <- make_edge(1,2,2, dat3)
dat3 <- make_edge(2,3,3, dat3)
dat3 <- make_edge(2,4,7, dat3)
dat3 <- make_edge(1,3,4, dat3)
dat3 <- make_edge(3,4,5, dat3)
dat3 <- make_edge(1,4,6, dat3)
res3 <- graph_mr(dat3)

res3

dat3$r
round(res3$b,1)
net3 <- graph.adjacency(t(dat3$r), weighted=TRUE, mode="directed")
E(net3)$width <- abs(E(net3)$weight)
plot(net3)

dat3$r[3,1] - (res3$b[3,1] - get_prods(get_paths(1, 3, 4), round(res3$b, 1)))
dat3$r[3,2] - (res3$b[3,2] - get_prods(get_paths(2, 3, 4), round(res3$b, 1)))
dat3$r[4,2] - (res3$b[4,2] - get_prods(get_paths(2, 4, 4), round(res3$b, 1)))
dat3$r[4,3] - (res3$b[4,3] - get_prods(get_paths(3, 4, 4), round(res3$b, 1)))
dat3$r[4,1] - (res3$b[4,1] - get_prods(get_paths(1, 4, 4), round(res3$b, 1)))


# Get all the direct effect products

# Problem: within the raw matrix there are compound effects, e.g. 1-2 is direct, 1-3 is direct plus indirect, so the 1-3 element is not the orthogonal effect. To get 1-4 (via 2 and 3) we need the direct and indirect paths for all intermediate nodes. 

# Solution: Work from the network rather than the matrix. Make the network sparse, identify all paths between 1 and 4

# Network deconvolution: http://www.nature.com/nbt/journal/v31/n8/full/nbt.2635.html
# https://github.com/gidonro/Network-Deconvolution/blob/master/ND.py



dat3$r %*% solve(diag(4) - dat3$r)

dat3$r %*% (diag(4) + dat3$r + dat3$r %*% dat3$r^2 + dat3$r %*% dat3$r^2 %*% dat3$r^2)

res3$b



# Find all the paths in a sparse adjacency matrix
# Work on these instead of the full matrix
#



dat4 <- init_dat(500000, 3)
dat4 <- make_edge(1,2,2, dat4)
dat4 <- make_edge(2,3,3, dat4)
dat4 <- make_edge(1,3,4, dat4)
res4 <- graph_mr(dat4)

res4

dat4$r

get_prods(get_paths(1, 3, 3), res4$b)



round(-solve(res1$b))


r <- matrix(0,5,5)



r <- matrix(0,5,5)
diag(r) <- 1
r[2,1] <- -1
r[3,1] <- -2
r[3,2] <- 3
r[4,3] <- 4
r[5,3] <- 5

plot_from_matrix(r)





r <- matrix(0,5,5)
diag(r) <- 1
r[1,2] <- 4
r[2,3] <- -2
r[3,4] <- 8
r[4,5] <- -3

plot_from_matrix(r)


effs <- c(1:5)
n <- 500000
g <- scale(matrix(rbinom(n*5, 2, 0.5), n))
a <- g[,1] + rnorm(n)
b <- g[,2] + rnorm(n)
c <- g[,3] + rnorm(n)
d <- g[,4] + rnorm(n)
e <- g[,5] + rnorm(n)
d <- data.frame(a, b, c, d, e)

d[,2] <- d[,2] + d[,1] * 

b <- b + a * -1
c <- c + a * -2
c <- c + b * 3
d <- d + c * 4
e <- e + c * 5



round(b)
round(solve(b))
r
table(r==-round(solve(b)))



# library(glasso)
# out <- glasso(b, rho=0.5)




library(igraph)

net3 <- graph_from_incidence_matrix(abs(round(solve(b))))

dev.new()
plot(net)
dev.new()
plot(net2)
dev.new()
plot(net3)













n <- 10000
nstep <- 100
ab <- 0.1
ba <- 0.2

a <- rnorm(n)
b <- rnorm(n)

for(i in 1:nstep)
{
	a <- a + b * ba
	b <- b + a * ab
}

cor(a,b)
mean(a)
mean(b)

plot(a,b)


```

