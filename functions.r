suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(gtools))

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

getFittedVals <- function(y, x)
{
	n <- length(x)
	bhat <- cov(y, x) / var(x)
	ahat <- mean(y) - bhat * mean(x)
	fitted <- ahat + x * bhat
	return(fitted)
}

tsls <- function(y, x, g)
{
	xhat <- getFittedVals(x, g)
	# res <- cov(y, xhat) / var(xhat)
	res <- fastAssoc(y, xhat)
	return(res)
}

init_dat <- function(n, p)
{
	r <- matrix(0, p, p)
	diag(r) <- 1
	g <- scale(matrix(rbinom(n*p, 2, 0.5), n))
	
	d <- list()
	for(i in 1:p)
	{
		d[[i]] <- g[,i] + rnorm(n)
	}
	d <- do.call(cbind, d)
	return(list(r=r, d=d, g=g))
}

make_edge <- function(i, j, eff, dat)
{
	dat$r[j, i] <- eff
	dat$d[,j] <- dat$d[,j] + dat$d[,i] * eff
	return(dat)
}

graph_mr <- function(dat)
{
	p <- ncol(dat$d)
	b <- matrix(1, p, p)
	se <- matrix(0, p, p)
	for(i in 1:p)
	{
		for(j in 1:p)
		{
			if(i != j)
			{
				message(i, j)
				a <- tsls(dat$d[,i], dat$d[,j], dat$g[,j])
				b[i,j] <- a$bhat
				se[i,j] <- a$se
			}
		}
	}
	return(list(b=b, se=se))
}

get_paths <- function(first, last, size)
{
	stopifnot(first <= size)
	stopifnot(last <= size)
	stopifnot(size > 2)
	a <- c(first, c(1:size)[-c(first, last)], last)
	combs <- 3:size
	l <- list()
	for(i in combs)
	{
		b <- permutations(size, i, a)
		b <- b[b[,1] == first & b[,i] == last, , drop=FALSE]
		b1 <- b[,-c(1, ncol(b)), drop=FALSE]
		index <- apply(b1, 1, function(x) all(diff(x) >= 1))
		l[[i-2]] <- b[index, , drop=FALSE]
	}
	return(l)
}


get_prods <- function(paths, mat)
{
	s <- 0
	for(i in 1:length(paths))
	{
		p <- paths[[i]]
		for(j in 1:nrow(p))
		{
			r <- p[j, ]
			l <- length(r) - 1
			out <- rep(0, l)
			for(k in 1:l)
			{
				out[k] <- mat[r[k], r[k+1]]
			}
			s <- s + prod(out)
		}
	}
	return(s)
}


inversion_method <- function(res)
{
	a <- -solve(res$b)
	diag(a) <- 1
	return(a)
}



mediation_method <- function(res)
{
	mat <- res$b
	n <- nrow(mat)
	mmat <- matrix(0, nrow(mat), ncol(mat))
	for(i in 1:n)
	{
		for(j in 1:n)
		{
			message(i, j, sep=" ")
			if(i == j)
			{
				mmat[i,j] <- 1
			} else {
				p <- get_paths(i, j, n)
				mmat[i,j] <- mat[i,j] - get_prods(p, mat)
			}
		}
	}
	return(mmat)
}


deconvolution_method <- function(res)
{
	mat <- res$b
	out <- mat %*% solve(diag(nrow(mat)) + mat)
	return(out)
}


plot_from_matrix <- function(mat, title="")
{
	diag(mat) <- 0
	net <- graph.adjacency(round(t(mat), 1), weighted=TRUE, mode="directed")
	# E(net)$width <- E(net)$weight
	plot(net, edge.label = E(net)$weight, main=title)
}


