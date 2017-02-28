# bmi
# chd
# t1d
# t2d
# urate
# ldl
# smoking

library(TwoSampleMR)
toggle_dev("test")

ao <- available_outcomes()

exposure_dat <- extract_instruments(c(2, 7, 23, 285, 1055, 300, 961))
table(exposure_dat$exposure)
outcome_dat <- extract_outcome_data(unique(exposure_dat$SNP), c(2, 7, 23, 285, 1055, 300, 961))

dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- subset(dat, exposure != "Cigarettes smoked per day || TAG || 2010 || Cigarettes per day")
dat <- subset(dat, outcome != "Cigarettes smoked per day || TAG || 2010")
table(dat$exposure, dat$outcome)

source("rucker.r")

out <- group_by(dat, exposure, outcome) %>%
do({
	x <- .
	x <- subset(x, mr_keep)
	message(x$exposure[1], " ", x$outcome[1])
	res <- rucker(x, Igxthresh = 0)
	dat <- data.frame(b=res$b, se=res$se, nsnp=res$nsnp, pval=res$pval, model=res$model)
	return(dat)
})
out$exposure <- sapply(strsplit(as.character(out$exposure), split="\\|"), function(x) x[1]) %>% gsub(" $", "", .)
out$outcome <- sapply(strsplit(as.character(out$outcome), split="\\|"), function(x) x[1]) %>% gsub(" $", "", .)

resr <- list()
resr$b <- spread(select(out, exposure, outcome, b), key=exposure, value=b) %>% subset(., select=-c(outcome)) %>% as.matrix()
resr$se <- spread(select(out, exposure, outcome, se), key=exposure, value=se) %>% subset(., select=-c(outcome)) %>% as.matrix()
resr$pval <- spread(select(out, exposure, outcome, pval), key=exposure, value=pval) %>% subset(., select=-c(outcome)) %>% as.matrix()
diag(resr$b) <- 1
diag(resr$se) <- 0
diag(resr$pval) <- 0
rownames(resr$b) <- colnames(resr$b)
rownames(resr$se) <- colnames(resr$se)
rownames(resr$pval) <- colnames(resr$pval)
resb <- bootstrap_graphs(resr)

save(exposure_dat, outcome_dat, dat, out, resr, resb, file="empirical_analysis.rdata")


par(mfrow=c(1,2))
pl1 <- resb$b
pl1[resb$pval > 0.01] <- 0
plot_from_matrix(pl1, "Deconvolved")

pl2 <- resr$b
pl2[resr$pval > 0.05] <- 0
plot_from_matrix(pl2, "Raw")

