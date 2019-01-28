# Causal network deconvolution

Simulations to evaluate how to deconvolve the matrix of indirect causal effects into one of direct causal effects. Using network deconvolution method.

To compile:

```
cd scripts/graphmr
Rscript -e "rmarkdown::render('graphmr.rmd', output_format='all')"
```

## To do

- Compare against ND of observational correlation matrix
- Use ND normalisation method, test on different levels of sparseness
- Look at graphical LASSO and use of shrinkage parameter
- Evaluate performance with
    - cycles in the graph
    - measurement error
    - unmeasured confounding / incomplete information
- Evaluate sensitivity to biases in the MR matrix
- Compare to multivariable MR

