
# Causal network deconvolution

The project is sub-divided into two parts:
1. Simulations
2. Application to real data




1. Simulation study

Questions to answer in this section:
• How bad is confounding for ND-correlation matrix?
• How bad is pleiotropy for ND-MR matrix?
• Sample size differences between ND-cor and ND-MR?
• How are cycles handled between the two methods?
• Missing nodes impact on ND-cor vs ND-MR
• Sparsity impact on ND-cor vs ND-MR
• Measurement error impact on ND-cor vs ND-MR

Decisions to make for analysis
• Number of iterations - look at the standard error of the points in the existing results - if it is very high then need more iterations
• Number of nodes range: 10, 100, 1000?
• Sample size range: 100 or 5000
• Confounding level: 0-100
• Pleiotropy level: 0-100
• Consider higher order confounding later on - e.g. confounding by cell type in bulk tissue analysis of gene expression levels
• Does the intersection between ND-cor and ND-MR improve AUC?



Simulations to evaluate how to deconvolve the matrix of indirect causal effects into one of direct causal effects. Using network deconvolution method.

## To do

- Compare against ND of observational correlation matrix
- Use ND normalisation method, test on different levels of sparseness
- Evaluate performance with
    - cycles in the graph
    - measurement error
    - unmeasured confounding / incomplete information
- Evaluate sensitivity to biases in the MR matrix
- Compare to multivariable MR

## Setup

Create a `config.json` file with relevant paths:

```json
{
    "genodir": "",
    "phendir": "",
    "igddir": ""
    "pipelinedir": "",
    "datadir": "",
    "resultsdir": "",
}
```

* `genodir` = path to ukbb bgen files
* `phendir` = path to phesant formatted phenotype files
* `igddir` = path to GWAS VCF files for ukb-b batch
* `pipelinedir` = path to original ukb-b gwas pipeline dir
* `datadir` = where to store generated data
* `resultsdir` = where to store generated results


## To run everything on bc4

Use Snakemake to orchestrate the analysis. This will submit the jobs to slurm. 
So when you're in a screen session, to run the snakemake process:

```
module add languages/anaconda3/5.2.0-tflow-1.11
source ~/.bash_profile
snakemake -prk \
-j 100 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```

