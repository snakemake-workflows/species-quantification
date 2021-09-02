# Snakemake workflow: Tumor-microbiome-calling

A Snakemake workflow for producing short (Illumina) and long (ONT) read samples given some bacteria species and mimics a tumor microbime environment by having a mixture of bacterial reads with human reads. This workflow also quantifies abundances of bacterial species found in the samples, using kraken2-bracken. The benchmarking results are compared with each other, outputting scatter plots for short and long read mixture samples.

As already stated above, this workflow contains two usages:

1. Benchmarking
2. Without benchmarking (only abundance quantification of given samples)

The required configuration of these usages can be seen in `config/README.md`.
