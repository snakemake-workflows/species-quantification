Tumor microbiome calling pipeline produces short (Illumina) and long (ONT) read samples given some bacteria species and mimics a tumor microbime environment by having a mixture of bacterial reads with human reads. This snakemake workflow also quantifies abundances of bacteria found in the samples, using kraken2-bracken and sourmash. The results are then compared with each other, outputting scatter plots of each method.

This workflow contains two usages:

1. Benchmarking: Please modify config to specify benchmarking.
2. Without benchmarking (only abundance quantification of given samples): Please modify config.yaml to specify benchmarking.

