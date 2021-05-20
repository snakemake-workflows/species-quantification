Tumor microbiome calling pipeline produces short (Illumina) and long (PacBio) read samples given some bacteria species and mimics a tumor microbime environment by having a mixture of bacterial reads with human reads. This snakemake workflow also quantifies abundances of bacteria found in the samples, using kraken2-bracken and sourmash. The results are then compared with each other, outputting scatter plots of each method.