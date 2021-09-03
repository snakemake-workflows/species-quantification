# General settings #

This workflow has to be configured prior to run. It quantifies the given fastq samples for the presence of species (e.g. bacteria, virus, archaea), using kraken2 and bracken in combination. One use case is that it can be used to analyse the microbiome of a tumor sample.

Optionally, it can also be used to do benchmarking of mixture samples that are generated within the workflow (for Illumina short reads and ONT (Oxford Nanopore Technologies) long reads at the same time). It starts with simulating short and long reads of desired bacterial species and mimics a biological environment that may consist of these species in a human sample. It results in the generation of scatter plots for the quantification of the presence of species in short and long read mixture samples, respectively. `benchmarking` should be set to `True` to perform benchmarking.

# Sample sheet - Abundance quantification #

This sheet **has to be** defined if the **purpose** is only to do **abundance quantification** of samples. 

Samples should be added to `config/samples.tsv`. All the columns `sample_name`, `fq1` and `fq2` should be defined.

Two things should be carried out depending on the samples.

If samples are paired-end:

* `fq1` and `fq2`should be defined accordingly.
* `paired` should be set to `True` in the `config/config.yaml`.

If samples are not paired-end:

* Only `fq1` column should be defined with the single-end fastq file.
* `paired` should be set to `False` in the `config/config.yaml`.

## Bacteria sheet - Benchmarking (Optional) ##

Although this functionality of the workflow is optional, if desired, the below configurations should be made prior to run.

The bacterial sheet **has to be** defined if the **purpose** is to do **benchmarking** of short reads (Illumina) and long reads (ONT).

After making sure that `benchmarking` is set to `True`, bacteria should be added to `config/bacteria.yaml`. All the columns should be defined.

* For `bacteria` column, any name to define bacteria can be selected (only important point is that it should **not** contain **whitespaces**.)
* For `fasta` column, the relative path belonging to reference fasta sequence of the bacterial species that are desired to be present in the mixture samples should be defined.
* For `bacterium_name` column, exact names of the bacterial species should be defined.

### Additional settings for the benchmarking ###

The following configurations can be made in the `config/config.yaml`.

* `number_of_samples`, `p` should be defined to have the desired number of mixture samples in the end and to select for which fractions (by number of reads) of bacterial species to add to mixture samples, respectively.
* `short_read_len` and `n_reads_per_seq` can be configured for Art simulator.
* `long_nreads` can be configured for Nanosim simulator.

