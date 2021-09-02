# General settings

This workflow has to be configured prior to run. It consists of two usages and the `config/config.yaml` should be configured firstly. `benchmarking` should be set to `True` if benchmarking is desired and to `False` if the purpose is only to perform abundance quantification given fastq reads.


# Unit sheet - Abundance quantification

This sheet **has to be** defined if the **purpose** is only to do **abundance quantification** of samples. 

Samples should be added to `config/samples.tsv`. All the columns `sample_name`, `fq1` and `fq2` should be defined.

Two things should be carried out depending on the samples.

If samples are paired-end:

* `fq1` and `fq`should be defined accordingly.
* `paired` should be set to `True` in the `config/config.yaml`.

If samples are not paired-end:

* Only `fq1` column should be defined with the single-end fastq file.
* `paired` should be set to `False` in the `config/config.yaml`.

# Bacteria sheet - Benchmarking

This sheet **has to be** defined if the **purpose** is to do **benchmarking** of short reads (Illumina) and long reads (ONT - Oxford Nanopore Technologies).

After making sure that `benchmarking` is set to `True`, bacteria should be added to `config/bacteria.yaml`. All the columns should be defined.

* For `bacteria` column, any name to define bacteria can be selected (only important point is that it should **not** contain **whitespaces**.)
* For `fasta` column, the relative path belonging to reference fasta sequence of the bacterial species that are desired to be present in the mixture samples should be defined.
* For `bacterium_name` column, exact names of the bacterial species should be defined.

# Additional settings for the benchmarking

The following configurations can be made in the `config/config.yaml'.

* `number_of_samples`, `p` should be defined to have the desired number of mixture samples in the end and to select for which fractions (by number of reads) of bacterial species to add to mixture samples, respectively.
* `short_read_len` and `n_reads_per_seq` can be configured for Art simulator.
* `long_nreads` can be configured for Nanosim simulator.

