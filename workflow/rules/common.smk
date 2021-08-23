import pandas as pd

configfile: "config/config.yaml"

#read sample sheet
samples = pd.read_csv(config["samples"], sep ="\t")
units = pd.read_csv(config["units"], sep ="\t", ).set_index(["sample_name", "unit_name"], drop=False)

#only if benchmarking is set to true
bacteria = pd.read_csv(config["bacteria"], sep ="\t"). set_index(["bacteria"], drop=False)

#input function to retrieve fastq samples
def get_fastq_input(wildcards):
   files = units.loc[wildcards.sample].loc[wildcards.unit, "fq"]
   return [files]

#input function to retrieve bacteria, only required when benchmarking set to TRUE
def get_bacteria_input(wildcards):
   files = bacteria.loc[wildcards.bac_ref, "fasta"]
   return [files]
