import pandas as pd

configfile: "config/config.yaml"

#read sample sheet
samples = pd.read_csv(config["samples"], sep ="\t")
units = pd.read_csv(config["units"], sep ="\t", ).set_index(["sample_name", "unit_name"], drop=False)

def get_fastq_input(wildcards):
   files = units.loc[wildcards.sample].loc[wildcards.unit, "fq"]
   return [files]

#only if benchmarking is set to true
bac = pd.read_csv(config["bacteria"], sep ="\t")

