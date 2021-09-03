import pandas as pd

# read sample sheet
samples = pd.read_csv(config["samples"], sep="\t").set_index(
    ["sample_name"], drop=False
)

# only if benchmarking is set to true
bacteria = pd.read_csv(config["bacteria"], sep="\t").set_index(["bacteria"], drop=False)

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    sample = samples.loc[wildcards.sample]

    if pd.isna(sample["fq2"]):
        return [sample["fq1"]]
    else:
        return [sample["fq1"], sample["fq2"]]


# input function to retrieve bacteria, only required when benchmarking set to TRUE
def get_bacteria_input(wildcards):
    files = bacteria.loc[wildcards.bac_ref, "fasta"]
    return [files]
