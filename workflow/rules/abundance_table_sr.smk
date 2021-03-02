
rule sourmash_table_sr:
    input:
        expand("results/sourmash/sr/lca-class/mixed_sample{n}_{p}_{e}.csv",
        n=range(1, config["number_of_samples"] + 1), 
        p=config["p"], e={"R1", "R2"})
    output:
        "results/final_abundance/sourmash_table.txt"
    log:
        notebook="logs/notebooks/abundance_table_sr.pynb"
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "../notebooks/abundance_table.r.ipynb"