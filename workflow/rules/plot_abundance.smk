rule compare_results:
    input:
        sourmash = "results/sourmash/sr/lca-class/Scaled_2000/Scaled_2000_mixed_sample{n}_{p}_R1.csv",
        kraken2 = "results/kraken2/sr/bacterial-db/evol1_Sample{n}_fraction{p}_bracken_species"
    params:
        n_samples = "{n}",
        fraction = "{p}"
    output:
        "results/final_abundance/second_try/sr/Sample{n}_fraction{p}_total_abundance.pdf"
    script:
        "../scripts/abundance_plot.R"