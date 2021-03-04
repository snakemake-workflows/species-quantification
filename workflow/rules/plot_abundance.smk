rule compare_results:
    input:
        sourmash = expand("results/sourmash/sr/lca-class/Scaled_2000/Scaled_2000_mixed_sample{n}_{p}_R1.csv",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"]),
        kraken2 = expand("results/kraken2/sr/bacterial-db/evol1_Sample{n}_fraction{p}_bracken_species",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"])
    output:
        "results/final_abundance/scatter_plot/sr/final_abundance_all_samples_coord_fixed.pdf"
    script:
        "../scripts/abundance_plot.R"

        