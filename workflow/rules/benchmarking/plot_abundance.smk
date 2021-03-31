rule compare_results_sr:
    input:
        sourmash = expand("results/sourmash/sr/lca-class/Scaled_2000_mixed_sample{n}_{p}_R1.csv",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"]),
        kraken2 = expand("results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"]),
        bracken = expand("results/bracken/sr/sb/evol1_Sample{n}_fraction{p}.bracken",
        n = range(1, config["number_of_samples"] + 1),
        p = config["p"])
    output:
        "results/final_abundance/scatter_plot/sr/sr_final_abundance_all_samples_coord_fixed.pdf",
        "results/final_abundance/scatter_plot/sr/sr_final_abundance_all_samples.csv"
    script:
        "../scripts/sr_abundance_plot.R"

rule compare_results_lr:
    input:
        sourmash = expand("results/sourmash/lr/lca-class/Scaled_2000_mixed_sample{n}_{p}_k21.csv",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"]),
        kraken2 = expand("results/kraken2/lr/sb/evol1_Sample{n}_fraction{p}",
        n = range(1, config["number_of_samples"] + 1),
        p =  config["p"]),
        bracken = expand("results/bracken/lr/sb/evol1_Sample{n}_fraction{p}.bracken",
        n = range(1, config["number_of_samples"] + 1),
        p = config["p"])
    output:
        "results/final_abundance/scatter_plot/lr/lr_final_abundance_all_samples_coord_fixed.pdf",
        "results/final_abundance/scatter_plot/lr/lr_final_abundance_all_samples.csv"
    script:
        "../scripts/lr_abundance_plot.R"