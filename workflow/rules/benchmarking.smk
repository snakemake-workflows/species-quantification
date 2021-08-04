rule get_hs_genome:
    output:
        "results/refs/hs_genome.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="102",
    log:
        "logs/ensembl/get_genome.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.72.0/bio/reference/ensembl-sequence"


rule art_sim_hs:
    input:
        "results/refs/hs_genome.fasta",
    output:
        multiext("results/art/hum/Sample{art_hum}", "1.fq", "2.fq"),
    log:
        "logs/art/hum/{art_hum}.log",
    threads: 2
    conda:
        "../envs/art.yaml"
    params:
        length=config["short_read_len"],
        f_cov=config["fold_coverage"],
    shell:
        "art_illumina -ss HS25 -i results/refs/hs_genome.fasta -p -l {params.length} -f {params.f_cov} -m 200 -s 10 -o"
        " results/art/hum/Sample{wildcards.art_hum} --noALN 2> {log}"


rule art_sim_bac:
    input:
        "resources/bac_refs/{bac_ref}genomic.fna",
    output:
        multiext("results/art/bac/{bac_ref}", "1.fq", "2.fq"),
    log:
        "logs/art/bac/{bac_ref}.log",
    threads: 2
    conda:
        "../envs/art.yaml"
    params:
        length=config["short_read_len"],
        f_cov=config["fold_coverage"],
    shell:
        "art_illumina -ss HS25 -i {input} -p -l {params.length} -f {params.f_cov} -m 200 -s 10 -o results/art/bac/{wildcards.bac_ref} --noALN 2> {log}"

rule nanosim_hs:
    input:
        ref="results/refs/hs_genome.fasta",
        model="resources/human_NA12878_DNA_FAB49712_guppy"
    output:
        "results/nanosim/hum/{n}_aligned_error_profile",
        "results/nanosim/hum/{n}_aligned_reads.fastq",
        "results/nanosim/hum/{n}_unaligned_reads.fastq",
    log:
        "logs/nanosim/hum/{n}.log",
    threads: 30
    conda:
        "../envs/nanosim.yaml"
    params:
        med_len=config["med_long_read_len"],  #add -med for median length
        nreads=config["long_nreads"],
    shell:
        "simulator.py genome -rg {input.ref} -c {input.model}/training -b guppy --num_threads {threads}"
        " --fastq -o results/nanosim/hum/{wildcards.n} -n {params.nreads} 2> {log}"


rule nanosim_bac_train:
    input:
        read="resources/ecoli_read/sra_data.fastq",
        r="resources/bac_refs/GCF_000008865.2_ASM886v2_chr_genomic.fna", #ecoli genome
    output:
        "results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_aligned_reads.pkl",
    log:
        "logs/nanosim_train/ecoli_train.log",
    threads: 30
    conda:	
        "../envs/nanosim.yaml"
    shell:
        "read_analysis.py genome -i {input.read} -rg {input.r} --num_threads {threads}"
        " -o results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2 2> {log}"


rule nanosim_bac_sim:
    input:
        r="resources/bac_refs/{bac_ref}_genomic.fna",
        model="results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_aligned_reads.pkl",
    output:
        "results/nanosim/bac/{bac_ref}_aligned_error_profile",
        "results/nanosim/bac/{bac_ref}_aligned_reads.fastq",
        "results/nanosim/bac/{bac_ref}_unaligned_reads.fastq",
    log:
        "logs/nanosim/bac/{bac_ref}.log",
    threads: 30
    conda:
        "../envs/nanosim.yaml"
    params:
        med_len=config["med_long_read_len"],
	nreads=config["long_nreads"]
    shell:
        "simulator.py genome -rg {input.r} -c results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2"
        " -b guppy --num_threads {threads} --fastq -o results/nanosim/bac/{wildcards.bac_ref} -dna_type circular -n {params.nreads} 2> {log}"
        #add -med for median length

# fractionation and concatenation
rule fraction_sr:
    input:
        bac_fq1="results/art/bac/{bac_ref}_1.fq",
        bac_fq2="results/art/bac/{bac_ref}_2.fq",
    output:
        out_fq1="results/fractions/short_reads/{bac_ref}_{p}_1.fastq",
        out_fq2="results/fractions/short_reads/{bac_ref}_{p}_2.fastq",
    log:
        "logs/seqtk/short_read/{bac_ref}_{p}.log",
    params:
        fraction="{p}",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.bac_fq2} {params} > {output.out_fq2} 2> {log}"


rule concat_fractions_sr:
    input:
        bac_fq1=expand(
            "results/fractions/short_reads/{bac_ref}_{{p}}_1.fastq", bac_ref=bac.bacteria
        ),
        bac_fq2=expand(
            "results/fractions/short_reads/{bac_ref}_{{p}}_2.fastq", bac_ref=bac.bacteria
        ),
        hum_fq1="results/art/hum/Sample{n}_1.fq",
        hum_fq2="results/art/hum/Sample{n}_2.fq",
    output:
        out_fq1="results/mixed_sr/mixed_{p}_Sample{n}_1.fastq",
        out_fq2="results/mixed_sr/mixed_{p}_Sample{n}_2.fastq",
    log:
        "logs/concatenation/sr/mixed_{p}_Sample{n}.log",
    shell:
        "cat {input.hum_fq1} {input.bac_fq1} > {output.out_fq1}; "
        "cat {input.hum_fq2} {input.bac_fq2} > {output.out_fq2}"


rule fraction_lr:
    input:
        bac_fq="results/nanosim/bac/{bac_ref}_aligned_reads.fastq",
    output:
        out_fq="results/fractions/long_reads/{bac_ref}_{p}.fastq",
    log:
        "logs/seqtk/long_read/{bac_ref}_{p}.log",
    params:
        fraction="{p}",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq} {params} > {output.out_fq} 2> {log}"


rule concat_fractions_lr:
    input:
        bac_fq=expand(
            "results/fractions/long_reads/{bac_ref}_{{p}}.fastq",
            bac_ref=bac.bacteria,
            p=config["p"],
        ),
        hum_fq="results/nanosim/hum/{n}_aligned_reads.fastq",
    output:
        out_fq="results/mixed_lr/mixed_{p}_Sample{n}.fastq",
    log:
        "logs/concatenation/sr/mixed_{p}_Sample{n}.log",
    shell:
        "cat {input.hum_fq} {input.bac_fq} > {output.out_fq}"

#rule kraken_build:
#        output:
#                "resources/kraken2-db/standard_db"
#        log:
#                "logs/kraken2-build/kraken_db.log"
#        threads: 20
#        params:
#                read_len = 100 #default value, please note that Bracken wasn't necessarily designed to run on nanopore data.
#        conda:
#                "../envs/kraken2.yaml"
#        cache: True
#        shell:
#                "kraken2-build --standard --threads {threads} --db {output} &&"
#                "bracken-build -d {output} -l {params.read_len}"
#
#
rule kraken2_sr:
	input:
		fq1 = "results/mixed_sr/mixed_{p}_Sample{n}_1.fastq",
		fq2 = "results/mixed_sr/mixed_{p}_Sample{n}_2.fastq",
		db = "resources/kraken2-db/standard_db"
	output:
		rep = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}",
		kraken = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}.kraken"
	log:
		"logs/kraken2/sr/sb/Sample{n}_{p}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input --report "
		" {output.rep} --paired {input.fq1} {input.fq2} > {output.kraken} 2> {log}"

rule bracken_sr:
	input:
		db = "resources/kraken2-db/standard_db",
		rep = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}"
	output:
		bracken = "results/bracken/sr/sb/evol1_Sample{n}_fraction{p}.bracken"
	log:
		"logs/bracken/sr/sb/Sample{n}_{p}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"		

rule sourmash_comp_sr:
	input:
		fq = "results/mixed_sr/mixed_{p}_Sample{n}_{e}.fastq"
	output:
		sig = "results/sourmash/sr/sig/Scaled_{s}_mixed_sample{n}_{p}_k51_R{e}.sig"
	log:
		"logs/sourmash-compute/sr/Scaled_{s}_Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	params:
		"{s}"
	shell:
		"sourmash compute --scaled {params} {input} -o {output} -k=51 2> {log}"


rule sourmash_lca_sr:
	input:
		sig = "results/sourmash/sr/sig/Scaled_{s}_mixed_sample{n}_{p}_k51_R{e}.sig",
		db = "resources/sourmash/genbank-k51.lca.json"
	output:
		sum = "results/sourmash/sr/lca-class/Scaled_{s}_mixed_sample{n}_{p}_R{e}.csv"
	log:
		"logs/sourmash-lca/sr/Scaled_{s}_Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"



rule kraken2_lr:
	input:
		fq = "results/mixed_lr/mixed_{p}_Sample{n}.fastq",
		db = "resources/kraken2-db/standard_db"
	output:
		rep = "results/kraken2/lr/sb/evol1_Sample{n}_fraction{p}",
		kraken = "results/kraken2/lr/sb/evol1_Sample{n}_fraction{p}.kraken"
	log:
		"logs/kraken2/lr/sb/Sample{n}_{p}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {input.fq} "
		" --report {output.rep} > {output.kraken} 2> {log}"

rule bracken_lr:
	input:
		db = "resources/kraken2-db/standard_db",
		rep = "results/kraken2/lr/sb/evol1_Sample{n}_fraction{p}"
	output:
		bracken = "results/bracken/lr/sb/evol1_Sample{n}_fraction{p}.bracken"
	log:
		"logs/bracken/lr/sb/Sample{n}_{p}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"		
 
rule sourmash_comp_lr:
	input:
		fq = "results/mixed_lr/mixed_{p}_Sample{n}.fastq"
	output:
		sig = "results/sourmash/lr/sig/Scaled_{s}_mixed_sample{n}_{p}_k21.sig"
	log:
		"logs/sourmash-compute/lr/Scaled_{s}_Sample{n}_{p}.log"
	conda:
		"../envs/sourmash.yaml"
	params:
		"{s}"
	shell:
		"sourmash compute --scaled {params} {input} -o {output} -k=21 2> {log}"

rule sourmash_lca_lr:
	input:
		sig = "results/sourmash/lr/sig/Scaled_{s}_mixed_sample{n}_{p}_k21.sig",
		db = "resources/sourmash/genbank-k21.lca.json"
	output:
		sum = "results/sourmash/lr/lca-class/Scaled_{s}_mixed_sample{n}_{p}_k21.csv"
	log:
		"logs/sourmash-lca/lr/Scaled_{s}_Sample{n}_{p}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"

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
