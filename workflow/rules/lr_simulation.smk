rule nanosim_hs:
    input:
        ref="results/refs/hs_genome.fasta",
        model=config["hum_pro"],
    output:
        "results/nanosim/hum/{n}_aligned_error_profile",
        "results/nanosim/hum/{n}_aligned_reads.fastq",
        "results/nanosim/hum/{n}_unaligned_reads.fastq",
    log:
        "logs/nanosim/hum/{n}.log",
    threads: 4
    conda:
        "../envs/nanosim.yaml"
    params:
        med_len=config["med_long_read_len"],  #add -med for median length
        nreads=config["long_nreads"],
    shell:
        "simulator.py genome -rg {input.ref} -c {input.model}/training -b albacore --num_threads {threads}"
        " --fastq -o results/nanosim/hum/{wildcards.n} -n {params.nreads} "


rule nanosim_bac_train:
    input:
        read=config["ecoli_read"] + "/sra_data.fastq",
        r="resources/bac_refs/{bac_ref}_genomic.fna",
    output:
        "results/nanosim_train/{bac_ref}/{bac_ref}_aligned_reads.pkl",
    log:
        "logs/nanosim_train/{bac_ref}_train.log",
    threads: 10
    params:
        nanosim=config["nanosim"],  #conda installation does not work for read_analysis.py
    conda:
        "../envs/nanosim_train.yaml"
    shell:
        "python {params.nanosim}/read_analysis.py genome -i {input.read} -rg {input.r} --num_threads {threads}"
        " -o results/nanosim_train/{wildcards.bac_ref}/{wildcards.bac_ref}"


# question: dna_type circular or linear????
# error: "Do not choose circular if there is more than one chromosome in the genome!" when the option is set to circular


rule nanosim_bac_sim:
    input:
        r="resources/bac_refs/{bac_ref}_genomic.fna",
        model="results/nanosim_train/{bac_ref}/{bac_ref}_aligned_reads.pkl",
    output:
        "results/nanosim/bac/{bac_ref}_aligned_error_profile",
        "results/nanosim/bac/{bac_ref}_aligned_reads.fastq",
        "results/nanosim/bac/{bac_ref}_unaligned_reads.fastq",
    log:
        "logs/nanosim/bac/{bac_ref}.log",
    threads: 10
    conda:
        "../envs/nanosim.yaml"
    params:
        med_len=config["med_long_read_len"],
        nreads=config["long_nreads"],
    shell:
        "simulator.py genome -rg {input.r} -c results/nanosim_train/{wildcards.bac_ref}/{wildcards.bac_ref}"
        " -b albacore --num_threads {threads} --fastq -o results/nanosim/bac/{wildcards.bac_ref} -n {params.nreads}"
        #add -med for median length
