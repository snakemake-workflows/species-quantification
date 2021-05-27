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
    threads: 30
    params:
        med_len=config["med_long_read_len"],  #add -med for median length
        nreads=config["long_nreads"],
    shell:
        "resources/NanoSim/src/simulator.py genome -rg {input.ref} -c {input.model}/training -b guppy --num_threads {threads}"
        " --fastq -o results/nanosim/hum/{wildcards.n} -n {params.nreads} 2> {log}"


rule nanosim_bac_train:
    input:
        read=config["ecoli_read"] + "/sra_data.fastq",
        r="resources/bac_refs/GCF_000008865.2_ASM886v2_chr_genomic.fna", #ecoli genome
    output:
        "results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_aligned_reads.pkl",
    log:
        "logs/nanosim_train/ecoli_train.log",
    threads: 30
    shell:
        "resources/NanoSim/src/read_analysis.py genome -i {input.read} -rg {input.r} --num_threads {threads}"
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
    params:
        med_len=config["med_long_read_len"],
	nreads=config["long_nreads"]
    shell:
        "resources/NanoSim/src/simulator.py genome -rg {input.r} -c results/nanosim_train/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2"
        " -b guppy --num_threads {threads} --fastq -o results/nanosim/bac/{wildcards.bac_ref} -dna_type circular -n {params.nreads} 2> {log}"
        #add -med for median length
