configfile: "config.yaml"


rule all:
    input:
        "refs/chrY.fasta",
        "refs/hs_genome.fasta",
        expand("art/hum/Sample{art_hum}_1.fq", art_hum=range(1,config['number_of_samples']+1)),
        expand("art/bac/{bac_ref}_1.fq", bac_ref=config['bac_ref']),
        expand("nanosim/hum/{n}_aligned_reads.fastq", n=range(1,config['number_of_samples']+1)),
        expand("nanosim_training/bac_model/{bac_ref}/{bac_ref}_aligned_reads.pkl",bac_ref=config['bac_ref']),
        expand("nanosim/bac/{bac_ref}_aligned_reads.fastq", bac_ref=config['bac_ref']),
        expand("fractions/{bac_ref}_{p}_{k}.fastq", bac_ref=config['bac_ref'] ,p=config['p'], k={1,2}),
        expand("mixed/mixed_{p}_Sample{n}_{k}.fastq", p=config['p'], n=range(1,config['number_of_samples']+1), k={1,2}),
        expand("fractions/long_read/{bac_ref}_{p}.fastq", bac_ref=config['bac_ref'], p=config['p']),
        expand("mixed/long_read/mixed_{p}_Sample{n}.fastq", p=config['p'], n=range(1,config['number_of_samples']+1))
        
        #expand("{assembly}/{assembly}.fa", assembly=config['bac']),
        #expand("nanosim/bac/simulated_sample{sample3}_aligned_reads.fastq", sample3=range(config['number_of_samples'])),
        #"nanosim/bac/simulated_sample0_aligned_reads.fastq",
        #expand("InSilicoSeq/output/simulated{sample3}_R1.fastq", sample3 = range(1,config['number_of_samples']+1))

rule get_chromosome:
    output:
        "refs/chrY.fasta"
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="102",
        chromosome="Y"
    log:
        "logs/ensembl/get_seq.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.70.0/bio/reference/ensembl-sequence"

rule get_hs_genome:
    output:
        "refs/hs_genome.fasta"
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="102"
    log:
        "logs/ensembl/get_genome.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.72.0/bio/reference/ensembl-sequence"

'''
rule genomepy:
    output:
        "{assembly}.fa", 
        "{assembly}.fa.fai"
    log:
        "logs/genomepy_{assembly}.log"
    params:
        provider="ncbi",  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
    cache: True  # mark as eligible for between workflow caching
    wrapper:
        "0.70.0/bio/genomepy"

'''
 #human sequence simulation 
rule art_sim_hs:
    input:
        "refs/chrY.fasta"
    output:
        multiext("art/hum/Sample{art_hum}", "1.fq", "2.fq")
    log:
        "art/logs/{art_hum}.log"
    threads: 2
    params:
        length = config['short_read_len'],
        f_cov = config['fold_coverage']
    shell:
        "art_illumina -ss HS25 -i refs/chrY.fasta -p -l {params.length} -f {params.f_cov} -m 200 -s 10 -o" 
        " art/hum/Sample{wildcards.art_hum} --noALN"
        # "&& cat art/hum/Sample{wildcards.sample}.temp | gzip > art/hum/Sample{wildcards.sample} "
        # "&& rm art/hum/Sample{wildcards.sample}.temp "
    #bigger fold of read coverage kills the mixing process

rule art_sim_bac:
    input:
        "refs/{bac_ref}genomic.fna"
    output:
        multiext("art/bac/{bac_ref}", "1.fq", "2.fq")
    log:
        "art/logs/{bac_ref}.log"
    threads: 2
    params:
        length = config['short_read_len'],
        f_cov = config['fold_coverage']
    shell:
        "art_illumina -ss HS25 -i {input} -p -l {params.length} -f {params.f_cov} -m 200 -s 10 -o art/bac/{wildcards.bac_ref} --noALN"

rule nanosim_hs:
    input:
        ref = "refs/chrY.fasta",
        model = config['hum_pro'] 
    output:
        "nanosim/hum/{n}_aligned_error_profile",
        "nanosim/hum/{n}_aligned_reads.fastq",
        "nanosim/hum/{n}_unaligned_reads.fastq"      
    log:
        "logs/nanosim/{n}.log"
    threads: 4
    conda:
        "nanosim_env.yaml"
    params:
        nanosim = config['nanosim'],
        min_reads = config['min_long_read_len'], #usage of these causes a bug
        max_reads = config['max_long_read_len']
    shell:
        "python {params.nanosim}/simulator.py genome -rg {input.ref} -c {input.model}/training -b albacore --num_threads {threads}"
        " --fastq -o nanosim/hum/{wildcards.n} "

#bacterial sequence simulation

######
''' simulation with abundances

rule nanosim_bac_train:
    input:
       read = config['ecoli_read'] + "/sra_data.fastq",
       g_list = config['bac_list'] + "/nanosim_metagenome_list_for_training.tsv"
    output:
        "nanosim_training/bac_model/e_coli"      
    log:
        "logs/nanosim/nanosim_bac_train.log"
    threads: 2
    conda:
        "nanosim_env.yaml"
    params:
        nanosim = config['nanosim']
    shell:
        "python {params.nanosim}/read_analysis.py metagenome -i {input.read} -gl {input.g_list} -o {output}"

rule nanosim_bac:
    input:
        genome_list = config['bac_list'] + "/nanosim_metagenome_list_for_sim.tsv",
        abun = config['bac_list'] + "/abundance_for_simulation_multi_sample.tsv",
        dna_list = config['bac_list'] + "/dna_type_list.tsv",
        model = model_dir2 + "/e_coli_model_profile"
    output:
        "nanosim/bac/simulated_sample0_aligned_error_profile",
        "nanosim/bac/simulated_sample0_aligned_reads.fastq",
        "nanosim/bac/simulated_sample0_unaligned_reads.fastq"      
    log:
        "logs/nanosim-bac/simulated.log"
    threads: 4
    conda:
        "nanosim_env.yaml"
    params:
        nanosim = config['nanosim']
    shell:
        "python {params.nanosim}/simulator.py metagenome -gl {input.genome_list} -a {input.abun} -dl {input.dna_list} -c {model_dir2}/e_coli -o nanosim/bac/simulated --num_threads 4 --fastq -b albacore"

#before this rule the genome files of each bacterial species should be concatenated to create one fasta file, cat *.fna > genomes.fasta  

rule iss_bac:
    input:
        genomes = "refs/genomes.fasta",
        #abun = "InSilicoSeq/abundance.txt" not use this time
    output:
        "InSilicoSeq/output/simulated{sample3}_R1.fastq",
        "InSilicoSeq/output/simulated{sample3}_R2.fastq"
    log: 
        "logs/insilicoseq/{sample3}.log"
    threads: 4
    conda: "test.yaml"
    shell:
        "iss generate -g {input.genomes} -m HiSeq -o InSilicoSeq/output/simulated{wildcards.sample3} --cpus {threads} --n_reads 1000"
'''
#######

rule nanosim_bac_train:
    input:
       read = config['ecoli_read'] + "/sra_data.fastq",
       r = "refs/{bac_ref}_genomic.fna"
    output:     
        "nanosim_training/bac_model/{bac_ref}/{bac_ref}_aligned_reads.pkl"
    log:
        "logs/nanosim/{bac_ref}_train.log"
    threads: 4
    conda:
        "nanosim_env.yaml"
    params:
        nanosim = config['nanosim']
    shell:
        "python {params.nanosim}/read_analysis.py genome -i {input.read} -rg {input.r} --num_threads {threads}"
        " -o nanosim_training/bac_model/{wildcards.bac_ref}/{wildcards.bac_ref}"


#question: dna_type circular or linear????
#error: "Do not choose circular if there is more than one chromosome in the genome!" when the option is set to circular

rule nanosim_bac_sim:
    input:
        r = "refs/{bac_ref}_genomic.fna",
        model = "nanosim_training/bac_model/{bac_ref}/{bac_ref}_aligned_reads.pkl"
    output:
        "nanosim/bac/{bac_ref}_aligned_error_profile",
        "nanosim/bac/{bac_ref}_aligned_reads.fastq",
        "nanosim/bac/{bac_ref}_unaligned_reads.fastq"      
    log:
        "logs/nanosim/{bac_ref}.log"
    threads: 4
    conda:
        "nanosim_env.yaml"
    params:
        nanosim = config['nanosim'],
        min_reads = config['min_long_read_len'],
        max_reads = config['max_long_read_len']
    shell:
        "python {params.nanosim}/simulator.py genome -rg {input.r} -c nanosim_training/bac_model/{wildcards.bac_ref}/{wildcards.bac_ref}"
        " -b albacore --num_threads {threads} --fastq -o nanosim/bac/{wildcards.bac_ref}"

#fractionation and concatenation
rule fraction_sr:
    input:
        bac_fq1 = "art/bac/{bac_ref}_1.fq",
        bac_fq2 = "art/bac/{bac_ref}_2.fq"
    output:
        out_fq1="fractions/{bac_ref}_{p}_1.fastq",
        out_fq2="fractions/{bac_ref}_{p}_2.fastq"
    log:
        "logs/seqtk/{bac_ref}_{p}.log"
    params:
        fraction = "{p}"       
    conda:
        "concat_env.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.bac_fq2} {params} > {output.out_fq2}"

rule concat_fractions_sr:
    input:
        bac_fq1 = expand("fractions/{bac_ref}_{{p}}_1.fastq", bac_ref=config['bac_ref']),
        bac_fq2 = expand("fractions/{bac_ref}_{{p}}_2.fastq", bac_ref=config['bac_ref']),
        hum_fq1 = "art/hum/Sample{n}_1.fq",
        hum_fq2 = "art/hum/Sample{n}_2.fq"
    output:
        out_fq1 = "mixed/mixed_{p}_Sample{n}_1.fastq",
        out_fq2 = "mixed/mixed_{p}_Sample{n}_2.fastq" 
    log:
        "logs/concat/mixed_{p}_Sample{n}.log"
    shell:
        "cat {input.hum_fq1} {input.bac_fq1} > {output.out_fq1}; "
        "cat {input.hum_fq2} {input.bac_fq2} > {output.out_fq2}"

rule fraction_lr:
    input:
        bac_fq = "nanosim/bac/{bac_ref}_unaligned_reads.fastq",
    output:
        out_fq="fractions/long_read/{bac_ref}_{p}.fastq",
    log:
        "logs/seqtk/lr/{bac_ref}_{p}.log"
    params:
        fraction = "{p}"       
    conda:
        "concat_env.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq} {params} > {output.out_fq}"

rule concat_fractions_lr:
    input:
        bac_fq = expand("fractions/long_read/{bac_ref}_{{p}}.fastq", bac_ref=config['bac_ref'], p=config['p']),
        hum_fq = "nanosim/hum/{n}_unaligned_reads.fastq",
    output:
        out_fq = "mixed/long_read/mixed_{p}_Sample{n}.fastq",
    log:
        "logs/concat/lr/mixed_{p}_Sample{n}.log"
    shell:
        "cat {input.hum_fq} {input.bac_fq} > {output.out_fq}"
