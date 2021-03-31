rule get_chromosome:
    output:
        "results/refs/chrY.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="102",
        chromosome="Y",
    log:
        "logs/ensembl/get_seq.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.70.0/bio/reference/ensembl-sequence"


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

