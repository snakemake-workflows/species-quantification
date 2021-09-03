
rule kraken2_build:
    output:
        db=directory("results/kraken2-db"),
        mock="results/kraken2-db/mock.txt",
    params:
        read_len=100,  #default value, please note that Bracken wasn't necessarily designed to run on nanopore data.
        dbtype=config["dbtype"],
    log:
        "logs/kraken2-build/kraken_db.log",
    conda:
        "../envs/kraken2.yaml"
    threads: 20
    priority: 1
    cache: True
    shell:
        "kraken2-build --download-taxonomy --skip-maps --db {output.db} && "
        "kraken2-build {params.dbtype} --threads {threads} --db {output.db} && kraken2-build --build --db {output.db} --threads {threads} && "
        "bracken-build -d {output.db} && touch {output.mock} && "
        "kraken2-build --clean --db {output.db}"


rule kraken2:
    input:
        fq=get_fastq_input,
        db="results/kraken2-db",
    output:
        rep="results/kraken2/{sample}/evol1_{sample}",
        kraken="results/kraken2/{sample}/evol1_{sample}.kraken",
    params:
        paired="--paired" if config["paired"] == True else "",
    log:
        "logs/kraken2/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    threads: 20
    shell:
        "kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {params.paired} {input.fq} "
        " --report {output.rep} > {output.kraken} 2> {log}"


rule bracken:
    input:
        db="results/kraken2-db",
        rep="results/kraken2/{sample}/evol1_{sample}",
    output:
        bracken="results/bracken/{sample}/evol1_{sample}.bracken",
    log:
        "logs/bracken/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    threads: 2
    shell:
        "bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"
