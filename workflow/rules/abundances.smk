
rule kraken2_build:
	output:
		db = directory("results/kraken2-db"),
                mock = "results/kraken2-db/mock.txt"
	log:
		"logs/kraken2-build/kraken_db.log"
	threads: 20
	params:
		read_len = 100, #default value, please note that Bracken wasn't necessarily designed to run on nanopore data.
                dbtype = config["dbtype"]
	conda:
		"../envs/kraken2.yaml"
	priority: 1
	cache: True
	shell:
		"kraken2-build --download-taxonomy --db {output.db} && " #only required to download test database
		"kraken2-build {params.dbtype} --threads {threads} --db {output.db} && kraken2-build --build --db {output.db} --threads {threads} && "
		"bracken-build -d {output.db} && touch {output.mock} && "
		"kraken2-build --clean --db {output.db}"
		
rule kraken2:
	input:
		fq = get_fastq_input,
		db = "results/kraken2-db",
	output:
		rep = "results/kraken2/{sample}/evol1_{sample}",
		kraken = "results/kraken2/{sample}/evol1_{sample}.kraken"
	log:
		"logs/kraken2/{sample}.log"
	threads: 20
	params:
		 paired = "--paired" if config["paired"] == True else ""
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {params.paired} {input.fq} "
		" --report {output.rep} > {output.kraken} 2> {log}"

rule bracken:
	input:
		db = "results/kraken2-db",
		rep = "results/kraken2/{sample}/evol1_{sample}"
	output:
		bracken = "results/bracken/{sample}/evol1_{sample}.bracken"
	log:
		"logs/bracken/{sample}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"
