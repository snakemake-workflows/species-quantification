rule kraken_build:
	output:
		db = "resources/kraken2-db",
                mock = "resources/kraken2-db/mock.txt"
	log:
		"logs/kraken2-build/kraken_db.log"
	threads: 20
	params:
		read_len = 100, #default value, please note that Bracken wasn't necessarily designed to run on nanopore data.
                dbtype = config["dbtype"]
	conda:
		"../envs/kraken2.yaml"
	cache: True
	shell:
		"kraken2-build --download-taxonomy --skip-maps --db {output.db} &&" #only required to download test database
		"kraken2-build {params.dbtype} --threads {threads} --db {output.db} && kraken2-build --clean "
		"bracken-build -d {output.db} -l {params.read_len} && touch {output.mock}"
		
rule kraken2:
	input:
		fq = get_fastq_input,
		db = "resources/kraken2-db",
	output:
		rep = "results/kraken2/{sample}/evol1_{sample}_{unit}",
		kraken = "results/kraken2/{sample}/evol1_{sample}_{unit}.kraken"
	log:
		"logs/kraken2/{sample}_{unit}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {input.fq} "
		" --report {output.rep} > {output.kraken} 2> {log}"

rule bracken:
	input:
		db = "resources/kraken2-db",
		rep = "results/kraken2/{sample}/evol1_{sample}_{unit}"
	output:
		bracken = "results/bracken/{sample}/evol1_{sample}_{unit}.bracken"
	log:
		"logs/bracken/{sample}_{unit}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"		
 
rule sourmash_comp:
	input:
		get_fastq_input,
	output:
		sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k21.sig"
	log:
		"logs/sourmash-compute/{sample}_{unit}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash compute --scaled 2000 {input} -o {output} -k=21 2> {log}"

rule sourmash_lca:
	input:
		sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k21.sig",
		db = "resources/sourmash/genbank-k21.lca.json"
	output:
		sum = "results/sourmash/lca-class/{sample}/{sample}_{unit}_k21.csv"
	log:
		"logs/sourmash-lca/{sample}_{unit}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"

