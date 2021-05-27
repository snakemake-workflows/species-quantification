rule kraken2:
	input:
		fq = "resources/data/{sample}.fastq",
		db = "resources/kraken2-bacteria/standard_db_LR"
	output:
		rep = "results/kraken2/evol1_{sample}",
		kraken = "results/kraken2/evol1_{sample}.kraken"
	log:
		"logs/kraken2/{sample}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {input.fq} "
		" --report {output.rep} > {output.kraken} 2> {log}"

rule bracken:
	input:
		db = "resources/kraken2-bacteria/standard_db_LR",
		rep = "results/kraken2/evol1_{sample}"
	output:
		bracken = "results/bracken/evol1_{sample}.bracken"
	log:
		"logs/bracken/{sample}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"		
 
rule sourmash_comp:
	input:
		fq = "resources/data/{sample}.fastq"
	output:
		sig = "results/sourmash/sig/{sample}_k21.sig"
	log:
		"logs/sourmash-compute/{sample}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash compute --scaled 2000 {input} -o {output} -k=21 2> {log}"

rule sourmash_lca:
	input:
		sig = "results/sourmash/sig/{sample}_k21.sig",
		db = "resources/sourmash/genbank-k21.lca.json"
	output:
		sum = "results/sourmash/lca-class/{sample}_k21.csv"
	log:
		"logs/sourmash-lca/{sample}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"

