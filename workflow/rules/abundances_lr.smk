rule kraken2_lr:
	input:
		fq = "results/mixed_lr/mixed_{p}_Sample{n}.fastq",
		db = "resources/minikraken2_v2_8GB_201904_UPDATE"
	output:
		rep = "results/kraken2/lr/evol1_Sample{n}_fraction{p}",
		kraken = "results/kraken2/lr/evol1_Sample{n}_fraction{p}.kraken"
	log:
		"logs/kraken2/lr/Sample{n}_{p}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input {input.fq} "
		" --report {output.rep} > {output.kraken}"

rule bracken_lr:
	input:
		db = "resources/minikraken2_v2_8GB_201904_UPDATE",
		rep = "results/kraken2/lr/evol1_Sample{n}_fraction{p}"
	output:
		bracken = "results/bracken/lr/evol1_Sample{n}_fraction{p}.bracken"
	log:
		"logs/bracken/lr/Sample{n}_{p}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output}"		

rule sourmash_comp_lr:
	input:
		fq = "results/mixed_lr/mixed_{p}_Sample{n}.fastq"
	output:
		sig = "results/sourmash/lr/sig/mixed_sample{n}_{p}_k51.sig"
	log:
		"logs/sourmash-compute/lr/Sample{n}_{p}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash compute --scaled 10000 {input} -o {output} -k=51"

rule sourmash_lca_lr:
	input:
		sig = "results/sourmash/lr/sig/mixed_sample{n}_{p}_k51.sig",
		db = "resources/sourmash/genbank-k51.lca.json"
	output:
		sum = "results/sourmash/lr/lca-class/mixed_sample{n}_{p}.csv"
	log:
		"logs/sourmash-lca/lr/Sample{n}_{p}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} > {output}"

