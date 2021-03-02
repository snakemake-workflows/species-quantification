rule kraken2_lr:
	input:
		fq = "results/mixed_lr/mixed_{p}_Sample{n}.fastq",
		db = "resources/kraken2-bacteria/bacterial-db"
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
		db = "resources/kraken2-bacteria/bacterial-db",
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
		sig = "results/sourmash/lr/sig/Scaled_{s}_mixed_sample{n}_{p}_k21.sig"
	log:
		"logs/sourmash-compute/lr/Scaled_{s}_Sample{n}_{p}.log"
	conda:
		"../envs/sourmash.yaml"
	params:
		"{s}"
	shell:
		"sourmash compute --scaled {params} {input} -o {output} -k=21"

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
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output}"

