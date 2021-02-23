rule kraken2_sr:
	input:
		fq1 = "results/mixed_sr/mixed_{p}_Sample{n}_1.fastq",
		fq2 = "results/mixed_sr/mixed_{p}_Sample{n}_2.fastq",
		db = "resources/minikraken2_v2_8GB_201904_UPDATE"
	output:
		rep = "results/kraken2/sr/evol1_Sample{n}_fraction{p}",
		kraken = "results/kraken2/sr/evol1_Sample{n}_fraction{p}.kraken"
	log:
		"logs/kraken2/sr/Sample{n}_{p}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input --report "
		" {output.rep} --paired {input.fq1} {input.fq2} > {output.kraken}"

rule bracken_sr:
	input:
		db = "resources/minikraken2_v2_8GB_201904_UPDATE",
		rep = "results/kraken2/sr/evol1_Sample{n}_fraction{p}"
	output:
		bracken = "results/bracken/sr/evol1_Sample{n}_fraction{p}.bracken"
	log:
		"logs/bracken/sr/Sample{n}_{p}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output}"		

rule sourmash_comp_sr:
	input:
		fq = "results/mixed_sr/mixed_{p}_Sample{n}_{e}.fastq"
	output:
		sig = "results/sourmash/sr/sig/mixed_sample{n}_{p}_k51_R{e}.sig"
	log:
		"logs/sourmash-compute/sr/Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash compute --scaled 10000 {input} -o {output} -k=51"


rule sourmash_lca_sr:
	input:
		sig = "results/sourmash/sr/sig/mixed_sample{n}_{p}_k51_R{e}.sig",
		db = "resources/sourmash/genbank-k51.lca.json"
	output:
		sum = "results/sourmash/sr/lca-class/mixed_sample{n}_{p}_R{e}.csv"
	log:
		"logs/sourmash-lca/sr/Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} > {output}"
