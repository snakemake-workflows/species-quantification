
rule kraken2_sr:
	input:
		fq1 = "results/mixed_sr/mixed_{p}_Sample{n}_1.fastq",
		fq2 = "results/mixed_sr/mixed_{p}_Sample{n}_2.fastq",
		db = "resources/kraken2-bacteria/standard_db"
	output:
		rep = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}",
		kraken = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}.kraken"
	log:
		"logs/kraken2/sr/sb/Sample{n}_{p}.log"
	threads: 20
	conda:
		"../envs/kraken2.yaml"
	shell:
		"kraken2 --use-names --threads {threads} --db {input.db} --fastq-input --report "
		" {output.rep} --paired {input.fq1} {input.fq2} > {output.kraken} 2> {log}"

rule bracken_sr:
	input:
		db = "resources/kraken2-bacteria/standard_db",
		rep = "results/kraken2/sr/sb/evol1_Sample{n}_fraction{p}"
	output:
		bracken = "results/bracken/sr/sb/evol1_Sample{n}_fraction{p}.bracken"
	log:
		"logs/bracken/sr/sb/Sample{n}_{p}.log"
	threads: 2
	conda:
		"../envs/kraken2.yaml"
	shell: 
		"bracken -d {input.db} -i {input.rep} -l S -o {output} 2> {log}"		

rule sourmash_comp_sr:
	input:
		fq = "results/mixed_sr/mixed_{p}_Sample{n}_{e}.fastq"
	output:
		sig = "results/sourmash/sr/sig/Scaled_{s}_mixed_sample{n}_{p}_k51_R{e}.sig"
	log:
		"logs/sourmash-compute/sr/Scaled_{s}_Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	params:
		"{s}"
	shell:
		"sourmash compute --scaled {params} {input} -o {output} -k=51 2> {log}"


rule sourmash_lca_sr:
	input:
		sig = "results/sourmash/sr/sig/Scaled_{s}_mixed_sample{n}_{p}_k51_R{e}.sig",
		db = "resources/sourmash/genbank-k51.lca.json"
	output:
		sum = "results/sourmash/sr/lca-class/Scaled_{s}_mixed_sample{n}_{p}_R{e}.csv"
	log:
		"logs/sourmash-lca/sr/Scaled_{s}_Sample{n}_{p}_{e}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"



