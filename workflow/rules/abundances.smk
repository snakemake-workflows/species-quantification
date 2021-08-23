
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
	priority: 2
	cache: True
	shell:
		"kraken2-build --download-taxonomy --skip-maps --db {output.db} && " #only required to download test database
		"kraken2-build {params.dbtype} --threads {threads} --db {output.db} && kraken2-build --build --db {output.db} --threads {threads} && "
		"bracken-build -d {output.db} && touch {output.mock}"
		#kraken2-build --clean --db {output.db
		
rule kraken2:
	input:
		fq = get_fastq_input,
		db = "results/kraken2-db",
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
		db = "results/kraken2-db",
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

rule sourmash_lca_db_k21:
	output:
		expand("results/sourmash_lca_db/{db}", db = config["sourmash_lca_name_k21"])
	params:
                name = config["sourmash_lca_name_k21"],
                link = config["sourmash_lca_link_k21"]
	priority: 1
	log:
		"logs/sourmash_lca_db/wget.log"
	shell:
		"cd results/sourmash_lca_db && wget {params.link} -O {params.name}.gz && gunzip {params.name}.gz"

rule sourmash_comp_k21:
	input:
		get_fastq_input,
	output:
		sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k21.sig"
	log:
		"logs/sourmash-compute-k21/{sample}_{unit}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash compute --scaled 2000 {input} -o {output} -k=21 2> {log}"

rule sourmash_lca_k21:
	input:
		sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k21.sig",
		db = expand("results/sourmash_lca_db/{db}", db = config["sourmash_lca_name_k21"]),
	output:
		sum = "results/sourmash/lca-class/{sample}/{sample}_{unit}_k21.csv"
	log:
		"logs/sourmash-lca-k21/{sample}_{unit}.log"
	conda:
		"../envs/sourmash.yaml"
	shell:
		"sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"

rule sourmash_comp_k51:
        input:
                get_fastq_input,
        output:
                sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k51.sig"
        log:
                "logs/sourmash-compute-k51/{sample}_{unit}.log"
        conda:
                "../envs/sourmash.yaml"
        shell:
                "sourmash compute --scaled 2000 {input} -o {output} -k=51 2> {log}"

rule sourmash_lca_k51:
        input:
                sig = "results/sourmash/sig/{sample}/{sample}_{unit}_k51.sig",
                db = expand("results/sourmash_lca_db/{db}", db = config["sourmash_lca_name_k51"]),
        output:
                sum = "results/sourmash/lca-class/{sample}/{sample}_{unit}_k51.csv"
        log:
                "logs/sourmash-lca-k51/{sample}_{unit}.log"
        conda:
                "../envs/sourmash.yaml"
        shell:
                "sourmash lca summarize --query {input.sig} --db {input.db} -o {output} 2> {log}"


