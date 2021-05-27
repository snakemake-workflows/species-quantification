
# fractionation and concatenation
rule fraction_sr:
    input:
        bac_fq1="results/art/bac/{bac_ref}_1.fq",
        bac_fq2="results/art/bac/{bac_ref}_2.fq",
    output:
        out_fq1="results/fractions/short_reads/{bac_ref}_{p}_1.fastq",
        out_fq2="results/fractions/short_reads/{bac_ref}_{p}_2.fastq",
    log:
        "logs/seqtk/short_read/{bac_ref}_{p}.log",
    params:
        fraction="{p}",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.bac_fq2} {params} > {output.out_fq2} 2> {log}"


rule concat_fractions_sr:
    input:
        bac_fq1=expand(
            "results/fractions/short_reads/{bac_ref}_{{p}}_1.fastq", bac_ref=bac.bacteria
        ),
        bac_fq2=expand(
            "results/fractions/short_reads/{bac_ref}_{{p}}_2.fastq", bac_ref=bac.bacteria
        ),
        hum_fq1="results/art/hum/Sample{n}_1.fq",
        hum_fq2="results/art/hum/Sample{n}_2.fq",
    output:
        out_fq1="results/mixed_sr/mixed_{p}_Sample{n}_1.fastq",
        out_fq2="results/mixed_sr/mixed_{p}_Sample{n}_2.fastq",
    log:
        "logs/concatenation/sr/mixed_{p}_Sample{n}.log",
    shell:
        "cat {input.hum_fq1} {input.bac_fq1} > {output.out_fq1}; "
        "cat {input.hum_fq2} {input.bac_fq2} > {output.out_fq2}"


rule fraction_lr:
    input:
        bac_fq="results/nanosim/bac/{bac_ref}_aligned_reads.fastq",
    output:
        out_fq="results/fractions/long_reads/{bac_ref}_{p}.fastq",
    log:
        "logs/seqtk/long_read/{bac_ref}_{p}.log",
    params:
        fraction="{p}",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk sample -s100 {input.bac_fq} {params} > {output.out_fq} 2> {log}"


rule concat_fractions_lr:
    input:
        bac_fq=expand(
            "results/fractions/long_reads/{bac_ref}_{{p}}.fastq",
            bac_ref=bac.bacteria,
            p=config["p"],
        ),
        hum_fq="results/nanosim/hum/{n}_aligned_reads.fastq",
    output:
        out_fq="results/mixed_lr/mixed_{p}_Sample{n}.fastq",
    log:
        "logs/concatenation/sr/mixed_{p}_Sample{n}.log",
    shell:
        "cat {input.hum_fq} {input.bac_fq} > {output.out_fq}"
