
rule sourmash_table_sr:
    input:
        "/home/uzuner/Desktop/sm-test/mixed_1000_Sample1.csv"
    output:
        "/home/uzuner/Desktop/sm-test/test.txt"
    log:
        notebook="logs/notebooks/sourmash_table-sr.ipynb"
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "notebooks/sourmash_table_sr.py.ipynb"