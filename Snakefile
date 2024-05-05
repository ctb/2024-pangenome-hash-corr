METAG_SAMPLES, = glob_wildcards('agatha-gather/{sample}.gather.csv')
print(f"loaded {len(METAG_SAMPLES)} samples.")

rule all:
    input:
        expand("agatha-fullgather/{sample}.gather.csv", sample=METAG_SAMPLES)

rule metag_extract:
    input:
        combined = "ihmp-k21.sigs.zip"
    output:
        "metag/{sample}.sig.zip",
    shell: """
        sourmash sig grep {wildcards.sample} {input.combined} -o {output}
    """

rule fullgather:
    input:
        fastgather = "agatha-gather/{sample}.gather.csv",
        db = "gtdb-rs214-agatha-k21.zip",
        metag = "metag/{sample}.sig.zip"
    output:
        csv = touch("agatha-fullgather/{sample}.gather.csv"),
        out = touch("agatha-fullgather/{sample}.gather.out"),
    shell: """
        sourmash gather {input.metag} {input.db} --picklist {input.fastgather}::prefetch -o {output.csv} > {output.out}
    """
        
