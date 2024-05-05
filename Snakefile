METAG_SAMPLES, = glob_wildcards('agatha-gather/{sample}.gather.csv')
print(f"loaded {len(METAG_SAMPLES)} samples.")

SCALED=10000

rule all:
    input:
        expand("agatha-fullgather/{sample}.gather.csv", sample=METAG_SAMPLES),

rule assoc:
    input:
        "agatha-genomes.10k.cmp.matrix.png",
        "agatha-genomes.10k.cmp.matrix.pdf",
        "agatha-ihmp.10k.cmp.matrix.png",
        "agatha-ihmp.10k.cmp.matrix.pdf",

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
        metag = "metag/{sample}.sig.zip",
    output:
        csv = touch("agatha-fullgather/{sample}.gather.csv"),
        out = touch("agatha-fullgather/{sample}.gather.out"),
    shell: """
        sourmash gather {input.metag} {input.db} --picklist {input.fastgather}::prefetch -o {output.csv} > {output.out}
    """
        
rule calc_genome_presence:
    input:
        pangenome_csv = "pangenome.agathobacter_faecis.csv",
        db = "gtdb-rs214-agatha-k21.zip"
    output: "agatha-genomes.10k.dump"
    shell: """
        ./calc-hash-presence.py {input.pangenome_csv} {input.db} -o {output} --scaled={SCALED}
    """

rule calc_metag_presence:
    input:
        pangenome_csv = "pangenome.agathobacter_faecis.csv",
        db = "ihmp-k21.sigs.zip",
    output: "agatha-ihmp.10k.dump"
    shell: """
        ./calc-hash-presence.py {input.pangenome_csv} {input.db} -o {output} --scaled={SCALED}
    """

rule make_assoc_matrix:
    input:
        "{sample}.dump",
    output:
        "{sample}.cmp",
    shell: """
        ./calc-hash-assoc.py {input} -o {output} --scaled=50000
    """

rule make_png_plot:
    input: "{cmp}"
    output: "{cmp}.matrix.png"
    shell: "sourmash plot {input}"

rule make_pdf_plot:
    input: "{cmp}"
    output: "{cmp}.matrix.pdf"
    shell: "sourmash plot {input} --pdf"
