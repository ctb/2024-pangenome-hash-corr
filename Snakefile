# metagenome/hash exploration.

# sourmash scripts fastmultigather ihmp-k21.sigs.zip gtdb-rs214-agatha-k21.zip  -k 21 -t 0
# sourmash scripts fastmultigather ../ihmp-k21.sigs.zip ../agatha-merged.sig.zip -k 21 -t 0


#METAG_SAMPLES, = glob_wildcards('agatha-gather/{sample}.gather.csv')
METAG_SAMPLES = [ x.strip() for x in open('METAG-AGATHA-5.list.txt') ]
print(f"loaded {len(METAG_SAMPLES)} samples.")

SCALED=10000

rule all:
    input:
        expand("agatha-fullgather/{sample}.gather.csv", sample=METAG_SAMPLES),
        expand("agatha-pangenome/{sample}.gather.csv", sample=METAG_SAMPLES),
        expand("agatha-pangenome/{sample}.agatha.sig.zip", sample=METAG_SAMPLES),

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
        
rule pang_gather:
    input:
        db = "agatha-merged.sig.zip",
        metag = "metag/{sample}.sig.zip",
    output:
        csv = touch("agatha-pangenome/{sample}.gather.csv"),
        out = touch("agatha-pangenome/{sample}.gather.out"),
    shell: """
        sourmash gather {input.metag} {input.db} -o {output.csv} > {output.out}
    """

rule pang_intersect:
    input:
        q = "agatha-merged.sig.zip",
        metag = "metag/{sample}.sig.zip",
    output:
        "agatha-pangenome/{sample}.agatha.sig.zip",
    shell: """
        sourmash sig intersect {input} -o - | sourmash sig rename - {wildcards.sample}.agatha -o {output}
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
        filter_list = "METAG-AGATHA-5.list.txt"
    output: protected("agatha-ihmp.10k.dump")
    shell: """
        ./calc-hash-presence.py {input.pangenome_csv} {input.db} -o {output} --scaled={SCALED} --filter-samples {input.filter_list}
    """

rule make_assoc_matrix_genomes:
    input:
        "agatha-genomes.10k.dump",
    output:
        "agatha-genomes.10k.cmp",
    shell: """
        ./hash-by-hash-assoc.py {input} -o {output} --scaled=10000 --min=2 # --pangenome 3
    """

rule make_assoc_matrix_ihmp:
    input:
        "agatha-ihmp.10k.dump",
    output:
        "agatha-ihmp.10k.cmp",
    shell: """
        ./hash-by-hash-assoc.py {input} -o {output} --scaled=10000 --min=10 --pangenome 12345
    """

rule make_png_plot:
    input: "{cmp}"
    output: "{cmp}.matrix.png"
    shell: "sourmash plot {input}"

rule make_pdf_plot:
    input: "{cmp}"
    output: "{cmp}.matrix.pdf"
    shell: "sourmash plot {input} --pdf"
