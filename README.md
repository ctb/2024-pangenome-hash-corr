# 2024-pangenome-hash-corr

## Calculating hash correlations across sketches (genomes or metagenomes).

This repository contains three primary command-line scripts:

1. `calc-hash-presence.py` - calculate hash presence/absence information for many samples.
2. `hash-by-hash-assoc.py` - using hash presence info, produce a square similarity matrix of hashval x hashval.
3. `hash-by-sample-assoc.py` - using hash presence info, produce a rectangular matrix of hashval x sample.

## Quickstart

Produce the presence/absence dump file by measuring the
presence/absence of the hashes from the given ranktable in the
provided sketches.

```
./calc-hash-presence.py ranktable.agathobacter_faecis.csv \
    gtdb-rs214-agatha-k21.zip -o agatha-genomes.1k.dump \
    --scaled=1000
```

Next, produce a square matrix of hash x hash correlations across samples:
```
./hash-by-hash-assoc.py agatha-genomes.1k.dump \
    -o agatha-genomes.10k.assoc --scaled=10000 --min-presence=2 \
    -C agatha-genomes.1k.assoc.categories.csv
sourmash scripts plot3 agatha-genomes.10k.assoc \
    agatha-genomes.10k.assoc.labels.csv -o agatha-genomes.10k.assoc.png \
    -C agatha-genomes.1k.assoc.categories.csv
```

These commands produce this plot:
![](example_output/agatha-genomes.10k.assoc.png)


Finally, produce a rectangular matrix showing hash x genome correlations:
```
./hash-by-sample.py agatha-genomes.1k.dump \
    -o agatha-genomes.10k.presence.csv \
    --categories-csv agatha-genomes.1k.presence.categories.csv
sourmash scripts clustermap1 agatha-genomes.10k.presence.csv -u presence \
    -C agatha-genomes.1k.presence.categories.csv \
    -o agatha-genomes.10k.presence.png
```

These commands produce this plot:

![](example_output/agatha-genomes.10k.presence.png)

## Docoumentation

### Running `calc-hash-presence.py`

Usage: 
```
./calc-hash-presence.py <ranktable_csv> <sample1> [<samples ...>] \
    -o <output>.dump
```
will calculate presence/absence info for the hash values in
`ranktable_csv` across `sample*` sketches, saving the info to
`<output>.dmp`.

`ranktable_csv` is in the format produced by the
`sourmash_plugin_pangenomics` command `pangenome_ranktable`.

Optional parameters:

* `-k`, `--ksize` - select k-mer size
* `filter-samples` - use only these samples
* `--moltype` (CTB: does not yet work)

