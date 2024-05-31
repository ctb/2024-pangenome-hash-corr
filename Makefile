.PHONY: test all

all: test

test:
	snakemake -s Snakefile.test -c 4

cleanall:
	snakemake -s Snakefile.test -c 4 --delete-all-output
	snakemake -s Snakefile.test -c 4
