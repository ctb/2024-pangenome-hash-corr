.PHONY: test all clean

all: test

test:
	snakemake -s Snakefile.test -c 4

clean:
	rm -f agatha-genomes.*.{csv,png,assoc,dump,dump.*} 
	snakemake -s Snakefile.test -c 4 --delete-all-output

cleanall: clean
	snakemake -s Snakefile.test -c 4
