#! /usr/bin/env python
import sys
import argparse
import sourmash
import csv
from collections import defaultdict
import numpy
import pickle


CENTRAL_CORE=1
EXTERNAL_CORE=2
SHELL=3
INNER_CLOUD=4
SURFACE_CLOUD=5


NAMES = { CENTRAL_CORE: 'central core',
          EXTERNAL_CORE: 'external core',
          SHELL: 'shell',
          INNER_CLOUD: 'inner cloud',
          SURFACE_CLOUD: 'surface cloud' }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('hashlist')
    p.add_argument('metagenomes')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-k', '--ksize', type=int, default=21)
    p.add_argument('--scaled', type=int, default=100000)
    p.add_argument('--filter-samples', default=None)
    args = p.parse_args()

    with open(args.hashlist, 'r', newline='') as fp:
        r = csv.DictReader(fp)

        classify_d = {}
        for row in r:
            hashval = int(row['hashval'])
            classify_as = int(row['pangenome_classification'])
            classify_d[hashval] = classify_as

    print(f"loaded {len(classify_d)} hashvals... downsampling soon.")

    hash_to_sample = defaultdict(set)

    idx = sourmash.load_file_as_index(args.metagenomes)
    idx = idx.select(ksize=args.ksize, scaled=args.scaled)

    print(f"found {len(idx)} metagenomes")

    query_minhash = next(iter(idx.signatures())).minhash.copy_and_clear()
    for hashval in classify_d:
        query_minhash.add_hash(hashval)
    query_minhash = query_minhash.downsample(scaled=args.scaled)

    hashes = list(sorted(set(query_minhash.hashes)))
    print(f"at downsampled scale={args.scaled}, {len(hashes)} hashes found.")

    filter_by_name = None
    if args.filter_samples:
        filter_by_name = set([ x.strip() for x in open(args.filter_samples) ])

    # calculate sample presence
    n_skipped = 0
    for n, metag_ss in enumerate(idx.signatures()):
        metag_name = metag_ss.name
        if n and n % 10 == 0:
            print('...', n)
        if filter_by_name is not None and metag_name not in filter_by_name:
            n_skipped += 1
            continue

#        if n > 100:
#            break
        metag_mh = metag_ss.minhash.downsample(scaled=args.scaled)
        if query_minhash.contained_by(metag_mh) > 0:
            metag_hashes = set(metag_mh.hashes)
            for hashval in hashes:
                if hashval in metag_hashes:
                    hash_to_sample[hashval].add(metag_name)

    to_save = (args.ksize, args.scaled, classify_d, hash_to_sample)
    with open(args.output, 'wb') as fp:
        pickle.dump(to_save, fp)

    print(f'skipped: {n_skipped}')


if __name__ == '__main__':
   sys.exit(main())
