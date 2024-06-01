#! /usr/bin/env python
"""
Calculate hash presence/absence information for many samples.
"""
import sys
import argparse
import sourmash
import csv
from collections import defaultdict
import numpy
import pickle

from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation, read_ranktable_csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('ranktable_csv')
    p.add_argument('metagenomes')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-k', '--ksize', type=int, default=21)
    p.add_argument('--scaled', type=int, default=100000)
    p.add_argument('--filter-samples', default=None)
    args = p.parse_args()

    classify_d = read_ranktable_csv(args.ranktable_csv)
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

    presence_info = HashPresenceInformation(ksize=args.ksize,
                                            scaled=args.scaled,
                                            moltype='DNA',
                                            classify_d=classify_d,
                                            hash_to_sample=hash_to_sample)

    presence_info.save_to_file(args.output)

    print(f'skipped: {n_skipped}')


if __name__ == '__main__':
   sys.exit(main())
