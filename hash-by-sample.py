#! /usr/bin/env python
"""
using hash presence info, produce a rectangular matrix of hashval x sample.
"""
import sys
import argparse
import sourmash
from sourmash import sourmash_args
import csv

from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('presence_pickle')
    p.add_argument('-o', '--output', required=True, help="output CSV")
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('-m', '--min-presence', type=int, default=5)
    p.add_argument('--pangenome-types', type=str, default=None)
    p.add_argument('-C', '--categories-csv', default=None)
    args = p.parse_args()

    presence_info = HashPresenceInformation.load_from_file(args.presence_pickle)
    print(f"loaded {len(presence_info.hash_to_sample)} hash to sample entries.")
    if args.scaled:
        presence_info = presence_info.downsample(args.scaled)
        print(f"downsampled to scaled={presence_info.scaled}; {len(presence_info.hash_to_sample)} hashes left.")

    if args.min_presence > 1:
        presence_info = presence_info.filter_by_min_samples(args.min_presence)
        print(f"filtered to min_presence={args.min_presence}; {len(presence_info.hash_to_sample)} hashes left.")

    # filter for pangenome_types
    if args.pangenome_types:
        typelist = list(map(int, list(args.pangenome_types)))
        presence_info = presence_info.filter_by_pangenome_type(typelist)

        print(f"After pangenome-hash-type filtering to {typelist}, {len(presence_info.hash_to_sample)} left.")

    classify_d = presence_info.classify_d
    hash_to_sample = presence_info.hash_to_sample

    # collect hashes
    hashes = list(sorted(hash_to_sample))

    # collect samples
    samples = set()
    for h in hashes:
        samples.update(hash_to_sample[h])
    samples = list(sorted(samples))
    sample_to_idx = {}
    for i, sample_name in enumerate(samples):
        sample_to_idx[sample_name] = i


    n_written = 0
    with open(args.output, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(['query_name', 'match_name', 'presence'])

        for hashval in hashes:
            presence_j = hash_to_sample[hashval]

            for sample_name in presence_j:
                w.writerow([sample_name, hashval, 1])
                n_written += 1

        print(f"wrote {n_written} entries to '{args.output}'")
        print(f"use 'sourmash scripts clustermap1' from betterplot to plot!")
        print(f"e.g. 'sourmash scripts clustermap1 {args.output} -o fig.png")

    if args.categories_csv:
        n_written = 0
        with open(args.categories_csv, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(['label','category'])

            for hashval in hashes:
                category = classify_d.get(hashval)
                if category in NAMES:
                    category_name = NAMES[category]
                    w.writerow([hashval, category_name])
                    n_written += 1
                else:
                    assert category is None

        print(f"{n_written} category entries written to '{args.categories_csv}'")


if __name__ == '__main__':
   sys.exit(main())
