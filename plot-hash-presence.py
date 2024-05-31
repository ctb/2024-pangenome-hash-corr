#! /usr/bin/env python
import sys
import argparse
import sourmash
from sourmash import sourmash_args
import csv
from collections import defaultdict
import numpy
import pickle
import seaborn as sns

from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('presence_pickle')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('--min-presence', type=int, default=5)
    p.add_argument('--pangenome-types', type=str, default=None)
    p.add_argument('--categories-csv', default=None)
    args = p.parse_args()

    presence_info = HashPresenceInformation.load_from_file(args.presence_pickle)
    ksize = presence_info.ksize
    scaled = presence_info.scaled
    if args.scaled is None:
        args.scaled = presence_info.scaled
    moltype = presence_info.moltype # @CTB
    classify_d = presence_info.classify_d
    hash_to_sample = presence_info.hash_to_sample

    print(f"loaded {len(hash_to_sample)} hash to sample entries.")

    if args.scaled > scaled:
        # downsample
        mh = sourmash.MinHash(n=0, ksize=ksize, scaled=args.scaled)
        for hashval in hash_to_sample:
            mh.add_hash(hashval)

        hashes = mh.hashes
        new_d = {}
        for hashval in hashes:
            new_d[hashval] = hash_to_sample[hashval]

        hash_to_sample = new_d
        print(f"after downsampling from {scaled} => {args.scaled}, {len(new_d)} left.")

        scaled = args.scaled
    elif args.scaled == scaled:
        pass
    else:
        assert 0, f"cannot downsample to {args.scaled}, lower than {scaled}"

    # filter hashes on presence
    new_d = {}
    for hashval, presence in hash_to_sample.items():
        if len(presence) >= args.min_presence:
            new_d[hashval] = presence

    hash_to_sample = new_d
    print(f"After presence-filtering to >= {args.min_presence}, {len(new_d)} left.")

    # filter for pangenome_types
    if args.pangenome_types:
        typelist = list(map(int, list(args.pangenome_types)))
        print(typelist)
        assert min(typelist) >= 1
        assert max(typelist) <= 5

        new_d = {}
        for hashval, presence in hash_to_sample.items():
            if classify_d.get(hashval) in typelist:
                new_d[hashval] = presence

        print(f"After pangenome-hash-type filtering to {typelist}, {len(new_d)} left.")

        hash_to_sample = new_d

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
