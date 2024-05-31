#! /usr/bin/env python
import sys
import argparse
import sourmash
from sourmash import sourmash_args
import csv
from collections import defaultdict
import numpy

from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('presence_pickle')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('--min-presence', type=int, default=5)
    p.add_argument('--pangenome-types', type=str, default=None)
    p.add_argument('--compare-csv', help="write compare CSV", default=None)
    p.add_argument('--categories-csv', help="write categories CSV",
                   default=None)
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

    hashes = list(sorted(hash_to_sample))

    pa = numpy.zeros((len(hashes), len(hashes)), dtype=float)

    for i in range(len(hashes)):
        hash_i = hashes[i]
        presence_i = hash_to_sample[hash_i]
        
        for j in range(i):
            hash_j = hashes[j]
            presence_j = hash_to_sample[hash_j]
            jaccard = len(presence_i.intersection(presence_j)) / \
                len(presence_i.union(presence_j))
            pa[i][j] = jaccard
            pa[j][i] = jaccard

        pa[i][i] = 1

    with open(args.output, 'wb') as fp:
        numpy.save(fp, pa)
    with open(args.output + '.labels.txt', 'wt') as fp:
        hashstr = [ NAMES[classify_d[x]] for x in hashes ]
        fp.write("\n".join(hashstr))
    with open(args.output + '.labels.csv', 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['sort_order', 'label', 'category'])
        for n, x in enumerate(hashes):
            w.writerow([n, x, NAMES[classify_d[x]]])

    if args.compare_csv:
        with sourmash_args.FileOutputCSV(args.compare_csv) as csv_fp:
            w = csv.writer(csv_fp)
            # hashes as column headers
            w.writerow([ str(h) for h in hashes ])

            for i in range(len(hashes)):
                y = []
                for j in range(len(hashes)):
                    y.append(str(pa[i][j]))
                w.writerow(y)

    if args.categories_csv:
        with sourmash_args.FileOutputCSV(args.categories_csv) as csv_fp:
            w = csv.writer(csv_fp)

            w.writerow(["labels", "category"])
            for hashval in hashes:
                w.writerow([hashval, NAMES[classify_d[hashval]]])


if __name__ == '__main__':
   sys.exit(main())
