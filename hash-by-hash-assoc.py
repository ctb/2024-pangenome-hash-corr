#! /usr/bin/env python
"""
using hash presence info, produce a square similarity matrix of hashval x hashval.
"""
import sys
import argparse
import sourmash
from sourmash import sourmash_args
import csv
import numpy

from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('presence_pickle')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('-m', '--min-presence', type=int, default=5)
    p.add_argument('--pangenome-types', type=str, default=None)
    p.add_argument('--compare-csv', help="write compare CSV", default=None)
    p.add_argument('-C', '--categories-csv', help="write categories CSV",
                   default=None)
    args = p.parse_args()

    presence_info = HashPresenceInformation.load_from_file(args.presence_pickle)
    print(f"loaded {len(presence_info.hash_to_sample)} hash to sample entries.")
    if args.scaled:
        presence_info = presence_info.downsample(args.scaled)
        print(f"downsampled to {presence_info.scaled}; {len(presence_info.hash_to_sample)} hashes left.")

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

    hashes, pa = presence_info.build_association_matrix()

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

    print(f"writing similarity matrix to '{args.output}'")
    with open(args.output, 'wb') as fp:
        numpy.save(fp, pa)
    # for 'sourmash plot'

    #with open(args.output + '.labels.txt', 'wt') as fp:
    #    hashstr = [ NAMES[classify_d[x]] for x in hashes ]
    #    fp.write("\n".join(hashstr))

    # for 'sourmash plot --labels-from'
    labels_to = args.output + '.labels.csv'
    with open(labels_to, 'w', newline='') as fp:
        print(f"writing labels_to file to '{labels_to}'")
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
        print(f"writing categories CSV to '{args.categories_csv}'")
        with sourmash_args.FileOutputCSV(args.categories_csv) as csv_fp:
            w = csv.writer(csv_fp)

            w.writerow(["label", "category"])
            for hashval in hashes:
                w.writerow([hashval, NAMES[classify_d[hashval]]])


if __name__ == '__main__':
   sys.exit(main())
