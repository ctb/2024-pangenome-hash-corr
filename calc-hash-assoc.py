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
    p.add_argument('presence_pickle')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('--min-presence', type=int, default=5)
    args = p.parse_args()

    with open(args.presence_pickle, 'rb') as fp:
        saved_info = pickle.load(fp)

    ksize, scaled, classify_d, hash_to_sample = saved_info
    if args.scaled is None:
        args.scaled = scaled

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

    print(f"After presence-filtering to >= {args.min_presence}, {len(new_d)} left.")

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


if __name__ == '__main__':
   sys.exit(main())