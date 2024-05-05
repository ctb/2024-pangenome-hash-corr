#! /usr/bin/env python
import sys
import argparse
import sourmash
import csv
from collections import defaultdict
import numpy


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
    p.add_argument('--min-presence', type=int, default=5)
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

    # calculate sample presence
    for n, metag_ss in enumerate(idx.signatures()):
        if n and n % 10 == 0:
            print('...', n)
#        if n > 100:
#            break
        metag_mh = metag_ss.minhash.downsample(scaled=args.scaled)
        metag_name = metag_ss.name
        if query_minhash.contained_by(metag_mh) > 0:
            metag_hashes = set(metag_mh.hashes)
            for hashval in hashes:
                if hashval in metag_hashes:
                    hash_to_sample[hashval].add(metag_name)
    
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
