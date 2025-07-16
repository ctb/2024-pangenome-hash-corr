#! /usr/bin/env python
"""
Calculate hash presence/absence information for many samples based on
hashes in a sig.
"""
import sys
import argparse
import sourmash
import csv
from collections import defaultdict
import numpy
import pickle

import sourmash_utils
from sourmash_plugin_pangenomics import NAMES
from hash_presence_lib import HashPresenceInformation, read_ranktable_csv



def main():
    p = argparse.ArgumentParser()
    p.add_argument('source_sketch')
    p.add_argument('sketches', nargs='+')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-C', '--category-out')
    sourmash_utils.add_standard_minhash_args(p)
    args = p.parse_args()

    hash_to_sample = defaultdict(set)
    samples = []

    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"selecting sketches: {select_mh}")

    # load the source sketch for the hashes for query
    print(f"loading source sketch from file '{args.source_sketch}'")
    query_ss = sourmash_utils.load_index_and_select(args.source_sketch,
                                                    select_mh)
    assert len(query_ss) == 1
    query_ss = list(query_ss.signatures())[0]
    query_minhash = query_ss.minhash
    query_minhash = query_minhash.downsample(scaled=args.scaled)
    query_hashes = set(query_minhash.hashes)

    print(f"loaded {len(query_hashes)} hashes at {args.scaled}.")

    # calculate sample presence
    n = 0
    for sketch_filename in args.sketches:
        idx = sourmash_utils.load_index_and_select(sketch_filename,
                                                   select_mh)
        for ss in idx.signatures():
            n += 1
            sig_name = ss.name
            if n and n % 10 == 0:
                print('...', n)

#        if n > 100:
#            break
            sig_mh = ss.minhash.downsample(scaled=args.scaled)
            if query_minhash.contained_by(sig_mh) > 0:
                samples.append(sig_name)
                sig_hashes = set(sig_mh.hashes) & query_hashes
                for hashval in sig_hashes:
                    hash_to_sample[hashval].add(sig_name)

    presence_info = HashPresenceInformation(ksize=args.ksize,
                                            scaled=args.scaled,
                                            moltype=select_mh.moltype,
                                            classify_d={},
                                            hash_to_sample=hash_to_sample)

    presence_info.save_to_file(args.output)

    if args.category_out:
        with open(args.category_out, 'wt') as fp:
            fp.write('label,category\n')
            for label in samples:
                fp.write(f'{label},default\n')


if __name__ == '__main__':
   sys.exit(main())
