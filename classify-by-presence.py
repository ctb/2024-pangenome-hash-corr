#! /usr/bin/env python
import sys
import argparse
import pickle
import csv


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
    args = p.parse_args()

    with open(args.presence_pickle, 'rb') as fp:
        saved_info = pickle.load(fp)

    ksize, scaled, classify_d, hash_to_sample = saved_info

    sample_presence = [ len(v) for v in hash_to_sample.values() ]
    max_sample_presence = max(sample_presence)

    print(f"max sample presence: {max_sample_presence}")

    central_core_threshold = 0.95 #0.95 is core , 90% is technically soft core
    external_core_threshold = 0.90
    shell_threshold = 0.10 #0.10
    inner_cloud_threshold = 0.01 # 0.0 is the full cloud, but trimming (0.001?) may be necessary to create the viz...?
    surface_cloud_threshold = 0.00

    central_core = []
    external_core = []
    shell = []
    inner_cloud = []
    surface_cloud = []

    xx = []
    for hashval, presence in hash_to_sample.items():
        f_present = len(presence) / max_sample_presence;
        if f_present >= central_core_threshold:
            tag = CENTRAL_CORE
        elif f_present >= external_core_threshold:
            tag = EXTERNAL_CORE
        elif f_present >= shell_threshold:
            tag = SHELL
        elif f_present >= inner_cloud_threshold:
            tag = INNER_CLOUD
        elif f_present >= surface_cloud_threshold:
            tag = SURFACE_CLOUD
        else:
            assert 0

        xx.append((hashval, tag))

    with open(args.output, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['hashval', 'pangenome_classification'])
        
        for hashval, tag in xx:
            w.writerow([hashval, tag])


if __name__ == '__main__':
   sys.exit(main())
