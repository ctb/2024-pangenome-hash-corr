#! /usr/bin/env python
"""
Cluster the output of 'hash-by-hash-assoc.py'
"""
import sys
import os
import argparse
import csv

import sourmash
from sourmash import sourmash_args
import numpy
import seaborn as sns
import sklearn.cluster
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.lines import Line2D


from sourmash_plugin_betterplot import load_categories_csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cmp_matrix')
    p.add_argument('cmp_labels_from')
    p.add_argument('-m', '--min-cluster-size', type=int, default=15,
                   help='hdbscan min_cluster_size parameter')
    p.add_argument('--output-tsne',
                   help='save (optional) tSNE plot to this file')
    p.add_argument(
        "--cluster-prefix",
        default=None,
        help="prefix to prepend to cluster names; default is cmp file",
    )
    # these are needed because we don't track them in the matrix :sob:
    p.add_argument('--output-ksize', default=31, type=int)
    p.add_argument('--output-moltype', default='DNA')
    p.add_argument('--output-scaled', default=1000, type=int)
    p.add_argument('--save-categories-csv')
    args = p.parse_args()

    # load matrix & turn into distance matrix
    cmp = numpy.load(args.cmp_matrix)
    dist = 1 - cmp

    # load matrix labels
    with sourmash_args.FileInputCSV(args.cmp_labels_from) as r:
        labelinfo = list(r)

    # cluster!
    print(f"clustering using hdbscan with min_cluster_size={args.min_cluster_size}")
    min_cluster_size=args.min_cluster_size

    hdbscan = sklearn.cluster.HDBSCAN(min_cluster_size=min_cluster_size)
    labels = hdbscan.fit_predict(dist)

    print(f'got {numpy.unique(labels).max()} clusters')

    ## pull out the clusters ;)
    hashinfo = list(labelinfo)

    hashinfo.sort(key=lambda row: int(row["sort_order"]))
    hashvals = [ int(row['label']) for row in hashinfo ]

    clusters_d = defaultdict(set)
    unclust = set()
    for hashval, cluster_num in zip(hashvals, labels):
        if cluster_num >= 0:
            clusters_d[cluster_num].add(hashval)
        else:
            unclust.add(hashval)

    # reorder by size:
    new_clusters_d = {}
    new_cluster_id = 0
    for n, v in enumerate(sorted(clusters_d.values(), key=lambda x: -len(x))):
        new_clusters_d[new_cluster_id] = v
        new_cluster_id += 1
    clusters_d = new_clusters_d

    print(f"5 largest clusters (of {len(clusters_d)}:")
    for n, (k, v) in enumerate(sorted(clusters_d.items())):
        print(f"\tcluster {k} has size {len(v)}")
        if n > 5:
            break

    # output clusters
    prefix = args.cluster_prefix or os.path.basename(args.cmp_matrix)
    print(f"outputting clusters with prefix '{prefix}'")

    assert args.output_moltype == 'DNA'
    mh_template = sourmash.MinHash(n=0, ksize=args.output_ksize,
                          scaled=args.output_scaled)

    cluster_n = 0
    for k, v in sorted(clusters_d.items()):
        hashvals = clusters_d[k]
        mh = mh_template.copy_and_clear()
        mh.add_many(hashvals)
        ss = sourmash.SourmashSignature(mh, name=f'cluster_{cluster_n}')

        outp = f"{prefix}.cluster_{cluster_n}.sig.zip"
        with sourmash_args.SaveSignaturesToLocation(outp) as save_sig:
            save_sig.add(ss)

        cluster_n += 1
    
    # plot tSNE?
    if args.output_tsne:
        print(f"running tSNE & saving to {args.output_tsne}")
        tsne = sklearn.manifold.TSNE(n_components=2, random_state=42, perplexity=50) # play with: perplexity
        tsne_coords = tsne.fit_transform(dist)

        palette = sns.color_palette('deep', numpy.unique(labels).max() + 1)
        cluster_colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]

        plt.scatter(tsne_coords[:, 0], tsne_coords[:, 1], color=cluster_colors)
        plt.xlabel("Dimension 1")
        plt.ylabel("Dimension 2")

        plt.savefig(args.output_tsne)

    # save categories file for clustermap1/plot3 plotting?
    if args.save_categories_csv:
        print(f"writing cluster categories CSV to '{args.save_categories_csv}'")
        with sourmash_args.FileOutputCSV(args.save_categories_csv) as csv_fp:
            w = csv.writer(csv_fp)

            w.writerow(["label", "category"])
            for k, v in clusters_d.items():
                name = f"cluster {k}"
                for hashval in v:
                    w.writerow([hashval, name])
            for hashval in unclust:
                w.writerow([hashval, "unclustered"])


if __name__ == '__main__':
    sys.exit(main())
