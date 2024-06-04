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


from hash_presence_lib import HashPresenceInformation

def main():
    p = argparse.ArgumentParser()
    p.add_argument('presence_pickle')
    p.add_argument('--min-cluster-size', type=int, default=15,
                   help='hdbscan min_cluster_size parameter')
    p.add_argument('--output-tsne-plot',
                   help='save (optional) tSNE plot to this file')
    p.add_argument('--output-assoc-plot',
                   help='save (optional) square hash association plot to this file')
    p.add_argument('--output-presence-plot',
                   help='save (optional) rectangular sample presence plot to this file')
    p.add_argument(
        "--cluster-prefix",
        default=None,
        help="prefix to prepend to cluster names; default is presence pickle file",
    )
    p.add_argument('--scaled', type=int, default=None)
    p.add_argument('--min-presence', type=int, default=5)
    p.add_argument('--save-categories-csv',
                   help="write categories CSV for clusters")
    args = p.parse_args()

    presence_info = HashPresenceInformation.load_from_file(args.presence_pickle)
    print(f"loaded {len(presence_info.hash_to_sample)} hash to sample entries.")
    if args.scaled:
        presence_info = presence_info.downsample(args.scaled)
        print(f"downsampled to {presence_info.scaled}; {len(presence_info.hash_to_sample)} hashes left.")

    if args.min_presence > 1:
        presence_info = presence_info.filter_by_min_samples(args.min_presence)
        print(f"filtered to min_presence={args.min_presence}; {len(presence_info.hash_to_sample)} hashes left.")

    # build similarity matrix
    hashvals, cmp = presence_info.build_association_matrix()

    # turn into distance matrix
    dist = 1 - cmp

    # cluster!
    print(f"clustering using hdbscan with min_cluster_size={args.min_cluster_size}")
    min_cluster_size=args.min_cluster_size

    hdbscan = sklearn.cluster.HDBSCAN(min_cluster_size=min_cluster_size)
    labels = hdbscan.fit_predict(dist)

    print(f'got {numpy.unique(labels).max()} clusters')

    ## pull out the clusters ;)
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
    prefix = args.cluster_prefix or os.path.basename(args.presence_pickle)
    print(f"outputting clusters with prefix '{prefix}'")

    mh_template = presence_info._make_minhash_obj()

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
    if args.output_tsne_plot:
        print(f"running tSNE & saving to {args.output_tsne_plot}")
        tsne = sklearn.manifold.TSNE(n_components=2, random_state=42, perplexity=50) # play with: perplexity
        tsne_coords = tsne.fit_transform(dist)

        palette = sns.color_palette('deep', numpy.unique(labels).max() + 1)
        cluster_colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]

        plt.scatter(tsne_coords[:, 0], tsne_coords[:, 1], color=cluster_colors)
        plt.xlabel("Dimension 1")
        plt.ylabel("Dimension 2")

        plt.savefig(args.output_tsne_plot)

    ## square association plot, our standard "cmp"

    if args.output_assoc_plot:
        print(f"running square association plot & saving to {args.output_assoc_plot}")

        palette = sns.color_palette('deep', numpy.unique(labels).max() + 1)
        # make uncluster => white
        cluster_colors = [palette[x] if x >= 0 else (1.0, 1.0, 1.0) for x in labels]

        fig = sns.clustermap(cmp, xticklabels=[], yticklabels=[], figsize=(8, 8), row_colors=cluster_colors) # , col_colors=category_colors)

        # create a custom legend of just the pangenome rank colors on the columns
        if 0:
            legend_elements = []
            for k, v in category_map.items():
                legend_elements.append(
                    Line2D([0], [0], color=v, label=k, marker="o", lw=0)
                )
                fig.ax_col_dendrogram.legend(handles=legend_elements)

        plt.savefig(args.output_assoc_plot)

    ## rectangular presence plot

    if args.output_presence_plot:
        print(f"running rectangular presence plot & saving to {args.output_presence_plot}")

        # make presence_mat!
        _, _, presence_mat = presence_info.build_presence_matrix()

        palette = sns.color_palette('deep', numpy.unique(labels).max() + 1)
        # make uncluster => white
        cluster_colors = [palette[x] if x >= 0 else (1.0, 1.0, 1.0) for x in labels]

        fig = sns.clustermap(presence_mat,
                             figsize=(11, 8),
                             xticklabels=[], yticklabels=[],
                             col_colors=cluster_colors,
                             cbar_pos=None)
        plt.savefig(args.output_presence_plot)

    # save categories file for clustermap1/plot3 plotting?
    if args.save_categories_csv:
        print(f"writing cluster categories CSV to '{args.save_categories_csv}'")
        with sourmash_args.FileOutputCSV(args.save_categories_csv) as csv_fp:
            w = csv.writer(csv_fp)

            w.writerow(["label", "category", "cluster"])
            for k, v in clusters_d.items():
                name = f"cluster {k}"
                for hashval in v:
                    w.writerow([hashval, name, k])
            for hashval in unclust:
                w.writerow([hashval, "unclustered", -1])


if __name__ == '__main__':
    sys.exit(main())
