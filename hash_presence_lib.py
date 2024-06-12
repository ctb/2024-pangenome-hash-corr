"""
Utilities for loading from & saving to files, as well as other common stuff.
"""
import pickle
import csv
import numpy

import sourmash
import sourmash_utils


class HashPresenceInformation:
    """
    Store a hash presence table & their associated rank classifications.
    """
    def __init__(self, *, ksize=21, scaled=1000, moltype='DNA',
                 classify_d=None, hash_to_sample=None):
        self.ksize = ksize
        self.scaled = scaled
        self.moltype = moltype
        self.classify_d = classify_d
        self.hash_to_sample = hash_to_sample

    def _make_minhash_obj(self):
        return sourmash_utils.FracMinHash(ksize=self.ksize,
                                          moltype=self.moltype,
                                          scaled=self.scaled)

    def downsample(self, new_scaled):
        "Downsample hashes to a new scaled value."
        if new_scaled < self.scaled:
            raise ValueError(f"cannot downsample to {new_scaled}: current scaled is {self.scaled}")

        # downsample
        mh = sourmash.MinHash(n=0, ksize=self.ksize, scaled=new_scaled)
        mh.add_many(self.hash_to_sample)

        # shift over the presence vectors
        hashes = mh.hashes
        hash_to_sample = self.hash_to_sample
        new_d = {}
        for hashval in hashes:
            new_d[hashval] = hash_to_sample[hashval]

        # CTB: could shift over the classify_d here too.
        return HashPresenceInformation(ksize=self.ksize,
                                       scaled=new_scaled,
                                       moltype=self.moltype,
                                       classify_d=self.classify_d,
                                       hash_to_sample=new_d)

    def filter_by_min_samples(self, min_presence):
        "Keep only hashes in a minimum of 'min_presence' samples."
        new_d = {}
        hash_to_sample = self.hash_to_sample
        for hashval, presence in hash_to_sample.items():
            if len(presence) >= min_presence:
                new_d[hashval] = presence

        # CTB: could filter classify_d here too.

        return HashPresenceInformation(ksize=self.ksize,
                                       scaled=self.scaled,
                                       moltype=self.moltype,
                                       classify_d=self.classify_d,
                                       hash_to_sample=new_d)

    def filter_by_pangenome_type(self, typelist):
        "Keep only hashes with specific pangenome ranks."
        assert min(typelist) >= 1
        assert max(typelist) <= 5

        new_d = {}
        hash_to_sample = self.hash_to_sample
        classify_d = self.classify_d
        for hashval, presence in hash_to_sample.items():
            if classify_d.get(hashval) in typelist:
                new_d[hashval] = presence

        # CTB: could filter classify_d here too.

        return HashPresenceInformation(ksize=self.ksize,
                                       scaled=self.scaled,
                                       moltype=self.moltype,
                                       classify_d=self.classify_d,
                                       hash_to_sample=new_d)

    def build_association_matrix(self):
        hash_to_sample = self.hash_to_sample

        hashvals = list(sorted(hash_to_sample))
        print(f"creating {len(hashvals)} by {len(hashvals)} array.")

        cmp = numpy.zeros((len(hashvals), len(hashvals)), dtype=float)

        for i in range(len(hashvals)):
            hash_i = hashvals[i]
            presence_i = hash_to_sample[hash_i]

            for j in range(i):
                hash_j = hashvals[j]
                presence_j = hash_to_sample[hash_j]
                jaccard = len(presence_i.intersection(presence_j)) / \
                    len(presence_i.union(presence_j))
                cmp[i][j] = jaccard
                cmp[j][i] = jaccard

            cmp[i][i] = 1

        return hashvals, cmp

    def build_presence_matrix(self):
        # get list of samples:
        all_samples = set()
        for k, vv in self.hash_to_sample.items():
            all_samples.update(vv)

        print(f"got {len(all_samples)} samples for presence plot.")

        sample_to_idx = {}
        for n, sample_name in enumerate(sorted(all_samples)):
            sample_to_idx[sample_name] = n

        hashval_to_idx = {}
        for n, hashval in enumerate(self.hash_to_sample):
            hashval_to_idx[hashval] = n

        print(f"creating presence matrix: {len(sample_to_idx)} x {len(hashval_to_idx)}")
        presence_mat = numpy.zeros((len(sample_to_idx), len(hashval_to_idx)))

        for hashval, sample_set in self.hash_to_sample.items():
            hashval_i = hashval_to_idx[hashval]
            for sample_name in sample_set:
                sample_j = sample_to_idx[sample_name]

                presence_mat[sample_j][hashval_i] = 1

        return sample_to_idx, hashval_to_idx, presence_mat
        
    def save_to_file(self, filename):
        "Save an object of this class to a file."
        with open(filename, 'wb') as fp:
            pickle.dump(self, fp)

    @classmethod
    def load_from_file(cls, filename):
        "Load an object of this class from a file."
        with open(filename, 'rb') as fp:
            obj = pickle.load(fp)
        assert isinstance(obj, cls)
        return obj


def read_ranktable_csv(filename):
    "Read a ranktable CSV."
    with open(filename, 'r', newline='') as fp:
        r = csv.DictReader(fp)

        classify_d = {}
        for row in r:
            hashval = int(row['hashval'])
            classify_as = int(row['pangenome_classification'])
            classify_d[hashval] = classify_as

        return classify_d
