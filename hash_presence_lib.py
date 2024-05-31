import pickle

class HashPresenceInformation:
    def __init__(self, *, ksize=21, scaled=1000, moltype='DNA',
                 classify_d=None, hash_to_sample=None):
        self.ksize = ksize
        self.scaled = scaled
        self.moltype = moltype
        self.classify_d = classify_d
        self.hash_to_sample = hash_to_sample
        
    def save_to_file(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(self, fp)

    @classmethod
    def load_from_file(cls, filename):
        with open(filename, 'rb') as fp:
            obj = pickle.load(fp)
        assert isinstance(obj, cls)
        return obj
