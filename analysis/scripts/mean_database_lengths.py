"""Retrieve mean sequence lengths of databases"""

import numpy as np

for db in ["grpI.db", "grpII.db", "RNaseP.db", "srp.db"]:
    lengths = []
    with open(f"analysis/rnastructure/data/{db}") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                header = line[1:].strip()
                seq = f.readline().strip()
                nat_ss = f.readline().strip()
                seqlen = len(seq)
                lengths.append(seqlen)
            line = f.readline()
    print(db, np.mean(lengths))