import random

def FilterDB(old_db, new_db, retain=0.1):
    """Create a new database with only a fraction of randomly selected database entries"""
    with open(f"{old_db}") as f:
        with open(f"{new_db}", "w") as f_new:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    if random.random() > retain:
                        f.readline()
                        f.readline()
                        line = f.readline()
                        continue
                    header = line.strip()
                    seq = f.readline().strip().upper()
                    native_ss = f.readline().strip()
                    seqlen = len(seq)
                    # if 40 <= seqlen <= 250:
                    if True:
                        # if native_ss.count("." * 30) == 0:
                        if True:
                            f_new.write(f"{header}\n{seq}\n{native_ss}\n")
                line = f.readline()

FilterDB("analysis/rnastructure/data/srp_middle.db", "analysis/rnastructure/data/srp_middle_retain0_2.db", retain=0.21)