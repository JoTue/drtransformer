"""Filter database based on sequence lengths"""

def FilterDB(old_db, new_db):
    with open(old_db) as f:
        with open(new_db, "w") as f_new:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    header = line.strip()
                    seq = f.readline().strip().upper()
                    native_ss = f.readline().strip()
                    seqlen = len(seq)
                    if 201 <= seqlen <= 300:
                        # if native_ss.count("." * 30) == 0:
                        if True:
                            f_new.write(f"{header}\n{seq}\n{native_ss}\n")
                line = f.readline()

def FilterDB_long(old_db, new_db):
    with open(f"analysis/rnastrand/data/{old_db}") as f:
        with open(f"analysis/rnastrand/data/{new_db}", "w") as f_new:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    header = line.strip()
                    seq = f.readline().strip().upper()
                    native_ss = f.readline().strip()
                    seqlen = len(seq)
                    if 251 <= seqlen <= 500:
                        # if native_ss.count("." * 30) == 0:
                        if True:
                            f_new.write(f"{header}\n{seq}\n{native_ss}\n")
                line = f.readline()

def FilterDB2(old_db, new_db): # sprinzl
    with open(f"analysis/rnastrand/data/{old_db}") as f:
        with open(f"analysis/rnastrand/data/{new_db}", "w") as f_new:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    header = line.strip()
                    seq = f.readline().strip().upper()
                    native_ss = f.readline().strip()
                    seqlen = len(seq)
                    if 40 <= seqlen <= 250:
                        if not "T" in seq:
                            f_new.write(f"{header}\n{seq}\n{native_ss}\n")
                line = f.readline()

FilterDB("analysis/rnastructure/data/srp.db", "analysis/rnastructure/data/srp_REMOVE.db")