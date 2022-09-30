"""Create a constraint-folded database. Not used in the project because it is included in analyze_db_cli.py."""

import subprocess
import json
import RNA
import os

def ConstraintFoldDB(filename):
    db_name = filename.split(".")[0]
    # change working directory
    dirpath = "analysis/rnastrand/"
    os.chdir(f"{dirpath}")

    # store db entries as dictionary
    db = {}
    with open(f"data/{filename}") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                header = line.strip() # different "header" than in analyze_db
                seq = f.readline().strip()
                native_ss = f.readline().strip()
                seqlen = len(seq)
                db[header] = {"seq": seq, "native_ss": native_ss, "seqlen": seqlen}
            line = f.readline()
    #print(db)

    out_f = open(f"data/constraint_folding/{db_name}_c.db", "w")
    for num, header in enumerate(db):
        print(f"{num+1}/{len(db)}, {header}")

        # RNAfold
        subprocess.run(f"RNAfold --noPS --canonicalBPonly --noLP -C > rnafold_c", text=True, input=f"{db[header]['seq']}\n{db[header]['native_ss']}\n@\n", shell=True)
        with open(f"rnafold_c") as f:
            lines = f.readlines()
            line_split = lines[-1].split()
            mfe_ss = line_split[0]
            mfe = float(line_split[-1].strip("()"))
            db[header]["mfe_c_ss"] = mfe_ss
            db[header]["mfe_c"] = mfe
        out_f.write(f"{header}\n{db[header]['seq']}\n{mfe_ss}\n")

        # RNAeval
        subprocess.run(f"RNAeval > rnaeval", text=True, input=f"{db[header]['seq']}\n{db[header]['native_ss']}\n@\n", shell=True)
        with open(f"rnaeval") as f:
            lines = f.readlines()
            native_e = float(lines[-1].split()[-1].strip("()"))
            db[header]["native_e"] = native_e

        # RNAdistance
        dist_nat_mfe_c = RNA.bp_distance(db[header]["native_ss"], db[header]["mfe_c_ss"])
        db[header]["dist_nat_mfe_c"] = dist_nat_mfe_c  

    out_f.close()

    # write dictionary to json file
    with open(f"data/constraint_folding/{db_name}_c_compare.json", "w") as f:
        print(json.dumps(db, indent=4), file=f)


ConstraintFoldDB("rnastrand.db")