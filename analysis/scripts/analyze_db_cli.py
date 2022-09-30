#!/usr/bin/env python

"""Command line script to analyze a database in dot-bracket format. 
It calculates the distances of the normal-DrTransformer/reversed-DrTransformer/RNAsubopt structures to the native structure."""

import os
import subprocess
import json
import RNA
import math
import sys
import argparse

def AnalyzeDB(filepath, dirpath, o_prune=0.05, t_fast=0.001):
    """Compute the distance to the native structure for the norma/reversed DrTransformer prediction and the thermodynamic equilibrium.
    Writes results to json-file."""
    # store db entries as dictionary
    db = {}
    with open(f"{filepath}") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                header = line[1:].strip()
                seq = f.readline().strip()
                nat_ss = f.readline().strip()
                seqlen = len(seq)
                db[header] = {"seq": seq, "nat_ss": nat_ss, "seqlen": seqlen}
            line = f.readline()
    #print(db)

    db_name = filepath.split("/")[-1].split(".")[0]
    # change working directory
    try:
        os.mkdir(f"{dirpath}/{db_name}")
    except FileExistsError:
        pass
    os.chdir(f"{dirpath}/{db_name}")

    constraint_f = open(f"{db_name}_c.db", "w")  # write constrained database to a file

    # analysis
    for num, header in enumerate(db):
        print(f"{num+1}/{len(db)}, {header}")

        # count N´s
        db[header]["N_count"] = db[header]['seq'].upper().count("N")

        # count >30bp unpaired stretches in native ss
        db[header]["long_unpaired"] = 0
        unpaired_stretches = sorted(db[header]['nat_ss'].replace("(","|").replace(")","|").split("|"), reverse=True)
        for unp_s in unpaired_stretches:
            if len(unp_s) >= 30:
                db[header]["long_unpaired"] += 1
            else:
                break

        # Constraint folding
        # RNAfold
        subprocess.run(f"RNAfold --noPS --canonicalBPonly -C > rnafold_c", text=True, input=f"{db[header]['seq']}\n{db[header]['nat_ss']}\n@\n", shell=True)
        with open(f"rnafold_c") as f:
            lines = f.readlines()
            line_split = lines[-1].split()
            natc_ss = line_split[0]
            natc = float(line_split[-1].strip("()"))
            db[header]["natc_ss"] = natc_ss  # natc: native-constrained
            db[header]["natc"] = natc
        constraint_f.write(f"{header}\n{db[header]['seq']}\n{natc_ss}\n")

        # RNAfold
        subprocess.run(f"echo {db[header]['seq']} | RNAfold --noPS > rnafold", shell=True)
        with open(f"rnafold") as f:
            lines = f.readlines()
            line_split = lines[-1].split()
            mfe_ss= line_split[0]
            mfe = float(line_split[-1].strip("()"))
            db[header]["mfe_ss"] = mfe_ss
            db[header]["mfe"] = mfe

        # RNAsubopt
        # Compute suboptimal structures in 2 kcal/mol range
        # Compute distance to part of the ensemble - alternative to computing distance only to MFE structure
        subprocess.run(f"echo {db[header]['seq']} | RNAsubopt -e 2 -s > rnasubopt", shell=True)
        with open(f"rnasubopt") as f:
            lines = f.readlines()
            subopt = [[line.split()[0], float(line.split()[1].strip("()"))] for line in lines[1:]]
        _kt = 1.987e-3 * (273.15 + 37) # kcal/mol, 37°C
        for i in range(len(subopt)):
            subopt[i].append(math.e**(-subopt[i][1] / _kt))  # format: [ss, energy, probability]
        sum_subopt = sum([el[2] for el in subopt])
        for el in subopt:
            el[2] = el[2]/sum_subopt
        dist_nat_subopt = 0
        dist_natc_subopt = 0
        for el in subopt:
            dist_nat_subopt += RNA.bp_distance(db[header]["nat_ss"], el[0]) * el[2]
            dist_natc_subopt += RNA.bp_distance(db[header]["natc_ss"], el[0]) * el[2]
        db[header]["dist_nat_subopt"] = dist_nat_subopt
        db[header]["dist_natc_subopt"] = dist_natc_subopt

        # RNAeval
        subprocess.run(f"RNAeval > rnaeval", text=True, input=f"{db[header]['seq']}\n{db[header]['nat_ss']}\n@\n", shell=True)
        with open(f"rnaeval") as f:
            lines = f.readlines()
            native_e = float(lines[-1].split()[-1].strip("()"))
            db[header]["native_e"] = native_e

        # RNAdistance
        dist_nat_mfe = RNA.bp_distance(db[header]["nat_ss"], db[header]["mfe_ss"])
        db[header]["dist_nat_mfe"] = dist_nat_mfe
        dist_natc_mfe = RNA.bp_distance(db[header]["natc_ss"], db[header]["mfe_ss"])
        db[header]["dist_natc_mfe"] = dist_natc_mfe
        dist_nat_natc = RNA.bp_distance(db[header]["nat_ss"], db[header]["natc_ss"])
        db[header]["dist_nat_natc"] = dist_nat_natc 

        # DrTransformer
        subprocess.run(f"echo {db[header]['seq']} | DrTransformer --name drtr_{header} --logfile --no-timecourse --o-prune {o_prune} --t-end 60 --t-fast {t_fast}", shell=True)
        subprocess.run(f"echo {db[header]['seq']} | DrTransformer --name drtr_{header}_rev --logfile --no-timecourse -rt --o-prune {o_prune} --t-end 60 --t-fast {t_fast}", shell=True)
        for mode in [0, 1]:
            with open(f"drtr_{header}{['', '_rev'][mode]}.log") as f:
                lines = f.readlines()
                d = {}
                for line in lines[::-1]:
                    words = line.split()
                    if words[0] != str(db[header]["seqlen"]):
                        break
                    d[words[2]] = [float(words[4]), float(words[6][:-1])] # format: {ss: [occ: after t-end, occ: equilibrium]}
                dist_nat_tend = 0
                dist_nat_equ = 0
                dist_natc_tend = 0
                dist_natc_equ = 0
                dist_mfe_tend = 0
                dist_subopt_tend = 0
                for el in d:
                    # dist to native ss
                    if d[el][0] + d[el][1] == 0: # continue if not populated
                        continue
                    dist_nat = RNA.bp_distance(db[header]["nat_ss"], el)
                    dist_nat_tend += dist_nat * d[el][0]
                    dist_nat_equ += dist_nat * d[el][1]
                    # dist to native-constrained ss
                    dist_natc = RNA.bp_distance(db[header]["natc_ss"], el)
                    dist_natc_tend += dist_natc * d[el][0]
                    dist_natc_equ += dist_natc * d[el][1]
                    # dist to mfe (after t-end)
                    dist_mfe = RNA.bp_distance(db[header]["mfe_ss"], el)
                    dist_mfe_tend += dist_mfe * d[el][0]
                    # dist to subopt (after t-end)
                    for el_subopt in subopt:
                        dist_subopt = RNA.bp_distance(el_subopt[0], el)
                        dist_subopt_tend += dist_subopt * d[el][0] * el_subopt[2]
                sum_tend = sum([d[el][0] for el in d])
                sum_equ = sum([d[el][1] for el in d])
                try:
                    dist_nat_tend /= sum_tend
                    dist_natc_tend /= sum_tend
                    dist_mfe_tend /= sum_tend
                    dist_subopt_tend /= sum_tend
                except ZeroDivisionError: # should not be possible
                    dist_nat_tend = float('nan')
                    dist_natc_tend = float('nan')
                    dist_mfe_tend = float('nan')
                    dist_subopt_tend = float('nan')
                try:
                    dist_nat_equ /= sum_equ
                    dist_natc_equ /= sum_equ
                except ZeroDivisionError: # is possible (edit: not anymore)
                    dist_nat_equ = float('nan')
                    dist_natc_equ = float('nan')
                db[header]["dist_nat_tend" + ['', '_rev'][mode]] = dist_nat_tend
                db[header]["dist_nat_equ" + ['', '_rev'][mode]] = dist_nat_equ
                db[header]["dist_natc_tend" + ['', '_rev'][mode]] = dist_natc_tend
                db[header]["dist_natc_equ" + ['', '_rev'][mode]] = dist_natc_equ
                db[header]["dist_mfe_tend" + ['', '_rev'][mode]] = dist_mfe_tend
                db[header]["dist_subopt_tend" + ['', '_rev'][mode]] = dist_subopt_tend
    constraint_f.close()
    # write dictionary to json file
    with open(f"{db_name}.json", "w") as f:
        print(json.dumps(db, indent=4), file=f)


def restricted_float(x):
    y = float(x)
    if y < 0.0 or y > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return y


def main():
    """ Analysis of a database in dot-bracket format.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Analysis of a database in dot-bracket format.')

    parser.add_argument("-i", "--input", metavar = '<str>',
            help = """Path of the input file""")
    parser.add_argument("-o", "--output", default = ".", metavar = '<str>',
            help = """Path of the output directory. Must be an existing directory""")
    parser.add_argument("--o-prune", type = restricted_float, default = 0.05, metavar = '<flt>',
            help = """DrTransformer: Occupancy threshold to prune structures from the 
            network. The structures with lowest occupancy are removed until
            at most o-prune occupancy has been removed from the total population.""")
    parser.add_argument("--t-fast", type = float, default = 0.001, metavar = '<flt>',
            help = """DrTransformer: Folding times faster than --t-fast are considered
            instantaneous. Structural transitions that are faster than
            --t-fast are considered part of the same macrostate. Directly
            translates to an energy barrier separating conformations using:
            dG = -RT*ln((1/t-fast)/k0). None: t-fast = 1/k_0 """)

    args = parser.parse_args()

    AnalyzeDB(args.input, args.output, args.o_prune, args.t_fast)
 
if __name__ == '__main__':
    main()
