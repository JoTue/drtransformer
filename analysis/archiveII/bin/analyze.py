#!/usr/bin/env python

"""Command line script to analyze a database in dot-bracket format. 
It calculates the distances of the normal-DrTransformer/reversed-DrTransformer/RNAsubopt structures to the native structure."""

import os
import subprocess
import json
from time import time
import RNA
import math
import sys
import argparse

def get_basepair_probabilities(seq):
    fc = RNA.fold_compound(seq)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()
    return basepair_probs

def get_basepairs(ss):
    basepairs = []
    stack = []
    for i, c in enumerate(ss):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if len(stack) == 0:
                raise IndexError("No matching closing parenthesis at: " + str(i))
            basepairs.append((stack.pop(), i))
    if len(stack) > 0:
        raise IndexError("No matching opening parenthesis at: " + str(stack.pop()))
    return basepairs

def get_distance(basepair_probs, basepair_list):
    d = 0
    for i in range(0, len(basepair_probs)-2):
        for j in range(i+1, len(basepair_probs[i])-1):
            if (i, j) in basepair_list:
                d += 1 - basepair_probs[i+1][j+1]
            else:
                d += basepair_probs[i+1][j+1]
    return d

def analyze(dataset, id1, id2, tmp_dir, t_list, t_ext):
    """Compute the distance to the native structure for the normal/reversed DrTransformer prediction and the thermodynamic equilibrium."""
    # # change working directory
    # os.chdir(f"results/{os.path.basename(dataset).split('.')[0]}")

    assert id1 <= id2
    with open(f"../../{dataset}") as f:
        i = 1
        line = f.readline()
        while line:
            if line.startswith(">"):
                if i < id1:
                    f.readline()
                    f.readline()
                    line = f.readline()
                    i += 1
                    continue
                # Parse
                db = {}
                header = line[1:].strip()
                seq = f.readline().strip()
                nat_ss = f.readline().strip()
                seqlen = len(seq)
                db[header] = {"seq": seq, "nat_ss": nat_ss, "seqlen": seqlen}

                # Analyze

                # Constraint folding
                # RNAfold
                subprocess.run(f"RNAfold --noPS --canonicalBPonly -C > {tmp_dir}/rnafold_c", text=True, input=f"{seq}\n{nat_ss}\n@\n", shell=True)
                with open(f"{tmp_dir}/rnafold_c") as f:
                    lines = f.readlines()
                    line_split = lines[-1].split()
                    natc_ss = line_split[0]
                    natc_e = float(line_split[-1].strip("()"))
                    db[header]["natc_ss"] = natc_ss  # natc: native-constrained
                    db[header]["natc_e"] = natc_e

                # RNAfold
                subprocess.run(f"echo {seq} | RNAfold --noPS > {tmp_dir}/rnafold", shell=True)
                with open(f"{tmp_dir}/rnafold") as f:
                    lines = f.readlines()
                    line_split = lines[-1].split()
                    mfe_ss= line_split[0]
                    mfe = float(line_split[-1].strip("()"))
                    db[header]["mfe_ss"] = mfe_ss
                    db[header]["mfe"] = mfe
                
                # RNAeval (get energy of native structure)
                subprocess.run(f"RNAeval > {tmp_dir}/rnaeval", text=True, input=f"{seq}\n{nat_ss}\n@\n", shell=True)
                with open(f"{tmp_dir}/rnaeval") as f:
                    lines = f.readlines()
                    nat_e = float(lines[-1].split()[-1].strip("()"))
                    db[header]["nat_e"] = nat_e
                
                # Compute distances of native structures to thermodynamic equilibirum
                basepair_probs = get_basepair_probabilities(seq)
                basepair_list_nat_ss = get_basepairs(nat_ss)
                basepair_list_natc_ss = get_basepairs(natc_ss)
                dist_nat_eq = get_distance(basepair_probs, basepair_list_nat_ss)
                dist_natc_eq = get_distance(basepair_probs, basepair_list_natc_ss)
                db[header]["dist_nat_eq"] = dist_nat_eq
                db[header]["dist_natc_eq"] = dist_natc_eq

                # RNAsubopt (only as comparison, dist_nat_eq/dist_natc_eq are more accurate)
                # Compute suboptimal structures in 2 kcal/mol range
                # Compute distance to part of the ensemble - alternative to computing distance only to MFE structure
                subprocess.run(f"echo {seq} | RNAsubopt -e 2 -s > {tmp_dir}/rnasubopt", shell=True)
                with open(f"{tmp_dir}/rnasubopt") as f:
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
                    dist_nat_subopt += RNA.bp_distance(nat_ss, el[0]) * el[2]
                    dist_natc_subopt += RNA.bp_distance(natc_ss, el[0]) * el[2]
                db[header]["dist_nat_subopt"] = dist_nat_subopt
                db[header]["dist_natc_subopt"] = dist_natc_subopt

                # RNAdistance
                dist_nat_mfe = RNA.bp_distance(nat_ss, mfe_ss)
                db[header]["dist_nat_mfe"] = dist_nat_mfe
                dist_natc_mfe = RNA.bp_distance(natc_ss, mfe_ss)
                db[header]["dist_natc_mfe"] = dist_natc_mfe
                dist_nat_natc = RNA.bp_distance(nat_ss, natc_ss)
                db[header]["dist_nat_natc"] = dist_nat_natc 

                # Parse DrTransformer output
                for mode in [0, 1]:
                    with open(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}{['', '_rev'][mode]}.log") as f:  # parse log-file for equilibrium occupancies (ignore t-end occupancies)
                        lines = f.readlines()
                        d = {}
                        for line in lines[::-1]:
                            words = line.split()
                            if words[0] != str(db[header]["seqlen"]):
                                break
                            d[words[2]] = [float(words[4]), float(words[6][:-1])] # format: {ss: [occ: after t-end, occ: equilibrium]}
                        dist_nat_dteq = 0
                        dist_natc_dteq = 0
                        dist_mfe_dteq = 0
                        dist_subopt_dteq = 0
                        dist_eq_dteq = 0
                        for el in d:
                            # dist to native ss
                            if d[el][0] + d[el][1] == 0: # continue if not populated
                                continue
                            dist_nat = RNA.bp_distance(db[header]["nat_ss"], el)
                            dist_nat_dteq += dist_nat * d[el][1]
                            
                            # dist to native-constrained ss
                            dist_natc = RNA.bp_distance(db[header]["natc_ss"], el)
                            dist_natc_dteq += dist_natc * d[el][1]
                            
                            # dist to mfe
                            dist_mfe = RNA.bp_distance(db[header]["mfe_ss"], el)
                            dist_mfe_dteq += dist_mfe * d[el][1]
                            
                            # dist to subopt
                            for el_subopt in subopt:
                                dist_subopt = RNA.bp_distance(el_subopt[0], el)
                                dist_subopt_dteq += dist_subopt * d[el][1] * el_subopt[2]

                            # dist to eq
                            basepair_list_el = get_basepairs(el)
                            dist_eq_dteq += get_distance(basepair_probs, basepair_list_el) * d[el][1]

                        sum_dteq = sum([d[el][1] for el in d])
                        try:
                            dist_nat_dteq /= sum_dteq
                            dist_natc_dteq /= sum_dteq
                            dist_mfe_dteq /= sum_dteq
                            dist_subopt_dteq /= sum_dteq
                            dist_eq_dteq /= sum_dteq
                        except ZeroDivisionError:
                            dist_nat_dteq = float('nan')
                            dist_natc_dteq = float('nan')
                            dist_mfe_dteq = float('nan')
                            dist_subopt_dteq = float('nan')
                            dist_eq_dteq = float('nan')

                        db[header]["dist_nat_dteq" + ['', '_rev'][mode]] = dist_nat_dteq
                        db[header]["dist_natc_dteq" + ['', '_rev'][mode]] = dist_natc_dteq
                        db[header]["dist_mfe_dteq" + ['', '_rev'][mode]] = dist_mfe_dteq
                        db[header]["dist_subopt_dteq" + ['', '_rev'][mode]] = dist_subopt_dteq
                        db[header]["dist_eq_dteq" + ['', '_rev'][mode]] = dist_eq_dteq
                    
                    with open(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}{['', '_rev'][mode]}.drf") as f:  # parse drf-file for occupancies at time points given by t-list
                        # get timepoints to look at
                        transcription_duration = (seqlen - 1) * t_ext
                        timepoints = [transcription_duration + t_ext]
                        timepoints.extend([time+transcription_duration for time in t_list])
                        # parse
                        lines = f.readlines()
                        d = {timepoint:{} for timepoint in timepoints}  # format: {timepoint: {structure_index: occupancy}}
                        d_ss = {}  # {index: [ss, dist_nat, dist_natc, dist_mfe, dist_eq]}
                        for line in lines[::-1]:
                            words = line.split()
                            index, timepoint, occupancy, ss = int(words[0]), float(words[1]), float(words[2]), words[3]
                            if timepoint not in timepoints:
                                break
                            if occupancy == 0:
                                continue
                            if int(words[0]) not in d_ss:
                                d_ss[index] = [ss, None, None, None, None]
                            d[timepoint][index] = occupancy
                        # calculate distances per ss
                        for index in d_ss:
                            ss = d_ss[index][0]
                            d_ss[index][1] = RNA.bp_distance(nat_ss, ss)
                            d_ss[index][2] = RNA.bp_distance(natc_ss, ss)
                            d_ss[index][3] = RNA.bp_distance(mfe_ss, ss)
                            d_ss[index][4] = get_distance(basepair_probs, get_basepairs(ss))
                        # calculate distances per timepoint
                        for timepoint in timepoints:
                            dist_nat = 0
                            dist_natc = 0
                            dist_mfe = 0
                            dist_eq = 0
                            for index in d[timepoint]:
                                occupancy = d[timepoint][index]
                                dist_nat += occupancy * d_ss[index][1]  # occ * dist
                                dist_natc += occupancy * d_ss[index][2]
                                dist_mfe += occupancy * d_ss[index][3]
                                dist_eq += occupancy * d_ss[index][4]
                            # store distances in dict
                            db[header][timepoint] = {}
                            db[header][timepoint]["dist_nat" + ['', '_rev'][mode]] = dist_nat
                            db[header][timepoint]["dist_natc" + ['', '_rev'][mode]] = dist_natc
                            db[header][timepoint]["dist_mfe" + ['', '_rev'][mode]] = dist_mfe
                            db[header][timepoint]["dist_eq" + ['', '_rev'][mode]] = dist_eq
                
                # write dictionary to json file (each database entry as separate json file, get concatenated later)
                with open(f"analysis_results/{os.path.basename(dataset).split('.')[0]}/{header}.json", "w") as json_f:
                    print(json.dumps(db, indent=4), file=json_f)

            line = f.readline()
            i += 1
            if i > id2:
                break


def main():
    """ Analysis of a database in dot-bracket format.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Analysis of a database in dot-bracket format.')

    parser.add_argument("dataset",
            help = """Path of the dataset file""")
    parser.add_argument("--id1", type = int,
            help = """Start-Index (1-based) of sequences to compute.""")
    parser.add_argument("--id2", type = int,
            help = """End-Index (1-based) of sequences to compute.""")
    parser.add_argument("--tmp",
            help = """Temporary directory.""")
    parser.add_argument("--t-list", nargs="+", default = ['1', '10', '60', '600', '3600'], metavar = '<str>',
            help = """List of times after transcription shown in drf-file.
            Highest time replaces end time (--t-end) in log-file.
            (--t-end and --t-log are ignored)""")
    parser.add_argument("--t-ext", type = float, default = 0.04, metavar = '<flt>',
            help = """Inverse of transcription rate, i.e. time per nucleotide extension
            [seconds per nucleotide].""")

    args = parser.parse_args()

    analyze(args.dataset, args.id1, args.id2, args.tmp, args.t_list, args.t_ext)
 
if __name__ == '__main__':
    main()
