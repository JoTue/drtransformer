#!/usr/bin/env python3

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
import tempfile
import shutil

def get_basepair_probabilities(seq):
    fc = RNA.fold_compound(seq)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()
    return basepair_probs

def get_basepairs(ss):
    basepairs = set()
    stack = []
    for i, c in enumerate(ss):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if len(stack) == 0:
                raise IndexError("No matching closing parenthesis at: " + str(i))
            basepairs.add((stack.pop(), i))
    if len(stack) > 0:
        raise IndexError("No matching opening parenthesis at: " + str(stack.pop()))
    return basepairs

def get_distance(basepair_probs, basepair_list):
    d = 0
    if basepair_list:  # dist of eq to structure
        for i in range(0, len(basepair_probs)-2):
            for j in range(i+1, len(basepair_probs[i])-1):
                if (i, j) in basepair_list:
                    d += 1 - basepair_probs[i+1][j+1]
                else:
                    d += basepair_probs[i+1][j+1]
    else: # dist of eq to itself
        for i in range(0, len(basepair_probs)-2):
            for j in range(i+1, len(basepair_probs[i])-1):
                d += 2* basepair_probs[i+1][j+1] * (1 - basepair_probs[i+1][j+1])
    return d

def metrics(n, ss1=None, ss2=None, bpp1=None, bpp2=None):
    # n=seqlen, 1=reference (native structure), 2=predicted, ss (secondary structure), bpp (base pair probabilities)
    if ss1 and ss2:  # comparing two structures
        bp1, bp2 = get_basepairs(ss1), get_basepairs(ss2)  # set of all basepairs as tuples of indices

        TP = len(bp1.intersection(bp2))
        FP = len(bp2.difference(bp1))
        FN = len(bp1.difference(bp2))
        TN = 0.5 * n * (n-1) - TP  # or: 0.5 * n * (n-1) - TP - FP - FN
    else:  # comparison involves equilibrium
        TP, FP, FN = 0, 0, 0
        if ss1 and bpp2:
            bp1 = get_basepairs(ss1)
            for i in range(0, len(bpp2)-2):
                for j in range(i+1, len(bpp2[i])-1):
                    if (i, j) in bp1:
                        TP += bpp2[i+1][j+1]
                        FN += 1 - bpp2[i+1][j+1]
                    else:
                        FP += bpp2[i+1][j+1]
        elif ss2 and bpp1:
            bp2 = get_basepairs(ss2)
            for i in range(0, len(bpp1)-2):
                for j in range(i+1, len(bpp1[i])-1):
                    if (i, j) in bp2:
                        TP += bpp1[i+1][j+1]
                        FP += 1 - bpp1[i+1][j+1]
                    else:
                        FN += bpp1[i+1][j+1]
        elif bpp1 and bpp2:  # TODO is this correct/valid??
            for i in range(0, len(bpp1)-2):
                for j in range(i+1, len(bpp1[i])-1):
                    TP += bpp1[i+1][j+1] * bpp2[i+1][j+1]
                    FP += (1 - bpp1[i+1][j+1]) * bpp2[i+1][j+1]
                    FN += bpp1[i+1][j+1] * (1 - bpp2[i+1][j+1])
        TN = 0.5 * n * (n-1) - TP  # or: 0.5 * n * (n-1) - TP - FP - FN

    try:
        TPR = TP / (TP + FN)    # True positive rate / sensitivity / sensibility
    except ZeroDivisionError:
        TPR = float('nan')
    try:
        PPV = TP / (TP + FP)    # positive predictive value / selectivity
    except ZeroDivisionError:
        PPV = float('nan')

    try:
        MCC = ((TP*TN) - (FP*FN)) / math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))  # alternative: with epsilon
    except ZeroDivisionError:
        MCC = float('nan')
    try:
        F1 = (2*PPV*TPR) / (PPV+TPR)
    except ZeroDivisionError:
        F1 = float('nan')

    return {"TP": TP, "FP": FP, "FN": FN, "TN": TN, "TPR": TPR, "PPV": PPV, "MCC": MCC, "F1": F1}

def analyze(dataset, id1, id_end, ids_per_job, tmp_dir, t_list, t_ext):
    """Compute the distance to the native structure for the normal/reversed DrTransformer prediction and the thermodynamic equilibrium."""
    # change working directory
    os.chdir(tmp_dir)
    tmp_dir = tempfile.mkdtemp()
    os.chdir("..")

    id2 = min([id1 + ids_per_job - 1, id_end])
    assert id1 <= id2
    missing_files = []
    missing_file_count = 0
    incomplete_seq_count = 0
    with open(f"{dataset}") as data_f:
        ind = 1
        line = data_f.readline()
        while line:
            if line.startswith(">"):
                if ind < id1:
                    data_f.readline()
                    data_f.readline()
                    line = data_f.readline()
                    ind += 1
                    if ind > id2:
                        break
                    continue

                # Parse
                db = {}
                header = line[1:].strip()
                seq = data_f.readline().strip()
                nat_ss = data_f.readline().strip()
                seqlen = len(seq)
                db[header] = {"seq": seq, "nat_ss": nat_ss, "seqlen": seqlen}
                db[header]["eq"] = {"nat": {}, "natc": {}, "mfe": {}, "eq": {}}
                db[header]["equilibrium"] = {"nat": {}, "natc": {}, "mfe": {}, "eq": {}, "nat_rev": {}, "natc_rev": {}, "mfe_rev": {}, "eq_rev": {}}

                if not os.path.isfile(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}.log") or not os.path.isfile(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}_rev.log"):
                    incomplete_seq_count += 1
                    if not os.path.isfile(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}.log"):
                        missing_files.append(f"drtr_{header}.log")
                        missing_file_count += 1
                    if not os.path.isfile(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}_rev.log"):
                        missing_files.append(f"drtr_{header}_rev.log")
                        missing_file_count += 1
                    print(f"Skipping: {ind} {header}")
                    line = data_f.readline()
                    ind += 1
                    if ind > id2:
                        break
                    continue
                print(f"{ind} {header}")

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

                # Compute distances to thermodynamic equilibirum
                basepair_probs = get_basepair_probabilities(seq)

                # basepair_list_nat_ss = get_basepairs(nat_ss)
                # basepair_list_natc_ss = get_basepairs(natc_ss)
                # metrics_eq2 = metrics(seqlen, bpp1=basepair_probs, ss2=el)
                # dist_nat_eq = get_distance(basepair_probs, basepair_list_nat_ss)
                # dist_natc_eq = get_distance(basepair_probs, basepair_list_natc_ss)
                # db[header]["bpdis_nat_eq"] = dist_nat_eq
                # db[header]["bpdis_natc_eq"] = dist_natc_eq

                # dist to native ss
                basepair_list = get_basepairs(nat_ss)
                dist_nat = get_distance(basepair_probs, basepair_list)
                metrics_nat = metrics(seqlen, ss1=nat_ss, bpp2=basepair_probs)
                metrics_nat["bpdis"] = dist_nat
                db[header]["eq"]["nat"] = metrics_nat

                # dist to native-constrained ss
                basepair_list = get_basepairs(natc_ss)
                dist_natc = get_distance(basepair_probs, basepair_list)
                metrics_natc = metrics(seqlen, ss1=natc_ss, bpp2=basepair_probs)
                metrics_natc["bpdis"] = dist_natc
                db[header]["eq"]["natc"] = metrics_natc

                # dist to mfe
                basepair_list = get_basepairs(mfe_ss)
                dist_mfe = get_distance(basepair_probs, basepair_list)
                metrics_mfe = metrics(seqlen, ss1=mfe_ss, bpp2=basepair_probs)
                metrics_mfe["bpdis"] = dist_mfe
                db[header]["eq"]["mfe"] = metrics_mfe

                # dist to eq (dist to itself)
                dist_eq = get_distance(basepair_probs, None)
                metrics_eq = metrics(seqlen, bpp1=basepair_probs, bpp2=basepair_probs)
                metrics_eq["bpdis"] = dist_eq
                db[header]["eq"]["eq"] = metrics_eq

                # RNAsubopt (only as comparison, dist_nat_eq/dist_natc_eq are more accurate)
                # Compute suboptimal structures in 2 kcal/mol range
                # Compute distance to part of the ensemble - alternative to computing distance only to MFE structure
                subprocess.run(f"echo {seq} | RNAsubopt -e 2 -s > {tmp_dir}/rnasubopt", shell=True)
                with open(f"{tmp_dir}/rnasubopt") as f:
                    lines = f.readlines()
                    subopt = [[line.split()[0], float(line.split()[1].strip("()"))] for line in lines[1:]]
                _kt = 1.987e-3 * (273.15 + 37) # kcal/mol, 37Â°C
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
                db[header]["bpdis_nat_subopt"] = dist_nat_subopt
                db[header]["bpdis_natc_subopt"] = dist_natc_subopt

                # RNAdistance
                dist_nat_mfe = RNA.bp_distance(nat_ss, mfe_ss)
                db[header]["bpdis_nat_mfe"] = dist_nat_mfe
                dist_natc_mfe = RNA.bp_distance(natc_ss, mfe_ss)
                db[header]["bpdis_natc_mfe"] = dist_natc_mfe
                dist_nat_natc = RNA.bp_distance(nat_ss, natc_ss)
                db[header]["bpdis_nat_natc"] = dist_nat_natc

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
                        dist_nat = 0
                        dist_natc = 0
                        dist_mfe = 0
                        dist_eq = 0
                        metrics_nat = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                        metrics_natc = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                        metrics_mfe = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                        metrics_eq = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                        for el in d:
                            # dist to native ss
                            if d[el][0] + d[el][1] == 0: # continue if not populated
                                continue
                            dist_nat = RNA.bp_distance(db[header]["nat_ss"], el)
                            dist_nat += dist_nat * d[el][1]
                            metrics_nat2 = metrics(seqlen, nat_ss, el)
                            metrics_nat = {metric: metrics_nat[metric] + d[el][1] * metrics_nat2[metric] for metric in metrics_nat}

                            # dist to native-constrained ss
                            dist_natc = RNA.bp_distance(db[header]["natc_ss"], el)
                            dist_natc += dist_natc * d[el][1]
                            metrics_natc2 = metrics(seqlen, natc_ss, el)
                            metrics_natc = {metric: metrics_natc[metric] + d[el][1] * metrics_natc2[metric] for metric in metrics_natc}

                            # dist to mfe
                            dist_mfe = RNA.bp_distance(db[header]["mfe_ss"], el)
                            dist_mfe += dist_mfe * d[el][1]
                            metrics_mfe2 = metrics(seqlen, mfe_ss, el)
                            metrics_mfe = {metric: metrics_mfe[metric] + d[el][1] * metrics_mfe2[metric] for metric in metrics_mfe}

                            # dist to eq
                            basepair_list_el = get_basepairs(el)
                            dist_eq += get_distance(basepair_probs, basepair_list_el) * d[el][1]
                            metrics_eq2 = metrics(seqlen, bpp1=basepair_probs, ss2=el)
                            metrics_eq = {metric: metrics_eq[metric] + d[el][1] * metrics_eq2[metric] for metric in metrics_eq}

                        # normalize
                        sum_occ = sum([d[el][1] for el in d])
                        try:
                            dist_nat /= sum_occ
                            dist_natc /= sum_occ
                            dist_mfe /= sum_occ
                            dist_eq /= sum_occ
                            metrics_nat = {metric: metrics_nat[metric] / sum_occ for metric in metrics_nat}
                            metrics_natc = {metric: metrics_natc[metric] / sum_occ for metric in metrics_natc}
                            metrics_mfe = {metric: metrics_mfe[metric] / sum_occ for metric in metrics_mfe}
                            metrics_eq = {metric: metrics_eq[metric] / sum_occ for metric in metrics_eq}
                        except ZeroDivisionError:
                            dist_nat = float('nan')
                            dist_natc = float('nan')
                            dist_mfe = float('nan')
                            dist_eq = float('nan')
                            metrics_nat = {metric: float('nan') for metric in metrics_nat}
                            metrics_natc = {metric: float('nan') for metric in metrics_natc}
                            metrics_mfe = {metric: float('nan') for metric in metrics_mfe}
                            metrics_eq = {metric: float('nan') for metric in metrics_eq}
                        # add bpdis to metrics dict
                        metrics_nat["bpdis"] = dist_nat
                        metrics_natc["bpdis"] = dist_natc
                        metrics_mfe["bpdis"] = dist_mfe
                        metrics_eq["bpdis"] = dist_eq
                        # store distances in dict
                        db[header]["equilibrium"]["nat" + ['', '_rev'][mode]] = metrics_nat
                        db[header]["equilibrium"]["natc" + ['', '_rev'][mode]] = metrics_natc
                        db[header]["equilibrium"]["mfe" + ['', '_rev'][mode]] = metrics_mfe
                        db[header]["equilibrium"]["eq" + ['', '_rev'][mode]] = metrics_eq

                    with open(f"drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}{['', '_rev'][mode]}.drf") as f:  # parse drf-file for occupancies at time points given by t-list
                        # get timepoints to look at
                        transcription_duration = (seqlen - 1) * t_ext
                        timepoints = [transcription_duration + t_ext]
                        timepoints.extend([float(time)+transcription_duration for time in t_list])
                        # parse
                        lines = f.readlines()
                        d = {timepoint:{} for timepoint in timepoints}  # format: {timepoint: {structure_index: occupancy}}
                        d_ss = {}  # {index: [ss, dist_nat, dist_natc, dist_mfe, dist_eq, metrics_nat, metrics_natc, metrics_mfe, metrics_eq]}
                        for line in lines[::-1]:
                            words = line.split()
                            index, timepoint, occupancy, ss = int(words[0]), float(words[1]), float(words[2]), words[3]
                            if timepoint not in timepoints:
                                break
                            if occupancy == 0:
                                continue
                            if int(words[0]) not in d_ss:
                                d_ss[index] = [ss, None, None, None, None, None, None, None, None]
                            d[timepoint][index] = occupancy
                        # check for timepoints missing in drf -> then output would be same as the one from timepoints[-1] (which is always present in output)
                        for timepoint in timepoints:
                            if not d[timepoint]:
                                d[timepoint] = d[timepoints[-1]]
                        # calculate distances per ss
                        for index in d_ss:
                            ss = d_ss[index][0]
                            d_ss[index][1] = RNA.bp_distance(nat_ss, ss)
                            d_ss[index][2] = RNA.bp_distance(natc_ss, ss)
                            d_ss[index][3] = RNA.bp_distance(mfe_ss, ss)
                            d_ss[index][4] = get_distance(basepair_probs, get_basepairs(ss))
                            d_ss[index][5] = metrics(seqlen, nat_ss, ss)
                            d_ss[index][6] = metrics(seqlen, natc_ss, ss)
                            d_ss[index][7] = metrics(seqlen, mfe_ss, ss)
                            d_ss[index][8] = metrics(seqlen, bpp1=basepair_probs, ss2=ss)
                        # calculate distances per timepoint
                        for timepoint in timepoints:
                            dist_nat = 0
                            dist_natc = 0
                            dist_mfe = 0
                            dist_eq = 0
                            metrics_nat = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                            metrics_natc = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                            metrics_mfe = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                            metrics_eq = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "TPR": 0, "PPV": 0, "MCC": 0, "F1": 0}
                            for index in d[timepoint]:
                                occupancy = d[timepoint][index]
                                dist_nat += occupancy * d_ss[index][1]  # occ * dist
                                dist_natc += occupancy * d_ss[index][2]
                                dist_mfe += occupancy * d_ss[index][3]
                                dist_eq += occupancy * d_ss[index][4]
                                metrics_nat = {metric: metrics_nat[metric] + occupancy * d_ss[index][5][metric] for metric in metrics_nat}
                                metrics_natc = {metric: metrics_natc[metric] + occupancy * d_ss[index][6][metric] for metric in metrics_natc}
                                metrics_mfe = {metric: metrics_mfe[metric] + occupancy * d_ss[index][7][metric] for metric in metrics_mfe}
                                metrics_eq = {metric: metrics_eq[metric] + occupancy * d_ss[index][8][metric] for metric in metrics_eq}
                            # normalize
                            sum_occ = sum([d[timepoint][index] for index in d[timepoint]])
                            try:
                                dist_nat /= sum_occ
                                dist_natc /= sum_occ
                                dist_mfe /= sum_occ
                                dist_eq /= sum_occ
                                metrics_nat = {metric: metrics_nat[metric] / sum_occ for metric in metrics_nat}
                                metrics_natc = {metric: metrics_natc[metric] / sum_occ for metric in metrics_natc}
                                metrics_mfe = {metric: metrics_mfe[metric] / sum_occ for metric in metrics_mfe}
                                metrics_eq = {metric: metrics_eq[metric] / sum_occ for metric in metrics_eq}
                            except ZeroDivisionError:
                                dist_nat = float('nan')
                                dist_natc = float('nan')
                                dist_mfe = float('nan')
                                dist_eq = float('nan')
                                metrics_nat = {metric: float('nan') for metric in metrics_nat}
                                metrics_natc = {metric: float('nan') for metric in metrics_natc}
                                metrics_mfe = {metric: float('nan') for metric in metrics_mfe}
                                metrics_eq = {metric: float('nan') for metric in metrics_eq}
                            # add bpdis to metrics dict
                            metrics_nat["bpdis"] = dist_nat
                            metrics_natc["bpdis"] = dist_natc
                            metrics_mfe["bpdis"] = dist_mfe
                            metrics_eq["bpdis"] = dist_eq
                            # store distances in dict
                            timepoint_print = round(timepoint - transcription_duration, 4)
                            if timepoint_print not in db[header]:
                                db[header][timepoint_print] = {"nat": {}, "natc": {}, "mfe": {}, "eq": {}, "nat_rev": {}, "natc_rev": {}, "mfe_rev": {}, "eq_rev": {}}
                            db[header][timepoint_print]["nat" + ['', '_rev'][mode]] = metrics_nat
                            db[header][timepoint_print]["natc" + ['', '_rev'][mode]] = metrics_natc
                            db[header][timepoint_print]["mfe" + ['', '_rev'][mode]] = metrics_mfe
                            db[header][timepoint_print]["eq" + ['', '_rev'][mode]] = metrics_eq

                # write dictionary to json file (each database entry as separate json file, get concatenated later)
                with open(f"analysis_results/{os.path.basename(dataset).split('.')[0]}/{header}.json", "w") as json_f:
                    print(json.dumps(db, indent=4), file=json_f)

            line = data_f.readline()
            ind += 1
            if ind > id2:
                break
    # cleanup
    shutil.rmtree(tmp_dir)

    # print missing files/counts
    print(f"{missing_files=}")
    print(f"{missing_file_count=}")
    print(f"{incomplete_seq_count=}")


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
    parser.add_argument("--id-end", type = int,
            help = """End-Index (1-based) of (total) sequences to compute.""")
    parser.add_argument("--ids-per-job", type = int,
            help = """Maximum number of IDs computed in this job.""")
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

    analyze(args.dataset, args.id1, args.id_end, args.ids_per_job, args.tmp, args.t_list, args.t_ext)

if __name__ == '__main__':
    main()
