#!/usr/bin/env python3

import argparse
import json
import statistics
import matplotlib.pyplot as plt
from scipy.stats import shapiro, kstest, ttest_ind, wilcoxon
import numpy as np
import plotly.graph_objects as go
import math
import scipy.stats as stat

TIMESTAMPS = ["0.04", "1.0", "10.0", "60.0", "600.0", "3600.0", "equilibrium"] #TODO: "eq"
STRUCTURES = ["nat", "natc", "mfe", "eq", "nat_rev", "natc_rev", "mfe_rev", "eq_rev"]
MEASUREMENTS = ["TP", "FP", "FN", "TN", "TPR", "PPV", "MCC", "F1", "bpdis"]

def stats(input_file):
    # convert json-file into dict
    with open(input_file) as f:
        db = json.load(f)

    d = {ts: {st: {ms: {"values": [], "mean": None, "stdev": None, "median": None} for ms in MEASUREMENTS} for st in STRUCTURES} for ts in TIMESTAMPS}  # stores statistics

    for header in db:
        for ts in TIMESTAMPS:
            for st in STRUCTURES:
                for ms in MEASUREMENTS:
                    d[ts][st][ms]["values"].append(float(db[header][ts][st][ms]))
    print("FINISHED writing to dict.")

    for ts in TIMESTAMPS:
        for st in STRUCTURES:
            for ms in MEASUREMENTS:
                d[ts][st][ms]["mean"] = statistics.mean(d[ts][st][ms]["values"])
                d[ts][st][ms]["stdev"] = statistics.stdev(d[ts][st][ms]["values"])
                d[ts][st][ms]["median"] = statistics.median(d[ts][st][ms]["values"])

    ts = "10.0"
    ms = "MCC"

    print(shapiro([x-y for x,y in zip(d[ts]["nat"][ms]["values"], d[ts]["nat_rev"][ms]["values"])]))

    print("two-sided")
    kstest_d = {"MCC": {}, "bpdis": {}}
    ks_array = np.zeros((len(TIMESTAMPS), 3), dtype=float)
    for i, ts in enumerate(TIMESTAMPS):
        for ms in ["MCC", "bpdis"]:
            out = wilcoxon(d[ts]["nat"][ms]["values"], d[ts]["nat_rev"][ms]["values"])
            kstest_d[ms][ts] = '{:.2e}'.format(out[1])
            if ms=="MCC":
                ks_array[i][0] = out[1]
    for ms in kstest_d:
        print(ms)
        for ts in kstest_d[ms]:
            print(f"{ts}\t{kstest_d[ms][ts]}")
    
    print("less")
    kstest_d = {"MCC": {}, "bpdis": {}}
    for i, ts in enumerate(TIMESTAMPS):
        for ms in ["MCC", "bpdis"]:
            out = wilcoxon(d[ts]["nat"][ms]["values"], d[ts]["nat_rev"][ms]["values"], alternative="less")
            kstest_d[ms][ts] = '{:.2e}'.format(out[1])
            if ms=="MCC":
                ks_array[i][1] = out[1]
    for ms in kstest_d:
        print(ms)
        for ts in kstest_d[ms]:
            print(f"{ts}\t{kstest_d[ms][ts]}")

    print("greater")
    kstest_d = {"MCC": {}, "bpdis": {}}
    for i, ts in enumerate(TIMESTAMPS):
        for ms in ["MCC", "bpdis"]:
            out = wilcoxon(d[ts]["nat"][ms]["values"], d[ts]["nat_rev"][ms]["values"], alternative="greater")
            kstest_d[ms][ts] = '{:.2e}'.format(out[1])
            if ms=="MCC":
                ks_array[i][2] = out[1]
    for ms in kstest_d:
        print(ms)
        for ts in kstest_d[ms]:
            print(f"{ts}\t{kstest_d[ms][ts]}")

    # np.set_printoptions(formatter={'float': '   {:.2e}   '.format})
    # print(ks_array)
    
    # plt.hist(t1, bins=50, alpha=0.5, label="normal")
    # plt.hist(t2, bins=50, alpha=0.5, label="reversed")
    # plt.xlabel(ms)
    # plt.ylabel("Frequency")
    # plt.legend()
    # plt.show()


def main():
    """Statistics on json-files.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Analysis of a database in dot-bracket format.')

    parser.add_argument("input", metavar = '<str>',
            help = """Path of the input file""")

    args = parser.parse_args()

    stats(args.input)
    # for inp in ["5s", "RNaseP", "srp", "tmRNA", "tRNA"]:
    #     print(inp)
    #     stats(f"analysis_results/{inp}/{inp}.JSON")
 
if __name__ == '__main__':
    main()
