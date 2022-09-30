#!/usr/bin/env python

import subprocess
import argparse
import os

def drtransformer(dataset, id, all, o_prune, t_list):
    # change working directory
    os.chdir(f"results/{dataset}")

    t_list_str = " ".join(t_list)
    if all:
        with open(f"{dataset}") as f:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    # Parse
                    header = line[1:].strip()
                    seq = f.readline().strip()
                    nat_ss = f.readline().strip()
                    
                    # DrTransformer
                    subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header} --logfile --o-prune {o_prune} --t-list {t_list_str}", shell=True)
                    subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header}_rev --logfile --o-prune {o_prune} --t-list {t_list_str} --reversed-transcription", shell=True)

                line = f.readline()

    elif id:
        i = 1
        with open(f"{dataset}") as f:
            line = f.readline()
            while line:
                if line.startswith(">"):
                    if i != id:
                        f.readline()
                        f.readline()
                        line = f.readline()
                        i += 1
                        continue
                    # Parse
                    header = line[1:].strip()
                    seq = f.readline().strip()
                    nat_ss = f.readline().strip()
                    
                    # DrTransformer
                    subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header} --logfile --o-prune {o_prune} --t-list {t_list_str}", shell=True)
                    subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header}_rev --logfile --o-prune {o_prune} --t-list {t_list_str} --reversed-transcription", shell=True)

                    break
    else:
        raise ValueError


def main():
    """ 
    Perform forward and reversed DrTransformer simulations.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Perform forward and reversed DrTransformer simulations.')

    parser.add_argument("dataset",
            help = """Path of the dataset file""")
    parser.add_argument("--id", type = int,
            help = """Index (1-based) of sequence to compute.""")
    parser.add_argument("--all", action = "store_true",
            help = """Compute whole database""")
    parser.add_argument("--o-prune", type = float, default = 0.05,
            help = """DrTransformer o-prune parameter""")
    parser.add_argument("--t-list", nargs="+", default = '1 60 600', metavar = '<str>',
            help = """List of times after transcription shown in drf-file. 
            Highest time replaces end time (--t-end) in log-file.
            (--t-end and --t-log are ignored)""")

    args = parser.parse_args()

    drtransformer(args.dataset, args.id, args.all, args.o_prune, args.t_list)
 
if __name__ == '__main__':
    main()
