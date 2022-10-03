#!/usr/bin/env python

import subprocess
import argparse
import os

def drtransformer(dataset, id1, id2, o_prune, t_ext, t_list):
    # change working directory
    os.chdir(f"drtr_tmp/{os.path.basename(dataset).split('.')[0]}")

    t_list_str = " ".join(t_list)

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
                header = line[1:].strip()
                seq = f.readline().strip()
                nat_ss = f.readline().strip()

                # DrTransformer
                subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header} --logfile --o-prune {o_prune} --t-ext {t_ext} --t-list {t_list_str} && mv drtr_{header}.log drtr_results/drtr_{header}.log", shell=True)
                subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header}_rev --logfile --o-prune {o_prune} --t-list {t_list_str} --reversed-transcription && mv drtr_{header}_rev.log drtr_results/drtr_{header}_rev.log", shell=True)
            line = f.readline()
            i += 1
            if i > id2:
                break


def main():
    """
    Perform forward and reversed DrTransformer simulations.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Perform forward and reversed DrTransformer simulations.')

    parser.add_argument("dataset",
            help = """Path of the dataset file""")
    parser.add_argument("--id1", type = int,
            help = """Start-Index (1-based) of sequences to compute.""")
    parser.add_argument("--id2", type = int,
            help = """End-Index (1-based) of sequences to compute.""")
    parser.add_argument("--o-prune", type = float, default = 0.05,
            help = """DrTransformer o-prune parameter""")
    parser.add_argument("--t-ext", type = float, default = 0.04, metavar = '<flt>',
            help = """Inverse of transcription rate, i.e. time per nucleotide extension
            [seconds per nucleotide].""")
    parser.add_argument("--t-list", nargs="+", default = ['1', '10', '60', '600', '3600'], metavar = '<str>',
            help = """List of times after transcription shown in drf-file.
            Highest time replaces end time (--t-end) in log-file.
            (--t-end and --t-log are ignored)""")

    args = parser.parse_args()

    drtransformer(args.dataset, args.id1, args.id2, args.o_prune, args.t_ext, args.t_list)

if __name__ == '__main__':
    main()