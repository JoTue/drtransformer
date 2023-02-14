#!/usr/bin/env python3

import RNA
import argparse
import os


def get_missing_ids(seq_file, results_dir, id1, id2):
    # missing_ids = []
    # get already existing result files
    files = os.listdir(results_dir)
    # parse seqfile
    with open(seq_file) as f:
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

                # check if results already exist
                fn_normal = f"drtr_{header}.log"
                fn_reversed = f"drtr_{header}_rev.log"
                if not fn_normal in files:
                    # missing_ids.append(f"{i}_n")
                    print(2*i - 1)
                if not fn_reversed in files:
                    # missing_ids.append(f"{i}_r")
                    print(2*i)
            line = f.readline()
            i += 1
            if i > id2:
                break
    # return missing_ids


def main():
    """Get IDs which have not yet been computed with DrTransformer (uneven=normal, even=reversed, eg. ID 2 -> 3=normal, 4=reversed).
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Get IDs which have not yet been computed with DrTransformer (uneven=normal, even=reversed, eg. ID 2 -> 3=normal, 4=reversed).')

    parser.add_argument("seq_file",
            help = """Path of the sequence file.""")
    parser.add_argument("results_dir",
            help = """Path of the results directory.""")
    parser.add_argument("--id1", type = int,
            help = """Start-Index (1-based) of sequences to compute.""")
    parser.add_argument("--id2", type = int,
            help = """End-Index (1-based) of sequences to compute.""")

    args = parser.parse_args()

    get_missing_ids(args.seq_file, args.results_dir, args.id1, args.id2)
 
if __name__ == '__main__':
    main()
