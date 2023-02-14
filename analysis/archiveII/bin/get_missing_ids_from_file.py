#!/usr/bin/env python3

import RNA
import argparse
import os

FILE_NAME = "missing_ids.txt"

def get_missing_ids(rna_class,  id1, id2):
    with open(FILE_NAME) as f:
        for line in f:
            if line.startswith(f"{rna_class}"):
                total_missing_ids = int(f.readline().strip())
                # assert that both id1 and id2 are in the correct range
                assert id1 <= total_missing_ids, f"{id1=}, {total_missing_ids=}"
                assert id2 <= total_missing_ids, f"{id2=}, {total_missing_ids=}"
                break

        for line in f:
            if line.startswith("Missing_ids_per_class:"):
                break
        for line in f:
            if line.startswith(f"{rna_class}"):
                i = 1
                while True:
                    if i < id1:
                        f.readline()
                        i += 1
                        continue
                    id = f.readline().strip()
                    print(id)
                    i += 1
                    if i > id2:
                        break


def main():
    """Get IDs which have not yet been computed with DrTransformer (uneven=normal, even=reversed, eg. ID 2 -> 3=normal, 4=reversed).
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Get IDs which have not yet been computed with DrTransformer (uneven=normal, even=reversed, eg. ID 2 -> 3=normal, 4=reversed).')

    parser.add_argument("rna_class",
            help = """Name of rna_class (e.g. 'tRNA').""")
    parser.add_argument("--id1", type = int,
            help = """Start-Index (1-based) of sequences to compute.""")
    parser.add_argument("--id2", type = int,
            help = """End-Index (1-based) of sequences to compute.""")

    args = parser.parse_args()

    get_missing_ids(args.rna_class, args.id1, args.id2)
 
if __name__ == '__main__':
    main()
