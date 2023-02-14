#!/usr/bin/env python

import subprocess
import argparse
import os

def drtransformer(dataset, id, id1, id_end, ids_per_job, o_prune, t_ext, t_list, print_only_id):
    # change working directory
    os.chdir(f"drtr_tmp/{os.path.basename(dataset).split('.')[0]}")

    t_list_str = " ".join(t_list)

    if not id:  # compute id1-id2, normal and reveresed at once
        id2 = min([id1 + ids_per_job - 1, id_end])
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
                    if os.path.isfile(f"../../drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}.log"):
                        print(f"Skipping id {i} (normal): drtr_{header} results already exist")
                    else:
                        if print_only_id:
                            print(f"{i=}\t{header}")
                        else:
                            subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header} --logfile --o-prune {o_prune} --t-ext {t_ext} --t-list {t_list_str} && mv drtr_{header}.* ../../drtr_results/{os.path.basename(dataset).split('.')[0]}", shell=True)
                    if os.path.isfile(f"../../drtr_results/{os.path.basename(dataset).split('.')[0]}/drtr_{header}_rev.log"):
                        print(f"Skipping id {i} (reversed): drtr_{header}_rev results already exist")
                    else:
                        if print_only_id:
                            print(f"{i=}\t{header}")
                        else:
                            subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header}_rev --logfile --o-prune {o_prune} --t-list {t_list_str} --reversed-transcription && mv drtr_{header}_rev.* ../../drtr_results/{os.path.basename(dataset).split('.')[0]}", shell=True)
                line = f.readline()
                i += 1
                if i > id2:
                    break
    else:  # compute only one id (uneven=normal, even=reversed, eg. ID 2 -> id 3=normal, 4=reversed)
        ID = (id+1) // 2
        with open(f"../../{dataset}") as f:
            i = 1
            line = f.readline()
            while line:
                if line.startswith(">"):
                    if i < ID:
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
                    if id%2!=0:
                        if print_only_id:
                            print(f"{ID=}\t{id=}\t{header}")
                        else:
                            subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header} --logfile --o-prune {o_prune} --t-ext {t_ext} --t-list {t_list_str} && mv drtr_{header}.* ../../drtr_results/{os.path.basename(dataset).split('.')[0]}", shell=True)
                    else:
                        if print_only_id:
                            print(f"{ID=}\t{id=}\t{header}")
                        else:
                            subprocess.run(f"echo {seq} | DrTransformer --name drtr_{header}_rev --logfile --o-prune {o_prune} --t-list {t_list_str} --reversed-transcription && mv drtr_{header}_rev.* ../../drtr_results/{os.path.basename(dataset).split('.')[0]}", shell=True)
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
    parser.add_argument("--id", type = int,
            help = """If specified, compute this id (uneven=normal, even=reversed, eg. ID 2 -> id 3=normal, 4=reversed)""")
    parser.add_argument("--id1", type = int,
            help = """Start-Index (1-based) of sequences to compute.""")
    parser.add_argument("--id-end", type = int,
            help = """End-Index (1-based) of (total) sequences to compute.""")
    parser.add_argument("--ids-per-job", type = int,
            help = """Maximum number of IDs computed in this job.""")
    parser.add_argument("--o-prune", type = float, default = 0.05,
            help = """DrTransformer o-prune parameter""")
    parser.add_argument("--t-ext", type = float, default = 0.04, metavar = '<flt>',
            help = """Inverse of transcription rate, i.e. time per nucleotide extension
            [seconds per nucleotide].""")
    parser.add_argument("--t-list", nargs="+", default = ['1', '10', '60', '600', '3600'], metavar = '<str>',
            help = """List of times after transcription shown in drf-file.
            Highest time replaces end time (--t-end) in log-file.
            (--t-end and --t-log are ignored)""")
    parser.add_argument("--print-only-id",  type = bool, default = False,
            help = """If True, no DrTransformer computation is done, only prints id´s that would be computed."""
    )

    args = parser.parse_args()

    drtransformer(args.dataset, args.id, args.id1, args.id_end, args.ids_per_job, args.o_prune, args.t_ext, args.t_list, args.print_only_id)

if __name__ == '__main__':
    main()
    