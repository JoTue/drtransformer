import json
import matplotlib.pyplot as plt
import os
import numpy as np

def LengthComparison(norm_by_length=False, mode=0, dirpath="analysis/slurm_out", native="natc", therm_equ="subopt", t="tend"):
    """Compare the folding of RNAs based on their sequence length.

    Input:
        norm_by_length: Boolean which determines whether the base pair distances are divided by the sequence length
        mode: 0 (compare normal/reversed DrTransformer), 1 (compare normal DrTransformer/thermodynamic equilibrium)
        dirpath: Path where all json-files are stored (e.g. "analysis/slurm_out").
        native: "nat" (native structure from original database) or "natc" (refolded native structure)
        therm_equ: "subopt" (structures in 2 kcal/mol range from MFE) or "mfe"
        t: "tend" or "equ" (timepoint which is used for the plots)
    """
    xs, ys = [], []
    short_xs, short_ys = [], []
    long_xs, long_ys = [], []
    for filename in ["5s", "grp_short", "RNaseP_short", "srp_middle", "srp_short", "tRNA"]:
        with open(f"{dirpath}/{filename}/{filename}.json") as f:
            db = json.load(f)
        for header in db:
            if mode == 0:
                x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
            elif mode == 1:
                x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
            if norm_by_length:
                x /= db[header]["seqlen"]
            y = db[header]["seqlen"]
            xs.append(x)
            ys.append(y)
            if filename in ["5s", "grp_short", "RNaseP_short", "srp_short", "tRNA"]:
                short_xs.append(x)
                short_ys.append(y)
            else:
                long_xs.append(x)
                long_ys.append(y)
    
    for subdir in ["telomerase", "tmRNA", "16s", "23s"]:
        for entry in os.scandir(f"{dirpath}/{subdir}/"):
            if entry.is_dir():
                with open(f"{dirpath}/{subdir}/{entry.name}/{entry.name}.json") as f:
                    db = json.load(f)
            else:
                continue
            for header in db:
                if mode == 0:
                    x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
                elif mode == 1:
                    x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
                if norm_by_length:
                    x /= db[header]["seqlen"]
                y = db[header]["seqlen"]
                xs.append(x)
                ys.append(y)
                long_xs.append(x)
                long_ys.append(y)

    x_mean = np.mean(xs)
    y_mean = np.mean(ys)
    print(f"{x_mean=}, {y_mean=}")
    short_x_mean = np.mean(short_xs)
    short_y_mean = np.mean(short_ys)
    print(f"{short_x_mean=}, {short_y_mean=}")
    long_x_mean = np.mean(long_xs)
    long_y_mean = np.mean(long_ys)
    print(f"{long_x_mean=}, {long_y_mean=}")
    
    # plot all
    plt.scatter(xs, ys, s=1)
    plt.scatter(x_mean, y_mean, c="red")
    plt.vlines(0, ymin=0, ymax=max(ys), color = 'black', ls="dashed", linewidth = 1)
    plt.title(f"Length comparison")
    plt.xlabel(f"{['normal - reversed', 'normal - equilibrium'][mode]}")
    plt.ylabel("Length")
    plt.savefig(f"analysis/length_comparison/length_comparison_{['nor-rev', 'nor-equ'][mode]}_{['unnormalized', 'normalized'][norm_by_length]}.png")
    plt.close()
    # plot short/long
    plt.scatter(short_xs, short_ys, s=1, label="a")
    plt.scatter(long_xs, long_ys, s=1, label="b")
    plt.scatter(short_x_mean, short_y_mean, c="darkblue", label="mean a")
    plt.scatter(long_x_mean, long_y_mean, c="orangered", label="mean b")
    plt.vlines(0, ymin=0, ymax=max(long_ys), color = 'black', ls="dashed", linewidth = 1)
    plt.title(f"Length comparison - colored by parameters")
    plt.xlabel(f"{['normal - reversed', 'normal - equilibrium'][mode]}")
    plt.ylabel("Length")
    plt.legend(title="Parameters:")
    plt.savefig(f"analysis/length_comparison/length_comparison_{['nor-rev', 'nor-equ'][mode]}_{['unnormalized', 'normalized'][norm_by_length]}_param.png")
    plt.close()

LengthComparison(norm_by_length=False, mode=0)
LengthComparison(norm_by_length=True, mode=0)
LengthComparison(norm_by_length=False, mode=1)
LengthComparison(norm_by_length=True, mode=1)
