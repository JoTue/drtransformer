"""
This script contains functions which create plots and rank the sequences from the json-files, created by analyze_db_cli.py.
"""

import json
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns

dirpath = "analysis/slurm_out"
filename = "srp_short"
subdir = "tmRNA"

def NormalSubopt_vs_NormalReversed(dirpath, filename, native="natc", therm_equ="subopt", t="tend", title=None):
    """Creates a 2D-plot which compares normal/reversed DrTransformer prediction and 
    normal DrTransformer prediction with thermodynamic equilibrium predictions.
    Input: 
        dirpath: Path where all json-files are stored (e.g. "analysis/slurm_out").
        filename: Name of RNA class (5s, grp_short, RNaseP_short, ...)
        native: "nat" (native structure from original database) or "natc" (refolded native structure)
        therm_equ: "subopt" (structures in 2 kcal/mol range from MFE) or "mfe"
        t: "tend" or "equ" (timepoint which is used for the plots)
        title: plot title
    Plot: Following distances to the native structure are compared:
        horizontal-axis: normal - reversed, vertical-axis: normal - equilibrium
    """
    with open(f"{dirpath}/{filename}/{filename}.json") as f:
        db = json.load(f)

    xs, ys, size = [], [], []
    # counts how often the data points fall into each quadrant
    quadrant_counts = [0, 0, 0, 0] # bottom-left, top-left, top-right, bottom-right
    for header in db:
        x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
        y = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
        
        # TODO remove
        # y = db[header][f"dist_{native}_{t}"] - db[header][f"dist_natc_equ"]
        # y = db[header][f"dist_{native}_{t}_rev"] - db[header][f"dist_{native}_{therm_equ}"]
        # x = db[header][f"dist_{native}_{t}"]
        # x = db[header][f"dist_{native}_{t}_rev"]
        # y = db[header][f"dist_{native}_{therm_equ}"]

        if x <= 0:
            if y <= 0:
                quadrant_counts[0] += 1
            else:
                quadrant_counts[1] += 1
        else:
            if y <= 0:
                quadrant_counts[3] += 1
            else:
                quadrant_counts[2] += 1
        xs.append(x)
        ys.append(y)
        size.append(db[header][f"dist_{native}_{t}"])
    
    x_mean = np.mean(xs)
    print(f"{filename}: {x_mean=}")
    y_mean = np.mean(ys)
    quadrant_frequencies = [c/len(db) for c in quadrant_counts]
    with open(f"{dirpath}/{filename}/quadrant_frequencies-{native}-{therm_equ}-{t}.txt", "w") as f:
        f.write(str(quadrant_frequencies))

    plt.scatter(xs, ys, s=5) # s=[500/x for x in size], alpha=0.5
    plt.scatter(x_mean, y_mean, c="red")
    #plt.plot([0, 80], [0, 80], c="black")
    plt.hlines(0, xmin=min(-5, min(xs)), xmax=max(5, max(xs)), color = 'black', ls="dashed", linewidth = 1)
    plt.vlines(0, ymin=min(-5, min(ys)), ymax=max(5, max(ys)), color = 'black', ls="dashed", linewidth = 1)
    if not title:
        plt.title(f"{filename}: normal - mfe/subopt vs. normal - reversed\n({native}, {therm_equ}, {t})")
    else:
        plt.title(title)
    plt.xlabel("normal - reversed")
    plt.ylabel("normal - equilibrium")
    plt.savefig(f"{dirpath}/{filename}/nor_sub_vs_nor_rev-{native}-{therm_equ}-{t}.png")
    plt.close()


# for telomerase, tmRNA
def NormalSubopt_vs_NormalReversed_Multi(dirpath, subdir, native="natc", therm_equ="subopt", t="tend", title=None):
    """Creates a 2D-plot which compares normal/reversed DrTransformer prediction and 
    normal DrTransformer prediction with thermodynamic equilibrium predictions.
    Considers multiple json-files which are stored in the input directory path (used for tmRNA, telomerase).
    Input: 
        dirpath: Path where all json-files are stored (e.g. "analysis/slurm_out").
        subdir: Name of RNA class (tmRNA, telomerase)
        native: "nat" (native structure from original database) or "natc" (refolded native structure)
        therm_equ: "subopt" (structures in 2 kcal/mol range from MFE) or "mfe"
        t: "tend" or "equ" (timepoint which is used for the plots)
        title: plot title
    Plot: Following distances to the native structure are compared:
        horizontal-axis: normal - reversed, vertical-axis: normal - equilibrium
    """
    xs, ys, size = [], [], []
    # counts how often the data points fall into each quadrant
    quadrant_counts = [0, 0, 0, 0] # bottom-left, top-left, top-right, bottom-right
    for entry in os.scandir(f"{dirpath}/{subdir}/"):
        if entry.is_dir():
            with open(f"{dirpath}/{subdir}/{entry.name}/{entry.name}.json") as f:
                db = json.load(f)
        else:
            continue

        for header in db:
            x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
            y = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
            # TODO remove
            # y = db[header][f"dist_{native}_{t}_rev"] - db[header][f"dist_{native}_{therm_equ}"]
            if x <= 0:
                if y <= 0:
                    quadrant_counts[0] += 1
                else:
                    quadrant_counts[1] += 1
            else:
                if y <= 0:
                    quadrant_counts[3] += 1
                else:
                    quadrant_counts[2] += 1
            xs.append(x)
            ys.append(y)
            size.append(db[header][f"dist_{native}_{t}"])

    x_mean = np.mean(xs)
    print(f"{subdir}: {x_mean=}")
    y_mean = np.mean(ys)
    quadrant_frequencies = [c/sum(quadrant_counts) for c in quadrant_counts]
    with open(f"{dirpath}/{subdir}/quadrant_frequencies-{native}-{therm_equ}-{t}.txt", "w") as f:
        f.write(str(quadrant_frequencies))

    plt.scatter(xs, ys, s=10) # s=[500/x for x in size], alpha=0.5
    plt.scatter(x_mean, y_mean, c="red")
    plt.hlines(0, xmin=min(-5, min(xs)), xmax=max(5, max(xs)), color = 'black', ls="dashed", linewidth = 1)
    plt.vlines(0, ymin=min(-5, min(ys)), ymax=max(5, max(ys)), color = 'black', ls="dashed", linewidth = 1)
    if not title:
        plt.title(f"{subdir}: normal - mfe/subopt vs. normal - reversed\n({native}, {therm_equ}, {t})")
    else:
        plt.title(title)
    plt.xlabel("normal - reversed")
    plt.ylabel("normal - equilibrium")
    plt.savefig(f"{dirpath}/{subdir}/nor_sub_vs_nor_rev-{native}-{therm_equ}-{t}.png")
    plt.close()


def NormalSubopt_vs_NormalReversed_16S(dirpath="analysis/slurm_out", subdir="16s", native="natc", therm_equ="subopt", t="tend", title=None):
    """Creates a 2D-plot which compares normal/reversed DrTransformer prediction and 
    normal DrTransformer prediction with thermodynamic equilibrium predictions.
    For 16S rRNA (each domain is plotted in a separate color).
    Input: 
        dirpath: Path where all json-files are stored (e.g. "analysis/slurm_out").
        subdir: Directory of 16S rRNA results.
        native: "nat" (native structure from original database) or "natc" (refolded native structure)
        therm_equ: "subopt" (structures in 2 kcal/mol range from MFE) or "mfe"
        t: "tend" or "equ" (timepoint which is used for the plots)
        title: plot title
    Plot: Following distances to the native structure are compared:
        horizontal-axis: normal - reversed, vertical-axis: normal - equilibrium
    """
    xs = {"domain1": [], "domain2": [], "domain3": [], "domain4": []}
    ys = {"domain1": [], "domain2": [], "domain3": [], "domain4": []}
    size = {"domain1": [], "domain2": [], "domain3": [], "domain4": []}

    for entry in os.scandir(f"{dirpath}/{subdir}/"):
        if entry.is_dir():
            with open(f"{dirpath}/{subdir}/{entry.name}/{entry.name}.json") as f:
                db = json.load(f)
        else:
            continue

        for header in db:
            domain = header.split("_")[-1]
            x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
            y = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
            # TODO remove
            # y = db[header][f"dist_{native}_{t}_rev"] - db[header][f"dist_{native}_{therm_equ}"]
            
            xs[domain].append(x)
            ys[domain].append(y)
            size[domain].append(db[header][f"dist_{native}_{t}"])

    x_mean = []
    for vals in xs.values():
        x_mean.extend(vals)
    x_mean = np.mean(x_mean)
    print(f"{subdir}: {x_mean=}")
    y_mean = []
    for vals in ys.values():
        y_mean.extend(vals)
    y_mean = np.mean(y_mean)

    for domain in ["domain1", "domain2", "domain3", "domain4"]:
        plt.scatter(xs[domain], ys[domain], s=10, label=domain) # s=[500/x for x in size], alpha=0.5
    plt.scatter(x_mean, y_mean, c="red")
    plt.hlines(y=0, xmin=min(-5, min([min(xs[domain]) for domain in xs])), xmax=max(5, max([max(xs[domain]) for domain in xs])), color = 'black', ls="dashed", linewidth = 1)
    plt.vlines(x=0, ymin=min(-5, min([min(ys[domain]) for domain in ys])), ymax=max(5, max([max(ys[domain]) for domain in ys])), color = 'black', ls="dashed", linewidth = 1)
    if not title:
        plt.title(f"16S: normal - mfe/subopt vs. normal - reversed\n({native}, {therm_equ}, {t})")
    else:
        plt.title(title)
    plt.xlabel("normal - reversed")
    plt.ylabel("normal - equilibrium")
    lgnd = plt.legend()
    for handle in lgnd.legendHandles:
        handle.set_sizes([20])
    plt.savefig(f"{dirpath}/{subdir}/nor_sub_vs_nor_rev-{native}-{therm_equ}-{t}.png")
    plt.close()


def NormalSubopt_vs_NormalReversed_23S(dirpath="analysis/slurm_out", subdir="23s", native="natc", therm_equ="subopt", t="tend", title=None):
    """Creates a 2D-plot which compares normal/reversed DrTransformer prediction and 
    normal DrTransformer prediction with thermodynamic equilibrium predictions.
    For 23S rRNA (each domain is plotted in a separate color).
    Input: 
        dirpath: Path where all json-files are stored (e.g. "analysis/slurm_out").
        subdir: Directory of 23S rRNA results.
        native: "nat" (native structure from original database) or "natc" (refolded native structure)
        therm_equ: "subopt" (structures in 2 kcal/mol range from MFE) or "mfe"
        t: "tend" or "equ" (timepoint which is used for the plots)
        title: plot title
    Plot: Following distances to the native structure are compared:
        horizontal-axis: normal - reversed, vertical-axis: normal - equilibrium
    """
    xs = {"domain1": [], "domain2": [], "domain3": [], "domain4": [], "domain5": [], "domain6": []}
    ys = {"domain1": [], "domain2": [], "domain3": [], "domain4": [], "domain5": [], "domain6": []}
    size = {"domain1": [], "domain2": [], "domain3": [], "domain4": [], "domain5": [], "domain6": []}

    for entry in os.scandir(f"{dirpath}/{subdir}/"):
        if entry.is_dir():
            with open(f"{dirpath}/{subdir}/{entry.name}/{entry.name}.json") as f:
                db = json.load(f)
        else:
            continue

        for header in db:
            domain = header.split("_")[-1]
            x = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{t}_rev"]
            y = db[header][f"dist_{native}_{t}"] - db[header][f"dist_{native}_{therm_equ}"]
            # TODO remove
            # y = db[header][f"dist_{native}_{t}_rev"] - db[header][f"dist_{native}_{therm_equ}"]
            
            xs[domain].append(x)
            ys[domain].append(y)
            size[domain].append(db[header][f"dist_{native}_{t}"])

    x_mean = []
    for vals in xs.values():
        x_mean.extend(vals)
    x_mean = np.mean(x_mean)
    print(f"{subdir}: {x_mean=}")
    y_mean = []
    for vals in ys.values():
        y_mean.extend(vals)
    y_mean = np.mean(y_mean)

    for domain in ["domain1", "domain2", "domain3", "domain4", "domain5", "domain6"]:
        plt.scatter(xs[domain], ys[domain], s=10, label=domain) # s=[500/x for x in size], alpha=0.5
    plt.hlines(y=0, xmin=min(-5, min([min(xs[domain]) for domain in xs])), xmax=max(5, max([max(xs[domain]) for domain in xs])), color = 'black', ls="dashed", linewidth = 1)
    plt.vlines(x=0, ymin=min(-5, min([min(ys[domain]) for domain in ys])), ymax=max(5, max([max(ys[domain]) for domain in ys])), color = 'black', ls="dashed", linewidth = 1)
    plt.scatter(x_mean, y_mean, c="red")
    if not title:
        plt.title(f"23S: normal - mfe/subopt vs. normal - reversed\n({native}, {therm_equ}, {t})")
    else:
        plt.title(title)
    plt.xlabel("normal - reversed")
    plt.ylabel("normal - equilibrium")
    lgnd = plt.legend()
    for handle in lgnd.legendHandles:
        handle.set_sizes([20])
    plt.savefig(f"{dirpath}/{subdir}/nor_sub_vs_nor_rev-{native}-{therm_equ}-{t}.png")
    plt.close()


def RankStructures(dirpath, filename):
    """Assign a score to each RNA.
    RNAs which fold better into the native structure in the normal transcription direction 
    compared to the reversed direction and thermodynamic equilibrium predictions, get higher scores.

    Score function: (rev-nor) + (sub-nor) - nor = rev + sub - 3*nor  
    """
    with open(f"{dirpath}/{filename}/{filename}.json") as f:
        db = json.load(f)

    score_dict = {}

    lengths = []
    with open(f"{dirpath}/{filename}/{filename}_scores.csv", "w") as f:
        for header in db:
            score = db[header]["dist_natc_tend_rev"] + db[header]["dist_natc_subopt"] - 3*db[header]["dist_natc_tend"]
            score_dict[header] = {"score": score, "dist_natc_tend": db[header]["dist_natc_tend"], "dist_natc_tend_rev": db[header]["dist_natc_tend_rev"], "dist_natc_subopt": db[header]["dist_natc_subopt"]}
            # length
            lengths.append(db[header]["seqlen"])

        sorted_scores = sorted(score_dict.items(), key=lambda x: x[1]["score"], reverse=True)
        
        print("header,score,dist_natc_tend,dist_natc_tend_rev,dist_natc_subopt", file=f)
        for item in sorted_scores:
            print(f"{item[0]},{item[1]['score']:.2f},{item[1]['dist_natc_tend']:.2f},{item[1]['dist_natc_tend_rev']:.2f},{item[1]['dist_natc_subopt']:.2f}", file=f)
    return np.mean(lengths)


def RankStructures_Multi(dirpath, subdir):
    """Assign a score to each RNA - considers multiple json-files which are stored in the input directory path (used for 16S, 23S, tmRNA, telomerase).
    RNAs which fold better into the native structure in the normal transcription direction 
    compared to the reversed direction and thermodynamic equilibrium predictions, get higher scores.

    Score function: (rev-nor) + (sub-nor) - nor = rev + sub - 3*nor  
    """
    score_dict = {}
    lengths = []
    with open(f"{dirpath}/{subdir}/{subdir}_scores.csv", "w") as f:
        print("header,score,dist_natc_tend,dist_natc_tend_rev,dist_natc_subopt", file=f)
        for entry in os.scandir(f"{dirpath}/{subdir}/"):
            if entry.is_dir():
                with open(f"{dirpath}/{subdir}/{entry.name}/{entry.name}.json") as json_f:
                    db = json.load(json_f)
            else:
                continue

            for header in db:
                score = db[header]["dist_natc_tend_rev"] + db[header]["dist_natc_subopt"] - 3*db[header]["dist_natc_tend"]
                score_dict[header] = {"score": score, "dist_natc_tend": db[header]["dist_natc_tend"], "dist_natc_tend_rev": db[header]["dist_natc_tend_rev"], "dist_natc_subopt": db[header]["dist_natc_subopt"]}
                # length
                lengths.append(db[header]["seqlen"])

        sorted_scores = sorted(score_dict.items(), key=lambda x: x[1]["score"], reverse=True)
            
        for item in sorted_scores:
            print(f"{item[0]},{item[1]['score']:.2f},{item[1]['dist_natc_tend']:.2f},{item[1]['dist_natc_tend_rev']:.2f},{item[1]['dist_natc_subopt']:.2f}", file=f)
    return np.mean(lengths)


# Function calls:

for filename in ["5s", "grp_short", "RNaseP_short", "srp_middle", "srp_short", "tRNA", "5s_param", "srp_short_param", "tRNA_param"]:
    NormalSubopt_vs_NormalReversed(dirpath, filename, "natc", "subopt", "tend")
    NormalSubopt_vs_NormalReversed(dirpath, filename, "nat", "subopt", "tend")
    NormalSubopt_vs_NormalReversed(dirpath, filename, "natc", "mfe", "tend")
    NormalSubopt_vs_NormalReversed(dirpath, filename, "natc", "subopt", "equ")
    mean_length = RankStructures(dirpath, filename)
    # print(filename, mean_length)

for subdir in ["telomerase", "tmRNA"]:
    NormalSubopt_vs_NormalReversed_Multi(dirpath, subdir, "natc", "subopt", "tend")
    NormalSubopt_vs_NormalReversed_Multi(dirpath, subdir, "nat", "subopt", "tend")
    NormalSubopt_vs_NormalReversed_Multi(dirpath, subdir, "natc", "mfe", "tend")
    NormalSubopt_vs_NormalReversed_Multi(dirpath, subdir, "natc", "subopt", "equ")
    mean_length = RankStructures_Multi(dirpath, subdir)
    # print(subdir, mean_length)

# 16S
NormalSubopt_vs_NormalReversed_16S("analysis/slurm_out", "16s", "natc", "subopt", "tend", title="16S rRNA")
NormalSubopt_vs_NormalReversed_16S("analysis/slurm_out", "16s", "nat", "subopt", "tend")
NormalSubopt_vs_NormalReversed_16S("analysis/slurm_out", "16s", "natc", "mfe", "tend")
NormalSubopt_vs_NormalReversed_16S("analysis/slurm_out", "16s", "natc", "subopt", "equ")
mean_length = RankStructures_Multi(dirpath, "16s")
# print("16s", mean_length)

# 23S
NormalSubopt_vs_NormalReversed_23S("analysis/slurm_out", "23s", "natc", "subopt", "tend", title="23S rRNA")
NormalSubopt_vs_NormalReversed_23S("analysis/slurm_out", "23s", "nat", "subopt", "tend")
NormalSubopt_vs_NormalReversed_23S("analysis/slurm_out", "23s", "natc", "mfe", "tend")
NormalSubopt_vs_NormalReversed_23S("analysis/slurm_out", "23s", "natc", "subopt", "equ")
mean_length = RankStructures_Multi(dirpath, "23s")
# print("23s", mean_length)

# redo plots with new titles
NormalSubopt_vs_NormalReversed(dirpath, "5s", title="5S rRNA")
NormalSubopt_vs_NormalReversed(dirpath, "grp_short", title="Group I intron (≤300bp)")
NormalSubopt_vs_NormalReversed(dirpath, "RNaseP_short", title="RNaseP (≤250bp)")
NormalSubopt_vs_NormalReversed(dirpath, "srp_middle", title="SRP (201-300bp)")
NormalSubopt_vs_NormalReversed(dirpath, "srp_short", title="SRP (≤200bp)")
NormalSubopt_vs_NormalReversed(dirpath, "tRNA", title="tRNA")
NormalSubopt_vs_NormalReversed(dirpath, "5s_param", title="5S rRNA\n(higher o-prune/t-fast)")
NormalSubopt_vs_NormalReversed(dirpath, "srp_short_param", title="SRP (≤200bp)\n(higher o-prune/t-fast)")
NormalSubopt_vs_NormalReversed(dirpath, "tRNA_param", title="tRNA\n(higher o-prune/t-fast)")
NormalSubopt_vs_NormalReversed_Multi(dirpath, "telomerase", title="Telomerase")
NormalSubopt_vs_NormalReversed_Multi(dirpath, "tmRNA", title="tmRNA")