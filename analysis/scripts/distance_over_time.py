import subprocess
import RNA
import os
import numpy as np
import matplotlib.pyplot as plt

def DistanceOverTimePlot(fasta, t8=1e10, title=None):
    """Takes one sequence-structure pair as an input and creates a plot of the distance to the native structure 
    of the normal/reversed transcription direction after transcription"""
    # parse fasta
    lines = fasta.split("\n")
    header = lines[0][1:]
    seq = lines[1]
    nat_ss = lines[2]
    seq_len = len(seq)

    # change working directory
    try:
        os.mkdir(f"analysis/distance_over_time/{header}")
    except FileExistsError:
        pass
    os.chdir(f"analysis/distance_over_time/{header}")

    # Constraint folding
    # RNAfold
    subprocess.run(f"RNAfold --noPS --canonicalBPonly -C > rnafold_c", text=True, input=f"{seq}\n{nat_ss}\n@\n", shell=True)
    with open(f"rnafold_c") as f:
        lines = f.readlines()
        line_split = lines[-1].split()
        natc_ss = line_split[0]
        # natc = float(line_split[-1].strip("()"))

    # DrTransformer
    subprocess.run(f"echo {seq} | DrTransformer --name drtr_normal --logfile --tmpdir tmp_normal --t-end {t8}", shell=True)
    subprocess.run(f"echo {seq} | DrTransformer --name drtr_rev --logfile --tmpdir tmp_rev -rt --t-end {t8}", shell=True)

    # plot .drf
    subprocess.run(f"cat drtr_normal.drf | DrPlotter --name drpl_normal --format png", shell=True)
    subprocess.run(f"cat drtr_rev.drf | DrPlotter --name drpl_rev --format png", shell=True)

    # retrieve p0 values
    # normal
    with open(f"tmp_normal/drtr_normal-{seq_len}_p0.txt") as f:
        line = f.readlines()[0]
        start_index = line.index("--p0")
        p0s_normal = line[start_index:].replace("--p0", "")
    # reversed
    with open(f"tmp_rev/drtr_rev-{seq_len}_p0.txt") as f:
        line = f.readlines()[0]
        start_index = line.index("--p0")
        p0s_rev = line[start_index:].replace("--p0", "")

    # DrSimulate
    with open("drsim_normal", "w") as outfile:
        subprocess.run(f"DrSimulate --t8 {t8} -r tmp_normal/drtr_normal-{seq_len}_rates.txt --transpose --p0 {p0s_normal}", shell=True, stdout=outfile)
    with open("drsim_rev", "w") as outfile:
        subprocess.run(f"DrSimulate --t8 {t8} -r tmp_rev/drtr_rev-{seq_len}_rates.txt --transpose --p0 {p0s_rev}", shell=True, stdout=outfile)

    # TODO: error if only one structure at transcription end

    # Get structures at end of transcription
    with open("drtr_normal.log") as f:
        lines = f.readlines()
        d_normal = {}
        for line in lines[::-1]:
            words = line.split()
            if words[0] != str(seq_len):
                break
            d_normal[int(words[1])] = words[2] # format: {index: ss}
    with open("drtr_rev.log") as f:
        lines = f.readlines()
        d_rev = {}
        for line in lines[::-1]:
            words = line.split()
            if words[0] != str(seq_len):
                break
            d_rev[int(words[1])] = words[2] # format: {index: ss}

    # Compute distances
    # normal
    dist_normal = {}
    for index in d_normal:
        dist = RNA.bp_distance(d_normal[index], natc_ss)
        dist_normal[index] = dist
    with open("drsim_normal") as f:
        lines = f.readlines()
    words = [line.split() for line in lines]
    arr_normal = np.array(words, dtype=float)
    for index in range(1, len(d_normal)+1):
        arr_normal[:, index] *= dist_normal[index]
    z = np.zeros((len(arr_normal),1), dtype=float)
    np.append(arr_normal, z, axis=1)
    for row_i in range(len(arr_normal)):
        arr_normal[row_i, -1] = np.sum(arr_normal[row_i, 1: -1])
    # rev
    dist_rev = {}
    for index in d_rev:
        dist = RNA.bp_distance(d_rev[index], natc_ss)
        dist_rev[index] = dist
    with open("drsim_rev") as f:
        lines = f.readlines()
    words = [line.split() for line in lines]
    arr_rev = np.array(words, dtype=float)
    for index in range(1, len(d_rev)+1):
        arr_rev[:, index] *= dist_rev[index]
    z = np.zeros((len(arr_rev),1), dtype=float)
    np.append(arr_rev, z, axis=1)
    for row_i in range(len(arr_rev)):
        arr_rev[row_i, -1] = np.sum(arr_rev[row_i, 1: -1])

    # plot
    plt.plot(arr_normal[:, 0], arr_normal[:, -1], label="normal")
    plt.plot(arr_rev[:, 0], arr_rev[:, -1], label="reversed")
    plt.xscale("log")
    plt.xlabel("time [sec]")
    plt.ylabel("Distance to native structure")
    if not title:
        plt.title(f"Distances to native structure over time\n{header}")
    else:
        plt.title(f"Distances to native structure over time\n{title}")
    plt.legend()
    plt.savefig("dot_plot")
    plt.close()



fasta = """>5s_Perinereis-brevicirris-1
GCCUACGGCCAUACUACGUUGAAAACACCGGUUCUCGUCUGAUCACCGAAGUUAAGCAACGUCGGGCCUGGUUAGUACUUGGAUGGGUGACCGCCUGGGAAUACCAGGUGCUGUAGGUUU
(((((((((....((((((((.....((((((............))))..))....)))))).)).((((((.....((.((.(((....))))).))....)))))).))))))))).."""

fasta2 = """>5s_Methanobrevibacter-ruminantium-1
UAGGUUUGGCGGUCAUAGCGAUGGGGUUACACCUGGUCUCGUUUCGAUCCCAGAAGUUAAGUCUUUUCGCGUUUUGUUUGUGUACUAUGGGUUCCGGUCUAUGGGAAUUUCAUUUAGCUGCCAGCUUUUU
.((((.(((((((....((((.((......((((((.............))))..)).....)).)).))....((........((((((((....))))))))......))....)))))))))))..."""

trna1 = """>tRNA_tdbR00000225-Synechococcus_elongatus_PCC_6301-269084-Leu-CAA
GGGCAAGUGGCGGAAUUGGUAGACGCAGCAGACUCAAAAUCUGCCGCUAGCGAUAGUGUGUGGGUUCGAGUCCCACCUUGCCCACCA
(((((((..(((...........))).(((((.......)))))...............(((((.......))))))))))))...."""

trna2 = """>tRNA_tdbR00000410-Candida_cylindracea-44322-Ser-GA
GGUGCAAUGGCCGAGUGGUUAAGGCGACGGAUUUGAACUCCGUUGGGAUUCUCCCUCGCAGGUUCGAAUCCUGUUUGCAUCGCCA
(((((((..(((..........)))((((((.......)))))).............(((((.......))))))))))))...."""

DistanceOverTimePlot(trna2, t8=1e8, title="tRNA Candida cylindracea")