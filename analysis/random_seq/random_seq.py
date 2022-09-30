import random
import subprocess
import matplotlib.pyplot as plt
import numpy as np

def GenerateRandomSequence(n, chars=["A", "C", "G", "U"]):
    s = ""
    for _ in range(n):
        s += random.choice(chars)
    return s

# n = 100
# s = GenerateRandomSequence(n)
# subprocess.run(f"echo {s} | DrTransformer --name random_{n} --outdir analysis/random_seq --logfile", shell=True)
# subprocess.run(f"echo {s} | DrTransformer --reversed-transcription --name random_{n}_rev --outdir analysis/random_seq --logfile", shell=True)

def Compare_53_35(iterations, seq_len):
    """Compare following characteristics of the normal/reversed mode on random sequences:
    minimum energy, total number of structures, structures at end of transcription"""
    energies = []
    end_structures = []
    max_id = []
    for i in range(iterations):
        print(i)
        energies.append([None, None]) # [normal, reversed]
        end_structures.append([None, None])
        max_id.append([None, None])
        s = GenerateRandomSequence(seq_len)
        # call DrTransformer and generate log files
        subprocess.run(f"echo {s} | DrTransformer --name compare_53_35 --outdir analysis/random_seq --t-end 60 --logfile", shell=True)
        subprocess.run(f"echo {s} | DrTransformer --reversed-transcription --name compare_53_35_rev --outdir analysis/random_seq --t-end 60 --logfile", shell=True)
        # parse log files
        with open("analysis/random_seq/compare_53_35.log") as f:
            current_min_energies = []          
            current_end_structures = 0
            current_max_id = 0
            lines = f.readlines()
            i = 0
            flag = False
            for line in lines:
                if not line or line.startswith("#"):
                    if flag:
                        break
                    continue
                flag = True
                words = line.split()
                # number of total structures
                if int(words[-1]) > current_max_id:
                    current_max_id = int(words[-1])
            for line in lines[::-1]:
                if not line:
                    continue
                words = line.split()
                if words[0] != str(seq_len):
                    break
                # min_energy, number of end structures
                current_min_energies.append(float(words[-7]))
                current_end_structures += 1
            energies[-1][0] = min(current_min_energies)
            end_structures[-1][0] = current_end_structures
            max_id[-1][0] = current_max_id
        with open("analysis/random_seq/compare_53_35_rev.log") as f:
            current_min_energies = []          
            current_end_structures = 0
            current_max_id = 0
            lines = f.readlines()
            i = 0
            flag = False
            for line in lines:
                if not line or line.startswith("#"):
                    if flag:
                        break
                    continue
                flag = True
                words = line.split()
                # number of total structures
                if int(words[-1]) > current_max_id:
                    current_max_id = int(words[-1])
            for line in lines[::-1]:
                if not line:
                    continue
                words = line.split()
                if words[0] != str(seq_len):
                    break
                # min_energy, number of end structures
                current_min_energies.append(float(words[-7]))
                current_end_structures += 1
            energies[-1][1] = min(current_min_energies)
            end_structures[-1][1] = current_end_structures
            max_id[-1][1] = current_max_id
    d_energies = [tup[0]-tup[1] for tup in energies]
    d_end_structures = [tup[0]-tup[1] for tup in end_structures]
    d_max_id = [tup[0]-tup[1] for tup in max_id]
    return energies, d_energies, end_structures, d_end_structures, max_id, d_max_id

def Compare_weighted_mean_energy(iterations, seq_len):
    """Plot the mean energy at the end of transcription of the normal/reversed mode on random sequences."""
    energies = []
    for i in range(iterations):
        print(i)
        energies.append([None, None]) # [normal, reversed]
        s = GenerateRandomSequence(seq_len)
        # call DrTransformer and generate log files
        subprocess.run(f"echo {s} | DrTransformer --name compare_53_35 --outdir analysis/random_seq --t-end 60 --logfile", shell=True)
        subprocess.run(f"echo {s} | DrTransformer --reversed-transcription --name compare_53_35_rev --outdir analysis/random_seq --t-end 60 --logfile", shell=True)
        # parse log files
        with open("analysis/random_seq/compare_53_35.log") as f:
            current_energy = 0  # weighted mean energy of ensemble at t-end
            total_occupancy = 0    
            lines = f.readlines()
            i = 0
            for line in lines[::-1]:
                if not line:
                    continue
                words = line.split()
                if words[0] != str(seq_len):
                    break
                # energy
                total_occupancy += float(words[4])
                current_energy += float(words[3]) * float(words[4])
            current_energy /= total_occupancy
            energies[-1][0] = current_energy
        with open("analysis/random_seq/compare_53_35_rev.log") as f:
            current_energy = 0
            total_occupancy = 0    
            lines = f.readlines()
            i = 0
            for line in lines[::-1]:
                if not line:
                    continue
                words = line.split()
                if words[0] != str(seq_len):
                    break
                # energy
                total_occupancy += float(words[4])
                current_energy += float(words[3]) * float(words[4])
            current_energy /= total_occupancy
            energies[-1][1] = current_energy
    d_energies = [tup[0]-tup[1] for tup in energies]
    
    # plot energy difference per sequence
    mi = np.floor(min(d_energies))
    ma = np.ceil(max(d_energies))
    y, x, _ = plt.hist(d_energies, bins=np.arange(mi-0.5, ma+1.5, 0.5))
    plt.vlines(x=np.mean(d_energies), ymin=0, ymax=y.max(), color = 'black', ls="dashed", linewidth = 1, label="mean energy difference")
    plt.title(f"Comparison of weighted mean energies per sequence\niterations={iterations}, seq_len={seq_len}")
    plt.xlabel("E(53) - E(35) [kcal/mol]")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(f"analysis/random_seq/weighted_mean_energy_per_seq_{iterations}_{seq_len}.png")
    plt.close()
    
    # plot energies
    bins = np.linspace(np.floor(min(min([x[0] for x in energies]), min([x[1] for x in energies]))), 
        np.ceil(max(max([x[0] for x in energies]), max([x[1] for x in energies])))+1, 20)
    plt.hist([x[0] for x in energies], bins, alpha=0.5, label="normal")
    plt.hist([x[1] for x in energies], bins, alpha=0.5, label="reversed")
    plt.title(f"Comparison of weighted mean energies\niterations={iterations}, seq_len={seq_len}")
    plt.xlabel("weighted mean energy [kcal/mol]")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(f"analysis/random_seq/weighted_mean_energy_{iterations}_{seq_len}.png")
    plt.close()

Compare_weighted_mean_energy(iterations=200, seq_len=120)

"""
iterations = 101
seq_len = 151
energies, d_energies, end_structures, d_end_structures, max_id, d_max_id = Compare_53_35(iterations, seq_len)
print(f"{energies=}\n{d_energies=}\n{end_structures=}\n{d_end_structures=}\n{max_id=}\n{d_max_id=}")

# min_energy plot
mi = np.floor(min(d_energies))
ma = np.ceil(max(d_energies))
plt.hist(d_energies, bins=np.arange(mi-0.5, ma+1.5, 1))
plt.title(f"Comparison of lowest energy at end of simulation\niterations={iterations}, seq_len={seq_len}")
plt.xlabel("dE = E(53) - E(35)")
plt.ylabel("Frequency")
plt.savefig(f"analysis/random_seq/min_energy_{iterations}_{seq_len}.png")
plt.close()

# number of end structures plot
mi = np.floor(min(d_end_structures))
ma = np.ceil(max(d_end_structures))
plt.hist(d_end_structures, bins=np.arange(mi-0.5, ma+1.5, 1))
plt.title(f"Comparison of number of structures at end of simulation\niterations={iterations}, seq_len={seq_len}")
plt.xlabel("d_end_structures = #end_structures(53) - #end_structures(35)")
plt.ylabel("Frequency")
plt.savefig(f"analysis/random_seq/end_structures_{iterations}_{seq_len}.png")
plt.close()

# number of total structures plot
# mi = np.floor(min(d_max_id))
# ma = np.ceil(max(d_max_id))
# m_abs = max(ma, -mi)
plt.hist(d_max_id)
plt.title(f"Comparison of total number of structures during simulation\niterations={iterations}, seq_len={seq_len}")
plt.xlabel("d_structures = #structures(53) - #structures(35)")
plt.ylabel("Frequency")
plt.savefig(f"analysis/random_seq/max_id_{iterations}_{seq_len}.png")
plt.close()
"""