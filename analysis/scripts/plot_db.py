"""
Old script to analyze the output json-file from analyze_db_cli.py - not used in the final project.
"""


import json
import matplotlib.pyplot as plt

dirpath = "analysis/slurm_out"
filename = "grp_short"
with open(f"{dirpath}/{filename}/{filename}.json") as f:
    db = json.load(f)

# write to analysis.tsv
with open(f"{dirpath}/{filename}/{filename}_analysis.tsv", "w") as f:
    print("All entries:\nheader\tlength\td_nat_sub\td_nat_tend\td_nat_tend_rev\td_sub_tend\td_sub_tend_rev", file=f)
    for header in db:
        print(f"{header}\t{db[header]['seqlen']}\t{db[header]['dist_nat_subopt']:.2f}\t{db[header]['dist_nat_tend']:.2f}\t{db[header]['dist_nat_tend_rev']:.2f}"
        f"\t{db[header]['dist_subopt_tend']:.2f}\t{db[header]['dist_subopt_tend_rev']:.2f}", file=f)
    c_normal = 0
    c_reversed = 0
    print("normal+5 < subopt:\nheader\tlength\td_nat_sub\td_nat_tend\td_nat_tend_rev\td_sub_tend\td_sub_tend_rev", file=f)
    for header in db:
        if db[header]['dist_nat_tend'] + 5 < db[header]['dist_nat_subopt']:
            print(f"{header}\t{db[header]['seqlen']}\t{db[header]['dist_nat_subopt']:.2f}\t{db[header]['dist_nat_tend']:.2f}\t{db[header]['dist_nat_tend_rev']:.2f}"
            f"\t{db[header]['dist_subopt_tend']:.2f}\t{db[header]['dist_subopt_tend_rev']:.2f}", file=f)
            c_normal += 1
        if db[header]['dist_nat_tend_rev'] + 5 < db[header]['dist_nat_subopt']:
            c_reversed += 1
    print(f"{c_normal=}, {c_reversed=}", file=f)
    
    c_normal = 0
    c_reversed = 0
    c_tie = 0
    print("normal+5 < reversed:\nheader\tlength\td_nat_sub\td_nat_tend\td_nat_tend_rev\td_sub_tend\td_sub_tend_rev", file=f)
    for header in db:
        if db[header]['dist_nat_tend'] + 5 < db[header]['dist_nat_tend_rev']:
            print(f"{header}\t{db[header]['seqlen']}\t{db[header]['dist_nat_subopt']:.2f}\t{db[header]['dist_nat_tend']:.2f}\t{db[header]['dist_nat_tend_rev']:.2f}"
            f"\t{db[header]['dist_subopt_tend']:.2f}\t{db[header]['dist_subopt_tend_rev']:.2f}", file=f)
            c_normal += 1
        elif db[header]['dist_nat_tend_rev'] + 5 < db[header]['dist_nat_tend']:
            c_reversed += 1
        else:
            c_tie += 1
    print(f"{c_normal=}, {c_reversed=}, {c_tie=}", file=f)

# plot dist_nat_subopt vs. dist_nat_mfe
dist_nat_mfe = []
dist_nat_subopt = []
for header in db:
    dist_nat_mfe.append(db[header]['dist_nat_mfe'])
    dist_nat_subopt.append(db[header]['dist_nat_subopt'])
plt.scatter(dist_nat_mfe, dist_nat_subopt, s=10)
plt.plot([0, 1000], [0, 1000], color = 'black', ls="dashed", linewidth = 1)
plt.xlim(0, max(max(dist_nat_mfe), max(dist_nat_subopt))+1)
plt.ylim(0, max(max(dist_nat_mfe), max(dist_nat_subopt))+1)
plt.title("dist_nat_subopt vs. dist_nat_mfe")
plt.xlabel("dist_nat_mfe")
plt.ylabel("dist_nat_subopt")
plt.savefig(f"{dirpath}/{filename}/nat_mfe_vs_subopt.png")
plt.close()

# plot dist_subopt_tend vs. dist_mfe_tend
dist_mfe_tend = []
dist_subopt_tend = []
for header in db:
    dist_mfe_tend.append(db[header]['dist_mfe_tend'])
    dist_subopt_tend.append(db[header]['dist_subopt_tend'])
plt.scatter(dist_mfe_tend, dist_subopt_tend, s=10)
plt.plot([0, 1000], [0, 1000], color = 'black', ls="dashed", linewidth = 1)
plt.xlim(0, max(max(dist_mfe_tend), max(dist_subopt_tend))+1)
plt.ylim(0, max(max(dist_mfe_tend), max(dist_subopt_tend))+1)
plt.title("dist_subopt_tend vs. dist_mfe_tend")
plt.xlabel("dist_mfe_tend")
plt.ylabel("dist_subopt_tend")
plt.savefig(f"{dirpath}/{filename}/tend_mfe_vs_subopt.png")
plt.close()

# plot dist_subopt_tend_rev vs. dist_mfe_tend_rev
dist_mfe_tend_rev = []
dist_subopt_tend_rev = []
for header in db:
    dist_mfe_tend_rev.append(db[header]['dist_mfe_tend_rev'])
    dist_subopt_tend_rev.append(db[header]['dist_subopt_tend_rev'])
plt.scatter(dist_mfe_tend_rev, dist_subopt_tend_rev, s=10)
plt.plot([0, 1000], [0, 1000], color = 'black', ls="dashed", linewidth = 1)
plt.xlim(0, max(max(dist_mfe_tend_rev), max(dist_subopt_tend_rev))+1)
plt.ylim(0, max(max(dist_mfe_tend_rev), max(dist_subopt_tend_rev))+1)
plt.title("dist_subopt_tend_rev vs. dist_mfe_tend_rev")
plt.xlabel("dist_mfe_tend_rev")
plt.ylabel("dist_subopt_tend_rev")
plt.savefig(f"{dirpath}/{filename}/tend_mfe_vs_subopt_rev.png")
plt.close()

# plot dist_nat_equ_rev vs. dist_nat_equ
dist_nat_equ = []
dist_nat_equ_rev = []
print("Differing Equilibria:\nheader\t\td_nat_equ\td_nat_equ_rev")
for header in db:
    dist_nat_equ.append(db[header]['dist_nat_equ'])
    dist_nat_equ_rev.append(db[header]['dist_nat_equ_rev'])
    if abs(db[header]['dist_nat_equ_rev'] - db[header]['dist_nat_equ']) > 5:
        print(f"{header}\t{db[header]['dist_nat_equ']:.2f}\t\t{db[header]['dist_nat_equ_rev']:.2f}")
plt.scatter(dist_nat_equ, dist_nat_equ_rev, s=10)
plt.plot([0, 1000], [0, 1000], color = 'black', ls="dashed", linewidth = 1)
plt.xlim(0, max(max(dist_nat_equ), max(dist_nat_equ_rev))+1)
plt.ylim(0, max(max(dist_nat_equ), max(dist_nat_equ_rev))+1)
plt.title("dist_nat_equ_rev vs. dist_nat_equ")
plt.xlabel("dist_nat_equ")
plt.ylabel("dist_nat_equ_rev")
plt.savefig(f"{dirpath}/{filename}/nat_equ_vs_nat_equ_rev.png")
plt.close()
