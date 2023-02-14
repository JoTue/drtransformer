#!/usr/bin/env python3

import subprocess

rna_classes = ["5s", "16s", "23s", "grpI", "grpII", "RNaseP", "srp", "telomerase", "tmRNA", "tRNA"]
file_name = "missing_ids.txt"
subprocess.run(f"> {file_name}", shell=True)  # create/overwrite output file

# write number of missing ids per class
subprocess.run(f"echo 'Number_of_missing_ids_per_class:' >> {file_name}", shell=True)
for rna_class in rna_classes:
    subprocess.run(f"echo '{rna_class}' >> {file_name}", shell=True)
    subprocess.run(f"bin/get_missing_ids.py data/{rna_class}.db drtr_results/{rna_class}/ --id1 1 --id2 10000 | wc -l >> {file_name}", shell=True)
    subprocess.run(f"echo '' >> {file_name}", shell=True)

# write number of missing ids per class
subprocess.run(f"echo 'Missing_ids_per_class:' >> {file_name}", shell=True)
for rna_class in rna_classes:
    subprocess.run(f"echo '{rna_class}' >> {file_name}", shell=True)
    subprocess.run(f"bin/get_missing_ids.py data/{rna_class}.db drtr_results/{rna_class}/ --id1 1 --id2 10000 >> {file_name}", shell=True)
    subprocess.run(f"echo '' >> {file_name}", shell=True)
