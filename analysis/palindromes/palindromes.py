"""Test normal/reversed mode on palindromic random sequences."""

import random
import subprocess

def GenerateRandomSequence(n, chars=["A", "C", "G", "U"]):
    s = ""
    for _ in range(n):
        s += random.choice(chars)
    return s

n = 40
s = GenerateRandomSequence(n)
p = s + s[::-1]

subprocess.run(f"echo {p} | DrTransformer --name pal_{n*2} --outdir analysis/palindromes --logfile", shell=True)
subprocess.run(f"echo {p} | DrTransformer --reversed-transcription --name pal_{n*2}_rev --outdir analysis/palindromes --logfile", shell=True)
