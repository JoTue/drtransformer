#!/usr/bin/env python3

import RNA
import argparse

def get_basepair_probabilities(seq):
    fc = RNA.fold_compound(seq)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()
    return basepair_probs

def get_basepairs(ss):
    basepairs = []
    stack = []
    for i, c in enumerate(ss):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if len(stack) == 0:
                raise IndexError("No matching closing parenthesis at: " + str(i))
            basepairs.append((stack.pop(), i))
    if len(stack) > 0:
        raise IndexError("No matching opening parenthesis at: " + str(stack.pop()))
    return basepairs

def get_distance(seq, ss):
    basepair_probs = get_basepair_probabilities(seq)
    basepair_list = get_basepairs(ss)
    d = 0
    for i in range(0, len(basepair_probs)-2):
        for j in range(i+1, len(basepair_probs[i])-1):
            if (i, j) in basepair_list:
                d += 1 - basepair_probs[i+1][j+1]
            else:
                d += basepair_probs[i+1][j+1]
    return d

# seq = "GGAUACGGCCAUACUGCGCAGAAAGCACCGCUUCCCAUCCGAACAGCGAAGUUAAGCUGCGCCAGGCGGUGUUAGUACUGGGGUGGGCGACCACCCGGGAAUCCACCGUGCCGUAUCCU"
# ss = "(((((((((....((((((((.......((((((......))..))))........)))))).))((((((......(((((.(((....)))))))).....)))))))))))))))."
# print(get_distance(seq, ss))

def main():
    """ Compute distance between structure and thermodynamic equilibrium by calculating pair probabilities.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'Compute distance between structure and thermodynamic equilibrium by calculating pair probabilities.')

    parser.add_argument("seq",
            help = """Primary sequence.""")
    parser.add_argument("ss",
            help = """Secondary structure in dot-bracket format.""")

    args = parser.parse_args()

    get_distance(args.seq, args.ss)
 
if __name__ == '__main__':
    main()
