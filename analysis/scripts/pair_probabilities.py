#!/usr/bin/env python3

import RNA

def get_basepair_probabilities(seq):
    fc = RNA.fold_compound(seq)
    (propensity, ensemble_energy) = fc.pf()
    basepair_probs = fc.bpp()
    return basepair_probs
    # for i in range(1, len(basepair_probs)):
    #     for j in range(i+1,len(basepair_probs[i])+1):
    #         print("pr(%d,%d) = %g" % (i, j, basepair_probs[i-1][j-1]))

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

seq = "GGAUACGGCCAUACUGCGCAGAAAGCACCGCUUCCCAUCCGAACAGCGAAGUUAAGCUGCGCCAGGCGGUGUUAGUACUGGGGUGGGCGACCACCCGGGAAUCCACCGUGCCGUAUCCU"
ss = "(((((((((....((((((((.......((((((......))..))))........)))))).))((((((......(((((.(((....)))))))).....)))))))))))))))."
print(get_distance(seq, ss))