#!/usr/bin/env python3

def read_fasta(filename):
    seqs = {}
    with open(filename) as f:
        name = ""
        seq = ""
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = seq
                name = line.strip()[1:].split()[0]
                seq = ""
            else:
                seq += line.strip()
        if name:
            seqs[name] = seq
    return seqs

print("PARSIMONIA ANALYSIS")
print("="*60)

seqs = read_fasta('alignment_6_mtdna.fasta')

sima = list(seqs.values())[0]
d3 = list(seqs.values())[1]
nean = list(seqs.values())[2]

diff_sd = sum(1 for i in range(min(len(sima), len(d3))) if sima[i] != d3[i])
diff_sn = sum(1 for i in range(min(len(sima), len(nean))) if sima[i] != nean[i])
diff_dn = sum(1 for i in range(min(len(d3), len(nean))) if d3[i] != nean[i])

print("DISTANCES:")
print("Sima-Denisova3:", diff_sd, "bp")
print("Sima-Neandertal:", diff_sn, "bp")
print("Denisova3-Neandertal:", diff_dn, "bp")
print("")

h1 = diff_sd + diff_sn + diff_dn
h2 = diff_sd + diff_dn

print("HYPOTHESIS 1 (sister):", h1, "steps")
print("HYPOTHESIS 2 (mother):", h2, "steps")
print("DIFFERENCE:", h1 - h2, "steps")
print("REDUCTION:", 100 * (h1 - h2) / h1, "%")
print("")

if h2 < h1:
    print("RESULT: Hypothesis 2 (mother) is MORE PARSIMONIOUS")
