codon_tbl = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
nt = 'A C G T'.split()
codons = [s+t+r for s in nt for t in nt for r in nt]

import sys
from mut_code import mut_code

for c in codons:
    sys.stdout.write("%s"%codon_tbl[c])
sys.stdout.write("\n")
for c in codons:
    sys.stdout.write("\"%s\","%c)
sys.stdout.write("\n")

idx = 8
for c1 in codons:
    for c2 in codons:
        c1c2 = c1+c2
        try:
            ns = mut_code[c1c2]
            #sys.stdout.write("%d,"%ns[1])
            print(c1c2,ns[0],ns[1])
        except KeyError:
            #sys.stdout.write("0,")
            print(c1c2,0,0)

        idx += 1

        #if idx%40==0:
        #    sys.stdout.write("\n")
        #    idx = 0
sys.stdout.write("\n")
