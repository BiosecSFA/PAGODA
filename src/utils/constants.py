from biotite.sequence import NucleotideSequence
import numpy as np

NUCLEOTIDES = 'ACGT'

codon_table = {}
for x in NUCLEOTIDES:
    for y in NUCLEOTIDES:
        for z in NUCLEOTIDES:
           codon_table[x + y + z] = str(NucleotideSequence(x + y + z).translate(complete=True))

inv_codon_table = {}
for key, value in codon_table.items():
    if value not in inv_codon_table.keys():
        inv_codon_table[value] = [key]
    else:
        inv_codon_table[value].append(key)


DEL_COST = -np.log(1/21/21)
MAX_SCORE_DIFF = 0.05
SCRATCH = '/Users/xu26/Work/overlapDP/scratch/'
# precompute translation of quartets to amino acids

