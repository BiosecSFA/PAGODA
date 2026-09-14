import Bio.Seq
from src.utils.constants import NUCLEOTIDES, codon_table
import numpy as np


def translate(kmer, frame, ckmer=None):
    if frame == 1:
        aa1p = codon_table[kmer[:3]]
        aa2p = codon_table[kmer[1:]]
    elif frame == 2:
        aa1p = codon_table[kmer[1:]]
        aa2p = codon_table[kmer[:3]]
    elif frame == -1:
        aa1p = codon_table[kmer[:3]]
        aa2p = codon_table[ckmer[1:]]
    elif frame == -2:
        aa1p = codon_table[kmer[1:]]
        aa2p = codon_table[ckmer[:3]]
    else:
        raise ValueError('Frame must be 1, 2, -1, or -2')
    return (aa1p, aa2p)


def reverse_complement(kmer):
    return Bio.Seq.reverse_complement(kmer)



# precompute translation of quartets to amino acids
def precompute_quartet_translation(frame):
    res = {}
    for start in NUCLEOTIDES:
        res[start] = {}
        for end in NUCLEOTIDES:
            res[start][end] = {}
            for n1 in NUCLEOTIDES:
                for n2 in NUCLEOTIDES:
                    kmer = f'{start}{n1}{n2}{end}'
                    aa1, aa2 = translate(kmer, frame, reverse_complement(kmer))
                    pair = aa1+aa2
                    if '*' not in pair:
                        if pair not in res[start][end].keys():
                            res[start][end][pair] = {}
                            res[start][end][pair]['count'] = 1
                            res[start][end][pair]['kmer'] = [kmer]
                        else:
                            res[start][end][pair]['count'] += 1
                            res[start][end][pair]['kmer'].append(kmer)
    return res




class Quartet:
    def __init__(self, start, end, aa1, aa2, prob1, prob2, frame):
        self.start = start
        self.end = end
        self.frame = frame
        self.aa1 = aa1
        self.aa2 = aa2
        self._max_kmer(prob1, prob2)

    def _max_kmer(self, prob1, prob2):
        self.max_score = -np.inf
        for x in NUCLEOTIDES:
            for y in NUCLEOTIDES:
                self.kmer = f'{self.start}{x}{y}{self.end}'
                if self.frame in [-1, -2]:
                    aa1p, aa2p = translate(self.kmer, self.frame, reverse_complement(self.kmer))
                else:
                    aa1p, aa2p = translate(self.kmer, self.frame)
                score1 = -prob1[aa1p]
                score2 = -prob2[aa2p]
                if score1 + score2 > self.max_score:
                    self.max_score = score1 + score2
                    self.max_kmer = self.kmer
                    self.aa1p = aa1p
                    self.aa2p = aa2p
        if self.max_score == -np.inf:
            pass
        else:
            self.kmer = self.max_kmer
            self.ckmer = reverse_complement(self.kmer)
            self.max_score = self.max_score


    def get_results(self):
        if self.max_score == -np.inf:
            return -np.inf, None, None, None
        else:
            return self.max_score, self.kmer, self.aa1p, self.aa2p

