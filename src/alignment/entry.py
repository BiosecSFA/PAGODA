from src.models.quartet import Quartet
import numpy as np
from src.utils.constants import NUCLEOTIDES, MAX_SCORE_DIFF


class MatrixEntry:
    """
    A class to represent an entry in the score matrix
    last_entry: list of tuples, each tuple is the index of the last entry
    state: str, 'd1' for deletion in sequence 1, 'd2' for deletion in sequence 2,
     'i1' for insertion in seq1, 'i2' for insertion in seq2, 'm' for match
    k: int, the index of the last nucleotide in the quartet
    init_score: list of floats, the initial score of the entry
    aa1: str, the amino acid in the first sequence
    aa2: str, the amino acid in the second sequence
    prob1: pandas Series, the probabilities of the amino acids in the first sequence
    prob2: pandas Series, the probabilities of the amino acids in the second sequence
    frame: int, the frame of the sequence
    """
    def __init__(self, last_entry, state=None, k=None, init_score=None,
                 aa1=None, aa2=None, prob1=None, prob2=None, frame=None):
        self.k = k
        self.state = state
        self.init_score = init_score
        self.last_entry = last_entry
        self.compute_score(aa1, aa2, prob1, prob2, frame)

    def compute_score(self, aa1, aa2, prob1, prob2, frame):
        last_entry = None
        if self.last_entry is None:
            # if the entry is the first entry
            self.score = self.init_score[0]
        elif self.state == 'd':
            # if the entry is a deletion
            self.score = np.max(self.init_score)
            self.aa1p = []
            self.aa2p = []
            if self.score > -np.inf:
                last_entry = [self.last_entry[i] for i, x in enumerate(self.init_score) if x == self.score]
                self.aa1p.append('')
                self.aa2p.append('')
        elif self.state in ['i1', 'i2', 'm']:
            scores = []
            kmers = []
            aa1ps = []
            aa2ps = []
            # import pandas as pd
            # temp1 = []
            # temp2 = []
            self.score = -np.inf
            for i, score in enumerate(self.init_score):
                quartet = Quartet(NUCLEOTIDES[self.last_entry[i][3]], NUCLEOTIDES[self.k],
                                  aa1, aa2, prob1, prob2, frame)
                score, kmer, aa1p, aa2p = quartet.get_results()
                # temp1.append(score)
                # if (aa1p is not None) and (aa2p is not None) and (kmer is not None):
                #     temp2.append(aa1+aa1p+':'+aa2+aa2p+':'+kmer)
                # else:
                #     temp2.append('')
                score = score + self.init_score[i]
                scores.append(score)
                kmers.append(kmer)
                aa1ps.append(aa1p)
                aa2ps.append(aa2p)
                if score > self.score:
                    self.score = score
                    self.quartets = [kmer]
                    self.aa1p = [aa1p]
                    self.aa2p = [aa2p]
                    last_entry = [self.last_entry[i]]
            # get second highest score
            other_scores = [x for x in scores if x != self.score]
            if self.score > -np.inf and len(other_scores)>0:
                second_score = np.max(other_scores)
            else:
                second_score = -np.inf
            if self.score - second_score < MAX_SCORE_DIFF and second_score > -np.inf:
                for i in range(len(scores)):
                    if scores[i] == second_score:
                    # when two scores are really similar add both (but keep the highest score)
                        self.quartets.append(kmers[i])
                        self.aa1p.append(aa1ps[i])
                        self.aa2p.append(aa2ps[i])
                        last_entry.append(self.last_entry[i])
        self.last_entry = last_entry


