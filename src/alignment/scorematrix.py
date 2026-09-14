import numpy as np
from src.alignment.entry import MatrixEntry
from src.utils.constants import NUCLEOTIDES, codon_table


class ScoreMatrix:
    def __init__(self, seq1, seq2, hmm1, hmm2, frame=1):
        self.seq1 = seq1
        self.seq2 = seq2
        self.hmm1 = hmm1
        self.hmm2 = hmm2
        self.frame = frame
        self._init_matrix()
        self._start_align()
        self._fill_matrix()
        self._find_max()

    def _find_max_helper(self, i, j, k, l):
        if self.matrix[i, j, k, l].score > self.max_score:
            self.max_score = self.matrix[i, j, k, l].score
            self.max_entry = self.matrix[i, j, k, l]
            self.max_idx = (i, j, k, l)

    def _find_max(self):
        self.max_score = -np.inf
        for k in range(5):
            for l in range(4):
                for i in range(len(self.seq1) // 2, len(self.seq1)):
                    self._find_max_helper(i, (len(self.seq2)), k, l)
                for j in range(len(self.seq2) // 2, len(self.seq2)):
                    self._find_max_helper((len(self.seq1)), j, k, l)

    # this needs to work with multiple last_entries
    # for each alternative entry, create a new dna, seq1_raw, seq2_raw, seq1_final, seq2_final
    def _trace_back(self, i, j, k, l, last_pos=False):
        if i == 0 or j == 0:
            return [''], [''], [''], [''], ['']
        elif k == 0:
            dna, raw2, final2 = [''], [''], ['']
            raw1 = ['-']
            final1 = ['-']
            if last_pos:
                dna = [NUCLEOTIDES[l]]
        elif k == 1:
            dna, raw1, final1 = [''], [''], ['']
            raw2 = ['-']
            final2 = ['-']
            if last_pos:
                dna = [NUCLEOTIDES[l]]
        else:
            dna, raw1, raw2, final1, final2 = [], [], [], [], []
            for i_q, quartet in enumerate(self.matrix[i, j, k, l].quartets):
                if last_pos:
                    dna.append(quartet)
                else:
                    dna.append(quartet[:3])
                if k == 4:
                    raw1.append(self.seq1[i - 1])
                    raw2.append(self.seq2[j - 1])
                elif k == 2:
                    raw1.append('+')
                    raw2.append(self.seq2[j - 1])
                elif k == 3:
                    raw1.append(self.seq1[i - 1])
                    raw2.append('+')
                final1.append(self.matrix[i, j, k, l].aa1p[i_q])
                final2.append(self.matrix[i, j, k, l].aa2p[i_q])

        new_dna, new_raw1, new_raw2, new_final1, new_final2 = [], [], [], [], []
        for entry in set(self.matrix[i, j, k, l].last_entry):
            last_dna, last_raw1, last_raw2, last_final1, last_final2 = self._trace_back(*entry)
            for x in range(len(dna)):
                for y in range(len(last_dna)):
                    if (last_final1[y] + final1[x] not in new_final1) and (
                            last_final2[y] + final2[x] not in new_final2):
                        new_dna.append(last_dna[y] + dna[x])
                        new_raw1.append(last_raw1[y] + raw1[x])
                        new_raw2.append(last_raw2[y] + raw2[x])
                        new_final1.append(last_final1[y] + final1[x])
                        new_final2.append(last_final2[y] + final2[x])
        return new_dna, new_raw1, new_raw2, new_final1, new_final2

    def _init_matrix(self):
        matrix = np.empty((len(self.seq1) + 1, len(self.seq2) + 1, 5, 4), dtype=MatrixEntry)
        self.matrix = matrix
        entry = MatrixEntry(last_entry=None, init_score=[0])
        entry.aa1p = ''
        entry.aa2p = ''
        self.matrix[0, 0, :, :] = entry

    def _start_align(self):
        invalid_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
        # for all deletions of seq1 before seq2 starts
        for i in range(len(self.seq1)):
            self.matrix[i + 1, 0, :, :] = invalid_entry
            self.matrix[i + 1, 0, 0, :] = self._delete(i + 1, 0, 0, 0, match_empty=True)

        # for all deletions of seq2 before seq1 starts
        for j in range(len(self.seq2)):
            self.matrix[0, j + 1, :, :] = invalid_entry
            self.matrix[0, j + 1, 1, :] = self._delete(0, j + 1, 1, 0, match_empty=True)

    def _fill_matrix(self):
        # for all other entries
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                for k in range(4):
                    self.matrix[i, j, 0, k] = self._delete(i, j, 0, k)  # delete in seq1
                    self.matrix[i, j, 1, k] = self._delete(i, j, 1, k)  # delete in seq2
                    self.matrix[i, j, 2, k] = self._insert(i, j, 0, k)  # insert in seq1
                    self.matrix[i, j, 3, k] = self._insert(i, j, 1, k)  # insert in seq2
                    self.matrix[i, j, 4, k] = self._match(i, j, k)  # match

    # s: s=0 means insert in seq1, s=1 means insert in seq2
    # k is the nucleotide index for the last nucleotide in the quartet
    def _insert(self, i, j, s, k):
        last_entry = []
        scores = []
        if s == 0:
            for x in [4, 2]:
                for y in range(4):
                    if self.matrix[i, j - 1, x, y].score == -np.inf:
                        continue
                    else:
                        entry = (i, j - 1, x, y)
                        last_entry.append(entry)
                        if x == 2:
                            scores.append(
                                self.matrix[entry].score - self.hmm1.trans_prob.loc[entry[0], 'i->i'] - \
                                self.hmm2.trans_prob.loc[entry[1], 'm->m'])
                        else:
                            scores.append(
                                self.matrix[entry].score - self.hmm1.trans_prob.loc[entry[0], 'm->i'] - \
                                self.hmm2.trans_prob.loc[entry[1], 'm->m'])
            if len(last_entry) > 0:
                new_entry = MatrixEntry(last_entry, state='i1', k=k, init_score=scores,
                                        aa1='+', aa2=self.seq2[j - 1],
                                        prob1=self.hmm1.insert_prob.loc[1],
                                        prob2=self.hmm2.match_prob.loc[j],
                                        frame=self.frame)
            else:
                new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
        elif s == 1:
            for x in [4, 3]:
                for y in range(4):
                    if self.matrix[i - 1, j, x, y].score == -np.inf:
                        continue
                    else:
                        entry = (i - 1, j, x, y)
                        last_entry.append(entry)
                        if x == 3:
                            scores.append(
                                self.matrix[entry].score - self.hmm2.trans_prob.loc[entry[1], 'i->i'] - \
                                self.hmm1.trans_prob.loc[entry[0], 'm->m'])
                        else:
                            scores.append(
                                self.matrix[entry].score - self.hmm2.trans_prob.loc[entry[1], 'm->i'] - \
                                self.hmm1.trans_prob.loc[entry[0], 'm->m'])
            if len(last_entry) > 0:
                new_entry = MatrixEntry(last_entry, state='i2', k=k, init_score=scores,
                                        aa1=self.seq1[i - 1], aa2='+',
                                        prob1=self.hmm1.match_prob.loc[i],
                                        prob2=self.hmm2.insert_prob.loc[j],
                                        frame=self.frame)
            else:
                new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
        else:
            raise ValueError('s must be 0 or 1')
        if new_entry.score == -np.inf:
            new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
            new_entry.score = -np.inf
        return new_entry

    def _delete(self, i, j, s, k, match_empty=False):
        # this option is to avoid penalizing deletions at the start of the sequence
        if match_empty:
            scores = np.array([0])
            if s == 0:
                last_entry = [(i - 1, j, 0, k)]
            elif s == 1:
                last_entry = [(i, j - 1, 1, k)]
            else:
                raise ValueError('s must be 0 or 1')
        else:
            last_entry = []
            scores = []
            if s == 0:
                for x in [0, 4]:
                    entry = (i - 1, j, x, k)
                    if self.matrix[entry].score == -np.inf:
                        continue
                    else:
                        last_entry.append(entry)
                        if x == 0:
                            scores.append(
                                self.matrix[entry].score - self.hmm1.trans_prob.loc[entry[0], 'd->d'])
                        else:
                            scores.append(
                                self.matrix[entry].score - self.hmm1.trans_prob.loc[entry[0], 'm->d'])
            elif s == 1:
                for x in [1, 4]:
                    entry = (i, j - 1, x, k)
                    if self.matrix[entry].score == -np.inf:
                        continue
                    else:
                        last_entry.append(entry)
                        if x == 1:
                            scores.append(
                                self.matrix[entry].score - self.hmm2.trans_prob.loc[entry[1], 'd->d'])
                        else:
                            scores.append(
                                self.matrix[entry].score - self.hmm2.trans_prob.loc[entry[1], 'm->d'])
            else:
                raise ValueError('s must be 0 or 1')
        if len(last_entry) == 0:
            new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
            new_entry.score = -np.inf
        else:
            new_entry = MatrixEntry(last_entry, state='d', k=k, init_score=scores)
        if new_entry.score == -np.inf:
            new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
        return new_entry

    def _match(self, i, j, k):
        last_entry = []
        scores = []
        for x in range(5):
            for y in range(4):
                if self.matrix[i - 1, j - 1, x, y].score == -np.inf:
                    continue
                elif i == 0 or j == 0:
                    continue
                else:
                    entry = (i - 1, j - 1, x, y)
                    last_entry.append((i - 1, j - 1, x, y))
                if x == 0:
                    scores.append(self.matrix[entry].score -
                                  self.hmm1.trans_prob.loc[entry[0], 'd->m'] - self.hmm2.trans_prob.loc[
                                      entry[1], 'm->m'])
                elif x == 1:
                    scores.append(self.matrix[entry].score -
                                  self.hmm2.trans_prob.loc[entry[1], 'd->m'] - self.hmm1.trans_prob.loc[
                                      entry[0], 'm->m'])
                elif x == 2:
                    scores.append(self.matrix[entry].score -
                                  self.hmm1.trans_prob.loc[entry[0], 'i->m'] - self.hmm2.trans_prob.loc[
                                      entry[1], 'm->m'])
                elif x == 3:
                    scores.append(self.matrix[entry].score -
                                  self.hmm2.trans_prob.loc[entry[1], 'i->m'] - self.hmm1.trans_prob.loc[
                                      entry[0], 'm->m'])
                elif x == 4:
                    scores.append(self.matrix[entry].score -
                                  self.hmm1.trans_prob.loc[entry[0], 'm->m'] - self.hmm2.trans_prob.loc[
                                      entry[1], 'm->m'])
        new_entry = MatrixEntry(last_entry, state='m', k=k, init_score=scores,
                                aa1=self.seq1[i - 1], aa2=self.seq2[j - 1],
                                prob1=self.hmm1.match_prob.loc[i], prob2=self.hmm2.match_prob.loc[j], frame=self.frame)
        if new_entry.score == -np.inf:
            new_entry = MatrixEntry(last_entry=None, init_score=[-np.inf])
            new_entry.score = -np.inf
        return new_entry
