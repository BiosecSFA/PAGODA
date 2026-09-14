import numpy as np
from biotite.sequence import NucleotideSequence
def find_max_helper(score_matrix, i, j, k, l, max_score):
    if (score_matrix.matrix[i, j, k, l].score > max_score) and (score_matrix.matrix[i,j,k,l].score < -0):
        max_score = score_matrix.matrix[i, j, k, l].score
        max_idx = (i, j, k, l)
    else:
        max_idx = score_matrix.max_idx
    return max_score, max_idx

def find_max(score_matrix, left=True, minlen=-1):
    max_score = -np.inf
    if minlen<0:
        minlen = min(len(score_matrix.seq1), len(score_matrix.seq2))//2
    for k in range(5):
        for l in range(4):
            if left:
                for i in range(minlen, len(score_matrix.seq1)):
                    max_score, score_matrix.max_idx = find_max_helper(score_matrix, i, (len(score_matrix.seq2)), k, l, max_score)
            else:
                for j in range(minlen, len(score_matrix.seq2)):
                    max_score, score_matrix.max_idx = find_max_helper(score_matrix, (len(score_matrix.seq1)), j, k, l,max_score)
    return score_matrix.max_idx, max_score


def correct_translation(final, final_withoutgaps):
    if final.replace('-', '') != final_withoutgaps:
        n_diff = 0
        j=0
        corrected = []
        for i, x in enumerate(final):
            if x == '-':
                corrected.append('-')
            else:
                corrected.append(final_withoutgaps[j])
                j += 1
                if corrected[-1] != final[i]:
                    n_diff += 1
        print(final, n_diff)
        return ''.join(corrected)
    else:
        return final

# dna = dna[0]
# idx = overlap_idx
# raw1 = raw1[0]
# raw2 = raw2[0]
# final1 = final1[0]
# final2 = final2[0]
def complete_sequence(dna, cds1, cds2, idx, raw1, raw2, final1, final2):
    assert len([x for x in raw1 if x!='-']) == len([x for x in raw2 if x!='-'])
    end1 = idx[0][0] * 3
    end2 = idx[0][1] * 3
    left1 = False
    right1 = False
    len_entangle1 = len([x for x in raw1 if x!='+'])
    len_entangle2 = len([x for x in raw2 if x!='+'])
    start1 = end1 - len_entangle1*3
    start2 = end2 - len_entangle2*3
    if start1 == 0 and start2 > 0:
        left = cds2[:start2]
    elif start1 > 0 and start2 == 0:
        left1 = True
        left = cds1[:start1]
    elif start1 == 0 and start2 == 0:
        left = ''
    else:
        # print error that the start of the sequence is not correct
        raise ValueError('error in starting position')
    if end1 == len(cds1) and end2 < len(cds2):
        right = cds2[end2:]
    elif end1 < len(cds1) and end2 == len(cds2):
        right = cds1[end1:]
        right1 = True
    elif end1 == len(cds1) and end2 == len(cds2):
        right = ''
    else:
        # print error that the end of the sequence is not correct
        raise ValueError('error in ending position')
    full_seq = ['a' + left.lower()[:-1], left.lower()][left1] + dna.upper() + [right.lower(), right.lower()[1:] + 'A'][right1]
    nt_seq1 = ['', left][left1] + dna[:-1] + ['', dna[-1:] + right[1:]][right1]
    recon_gene1 = NucleotideSequence(nt_seq1).translate(
        complete=True)
    # occasionally the left part of gene 2 + dna[0] could constitute a stop codon... no good way of avoiding it because that part of the sequence is not entangled
    # just filter those sequences out
    nt_seq2 = [left[:-1] + dna[0], ''][left1] + dna[1:] + [right, ''][right1]
    recon_gene2 = NucleotideSequence(nt_seq2).translate(
        complete=True)
    final1_withoutgaps = str(NucleotideSequence(dna[:-1]).translate(complete=True))
    final2_withoutgaps = str(NucleotideSequence(dna[1:]).translate(complete=True))
    # correct final1 sequence with final2_withoutgaps
    final1 = correct_translation(final1, final1_withoutgaps)
    final2 = correct_translation(final2, final2_withoutgaps)
    assert final2.replace('-', '') == final2_withoutgaps
    assert final1.replace('-', '') == final1_withoutgaps
    no_insert1 = ''.join([x for i,x in enumerate(final1) if raw1[i]!='+'])
    no_insert2 = ''.join(x for i,x in enumerate(final2) if raw2[i]!='+')
    aligned1 = str(recon_gene1).replace(final1_withoutgaps, no_insert1)
    aligned2 = str(recon_gene2).replace(final2_withoutgaps, no_insert2)
    return full_seq, nt_seq1, nt_seq2, len(dna)//3, aligned1, aligned2
