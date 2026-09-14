from potts import Potts
from potts.potts import get_subset_potts_optimized
from pycameox.ilp import CODON_TABLE as codon_table
from pycameox.ilp import AA_TO_I as aa_to_i
from itertools import product
from src.alignment.traceback import complete_sequence
from src.utils.utils import substring_coord
from src.visualization.alignment_viz import align_entangled_seq, make_entanglement_viz, update_final
from src.optimization.greedy import greedy_optimize
from src.utils.utils import compute_potts_score

import torch
aa_to_i['*'] = 21
from src.utils.constants import codon_table

tetra_nts = list(map(''.join, product('ACGT', repeat=4)))
i_to_tetra_nt = {i: tetra_nt for i, tetra_nt in enumerate(tetra_nts)}
tetra_nt_to_i = {tetra_nt: i for i, tetra_nt in i_to_tetra_nt.items()}
from pycameox.ilp import AA_TO_I as aa_to_i

def to_model_without_gaps(model, aa_seq):
    large_non_gap_idxs = [i for i, c in enumerate(aa_seq) if c != '-']
    L, A = model.L, model.A
    h = model.h.reshape(L, A)
    W = model.W.weight
    h_new, W_new = get_subset_potts_optimized(h, W, large_non_gap_idxs, [aa_to_i[aa] for aa in aa_seq])
    new_L = h_new.shape[0]
    h_new_no_gap = h_new[:, :21]
    W_new_no_gap = W_new.reshape(new_L, A, new_L, A)[:, :21, :, :21].reshape(new_L * 21, new_L * 21)
    model_no_gaps = Potts(h=h_new_no_gap, W=W_new_no_gap)
    return model_no_gaps



def from_no_gaps_to_with_insertions(model_no_gaps, nt_seq, aa_seq_no_gaps):
    codons_without_gaps = [nt_seq[i * 3:(i + 1) * 3].upper() for i in range(len(nt_seq) // 3)]
    insertions = []
    codon_idx = 0
    for i, aa in enumerate(aa_seq_no_gaps):
        codon = codons_without_gaps[codon_idx].upper()
        while codon_table[codon] != aa_seq_no_gaps[i]:
            # codon part of insertion
            insertions.append((codon_idx, codon))
            codon_idx += 1
            codon = codons_without_gaps[codon_idx].upper()
        codon_idx += 1

    L, A = model_no_gaps.L, model_no_gaps.A
    model_no_gaps_with_ins = Potts(h=model_no_gaps.h.reshape(L, A).detach().clone(),
                                   W=model_no_gaps.W.weight.detach().clone())
    for idx, codon in insertions:
        L, A = model_no_gaps_with_ins.L, model_no_gaps_with_ins.A
        h_c = model_no_gaps_with_ins.h.reshape(L, A).detach().clone()
        W_c = model_no_gaps_with_ins.W.weight.reshape(L, A, L, A).transpose(1, 2).detach().clone()
        new_h = torch.concat((h_c[:idx], torch.zeros(1, A), h_c[idx:]), dim=0, )
        new_W = torch.zeros(L + 1, L + 1, A, A)
        new_W[:idx, :idx] = W_c[:idx, :idx].clone()
        new_W[idx + 1:, idx + 1:] = W_c[idx:, idx:].clone()
        new_W[:idx, idx + 1:] = W_c[:idx, idx:].clone()
        new_W[idx + 1:, :idx] = W_c[idx:, :idx].clone()
        model_no_gaps_with_ins = Potts(h=new_h, W=new_W.transpose(1, 2).reshape((L + 1) * A, (L + 1) * A))
    return model_no_gaps_with_ins





def get_insertions(nt_seq, aa_seq_no_gaps):
    codons_without_gaps = [nt_seq[i * 3:(i + 1) * 3].upper() for i in range(len(nt_seq) // 3)]
    insertions = []
    codon_idx = 0
    for i, aa in enumerate(aa_seq_no_gaps):
        codon = codons_without_gaps[codon_idx].upper()
        while codon_table[codon] != aa_seq_no_gaps[i]:
            # codon part of insertion
            insertions.append((codon_idx, codon))
            codon_idx += 1
            codon = codons_without_gaps[codon_idx].upper()
        codon_idx += 1
    return insertions

def add_stop_codon_to_model(model):
    L, A = model.L, model.A
    h = model.h.reshape(L, A)
    W = model.W.weight.reshape(L, A, L, A).transpose(1,2)
    new_h = torch.zeros(L+1, A+1)
    new_h[:L, :A] = h
    new_W = torch.zeros(L+1, L+1, A+1, A+1)
    new_W[:L, :L, :A, :A] = W
    for aa in aa_to_i.keys():
        if aa != "*":
            new_h[L, aa_to_i[aa]] += 10000.
    for pos in range(L):
        new_h[pos, A] += 10000.
    model_with_stop = Potts(h=new_h, W=new_W.transpose(2,1).reshape((L+1)*(A+1), (L+1)*(A+1)))
    return model_with_stop

def add_start_codon_to_model(model):
    L, A = model.L, model.A
    h = model.h.reshape(L, A).detach().clone()
    for aa in aa_to_i.keys():
        if aa != "M":
            h[0, aa_to_i[aa]] += 10000.
    model_with_start = Potts(h=h.detach().clone(), W=model.W.weight.detach().clone())
    return model_with_start

def get_adjusted_model(nt_seq, aa_seq_no_ins_aligned, model):
    if aa_seq_no_ins_aligned.endswith("*"):
        aa_seq = aa_seq_no_ins_aligned[:-1]
    else:
        aa_seq = aa_seq_no_ins_aligned
    if '-' not in aa_seq:
        model_no_gaps = model
    else:
        model_no_gaps = to_model_without_gaps(model, aa_seq)
    aa_seq_no_ins_no_gaps = aa_seq_no_ins_aligned.replace('-', '')
    insertions = get_insertions(nt_seq, aa_seq_no_ins_no_gaps)
    model_no_gaps_with_ins = from_no_gaps_to_with_insertions(model_no_gaps, nt_seq, aa_seq_no_ins_no_gaps)
    model_no_gaps_with_ins = add_stop_codon_to_model(model_no_gaps_with_ins)
    model_no_gaps_with_ins = add_start_codon_to_model(model_no_gaps_with_ins)
    return model_no_gaps_with_ins, insertions




def get_entanglement_seq(dna, raw1, raw2, final1, final2, overlap_idx, minlen,order, cds1, cds2, potts1, potts2, gene1, gene2, prot1, prot2, id_gen):
    # if minlen=='full':
    #     plot_prefix = 'full'
    # else:
    #     plot_prefix = f'minlen{minlen}order{order}'
    res = []
    for i in range(len(raw1)):
        full_seq, nt_seq1, nt_seq2, overlaplen, aligned1, aligned2 = complete_sequence(
            dna[i], cds1, cds2, overlap_idx, raw1[i], raw2[i], final1[i], final2[i])
        potts1_adjusted, insertions1 = get_adjusted_model(nt_seq1, aligned1, potts1)
        potts2_adjusted, insertions2 = get_adjusted_model(nt_seq2, aligned2, potts2)
        coord1 = substring_coord(nt_seq1.upper(), full_seq.upper())
        coord2 = substring_coord(nt_seq2.upper(), full_seq.upper())
        codon1_updated, codon2_updated, initial_potts1, initial_potts2, updated_potts1, updated_potts2 = \
            greedy_optimize(full_seq, nt_seq1, nt_seq2, potts1_adjusted, potts2_adjusted)
        aligned_seqs, alignment = align_entangled_seq(full_seq, nt_seq1, nt_seq2, final1[i], final2[i], raw1[i],
                                                      raw2[i], prot1, prot2)
        # highlight, gap1, mismatch1, gap2, mismatch2 = make_entanglement_viz(f'data/msa/{gene1}.msa',
        #                                                                     f'data/msa/{gene2}.msa',
        #                                                                     alignment, f'{gene1}_{gene2}/{plot_prefix}_{i}',
        #                                                                     plot_size=(len(final1[i]), 10))
        updated_final1 = update_final(final1[i], codon1_updated)
        updated_final2 = update_final(final2[i], codon2_updated)
        full_seq, nt_seq1, nt_seq2, _, updated_aligned1, updated_aligned2 = complete_sequence(
            ''.join(codon1_updated) + codon2_updated[-1][-1], cds1, cds2, overlap_idx, raw1[i], raw2[i], updated_final1,
            updated_final2)
        aligned_seqs, alignment = align_entangled_seq(full_seq, nt_seq1, nt_seq2, updated_final1, updated_final2,
                                                      raw1[i], raw2[i], prot1, prot2)
        # highlight_opt, gap1, mismatch1_opt, gap2, mismatch2_opt = make_entanglement_viz(
        #     f'data/msa/{gene1}.msa', f'data/msa/{gene2}.msa', alignment,
        #     f'{gene1}_{gene2}/{plot_prefix}_optimized_{i}', plot_size=(len(final1[i]), 10))
        seqid = id_gen.get_next_id()
        # this is somewhat unavoidable (see complete_sequence comments)
        if (aligned1[-1] != '*' or aligned2[-1] != '*' or
                updated_aligned1[-1] != '*' or updated_aligned2[-1] != '*' or
                '*' in aligned1[:-1] or '*' in aligned2[:-1] or
                '*' in updated_aligned1[:-1] or '*' in updated_aligned2[:-1]):
            continue
        else:
            score1 = compute_potts_score([aa_to_i[aa] for aa in aligned1.replace('*', '')], potts1, E=False, encode=False)
            score2 = compute_potts_score([aa_to_i[aa] for aa in aligned2.replace('*', '')], potts2, E=False, encode=False)
            score3 = compute_potts_score([aa_to_i[aa] for aa in updated_aligned1.replace('*', '')], potts1, E=False,
                                         encode=False)
            score4 = compute_potts_score([aa_to_i[aa] for aa in updated_aligned2.replace('*', '')], potts2, E=False,
                                         encode=False)
            res.append(f'{seqid},{gene1},{gene2},{minlen},{order},{full_seq},{coord1[0]},{coord1[1]},'
                    f'{coord2[0]},{coord2[1]},{overlaplen},{raw1[i]},{raw2[i]},{final1[i]},{final2[i]},'
                    f'{updated_final1},{updated_final2},{score1},{score2},{score3},{score4}\n')
    return res

