import numpy as np
from potts.potts import get_subset_potts
from src.utils.constants import codon_table
from pycameox.ilp import AA_TO_I as aa_to_i
from src.utils.utils import compute_potts_score, substring_coord
aa_to_i['*'] = 21


def get_entangled_model(full_seq, nt_seq1, nt_seq2, potts1_adjusted, potts2_adjusted):
    # the output seq encoding has both insertions and deletions
    codon_without_gaps1 = [nt_seq1[i * 3:(i + 1) * 3] for i in range(len(nt_seq1) // 3)]
    codon_without_gaps2 = [nt_seq2[i * 3:(i + 1) * 3] for i in range(len(nt_seq2) // 3)]
    # check that both genes stop with a stop codon
    assert codon_without_gaps1[-1] in ['TAA', 'TAG', 'TGA']
    assert codon_without_gaps2[-1] in ['TAA', 'TAG', 'TGA']
    # if gene 1 is on the left and gene 1 ends first
    coord1 = substring_coord(nt_seq1.upper(), full_seq.upper())
    coord2 = substring_coord(nt_seq2.upper(), full_seq.upper())
    start_entangled = max(coord1[0], coord2[0])
    end_entangled = min(coord1[1], coord2[1])
    entangled_idx1 = [i-coord1[0]//3 for i in range(coord1[0] // 3, coord1[1] // 3 ) if start_entangled//3 <= i < end_entangled//3]
    entangled_idx2 = [i-coord2[0]//3 for i in range(coord2[0] // 3, coord2[1] // 3 ) if start_entangled//3 <= i < end_entangled//3]
    codon1 = [codon_without_gaps1[i] for i in entangled_idx1]
    codon2 = [codon_without_gaps2[i] for i in entangled_idx2]
    assert ''.join(codon1)[1:] == ''.join(codon2)[:-1]
    assert len(codon1) == len(codon2)
    aa_seq_enc1 = [aa_to_i[codon_table[c.upper()]] for c in codon_without_gaps1]
    aa_seq_enc2 = [aa_to_i[codon_table[c.upper()]] for c in codon_without_gaps2]
    assert len(aa_seq_enc1) == potts1_adjusted.L
    assert len(aa_seq_enc2) == potts2_adjusted.L
    if len(aa_seq_enc1)>len(entangled_idx1):
        model_entangled1 = get_subset_potts(potts1_adjusted, entangled_idx1, aa_seq_enc1)
    else:
        model_entangled1 = potts1_adjusted
    if len(aa_seq_enc2)>len(entangled_idx2):
        model_entangled2 = get_subset_potts(potts2_adjusted, entangled_idx2, aa_seq_enc2)
    else:
        model_entangled2 = potts2_adjusted
    return model_entangled1, model_entangled2, codon1, codon2

# there seems to be some issue with subsetting potts model that caused the same position to be flipped back and forth
# prevent this temporarily by preventing the same position to be flipped
def greedy_optimize(full_seq, nt_seq1, nt_seq2, potts1_adjusted, potts2_adjusted, verbose_level=1):
    model_entangled1, model_entangled2, codon1, codon2 = get_entangled_model(full_seq, nt_seq1, nt_seq2, potts1_adjusted, potts2_adjusted)
    initial_potts1 = compute_potts_score([aa_to_i[codon_table[c]] for c in codon1], model_entangled1, E=False, encode=False)
    initial_potts2 = compute_potts_score([aa_to_i[codon_table[c]] for c in codon2], model_entangled2, E=False, encode=False)
    all_changes = [x for x in codon_table.keys()]
    assert len(codon1) == len(codon2)
    assert len(codon1) == model_entangled1.L
    assert len(codon2) == model_entangled2.L
    assert model_entangled1.A == 22
    assert model_entangled2.A == 22
    codon1_updated = codon1.copy()
    codon2_updated = codon2.copy()
    n_changed = 1
    it = 0
    last_pos = -1
    # count of number of cycles the potts score didn't improve
    not_improved = 0
    while n_changed > 0 and not_improved <= 5:
        updated_potts1 = compute_potts_score([aa_to_i[codon_table[c]] for c in codon1_updated], model_entangled1,
                                             E=False, encode=False)
        updated_potts2 = compute_potts_score([aa_to_i[codon_table[c]] for c in codon2_updated], model_entangled2,
                                             E=False, encode=False)
        if it!=0:
            if verbose_level >=1:
                print(f"{n_changed}% of all changes accepted")
            print(
                f"updated potts1: {updated_potts1}, updated potts2: {updated_potts2}, sum: {updated_potts1 + updated_potts2}")
            if updated_potts1 + updated_potts2 <= last_score:
                not_improved += 1
        else:
            last_score = initial_potts1 + initial_potts2
            print(f"initial potts1: {initial_potts1}, initial potts2: {initial_potts2}, sum: {initial_potts1+initial_potts2}")
        n_changed = 0
        for it in range(200):
            mut_seq = np.random.randint(0, 2)
            pos = np.random.randint(0, len(codon1_updated))
            if mut_seq == 0: # change a codon in the large gene
                x1 = [aa_to_i[codon_table[c]] for c in all_changes]
                y1_codons = [c[1:] + codon2_updated[pos][2] for c in all_changes]
                y1 = [aa_to_i[codon_table[c]] for c in y1_codons]
                if pos==0:
                    y0 = [aa_to_i['-'] for c in all_changes]
                else:
                    y0_codons = [codon2_updated[pos-1][:2] + c[0] for c in all_changes]
                    y0 = [aa_to_i[codon_table[c]] for c in y0_codons]
                options = [(x1[i], y0[i], y1[i], all_changes[i]) for i in range(len(x1))]
                # remove duplicates in options
                options = list(set(options))
            elif mut_seq == 1: # change a codon in the small gene
                y1 = [aa_to_i[codon_table[c]] for c in all_changes]
                x1_codon = [codon1_updated[pos][0] + c[:2] for c in all_changes]
                x1 = [aa_to_i[codon_table[c]] for c in x1_codon]
                if pos == len(codon1_updated)-1:
                    x2 = [aa_to_i['-'] for c in all_changes]
                else:
                    x2_codon = [c[2] + codon1_updated[pos+1][1:] for c in all_changes]
                    x2 = [aa_to_i[codon_table[c]] for c in x2_codon]
                options = [(x1[i], x2[i], y1[i], all_changes[i]) for i in range(len(x1))]
                # remove duplicates in options
                options = list(set(options))
            if len(options) == 0:
                continue
            aa_enc1 = [aa_to_i[codon_table[c]] for c in codon1_updated]
            aa_enc2 = [aa_to_i[codon_table[c]] for c in codon2_updated]
            if mut_seq == 0:
                if pos!=0:
                    model_mut1 = get_subset_potts(model_entangled1, [pos], aa_enc1)
                    model_mut2 = get_subset_potts(model_entangled2, [pos-1, pos], aa_enc2)
                    mut_psl = [
                        compute_potts_score([a], model_mut1, E=False, encode=False) + \
                        compute_potts_score([b, c], model_mut2, E=False, encode=False)
                    for a, b, c, d in options]
                else:
                    model_mut1 = get_subset_potts(model_entangled1, [pos], aa_enc1)
                    model_mut2 = get_subset_potts(model_entangled2, [pos], aa_enc2)
                    mut_psl = [
                        compute_potts_score([a], model_mut1, E=False, encode=False) + \
                        compute_potts_score([c], model_mut2, E=False, encode=False)
                        for a,b,c,d in options]
            else:
                if pos<len(codon1_updated)-1:
                    model_mut1 = get_subset_potts(model_entangled1, [pos, pos+1], aa_enc1)
                    model_mut2 = get_subset_potts(model_entangled2, [pos], aa_enc2)
                    mut_psl = [
                        compute_potts_score([a, b], model_mut1, E=False, encode=False) + \
                        compute_potts_score([c], model_mut2, E=False, encode=False)
                        for a,b,c,d in options]
                else:
                    model_mut1 = get_subset_potts(model_entangled1, [pos], aa_enc1)
                    model_mut2 = get_subset_potts(model_entangled2, [pos], aa_enc2)
                    mut_psl = [
                        compute_potts_score([a], model_mut1, E=False, encode=False) + \
                        compute_potts_score([c], model_mut2, E=False, encode=False)
                        for a,b,c,d in options]
            max_option = options[np.argmax(mut_psl)]
            max_psl = mut_psl[np.argmax(mut_psl)]
            if mut_seq == 0:
                if pos!=0:
                    orig_psl = compute_potts_score([aa_enc1[pos]], model_mut1, E=False, encode=False) + \
                        compute_potts_score([aa_enc2[pos-1], aa_enc2[pos]], model_mut2, E=False, encode=False)
                else:
                    orig_psl = compute_potts_score([aa_enc1[pos]], model_mut1, E=False, encode=False) + \
                        compute_potts_score([aa_enc2[pos]], model_mut2, E=False, encode=False)
                if max_psl > orig_psl and pos!=last_pos:
                    codon1_updated[pos] = max_option[3]
                    if pos!=0:
                        codon2_updated[pos-1] = codon2_updated[pos-1][:2] + max_option[3][0]
                    codon2_updated[pos] = max_option[3][1:] + codon2_updated[pos][2]
                    if verbose_level==2:
                        print(f'change at seq {mut_seq+1} pos {pos} accepted')
                    last_pos = pos
                    n_changed += 1
            else:
                if pos < len(codon1_updated)-1:
                    orig_psl = compute_potts_score([aa_enc1[pos], aa_enc1[pos+1]], model_mut1, E=False, encode=False) + \
                        compute_potts_score([aa_enc2[pos]], model_mut2, E=False, encode=False)
                else:
                    orig_psl = compute_potts_score([aa_enc1[pos]], model_mut1, E=False, encode=False) + \
                        compute_potts_score([aa_enc2[pos]], model_mut2, E=False, encode=False)
                if max_psl > orig_psl and pos!=last_pos:
                    codon2_updated[pos] = max_option[3]
                    codon1_updated[pos] = codon1_updated[pos][0] + max_option[3][:2]
                    if pos<len(codon1_updated)-2:
                        codon1_updated[pos+1] = max_option[3][2] + codon1_updated[pos+1][1:]
                    if verbose_level==2:
                        print(f'change at seq {mut_seq+1} pos {pos} accepted')
                    n_changed += 1
                    last_pos = pos
    return codon1_updated, codon2_updated, initial_potts1, initial_potts2, updated_potts1, updated_potts2

