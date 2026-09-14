from Bio import AlignIO
from biotite.sequence.align import Alignment
from biotite.sequence import ProteinSequence
from matplotlib import pyplot as plt
from biotite.sequence import SequenceProfile
from biotite.sequence import graphics
import numpy as np
from src.utils.utils import substring_coord
from src.utils.constants import codon_table
from copy import deepcopy

def check_sequence(seq):
    special_symbols = ['B', 'Z', 'J', 'X', 'U', 'O']
    for symbol in special_symbols:
        if symbol in seq:
            return False
    return True

def get_msa_alignement(msa_file, subset_pos = None, format='stockholm'):
    """
    Read in a stockholm format MSA file and return the biotite alignment object
    :param msa_file:
    :param subset_pos:
    :return:
    """
    # alignment = AlignIO.read(msa_file, format)
    alignments = list(AlignIO.parse(msa_file, 'stockholm'))
    alignment = alignments[0]
    alignment_len = alignment.get_alignment_length()
    if format=='fasta':
        query_positions = [i for i in range(alignment_len)]
    else:
        query_positions = [i for i in range(alignment_len) if
                           alignment.column_annotations['reference_annotation'][i] != '.']
    alignment_array = np.array(alignment)
    if subset_pos is not None:
        query_positions = query_positions[subset_pos]
    filtered_alignment = alignment_array[:, query_positions]
    alignment_seqs = [''.join(x) + '*' for x in filtered_alignment if check_sequence(''.join(x))]
    trace = Alignment.trace_from_strings(alignment_seqs)
    alignment = Alignment([ProteinSequence(x.replace('-', '').replace('X', '-').replace('Z', '-').replace('B', '-').replace('U', '-')) for x in alignment_seqs], trace, 0)
    return alignment


def plot_gene_profile(alignment):
    profile = SequenceProfile.from_alignment(alignment)
    fig = plt.figure(figsize=(50, 4))
    ax = fig.add_subplot(111)
    graphics.plot_sequence_logo(ax, profile, scheme="rainbow")
    ax.axis("off")
    fig.tight_layout()
    plt.savefig(f'figures/temp.png')


def subset_alignment(alignment, idx):
    """
    Get a subset of the alignment
    :param alignment:
    :param idx: the index in the alignment to keep
    :return:
    """
    subset = Alignment([], alignment.trace[idx, :], 0)
    non_empty = []
    for i in range(len(alignment.sequences)):
        t = subset.trace[:, i]
        if len(t[t != -1]) == 0:
            continue
        seq = alignment.sequences[i]
        new_seq = [seq[i] for i in t if i != -1]
        subset.sequences.append(ProteinSequence(''.join(new_seq)))
        # convert non-empty trace to be consecutive integers starting from 0
        non_empty.append(i)
        new_idx = 0
        for j in range(len(t)):
            if t[j] != -1:
                t[j] = new_idx
                new_idx += 1
        subset.trace[:, i] = t
    subset.trace = subset.trace[:, non_empty]
    return subset

def add_msa_stop_codon(msa):
    """
    Add a stop codon to the msa
    :param msa:
    :return:
    """
    for i in range(len(msa.sequences)):
        msa.sequences[i] = ProteinSequence(str(msa.sequences[i]) + '*')
    msa.trace = np.vstack((msa.trace, msa.trace.max(axis=0)+1))
    return msa


from difflib import SequenceMatcher

def pad_aligned(a, b):
    matcher = SequenceMatcher(None, a, b)
    matches = matcher.get_matching_blocks()

    a_pos = b_pos = 0
    padded_b = []

    for match in matches:
        i, j, size = match

        # Add padding for unmatched region
        unmatched_a_len = i - a_pos
        unmatched_b_len = j - b_pos

        # If there's a gap in A but not in B, pad B
        if unmatched_a_len > unmatched_b_len:
            gap = b[b_pos:j] + '+' * (unmatched_a_len - unmatched_b_len)
        else:
            gap = b[b_pos:j]

        padded_b.append(gap)
        # Add the matching block
        padded_b.append(b[j:j+size])

        # Update positions
        a_pos = i + size
        b_pos = j + size

    padded_b_str = ''.join(padded_b)

    # Final padding if A is still longer
    if len(padded_b_str) < len(a):
        padded_b_str += '+' * (len(a) - len(padded_b_str))
    # add gaps from B back to A
    j=0
    padded_a_str = ''
    for i in range(len(padded_b_str)):
        if padded_b_str[i] == '-':
            padded_a_str += '-'
        else:
            assert a[j] == padded_b_str[i] or padded_b_str[i] == '+', f"Mismatch at position {i}: {a[i]} != {padded_b_str[i]}"
            padded_a_str += a[j]
            j += 1

    assert len(padded_a_str) == len(padded_b_str), f"Length mismatch: {len(a)} != {len(padded_b_str)}"
    return padded_b_str


def align_entangled_seq_old(full_seq, start, frame, seq1, seq2, seqaligned1, seqaligned2, wt1, wt2):
    """
    Align the entangled sequences
    :param full_seq: the full sequence of the entangled region
    :param start: the start position of the entangled region in the full sequence
    :param end: the end position of the entangled region in the full sequence
    :param seq1: the first sequence to align
    :param seq2: the second sequence to align
    :param seqaligned1: the aligned first sequence
    :param seqaligned2: the aligned second sequence
    :return: an alignment with the original sequences at top and bottom, and the entangled sequences in between
    """
    pos1_ori = 0
    pos2_ori = 0
    pos1_ent = 0
    pos2_ent = 0
    ngaps1 = 0
    ngaps2 = 0
    aligned_seqs = []
    trace = []
    if frame == 0:
        start = start + 1
        if full_seq[start:start+3] != 'ATG':
            start = start -1
    first_codon = full_seq[start:start+3]
    assert first_codon == 'ATG'
    seqaligned1 = pad_aligned(seq1, seqaligned1)
    seqaligned2 = pad_aligned(seq2, seqaligned2)
    if frame == 0:
        # first gene in frame 0
        coord1 = [0, len(seq1)*3]
        coord2 = [start, len(full_seq)]
    else:
        # first gene in frame 1
        coord1 = [1, len(seq1)*3+1]
        coord2 = [start, len(full_seq)-1]
    if frame == 0:
        a = 0
    else:
        a = 1
    i = -1
    while i < len(full_seq)//3:
        i += 1
        # gene1 only region
        if (i*3 >= coord1[0] and i*3+1 < coord2[0]):
            aligned_seqs.append([wt1[pos1_ori],'-', '-', '-'])
            trace.append([pos1_ori, -1, -1, -1])
            pos1_ori += 1
            pos1_ent += 1
        # gene2 only region
        elif (i*3 >= coord1[1] and i*3 < coord2[1]):
            aligned_seqs.append(['-','-', '-', wt2[pos2_ori]])
            trace.append([-1, -1, -1, pos2_ori])
            pos2_ori += 1
            pos2_ent += 1
        # gene 1 and gene 2 overlap
        elif (i*3 < coord1[1] and i*3 >= coord2[0]):
            # deletion in seq1
            if seqaligned1[pos1_ent] == '-':
                aligned_seqs.append([wt1[pos1_ori], '-', '-', '-'])
                trace.append([pos1_ori, -1, -1, -1])
                pos1_ent += 1
                pos1_ori += 1
                ngaps1 += 1
                i = i-1
            # deletion in seq2
            elif seqaligned2[pos2_ent] == '-':
                aligned_seqs.append(['-', '-', '-', wt2[pos2_ori]])
                trace.append([-1, -1, -1, pos2_ori])
                pos2_ori += 1
                pos2_ent += 1
                ngaps2 += 1
                i = i-1
            #insertion in gene1
            elif seqaligned1[pos1_ent] == '+':
                aligned_seqs.append(['-', seq1[pos1_ent],
                                     seqaligned2[pos2_ent], wt2[pos2_ori]])
                trace.append([-1, pos1_ent-ngaps1, pos2_ent-ngaps2, pos2_ori])
                pos1_ent += 1
                pos2_ent += 1
                pos2_ori += 1
            # insertions in seq2
            elif seqaligned2[pos2_ent] == '+':
                aligned_seqs.append([wt1[pos1_ori],seqaligned1[pos1_ent],
                                     seqaligned2[pos2_ent], '-'])
                trace.append([pos1_ori, pos1_ent-ngaps1, pos2_ent-ngaps2, -1])
                pos1_ent += 1
                pos1_ori += 1
                pos2_ent += 1
            # match
            else:
                aligned_seqs.append([wt1[pos1_ori], seqaligned1[pos1_ent],
                                     seqaligned2[pos2_ent], wt2[pos2_ori]])
                trace.append([pos1_ori, pos1_ent-ngaps1, pos2_ent-ngaps2, pos2_ori])
                pos1_ori += 1
                pos1_ent += 1
                pos2_ori += 1
                pos2_ent += 1

    aligned_seqs = np.array(aligned_seqs).T
    aligned_seqs_nogaps = [''.join(x).replace('-', '').replace('+','') for x in aligned_seqs]
    # modify position to be relative to the entangled region (only for gene1)
    ent_start = min([x[1] for x in trace if x[1] != -1])
    for i in range(len(trace)):
        if trace[i][1] != -1:
            trace[i][1] -= ent_start
    alignment = Alignment([ProteinSequence(x) for x in aligned_seqs_nogaps], np.array(trace), 0)
    return aligned_seqs, alignment


def align_entangled_seq(full_seq, nt_seq1, nt_seq2, final1, final2, raw1, raw2, seq1, seq2):
    """
    Align the entangled sequences
    :param final1:
    :param final2:
    :param raw1:
    :param raw2:
    :param seq1:
    :param seq2:
    :return: an alignment with the original sequences at top and bottom, and the entangled sequences in between
    """
    coord1 = substring_coord(nt_seq1, full_seq.upper())
    coord2 = substring_coord(nt_seq2, full_seq.upper())
    pos1_ori = 0
    pos2_ori = 0
    pos1_ent = 0
    pos2_ent = 0
    ngaps1 = 0
    ngaps2 = 0
    aligned_seqs = []
    trace = []
    i=-1
    while i < len(full_seq)//3:
        i += 1
        # gene1 only region
        if (i*3 >= coord1[0] and i*3 < coord1[1]) and (i*3+1 < coord2[0] or i*3+1 >= coord2[1]):
            # do something
            aligned_seqs.append([seq1[pos1_ori],'-', '-', '-'])
            trace.append([pos1_ori, -1, -1, -1])
            pos1_ori += 1
        # gene2 only region
        elif (i*3+1 >= coord2[0] and i*3+1 < coord2[1]) and (i*3 < coord1[0] or i*3 >= coord1[1]):
            aligned_seqs.append(['-','-', '-', seq2[pos2_ori]])
            trace.append([-1, -1, -1, pos2_ori])
            pos2_ori += 1
            # do something
        # gene 1 and gene 2 overlap
        elif (i*3 >= coord1[0] and i*3 < coord1[1]) and (i*3+1 >= coord2[0] and i*3+1 < coord2[1]):
            #insertion in gene1
            if raw1[pos1_ent] == '+':
                aligned_seqs.append(['-', final1[pos1_ent],
                                     final2[pos2_ent], raw2[pos2_ent]])
                trace.append([-1, pos1_ent-ngaps1, pos2_ent-ngaps2, pos2_ori])
                pos1_ent += 1
                pos2_ent += 1
                pos2_ori += 1
            # insertions in seq2
            elif raw2[pos2_ent] == '+':
                aligned_seqs.append([raw1[pos1_ent], final1[pos1_ent],
                                     final2[pos2_ent], '-'])
                trace.append([pos1_ori, pos1_ent-ngaps1, pos2_ent-ngaps2, -1])
                pos1_ent += 1
                pos1_ori += 1
                pos2_ent += 1
            # deletion in seq1
            elif raw1[pos1_ent] == '-':
                aligned_seqs.append([seq1[pos1_ori], '-', '-', '-'])
                trace.append([pos1_ori, -1, -1, -1])
                pos1_ori += 1
                pos1_ent += 1
                i -= 1
                ngaps1 += 1
            # deletion in seq2
            elif raw2[pos2_ent] == '-':
                aligned_seqs.append(['-', '-', '-', seq2[pos2_ori]])
                trace.append([-1, -1, -1, pos2_ori])
                pos2_ori += 1
                pos2_ent += 1
                i -= 1
                ngaps2 += 1
            # match
            else:
                aligned_seqs.append([raw1[pos1_ent], final1[pos1_ent],
                                     final2[pos2_ent], raw2[pos2_ent]])
                trace.append([pos1_ori, pos1_ent-ngaps1, pos2_ent-ngaps2, pos2_ori])
                pos1_ori += 1
                pos1_ent += 1
                pos2_ori += 1
                pos2_ent += 1
    aligned_seqs = np.array(aligned_seqs).T
    aligned_seqs_nogaps = [''.join(x).replace('-', '') for x in aligned_seqs]
    alignment = Alignment([ProteinSequence(x) for x in aligned_seqs_nogaps], np.array(trace), 0)
    return aligned_seqs, alignment


def get_entangled_alignment(alignment, highlight=None):
    """
    Get only the entangled part of the alignment for visualization
    :param alignment:
    :return:
    """
    idx = []
    for i in range(alignment.trace.shape[0]):
        if sum(alignment.trace[i, :]==-1) < 3:
            idx.append(i)
    idx = [i for i in range(min(idx), max(idx)+1)]
    if highlight is not None:
        new_highlight = []
        for i in highlight:
            new_highlight.append(idx.index(i))
        return subset_alignment(alignment, idx), new_highlight
    else:
        return subset_alignment(alignment, idx)

def add_gaps(alignment, gapped_seq):
    """
    Add gaps to the alignment
    :param alignment:
    :param gapped_seq:
    :return:
    """
    alignment_pos = 0
    new_trace = []
    assert len(gapped_seq.replace('-', '')) == alignment.trace.shape[0]
    for i in range(len(gapped_seq)):
        if gapped_seq[i] == '-':
            new_trace.append([-1]*alignment.trace.shape[1])
        else:
            new_trace.append(list(alignment.trace[alignment_pos, :]))
            alignment_pos += 1
    alignment.trace = np.array(new_trace)
    return alignment


def get_profile_entropy(profile, pos):
    """
    Get the entropy of a position in the profile
    :param profile:
    :param pos:
    :return:
    """
    freq = profile.symbols[pos, :]
    freq = freq / np.sum(freq)
    no_zeros = freq != 0
    pre_entropies = np.zeros(freq.shape)
    pre_entropies[no_zeros] = freq[no_zeros] * np.log2(freq[no_zeros])
    entropy = -np.sum(pre_entropies)
    return entropy, freq


def get_gapped_msa_profile(msa, alignment, orig_trace_i, ent_trace_i):
    """
    Get the profile of the MSA that is aligned to the original sequence
    Also provides the locations of the gaps and the conserved sites that are changed in the entangled sequences
    the index is in the coordinate of the msa (query sequence)
    :param msa_file:
    :param alignment:
    :param orig_trace_i:
    :param ent_trace_i:
    :return:
    """
    msa = add_msa_stop_codon(msa)
    profile = SequenceProfile.from_alignment(msa)
    symbol_to_i = {symbol: i for i, symbol in enumerate(profile.alphabet)}
    # add gaps to to msa where there is insertion
    idx = []
    insertion_sequence = []
    gap_idx = []
    low_entropy_mismatch = []
    orig_seq = alignment.sequences[orig_trace_i]
    ent_seq = alignment.sequences[ent_trace_i]
    start_idx = 0
    nins=0
    while alignment.trace[start_idx, orig_trace_i] == -1:
        start_idx += 1
    highlight_pos = []
    for i in range(start_idx,alignment.trace.shape[0]):
        # MATCH
        if alignment.trace[i, orig_trace_i] != -1 and alignment.trace[i, ent_trace_i] != -1:
            idx.append(alignment.trace[i, orig_trace_i])
            insertion_sequence.append(orig_seq[alignment.trace[i, orig_trace_i]])
            if orig_seq[alignment.trace[i, orig_trace_i]] != ent_seq[alignment.trace[i, ent_trace_i]]:
                entropy, freq = get_profile_entropy(profile, alignment.trace[i, orig_trace_i])
                if entropy < 2.2 and freq[symbol_to_i[ent_seq[alignment.trace[i, ent_trace_i]]]] < 0.05:
                    low_entropy_mismatch.append(alignment.trace[i, orig_trace_i])
                    highlight_pos.append(i)
        # deletion in the entangled partner sequence
        elif alignment.trace[i, orig_trace_i] == -1 and alignment.trace[i, ent_trace_i] == -1:
            # only if the gene has not ended
            if insertion_sequence[-1]!='*':
                insertion_sequence.append('-')
                nins += 1
        # DELETION and before entanglement starts or after entanglement ends
        elif alignment.trace[i, ent_trace_i] == -1:
            # if the gap happens after the first match
            if len(idx)>0:
                insertion_sequence.append(orig_seq[alignment.trace[i, orig_trace_i]])
                gap_idx.append(alignment.trace[i, orig_trace_i])
        # INSERTION
        elif alignment.trace[i, orig_trace_i] == -1:
            insertion_sequence.append('-')
            nins += 1
            # if the entropy in the msa is low and the alignment aa is different
    # remove insertion sequences added after the last match
    gap_idx = [x for x in gap_idx if (x>min(idx) and x<max(idx))]
    insertion_sequence = insertion_sequence[:(len(idx) + len(gap_idx) + nins)]
    subset = subset_alignment(msa, np.arange(min(idx), max(idx)+1))
    gapped_subset = add_gaps(subset,  ''.join(insertion_sequence))
    return SequenceProfile.from_alignment(gapped_subset), gap_idx, low_entropy_mismatch, highlight_pos


def make_entanglement_viz(msa_file1, msa_file2, alignment, plotname,
                          plot_size=(40, 10), make_plot=True):
    """
    Make a visualization of the entangled alignment
    :param msa_file1:
    :param msa_file2:
    :param alignment:
    :param plotname:
    :param plot_size:
    :return:
    """
    fig = plt.figure(figsize=plot_size)
    ax = fig.add_subplot(311)
    profile1, gap1, mismatch1, highlight1 = get_gapped_msa_profile(msa_file1, alignment, 0, 1)
    profile2, gap2, mismatch2, highlight2 = get_gapped_msa_profile(msa_file2, alignment, 3, 2)
    highlight = list(set(highlight1 + highlight2))
    entangled_alignment, highlight = get_entangled_alignment(alignment, highlight)
    graphics.plot_sequence_logo(ax, profile1, scheme="rainbow")
    ax.axis("off")
    fig.tight_layout()
    ax = fig.add_subplot(312)
    graphics.plot_alignment_type_based(
        ax,
        entangled_alignment, symbols_per_line=entangled_alignment.trace.shape[0],
        show_numbers=False, show_line_position=True)
    for text_obj in ax.texts:
        pos = int(text_obj.get_position()[0])
        # Only consider sequence characters (not position numbers)
        if pos in highlight:
            text_obj.set_fontweight('bold')
            text_obj.set_fontsize(20)
            text_obj.set_color('black')
        else:
            text_obj.set_fontsize(18)
            text_obj.set_color('white')
    ax.axis("off")
    fig.tight_layout()
    ax = fig.add_subplot(313)
    graphics.plot_sequence_logo(ax, profile2, scheme="rainbow")
    ax.axis("off")
    fig.tight_layout()
    plt.savefig(f'figures/{plotname}.png')
    plt.close()
    return len(highlight), len(mismatch1), len(gap1), len(mismatch2), len(gap2)


def make_entanglement_viz_old(alignment, plotname, profile1, highlight1,profile2, highlight2,
                          plot_size=(40, 10)):
    """
    Make a visualization of the entangled alignment
    :param msa_file1:
    :param msa_file2:
    :param alignment:
    :param plotname:
    :param plot_size:
    :return:
    """
    fig = plt.figure(figsize=plot_size)
    ax = fig.add_subplot(311)
    highlight = list(set(highlight1 + highlight2))
    entangled_alignment, highlight = get_entangled_alignment(alignment, highlight)
    entangled_len = entangled_alignment.trace.shape[0]
    new_profile1 = SequenceProfile(profile1.symbols[:entangled_len, :], profile1.gaps[:entangled_len],
                                   profile1.alphabet)
    graphics.plot_sequence_logo(ax, new_profile1, scheme="rainbow")
    ax.axis("off")
    fig.tight_layout()
    ax = fig.add_subplot(312)
    graphics.plot_alignment_type_based(
        ax,
        entangled_alignment, symbols_per_line=entangled_alignment.trace.shape[0],
        show_numbers=False, show_line_position=True)
    for text_obj in ax.texts:
        pos = int(text_obj.get_position()[0])
        # Only consider sequence characters (not position numbers)
        if pos in highlight:
            text_obj.set_fontweight('bold')
            text_obj.set_fontsize(20)
            text_obj.set_color('black')
        else:
            text_obj.set_fontsize(18)
            text_obj.set_color('white')
    ax.axis("off")
    fig.tight_layout()
    ax = fig.add_subplot(313)
    graphics.plot_sequence_logo(ax, profile2, scheme="rainbow")
    ax.axis("off")
    fig.tight_layout()
    plt.savefig(f'figures/{plotname}.png')
    plt.close()
    return len(highlight)


def update_final(final, updated_codons):
    codons = deepcopy(updated_codons)
    updated_final = []
    for i in range(len(final)):
        if final[i] == '-':
            updated_final.append('-')
        else:
            updated_final.append(codon_table[codons.pop(0)])
    return ''.join(updated_final)

# msa = get_msa_alignement('datasets/hicA/hicA.msa', subset_pos=None, format='stockholm')
