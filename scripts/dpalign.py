#!/usr/bin/env python3
"""
Overlapping Gene Design using Dynamic Programming

This script designs overlapping gene pairs by finding optimal alignments
using profile HMMs and Potts models, then optimizing codons for both reading frames.
"""

import argparse
import sys
import os


def main():
    """Main function - imports are inside to speed up --help"""
    # Heavy imports only when actually running
    from src.alignment.scorematrix import ScoreMatrix
    from src.alignment.traceback import find_max
    from src.optimization.entangle_utils import get_entanglement_seq
    from src.utils.utils import read_dna, UniqueIdGenerator
    from biotite.sequence import NucleotideSequence
    from src.models.hmm import HMM, overweight_wt
    from pycameox.ilp import load_ccmpred
    import pickle
    
    sys.setrecursionlimit(10000)
    
    # Get arguments from global parser
    args = parse_arguments()
    gene1 = args.gene1
    gene2 = args.gene2
    rep = args.rep
    p_del = args.del_penalty
    p_ins = args.ins_penalty
    w_match = args.match_weight
    data_dir = args.data_dir
    output_dir = args.output_dir
    min_overlap = args.min_overlap
    max_overlap = args.max_overlap
    
    # Load sequences
    cds1 = read_dna(f'{data_dir}/split_cds/{gene1}.fasta')
    cds2 = read_dna(f'{data_dir}/split_cds/{gene2}.fasta')
    prot1 = str(NucleotideSequence(cds1).translate(complete=True))
    prot2 = str(NucleotideSequence(cds2).translate(complete=True))
    
    # Load Potts models and HMMs
    potts1 = load_ccmpred(f'{data_dir}/raw/{gene1}.raw')
    potts2 = load_ccmpred(f'{data_dir}/raw/{gene2}.raw')
    hmm1 = HMM(f"{data_dir}/hmm/{gene1}.hmm")
    hmm2 = HMM(f"{data_dir}/hmm/{gene2}.hmm")
    
    # Add stop codons if missing
    if not prot1.endswith('*'):
        cds1 = cds1 + 'TAG'
    
    if not prot2.endswith('*'):
        cds2 = cds2 + 'TAG'
    
    # Load or create score matrix
    score_matrix_dir = f'{data_dir}/score_matrix/'
    os.makedirs(score_matrix_dir, exist_ok=True)
    pickle_file = f'{score_matrix_dir}{gene1}_{gene2}_del{p_del}_ins{p_ins}_{w_match}.pkl'
    
    if os.path.exists(pickle_file):
        score_matrix = pickle.load(open(pickle_file, 'rb'))
    else:
        # Adjust HMM parameters
        hmm1.trans_prob['m->d'] += p_del
        hmm2.trans_prob['m->d'] += p_del
        hmm1.trans_prob['d->d'] += p_del
        hmm2.trans_prob['d->d'] += p_del
        hmm1.trans_prob['m->i'] += p_ins
        hmm2.trans_prob['m->i'] += p_ins
        hmm1.trans_prob['i->i'] += p_ins
        hmm2.trans_prob['i->i'] += p_ins
        hmm1.match_prob = overweight_wt(hmm1.match_prob, prot1, w_match)
        hmm2.match_prob = overweight_wt(hmm2.match_prob, prot2, w_match)
        
        score_matrix = ScoreMatrix(prot1, prot2, hmm1, hmm2)
        with open(pickle_file, 'wb') as file:
            pickle.dump(score_matrix, file)
    
    # Create output file
    os.makedirs(output_dir, exist_ok=True)
    output = open(f'{output_dir}/{gene1}_{gene2}_del{p_del}_ins{p_ins}_{w_match}_rep{rep}.csv', 'w')
    
    # Write CSV header
    output.write('seqid,gene1,gene2,minlen,order,full_seq,start1,end1,start2,end2,overlaplen,raw1,raw2,init1,init2,optim1,'
                 'optim2,init1_score,init2_score,optim1_score,optim2_score\n')
    
    # Generate unique IDs
    id_gen = UniqueIdGenerator(prefix=f'{gene1}_{gene2}_', digits=5)
    
    if len(prot1) > len(prot2):
        full_left = True
    else:
        full_left = False
    
    # Design overlapping sequences for both gene orders
    for left in [True, False]:
        if left:
            order = '21'
            max_len = min(min(len(prot1)-2, len(prot2)-1), max_overlap)
        else:
            order = '12'
            max_len = min(min(len(prot1)-1, len(prot2)-1), max_overlap)

        for minlen in range(min_overlap, min(max_len, max_overlap)):
            overlap_idx = find_max(score_matrix, left, minlen)
            dna, raw1, raw2, final1, final2 = score_matrix._trace_back(*overlap_idx[0], last_pos=True)
            res = get_entanglement_seq(dna, raw1, raw2, final1, final2, overlap_idx, minlen, order, 
                                      cds1, cds2, potts1, potts2, gene1, gene2, prot1, prot2, id_gen)
            if len(res) > 0:
                for x in res:
                    output.write(x)
                output.flush()
            else:
                print('no results for', minlen, order)
    
    output.close()


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Design overlapping gene pairs using dynamic programming and Potts models',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s trpC proB 1
  %(prog)s gene1 gene2 1 --del_penalty 5 --ins_penalty 2
  %(prog)s hisF argH 1 --match_weight 0.5

Output:
  CSV file in specified output directory (default: results/) with designed overlapping sequences
        '''
    )
    
    parser.add_argument('gene1', type=str, 
                       help='Name of the first gene')
    parser.add_argument('gene2', type=str, 
                       help='Name of the second gene')
    parser.add_argument('rep', type=int, 
                       help='Replicate number for this gene pair')
    parser.add_argument('--del_penalty', type=float, default=4.0,
                       help='Deletion penalty (higher = fewer deletions, default: 4.0)')
    parser.add_argument('--ins_penalty', type=float, default=1.0,
                       help='Insertion penalty (higher = fewer insertions, default: 1.0)')
    parser.add_argument('--match_weight', type=float, default=0.3,
                       help='Wild-type sequence preference weight (default: 0.3)')
    parser.add_argument('--data_dir', type=str, default='data',
                       help='Data directory containing split_cds/, raw/, hmm/ subdirectories (default: data)')
    parser.add_argument('--output_dir', type=str, default='results',
                       help='Output directory for results (default: results)')
    parser.add_argument('--min_overlap', type=int, default=10,
                       help='Minimum overlap length in amino acids (default: 10)')
    parser.add_argument('--max_overlap', type=int, default=100,
                       help='Maximum overlap length in amino acids (default: 100)')

    return parser.parse_args()


if __name__ == "__main__":
    main()