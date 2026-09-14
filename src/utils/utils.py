from biotite.sequence.io.fasta import FastaFile
from biotite.sequence.io.fasta import FastaFile
from pycameox.ilp import AA_TO_I as aa_to_i
import torch.nn.functional as F
import torch
import os

import uuid


class UniqueIdGenerator:
    def __init__(self, prefix="", digits=5, start=0):
        self.prefix = prefix
        self.counter = start
        self.digits = digits

    def get_next_id(self):
        """Generate a sequential numeric ID"""
        self.counter += 1
        formatted_counter = str(self.counter).zfill(self.digits)
        return f"{self.prefix}{formatted_counter}"

    def get_uuid(self):
        """Generate a UUID (universally unique identifier)"""
        return f"{self.prefix}{uuid.uuid4()}"


def read_dna(file, single=True):
    seqs = FastaFile.read(file)
    genes = {}
    for genename, seq in seqs.items():
        genes[genename] = seq
    if single:
        return seq
    else:
        return genes

def check_match_prob(raw, final, hmm, prot, start_pos):
    for i in range(len(final[0])):
        if raw[0][i]!=final[0][i] or raw[0][i] != prot[i+start_pos]:
            print(i, raw[0][i], final[0][i], prot[i+start_pos],
                  hmm.match_prob.loc[i+start_pos+1, raw[0][i]],
                  hmm.match_prob.loc[i+start_pos+1, final[0][i]],
                  hmm.match_prob.loc[i+start_pos+1, prot[i+start_pos]])

def final_with_gaps(final, final_withoutgaps):
    if final.replace('-', '') != final_withoutgaps:
        n_diff = 0
        j = 0
        corrected = []
        for i, x in enumerate(final):
            if x == '-':
                corrected.append('-')
            else:
                corrected.append(final_withoutgaps[j])
                j += 1
                if corrected[-1] != final[i]:
                    n_diff += 1
        return ''.join(corrected)
    else:
        return final

def compute_potts_score(seq, model, E=True, encode=True, beta=1):
    if encode:
        seq_enc = [aa_to_i[aa] for aa in seq]
    else:
        seq_enc = seq
    assert len(seq) == model.L
    # E: closer to 0 is worse
    # pLL: closer to 0 is better
    if E:
        res = model(F.one_hot(torch.tensor(seq_enc), num_classes=model.A).float().unsqueeze(0), beta=beta).item()
    else:
        res = model.pseudolikelihood(
            F.one_hot(torch.tensor(seq_enc), num_classes=model.A).float().unsqueeze(0)).item()
    return res

def substring_coord(substring, string):
    start = string.find(substring)
    return start, start + len(substring)


####################### indel analysis #######################
def compare_seq(seq1, seq2):
    diff = 0
    for i, j in zip(seq1, seq2):
        if i != j and i != '-' and i != '+':
            diff += 1
    return diff

def p_change(row):
    ndel1 = row['raw1'].count('-')
    nins1 = row['raw1'].count('+')
    ndel2 = row['raw2'].count('-')
    nins2 = row['raw2'].count('+')
    # row['init_score'] = row['init1_score'] + row['init2_score']
    # row['optim_score'] = row['optim1_score'] + row['optim2_score']
    # row['pdel1'] = row['ndel1'] / row['overlaplen']
    # row['pins1'] = row['nins1'] / row['overlaplen']
    # row['pdel2'] = row['ndel2'] / row['overlaplen']
    # row['pins2'] = row['nins2'] / row['overlaplen']
    row['p_del'] = (ndel1 + ndel2) / row['overlaplen'] / 2
    row['p_ins'] = (nins1 + nins2) / row['overlaplen'] / 2
    row['mut1'] = compare_seq(row['raw1'], row['optim1'])
    row['mut2'] = compare_seq(row['raw2'], row['optim2'])
    row['mut'] = row['mut1'] + row['mut2']
    row['p_mut'] = row['mut'] / row['overlaplen'] / 2
    return row

###################### Correlation analysis ######################

import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence


def fit_linear(X, y):
    # Fit model
    X = sm.add_constant(X)  # Add intercept
    model = sm.OLS(y, X).fit()
    influence = OLSInfluence(model)
    # Studentized residuals
    studentized_resid = influence.resid_studentized_external
    # Cook's distance
    cooks_d = influence.cooks_distance[0]
    # Identify outliers
    outliers = (abs(studentized_resid) > 3) | (cooks_d > 1)
    res = [sum(outliers)]
    for x in X.columns:
        res.append(model.params[x])
    res += [model.rsquared, model.f_pvalue, X.shape[0]]
    return res

################### check if partial alignment was successfully run ###################
import pandas as pd

# Get the path to gene_info_final_table.tsv relative to this file
_UTILS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_UTILS_DIR))
_GENE_INFO_PATH = os.path.join(_REPO_ROOT, 'data', 'gene_info_final_table.tsv')

# Lazy-load gene_meta only when needed (not at import time)
_gene_meta = None

def _get_gene_meta():
    """Lazy-load gene metadata table."""
    global _gene_meta
    if _gene_meta is None:
        if os.path.exists(_GENE_INFO_PATH):
            _gene_meta = pd.read_csv(_GENE_INFO_PATH, sep='\t')
        else:
            # Return empty DataFrame if table doesn't exist
            _gene_meta = pd.DataFrame(columns=['gene', 'wt_length'])
    return _gene_meta

def is_finished(df):
    gene_meta = _get_gene_meta()
    if gene_meta.empty:
        # If no metadata table, can't determine if finished
        return False
    l1 = gene_meta.loc[gene_meta['gene'] == df.iloc[0, 1], 'wt_length'].values[0]
    l2 = gene_meta.loc[gene_meta['gene'] == df.iloc[0, 2], 'wt_length'].values[0]
    minlen = min([99, l1, l2])
    df.dropna(subset=['order'], inplace=True)
    if df.shape[0] == 0:
        return False
    df['minlen'] = pd.to_numeric(df['minlen'], errors='coerce')
    minlen_by_order = df.groupby('order')['minlen'].max()
    if 12 in minlen_by_order.index and 21 in minlen_by_order.index:
        if minlen_by_order[12] == minlen and minlen_by_order[21] == minlen:
            return True
        else:
            return False
    else:
        return False
