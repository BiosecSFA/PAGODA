import biotite.database.entrez as entrez
import numpy as np
import pandas as pd
import subprocess
import os

class HMM:
    def __init__(self, filename, description=None):
        self.filename = filename
        self.description = description
        self.trans_prob, self.match_prob, self.insert_prob = self.read_hmm_file()
        self.add_stop_codon()

    def subset(self, n, start=True):
        if start:
            self.trans_prob = self.trans_prob.loc[:n, :]
            self.match_prob = self.match_prob.loc[:n, :]
            self.insert_prob = self.insert_prob.loc[:n, :]
        else:
            npos = self.match_prob.shape[0]
            self.trans_prob = self.trans_prob.loc[npos-n:, :].reset_index(drop=True)
            self.match_prob = self.match_prob.loc[npos-n:, :].reset_index(drop=True)
            self.match_prob.index = [x+1 for x in self.match_prob.index]
            self.insert_prob = self.insert_prob.loc[npos-n:, :].reset_index(drop=True)
            self.insert_prob.index = [x+1 for x in self.insert_prob.index]

    def add_stop_codon(self):
        self.match_prob.loc[self.match_prob.index[-1] + 1, :] = np.inf
        self.insert_prob.loc[self.insert_prob.index[-1] + 1, :] = self.insert_prob.loc[self.insert_prob.index[-1], :]
        self.trans_prob.loc[self.trans_prob.index[-1] + 1, :] = np.inf
        self.match_prob['*'] = np.inf
        self.match_prob.loc[self.match_prob.index[-1], '*'] = 0
        self.insert_prob['*'] = np.inf
        self.insert_prob.loc[self.insert_prob.index[-1], '*'] = 0
        self.trans_prob.loc[self.trans_prob.index[-1], ['m->m', 'i->m', 'd->m']] = 0

    def read_hmm_file(self):
        start = False
        data = []
        with open(self.filename, 'r') as file:
            for line in file:
                if line.startswith('HMMER'):
                    continue
                elif line.startswith('HMM'):
                    start = True
                    data.append(line.strip().split())
                elif line.startswith('//') and start:
                    break
                elif start:
                    data.append(line.strip().split())

            trans_prob = []
            match_prob = []
            match_index = []
            ins_prob = []
            for i in range(len(data) // 3):
                trans_prob.append(data[3 * i + 1])
                match_index.append(data[3 * i + 2][0])
                match_prob.append(data[3 * i + 2][1:21])
                ins_prob.append(data[3 * i + 3])

            match_prob = pd.DataFrame(match_prob[1:], index=[x + 1 for x in range(len(match_prob) - 1)],
                                      columns=data[0][1:])
            trans_prob = pd.DataFrame(trans_prob[1:], index=[x for x in range(len(trans_prob) - 1)],
                                      columns=trans_prob[0])
            ins_prob = pd.DataFrame(ins_prob[1:], index=[x + 1 for x in range(len(ins_prob) - 1)],
                                    columns=data[0][1:])
            trans_prob = trans_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            trans_prob.loc[0] = [0, np.inf, np.inf, 0, np.inf, 0, 0]
            trans_prob.loc[len(trans_prob)] = [0, np.inf, np.inf, 0, np.inf, 0, np.inf]
            match_prob = match_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            ins_prob = ins_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            match_prob.loc[1] = np.inf
            match_prob.loc[1,'M'] = 0
            return (trans_prob, match_prob, ins_prob)


def get_hmm(entrez_search, skip_search=False, gene_ids=None):
    if not skip_search:
        id = entrez.search(entrez_search, 'protein')[0]
        print(id)
        filename = entrez.fetch(id, suffix='fa', target_path='data', db_name='protein', ret_type='fasta')
        entrez.fetch(id, suffix='cds.fa', target_path='data', db_name='protein', ret_type='fasta_cds_na')
    else:
        id = gene_ids.loc[entrez_search, 'file_name']
        filename = f'data/{id}.fa'
    # if data/id.hmm exists, load it
    if not os.path.exists(f'data/{id}.hmm'):
        subprocess.run(['phmmer', '-E', '1e-10',
                        '-A', f"{filename.replace('fa','msa')}",
                        filename, 'data/uniref50.fasta'])
        # build the HMM model from the MSA
        subprocess.run(['hmmbuild', "--hand", f"{filename.replace('fa','hmm')}",
                        f"{filename.replace('fa','msa')}"])
    hmm = HMM(f"{filename.replace('fa','hmm')}")
    return id, hmm


def overweight_wt(match_prob, seq, weight):
    adjusted_probs = match_prob.copy()
    for i in range(1, len(seq) + 1):
        prob = match_prob.loc[i].apply(lambda x: np.exp(-x) * (1-weight))
        aa = seq[i - 1]
        prob[aa] = prob[aa] + weight
        prob = prob/prob.sum()
        adjusted_probs.loc[i] = prob.apply(lambda x: -np.log(x) if x!=0 else np.inf)
    return adjusted_probs

