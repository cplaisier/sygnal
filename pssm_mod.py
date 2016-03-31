"""pssm_mod.py - PSSM functionality
"""
import json


class PSSM:
    def __init__(self, name, nsites, evalue, matrix, genes, de_novo_method='meme'):
        self.name = name
        self.nsites = int(nsites)  # remove the cast after sucessful loading
        self.evalue = evalue
        self.matrix = matrix
        self.genes = genes
        self.de_novo_method = de_novo_method

        self.matches = []
        self.expanded_matches = []
        self.correlatedMatches = {}
        self.permuted_pvalue = None

    def num_genes(self):
        return len(self.genes)

    def add_match(self, factor, confidence):
        self.matches.append({'factor': factor, 'confidence':confidence})


def make_pssm_json(pssm_json):
    return PSSM(pssm_json['name'], pssm_json['nsites'], pssm_json['evalue'],
                pssm_json['matrix'], pssm_json['genes'])


def load_pssms_json(path):
    with open(path, 'r') as infile:
        pssms_json = json.load(infile)
        return {'name': make_pssm_json(pssm_json)
                for name, pssm_json in pssms_json.items()}


def consensus_motif(pssm, lim1=0.6, lim2=0.8, three=0):
    consensus = ''
    for i in range(len(pssm.matrix)):
        consensus += __col_consensus(pssm.matrix, i, lim1, lim2, three)
    return consensus


def __col_consensus(matrix, i, lim1, lim2, three):
    two_base_l = ['Y','R','W','S','K','M']
    three_base_l = ['V','H','D','B']
    conLet = 'N'
    if matrix[i][0] >= lim1:
        conLet = 'A'
    elif matrix[i][1] >= lim1:
        conLet = 'C'
    elif matrix[i][2] >= lim1:
        conLet = 'G'
    elif matrix[i][3] >= lim1:
        conLet = 'T'
    else:
        two_base_c = [matrix[i][1] + matrix[i][3],
                      matrix[i][0] + matrix[i][2],
                      matrix[i][0] + matrix[i][3],
                      matrix[i][1] + matrix[i][2],
                      matrix[i][2] + matrix[i][3],
                      matrix[i][0] + matrix[i][1]]

        three_base_c = [matrix[i][0] + matrix[i][1] + matrix[i][2],
                        matrix[i][0] + matrix[i][1] + matrix[i][3],
                        matrix[i][0] + matrix[i][2] + matrix[i][3],
                        matrix[i][1] + matrix[i][2] + matrix[i][3]]

        pMax = 0
        for k in range(6):
            if two_base_c[k] > pMax:
                pMax = two_base_c[k]
                conLet = two_base_l[k]

        if not pMax > lim2 and three == 1:
            for k in range(4):
                if three_base_c[k] > pMax:
                    pMax = three_base_c[k]
                    conLet = three_base_l[k]
        if not pMax > lim2:
            conLet = 'N'
    return conLet


def pad(s):
    times = 8 - len(s) if len(s) < 8 else 0
    return s + '0' * times


def to_meme_format(pssm):
    """Returns a meme 3 formatted string (letter-probability matrix)"""
    result = 'MOTIF ' + pssm.name + '\n'
    result += 'BL   MOTIF ' + pssm.name + ' width=0 seqs=0\n'
    if pssm.de_novo_method == 'meme':
        nsites = pssm.nsites
        evalue = pssm.evalue
    elif pssm.de_novo_method == 'weeder':
        nsites = len(pssm.nsites)
        evalue = 0.05

    result += 'letter-probability matrix: alength= 4 w= '+str(len(pssm.matrix))+' nsites= '+str(nsites)+' E= '+str(evalue)
    for i in pssm.matrix:
        result += ('\n '+ pad(str(round(float(i[0]),6))) + '  ' + pad(str(round(float(i[1]),6))) +
                   '  ' + pad(str(round(float(i[2]),6))) + '  ' + pad(str(round(float(i[3]),6))))
    return result
