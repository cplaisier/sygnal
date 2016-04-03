#################################################################
# @Program: pssm.py                                             #
# @Version: 1                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  12/4/2009                      #
#################################################################

from math import log
from numpy import array, float64, log10


class pssm:
    """
    A class designed to hold a position specific scoring matrix
    and be able to output this in many different formats
    """
    def __init__(self, name, nsites, eValue, matrix, genes, de_novo_method='meme'):
        self.de_novo_method = de_novo_method
        # each entry in matches should be a dictionary of
        # {'factor':<factor_name>,'confidence':<confidence_value>}
        self.matches = []  
        self.permutedPValue = None
        self.name = name
        self.nsites = nsites
        self.eValue = eValue
        self.matrix = matrix
        self.genes = genes

    def num_genes(self):
        return len(self.genes)

    def add_match(self, factor, confidence):
        if not hasattr(self,'matches'):
            self.matches = []
        self.matches.append({'factor':factor, 'confidence':confidence})

    def get_matches(self):
        if hasattr(self, 'matches') and not self.matches==[]:
            return self.matches
        else:
            return None
    
    def add_expanded_match(self, factor, seedFactor):
        if not hasattr(self,'expandedMatches'):
            self.expandedMatches = []
        self.expandedMatches.append({'factor':factor, 'seedFactor':seedFactor})

    def get_expanded_matches(self):
        if hasattr(self, 'expandedMatches') and not self.expandedMatches==[]:
            return self.expandedMatches
        else:
            return None

    def add_correlated_match(self, subset, factor, rho, pValue):
        if not hasattr(self,'correlatedMatches'):
            self.correlatedMatches = {}
        if not subset in self.correlatedMatches:
            self.correlatedMatches[subset] = []
        self.correlatedMatches[subset].append({'factor':factor, 'rho':rho, 'pValue':pValue})

    def get_correlated_matches(self, subset):
        if hasattr(self, 'correlatedMatches') and subset in self.correlatedMatches and not self.correlatedMatches[subset]==[]:
            return self.correlatedMatches[subset]
        else:
            return None

    def set_permuted_pvalue(self, permutedPValue):
        self.permutedPValue = permutedPValue

    def get_permuted_pvalue(self):
        if hasattr(self,'permutedPValue'):
            return self.permutedPValue
        else:
            return None


def pad(s):
    """Pads the meme nucleotide frequencies with zeros"""
    times = 8 - len(s) if len(s) < 8 else 0
    return s + '0' * times


def to_meme_str(pssm):
    """Returns a meme 3 formatted string (letter-probability matrix)"""
    result = 'MOTIF ' + pssm.name + '\n'
    result += 'BL   MOTIF ' + pssm.name + ' width=0 seqs=0\n'
    if pssm.de_novo_method == 'meme':
        nsites = pssm.nsites
        evalue = pssm.eValue
    elif pssm.de_novo_method == 'weeder':
        nsites = len(pssm.nsites)
        evalue = 0.05

    result += 'letter-probability matrix: alength= 4 w= '+str(len(pssm.matrix))+' nsites= '+str(nsites)+' E= '+str(evalue)
    for i in pssm.matrix:
        result += ('\n '+ pad(str(round(float(i[0]),6))) + '  ' + pad(str(round(float(i[1]),6))) +
                   '  ' + pad(str(round(float(i[2]),6))) + '  ' + pad(str(round(float(i[3]),6))))
    return result


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
