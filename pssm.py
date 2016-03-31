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
    # Initialize the pssm
    def __init__(self, name, nsites, eValue, matrix, genes, de_novo_method='meme'):
        self.de_novo_method = de_novo_method
        self.matches = [] # each entry should be a dictionary of {'factor':<factor_name>,'confidence':<confidence_value>}
        self.permutedPValue = None
        self.name = name
        self.nsites = nsites
        self.eValue = eValue
        self.matrix = matrix
        self.genes = genes

    # Returns the name of the PSSM
    def getName(self):
        return self.name

    # Sets the name of the PSSM
    def setName(self,name):
        self.name = name

    # Returns the name of the PSSM
    def getMethod(self):
        return self.de_novo_method
    
    # Returns the name of the PSSM
    def setMethod(self, de_novo_method):
        self.de_novo_method = de_novo_method

    # Returns the E-value of the PSSM
    def getEValue(self):
        return self.eValue

    # Returns the number of sites for the PSSM
    def getNSites(self):
        return self.nsites

    # Returns the number of genes of the PSSM
    def getNumGenes(self):
        return len(self.genes)

    # Returns the genes of the PSSM
    def getGenes(self):
        return self.genes

    # Returns the matrix
    def getMatrix(self):
        return self.matrix

    # Pads the meme nucleotide frequencies with zeros
    def padMe(self,str1):
        if len(str1)<8:
            for i in range(8-len(str1)):
                str1 += '0'
        return str1

    # Retunrs a log-odds value
    def logOdds(self, p, f):
        p = float(p)
        f = float(f)
        if p==float(0):
            v1 = str(int(round(log(float(1E-300)/f,2),0)))
        else:
            v1 = str(int(round(log(p/f,2),0)))
        if len(v1)<6:
            for i in range(6-len(v1)):
                v1 = ' ' + v1
        return v1

    # Returns a meme 3 formatted string (letter-probability matrix)
    def getMemeFormatted(self,atFreq=0.25,cgFreq=0.25):
        memeFormatted = 'MOTIF '+self.name+'\n'
        memeFormatted += 'BL   MOTIF '+self.name+' width=0 seqs=0\n'
        if self.de_novo_method=='meme':
            nsites = self.nsites
            eValue = self.eValue
        elif self.de_novo_method=='weeder':
            nsites = len(self.nsites)
            eValue = 0.05
        memeFormatted += 'letter-probability matrix: alength= 4 w= '+str(len(self.matrix))+' nsites= '+str(nsites)+' E= '+str(eValue)
        for i in self.matrix:
            memeFormatted += '\n '+self.padMe(str(round(float(i[0]),6)))+'  '+self.padMe(str(round(float(i[1]),6)))+'  '+self.padMe(str(round(float(i[2]),6)))+'  '+self.padMe(str(round(float(i[3]),6)))
        return memeFormatted

    # Returns a mast formatted string (log-odds matrix)
    def getMastFormatted(self,atFreq=0.25,cgFreq=0.25):
        mastFormatted = 'log-odds matrix: alength= 4 w= '+str(len(self.matrix))
        for i in self.matrix:
            mastFormatted += '\n '+self.logOdds(i[0],atFreq)+'  '+self.logOdds(i[1],cgFreq)+'  '+self.logOdds(i[2],cgFreq)+'  '+self.logOdds(i[3],atFreq)
        return mastFormatted

   # Returns the consensus word for a motif
    def getConsensusMotif(self, lim1=0.6, lim2=0.8, three=0):
        consensus = ''
        for i in range(len(self.matrix)):
            consensus += self.colConsensus(self.matrix, i, lim1, lim2, three)
        return consensus

    # Returns the consensus letter for a motif column
    def colConsensus(self, pssm, i, lim1, lim2, three):
        two_base_l = ['Y','R','W','S','K','M']
        three_base_l = ['V','H','D','B']
        conLet = 'N'
        if float(pssm[i][0])>=lim1:
            conLet = 'A'
        elif float(pssm[i][1])>=lim1:
            conLet = 'C'
        elif float(pssm[i][2])>=lim1:
            conLet = 'G'
        elif float(pssm[i][3])>=lim1:
            conLet = 'T'
        else:
            two_base_c = [float(pssm[i][1])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][2]), float(pssm[i][0])+float(pssm[i][3]), float(pssm[i][1])+float(pssm[i][2]), float(pssm[i][2])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][1])]
            three_base_c = [float(pssm[i][0])+float(pssm[i][1])+float(pssm[i][2]), float(pssm[i][0])+float(pssm[i][1])+float(pssm[i][3]), float(pssm[i][0])+float(pssm[i][2])+float(pssm[i][3]), float(pssm[i][1])+float(pssm[i][2])+float(pssm[i][3])]
            pMax = 0
            for k in range(0,6):
                if two_base_c[k] > pMax:
                    pMax = two_base_c[k]
                    conLet = two_base_l[k]
            if not pMax>lim2 and three==1:
                for k in range(0,4):
                    if three_base_c[k] > pMax:
                        pMax = three_base_c[k]
                        conLet = three_base_l[k]
            if not pMax>lim2:
                conLet = 'N'
        return conLet

    # Plot a PSSM using weblogo
    def plot(self, fileName):
        dist = numpy.array( self.getMatrix(), numpy.float64 ) 
        data = LogoData.from_counts(corebio.seq.unambiguous_dna_alphabet, dist*100)
        options = LogoOptions()
        options.color_scheme = colorscheme.nucleotide
        format = LogoFormat(data, options)
        fout = open(fileName, 'w')
        png_formatter(data, format, fout)
        fout.close()

    # Add a match to a pssm
    def addMatch(self, factor, confidence):
        if not hasattr(self,'matches'):
            self.matches = []
        self.matches.append({'factor':factor, 'confidence':confidence})

    # Add a match to a pssm
    def getMatches(self):
        if hasattr(self, 'matches') and not self.matches==[]:
            return self.matches
        else:
            return None
    
    # Add an expanded match to a pssm
    def addExpandedMatch(self, factor, seedFactor):
        if not hasattr(self,'expandedMatches'):
            self.expandedMatches = []
        self.expandedMatches.append({'factor':factor, 'seedFactor':seedFactor})

    # Get expanded matches for a pssm
    def getExpandedMatches(self):
        if hasattr(self, 'expandedMatches') and not self.expandedMatches==[]:
            return self.expandedMatches
        else:
            return None

    # Add a correlated match to a pssm
    def addCorrelatedMatch(self, subset, factor, rho, pValue):
        if not hasattr(self,'correlatedMatches'):
            self.correlatedMatches = {}
        if not subset in self.correlatedMatches:
            self.correlatedMatches[subset] = []
        self.correlatedMatches[subset].append({'factor':factor, 'rho':rho, 'pValue':pValue})

    # Get correlated matches for a pssm
    def getCorrelatedMatches(self, subset):
        if hasattr(self, 'correlatedMatches') and subset in self.correlatedMatches and not self.correlatedMatches[subset]==[]:
            return self.correlatedMatches[subset]
        else:
            return None

    # Add a match to a pssm
    def setPermutedPValue(self, permutedPValue):
        self.permutedPValue = permutedPValue

    # Add a match to a pssm
    def getPermutedPValue(self):
        if hasattr(self,'permutedPValue'):
            return self.permutedPValue
        else:
            return None

