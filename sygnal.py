#!/usr/bin/env python
#################################################################
# @Program: offYerBackV2.py                                     #
# @Version: 2                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  5/21/2012                      #
#################################################################

from math import log10
import cPickle, os, re, sys
from sys import stdout
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
import subprocess
from shutil import rmtree
from copy import deepcopy

import rpy2.robjects as robj
from rpy2.robjects import FloatVector


# Custom offYerBack libraries
from cMonkeyWrapper import cMonkeyWrapper
from pssm import pssm
from miRvestigator import miRvestigator
from tomtom import tomtom
# Functions needed to run this script
from utils import *
from sys import stdout, exit
#from weeder import *
import gzip


#################################################################
## Parameters                                                  ##
#################################################################

# For MEME analysis
bgFile = 'seqs/bgFile.meme'
nMotifs = 2
regions = [ 'upstream' ]
motifWidth = { 'upstream': [6, 12] }
revComp = { 'upstream': True }
# Parameters for filtering results
maxEValue = 10
leo_nb_AtoB = 0.5 # Equates to ~3 times better model fit than next best model
mlogp_M_AtoB = 0.05
rhoCut = 0.3
pVCut = 0.05
percTargets = 0.1
postProcessedFile = 'postProcessed_gbmTCGA_pita.csv'
fpcFile = 'biclusterFirstPrincComponents.csv'
subsets = ['all'] # Might be nice to include subtypes
subsetsPos = { 'all': [0,422] } # Might be nice to include subtypes
randPssmsDir = 'randPSSMs'

SYNONYM_PATH = '../synonymThesaurus.csv.gz'
MIRNA_FASTA_PATH = 'miRNA/hsa.mature.fa'  # Need to update this to the latest database
RATIOS_PATH = '../gbmTCGA_exprMat_medianFiltered.tsv'
ALL_RATIOS_PATH = 'expression/gbmTCGA_exprMat.tsv'
CMONKEY2_RUNDB = '../out/cmonkey_run.db'

#################################################################
## Functions                                                   ##
#################################################################

# Run meme and get the output into PSSMs
def meme(num, seqFile=None, bgFile=None, nMotifs=1, minMotifWidth=6, maxMotifWidth=12, revComp=True, seed=None):
    # Arguments for meme
    memeArgs = str(seqFile)+' -bfile '+ str(bgFile) +' -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw ' + str(minMotifWidth) + ' -maxw ' + str(maxMotifWidth) + ' -nmotifs ' + str(nMotifs)
    if revComp==True:
        memeArgs += ' -revcomp'
    if not seed==None:
        memeArgs += ' -cons ' + str(seed)
    print memeArgs
    memeProc = Popen("meme " + memeArgs, shell=True,stdout=PIPE)
    output = memeProc.communicate()[0].split('\n')

    PSSMs = []
    # Now iterate through output and save data
    for i in range(len(output)):
        splitUp1 = output[i].strip().split(' ')
        if splitUp1[0]=='Motif' and splitUp1[2]=='position-specific' and splitUp1[3]=='probability':
            i += 2 # Skip the separator line, go to the summary line
            splitUp = output[i].strip().split(' ')
            width = int(splitUp[5])
            sites = splitUp[7]
            eValue = splitUp[9]
            matrix = []
            for j in range(width):
                i += 1
                matrix += [[float(let) for let in output[i].strip().split(' ') if let]]
            PSSMs.append(pssm(biclusterName=str((seqFile.split('_')[1]).split('.')[0])+'_motif'+str(splitUp1[1])+'_meme',nsites=sites,eValue=eValue,pssm=matrix,genes=[], de_novo_method='meme'))
    clusterMemeRuns[num] = PSSMs

# Wrapper function to run the meme function using a multiprocessing pool
def runMeme(i):
    meme(i,seqFile=clusterFileNames[i],bgFile=memeVars['bgFile'],nMotifs=memeVars['nMotifs'],minMotifWidth=memeVars['minMotifWidth'], maxMotifWidth=memeVars['maxMotifWidth'], revComp=memeVars['revComp'])


# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
def weeder(bicluster, seqFile=None, bgFile='MM', size='small', enriched='T50', revComp=False):
    print seqFile

    # First run weederTFBS
    weederArgs = str(seqFile)+' '+str(bgFile)+' '+str(size)+' '+str(enriched)
    if revComp==True:
        weederArgs += ' S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weederlauncher " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()

    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqFile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads will be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    if size=='medium':
        # Scroll to where the 10bp reads wll be
        while 1:
            outLine = outLines.pop(0)
            if not outLine.find('Searching for motifs of length 10') == -1:
                break

        # Get top hit of 10bp look for "1)"
        while 1:
            outLine = outLines.pop(0)
            if not outLine.find('1) ') == -1:
                break
        hitBp[10] = outLine.strip().split(' ')[1:]

    # Scroll to where the 10bp reads will be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break

    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    motifId = 1
    biclusterId = str((seqFile.split('_')[1]).split('.')[0])
    while 1:
        if len(outLines)<=1:
            break
        if revComp==True:
            name = outLines.pop(0).strip() # Get match
            name += '_'+outLines.pop(0).strip()
        if revComp==False:
            name = outLines.pop(0).strip() # Get match
        if not name.find('(not highest-ranking)') == -1:
            break
        # Get redundant motifs
        outLines.pop(0)
        redMotifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [i for i in col.strip().split('\t')[1].split(' ') if i]
            colSum = 0
            for i in nums:
                colSum += int(i.strip())
            matrix += [[ float(nums[0])/float(colSum), float(nums[1])/float(colSum), float(nums[2])/float(colSum), float(nums[3])/float(colSum)]]
            col = outLines.pop(0)
        PSSMs += [pssm(biclusterName=str(biclusterId)+'_motif'+str(motifId)+'_weeder',nsites=instances,eValue=hitBp[len(matrix)][1],pssm=matrix,genes=redMotifs, de_novo_method='weeder')]
        motifId += 1
    weederResults[bicluster] = PSSMs

# Wrapper function to run weeder using a multiprocessing pool
def runWeeder(i):
    weeder(i,seqFile=clusterFileNames[i],bgFile=weederVars['bgFile'], size=weederVars['size'], enriched=weederVars['enriched'], revComp=weederVars['revComp'])


# Run TomTom on the files
def TomTom(num, distMeth='ed', qThresh='1', minOverlap=6):
    # Arguments for tomtom
    tomtomArgs = ' -dist '+str(distMeth)+' -o tmp/tomtom_out -text -thresh '+str(qThresh)+' -min-overlap '+str(minOverlap)+' -verbosity 1 tmp/query'+str(num)+'.tomtom tmp/target'+str(num)+'.tomtom'
    print tomtomArgs
    with open('tmp/stderr_'+str(num)+'.out','w') as errOut:
        tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE, stderr=errOut)
        with open('tmp/tomtom_out/tomtom'+str(num)+'.out', 'w') as outputFile:
            output = tomtomProc.communicate()[0]
            outputFile.write(output)


# Wrapper function to run TomTom using multiprocessing pool
def runTomTom(i):
    TomTom(i, distMeth='ed', qThresh='1', minOverlap=6) #blic5


def phyper(q, m, n, k, lower_tail=False):
    """calls the R function phyper, input values are lists and returns a list"""
    r_phyper = robj.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return [f for f in 
            r_phyper(FloatVector(q), FloatVector(m), FloatVector(n), FloatVector(k), **kwargs)]


# Get a correlation p-value from R
def correlation(a1, a2):
    """
    Calculate the correlation coefficient and p-value between two variables.
    Input: Two arrays of float or integers.
    Returns: Corrleation coefficient and p-value.
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    runMe.append('c1 = cor.test(c('+','.join([str(i) for i in a1])+'),c('+','.join([str(i) for i in a2])+'))')
    runMe.append('c1$estimate')
    runMe.append('c1$p.value')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    rho = float(splitUp[1])
    pValue = float((splitUp[2].split(' '))[1])
    return [rho, pValue]

# Compute survival p-value from R
def survival(survival, dead, pc1, age):
    """
    Calculate the survival correlation coefficient and p-value between two variables.
    Input: Four arrays of float or integers.
    Returns:
    """
    # Fire up R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    # Make the data into an R matrix
    runMe.append('library(survival)')
    runMe.append('s1 = c('+','.join([str(i) for i in survival])+')')
    runMe.append('dead1 = c('+','.join(['\''+str(i)+'\'' for i in dead])+')')
    runMe.append('pc1 = c('+','.join([str(i) for i in pc1])+')')
    runMe.append('age1 = c('+','.join([str(i) for i in age])+')')
    runMe.append('scph1 = summary(coxph(Surv(s1,dead1==\'DEAD\') ~ pc1))')
    runMe.append('scph2 = summary(coxph(Surv(s1,dead1==\'DEAD\') ~ pc1 + age1))')
    runMe.append('scph1$coef[1,4]')
    runMe.append('scph1$coef[1,5]')
    runMe.append('scph2$coef[1,4]')
    runMe.append('scph2$coef[1,5]')
    runMe = '\n'.join(runMe)+'\n'
    out = rProc.communicate(runMe)
    # Process output
    splitUp = out[0].strip().split('\n')
    z1 = float((splitUp[0].split(' '))[1])
    pValue1 = float((splitUp[1].split(' '))[1])
    z2 = float((splitUp[2].split(' '))[1])
    pValue2 = float((splitUp[3].split(' '))[1])
    return [[z1, pValue1], [z2, pValue2]]

# Compare miRNA names
def compareMiRNANames(a, b):
    """
    Match miRNA names.
    Input: Two miRNA names to be compared.
    Returns:
    """
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
            re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(b+'[a-oq-z]?(-\d)?(-5p)?$')
        if re1.match(a):
            return 1
    return 0


def sygnal_init():
    """Preparing directories for output"""
    if not os.path.exists('output'):
        os.makedirs('output')


def read_synonyms():
    """ Load synonym thesaurus to get UCSC ID to Entrez ID"""
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    id2entrez = {}
    entrez2id = {}

    with gzip.open(SYNONYM_PATH, 'r') as infile:
        for line in infile:
            splitUp = line.strip().split(',')
            id = splitUp[0]
            entrez = [i for i in splitUp[1].split(';') if is_number(i)]
            if len(entrez)==1:
                id2entrez[id] = entrez[0]
                if not entrez[0] in entrez2id:
                    entrez2id[entrez[0]] = [id]
                else:
                    entrez2id[entrez[0]].append(id)
    return id2entrez, entrez2id


def miRNA_mappings():
    """Create a dictionary to convert the miRNAs to there respective ids"""
    miRNAIDs = {}
    miRNAIDs_rev = {}

    with open(MIRNA_FASTA_PATH, 'r') as infile:
        for line in infile:
            splitUp = line.split(' ')
            if not splitUp[1] in miRNAIDs_rev:
                miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()

            if not splitUp[0].lower() in miRNAIDs:
                miRNAIDs[splitUp[0].lower()] = splitUp[1]
            else:
                print 'Uh oh!', splitUp
    return miRNAIDs, miRNAIDs_rev


def read_cmonkey_run(path):
    """ Load cMonkey Object - turns cMonkey data into objects"""
    output_path = 'output/c1.pkl'

    # If this is the first time then load from the RData file
    if not os.path.exists(output_path):
        c1 = cMonkeyWrapper(path, meme_upstream=False, weeder_upstream=False,
                            weeder_3pUTR=False, tfbs_db=False,
                            pita_3pUTR=False,
                            targetscan_3pUTR=False) #, geneConv=entrez2id)

        with open(output_path, 'wb') as pklFile:
            cPickle.dump(c1, pklFile)

    # Otherwise load the dumped pickle file if it exists
    else:
        print 'Loading precached cMonkey object...'
        with open(output_path, 'rb') as pklFile:
            c1 = cPickle.load(pklFile)
        print 'Done.\n'

    return c1


# Initialize sygnal output directory and conversion dictionaries
sygnal_init()
id2entrez, entrez2id = read_synonyms()
miRNAIDs, miRNAIDs_rev = miRNA_mappings()

clusterFileNames = {}

if not os.path.exists('output/c1_all.pkl'):
    c1 = read_cmonkey_run(CMONKEY2_RUNDB)
    
    #################################################################
    ## Fill in the missing parts                                   ##
    #################################################################
     #  Check to see if all parts are there:                       #
     #   A. Upstream motifs (MEME)                                 #
     #   B. Upstream motif (Weeder)                                #
     #   C. 3' UTR Weeder-miRvestigator (Weeder)                   #
     #   D. 3' UTR PITA (Set Enrichment)                           #
     #   E. 3' UTR TargetScan (Set Enrichment)                     #
     ###############################################################

    ## A. Upstream motifs (MEME) ##
    # If MEME hasn't been run on the biclusters upstream sequences then do so
    if not c1.meme_upstream:
        print 'Running MEME on Upstreams:'
        # Use already run MEME results if they exist
        if not os.path.exists('output/meme_upstream.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/meme/fasta'):
                os.makedirs('tmp/meme/fasta')

            # Run MEME on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for MEME Upstream run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/meme/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqsUpstream(b1)
                if len(seqs) > 0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    with open(clusterFileName, 'w') as fastaFile:
                        fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

            # Where all the results will be stored
            clusterMemeRuns = mgr.dict()

            # Parameters to use for running MEME
            memeVars = mgr.dict( { 'bgFile': bgFile, 'nMotifs': nMotifs, 'minMotifWidth': motifWidth['upstream'][0], 'maxMotifWidth': motifWidth['upstream'][1], 'revComp': revComp['upstream'] } )

            # Then run MEME using all cores available
            print 'Running MEME on Upstream sequences...'
            cpus = cpu_count()
            print 'There are %d CPUs available.' % cpus
            pool = Pool(processes=cpus)
            pool.map(runMeme,[i for i in o1])
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            with open('output/meme_upstream.pkl','wb') as pklFile:
                cPickle.dump(deepcopy(clusterMemeRuns), pklFile)
        else:
            print 'Loading from precached object...'
            with open('output/meme_upstream.pkl','rb') as pklFile:
                clusterMemeRuns = cPickle.load(pklFile)

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in clusterMemeRuns.keys():
            for pssm1 in clusterMemeRuns[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('meme')
                b1.addPssmUpstream(pssm1)

        print 'Done with MEMEing.\n'

        # MEME upstream has been run on cMonkey run
        c1.meme_upstream = True

    ## B. Upstream motifs (Weeder) ##
    # If Weeder hasn't been run on the biclusters upstream sequences then do so
    if not c1.weeder_upstream:
        print 'Running Weeder on Upstreams:'
        # If this has been run previously just load it up
        if not os.path.exists('output/weeder_upstream.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/weeder/fasta'):
                os.makedirs('tmp/weeder/fasta')

            # Run Weeder on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for Upstream Weeder re-run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/weeder/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqsUpstream(b1)
                if len(seqs)>0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    with open(clusterFileName, 'w') as fastaFile:
                        fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

            # Where all the results will be stored
            weederResults = mgr.dict()

            # Parameters to use for running Weeder
            # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
            weederVars = mgr.dict( { 'bgFile': 'HS', 'size': 'small', 'enriched': 'T50', 'revComp': True } )

            # Run this using all cores available
            print 'Running Weeder...'
            cpus = cpu_count()
            print 'There are %d CPUs available.' % cpus
            pool = Pool(processes=cpus)
            pool.map(runWeeder, [i for i in o1])
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            with open('output/weeder_upstream.pkl','wb') as pklFile:
                cPickle.dump(deepcopy(weederResults), pklFile)
        else:
            print 'Loading from precached object...'
            with open('output/weeder_upstream.pkl','rb') as pklFile:
                weederResults = cPickle.load(pklFile)

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in weederResults.keys():
            for pssm1 in weederResults[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('weeder')
                b1.addPssmUpstream(pssm1)
        print 'Done with Weedering.\n'

        # Weeder upstream has been run on cMonkey run
        c1.weeder_upstream = True

    # C. 3' UTR Weeder-miRvestigator (Weeder)
    # If Weeder hasn't been run on the biclusters 3' UTR sequences then do so
    if not c1.weeder_3pUTR:
        print 'Running Weeder on 3\' UTRs:'
        # If this has been run previously just load it up
        if not os.path.exists('output/weeder_3pUTR.pkl'):
            # Make needed directories
            if os.path.exists('tmp'):
                rmtree('tmp')
            if not os.path.exists('tmp/weeder/fasta'):
                os.makedirs('tmp/weeder/fasta')

            # Run MEME on all biclusters
            mgr = Manager()
            clusterFileNames = mgr.dict()
            o1 = []

            # First make fasta files for all biclusters
            print 'Making Files for Upstream Weeder re-run...'
            for b1 in c1.getBiclusters():
                clusterFileName = 'tmp/weeder/fasta/bicluster_'+str(b1)+'.fasta'
                seqs = c1.getBiclusterSeqs3pUTR(b1)
                if len(seqs)>0:
                    o1.append(b1)
                    clusterFileNames[b1] = clusterFileName
                    with open(clusterFileName, 'w') as fastaFile:
                        fastaFile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

            # Where all the results will be stored
            weederResults = mgr.dict()

            # Parameters to use for running Weeder
            # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
            weederVars = mgr.dict( { 'bgFile': 'HS3P', 'size': 'small', 'enriched': 'T50', 'revComp': False } )

            # Run this using all cores available
            print 'Running Weeder...'
            cpus = cpu_count()
            print 'There are', cpus,'CPUs avialable.'
            pool = Pool(processes=cpus)
            pool.map(runWeeder,[i for i in o1])
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            with open('output/weeder_3pUTR.pkl','wb') as pklFile:
                cPickle.dump(deepcopy(weederResults),pklFile)
        else:
            print 'Loading from precached object...'
            with open('output/weeder_3pUTR.pkl','rb') as pklFile:
                weederResults = cPickle.load(pklFile)

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i in weederResults.keys():
            for pssm1 in weederResults[i]:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('weeder')
                b1.addPssm3pUTR(pssm1)
        print 'Done with 3\'UTR Weedering.\n'

        # Weeder 3pUTR has been run on cMonkey run
        c1.weeder_3pUTR = True

    # D. Upstream TFBS DB enrichment Analysis
    # If tfbs_db enrichment hasn't been calculated for the biclusters then do so
    if not c1.tfbs_db:
        print 'Running TFBS_DB on Biclusters:'
        # If this has been run previously just load it up
        if not os.path.exists('output/tfbs_db.pkl'):
            # Get ready for multiprocessor goodness
            mgr = Manager()

            # Get a list of all genes in the biclusters
            print 'Get a list of all genes in run...'
            tmpDict = c1.getBiclusters()
            genesInBiclusters = []
            for bicluster in tmpDict:
                genes = tmpDict[bicluster].getGenes()
                for gene in genes:
                    if not gene in genesInBiclusters:
                        genesInBiclusters.append(gene)
            biclusters = mgr.dict(tmpDict)
            del tmpDict

            # Load up TFBS_DB TF ids
            if not os.path.exists('TF/tfbs_db.pkl'):
                print 'Loading TFBS_DB predictions...'
                tmpList = []
                tmpDict = {}
                with gzip.open('TF/tfbsDb_5000_gs.csv.gz','r') as infile:
                    inLines = [i.strip().split(',') for i in infile.readlines() if i.strip()]

                for line in inLines:
                    if line[1] in genesInBiclusters:
                        if not line[1] in tmpList:
                            tmpList.append(line[1])
                        if not line[0] in tmpDict:
                            tmpDict[line[0]] = []
                        tmpDict[line[0]].append(line[1])

                with open('TF/tfbs_db.pkl','wb') as pklFile:
                    cPickle.dump(tmpDict, pklFile)
                    cPickle.dump(tmpList, pklFile)

            # Otherwise load the dumped pickle file if it exists
            else:
                print 'Loading pickled TFBS_DB predictions...'
                with open('TF/tfbs_db.pkl','rb') as pklFile:
                    tmpDict = cPickle.load(pklFile)
                    tmpList = cPickle.load(pklFile)

            # Setup for analysis
            predDict = mgr.dict(tmpDict)
            pred_totalTargets = mgr.list()
            pred_totalTargets.append(set(tmpList))
            del tmpDict
            del tmpList
            print 'TFBS_DB has %d TFs.' % len(predDict.keys())

            def clusterHypergeo_tfbs(biclustId, db = predDict, allGenes = pred_totalTargets[0]):
                # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
                # Take gene list and compute overlap with each miRNA
                genes = allGenes.intersection(biclusters[biclustId].getGenes())
                writeMe = []
                keys1 = db.keys()
                m1s = []
                q = []
                m = []
                n = []
                k = []
                for m1 in keys1:
                    m1s.append(m1)
                    miRNAGenes = allGenes.intersection(db[m1])
                    q.append(len(set(miRNAGenes).intersection(genes)))
                    m.append(len(miRNAGenes))
                    n.append(len(allGenes)-len(miRNAGenes))
                    k.append(len(genes))
                results = phyper(q,m,n,k)
                min_miRNA = []
                perc_targets = []
                min_pValue = float(1)
                for i in range(1,len(results)):
                    if float(results[i]) <= float(0.05)/float(674) and not q[i]==0 and float(q[i])/float(k[i]) >= 0.1:
                        if min_miRNA==[] or float(results[i]) < min_pValue:
                            min_miRNA = [i]
                            perc_targets = [float(q[i])/float(k[i])]
                            min_pValue = float(results[i])
                        elif float(results[i])==min_pValue:
                            min_miRNA.append(i)
                            perc_targets.append(float(q[i])/float(k[i]))
                print 'Bicluster #', biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA])
                return [biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA]), ' '.join([str(j) for j in perc_targets]), min_pValue]

            # Set allGenes as intersect of pred_totalTargets
            # allGenes = pred_totalTargets[0].intersection(genesInBiclusters)
            allGenes = pred_totalTargets[0]

            # Run this using all cores available
            print 'Running TFBS_DB enrichment analyses...'
            cpus = cpu_count()
            biclustIds = biclusters.keys()
            pool = Pool(processes=cpus)
            res1 = pool.map(clusterHypergeo_tfbs, biclustIds)
            pool.close()
            pool.join()

            # Dump TFBS_DB enrichment results as a pickle file
            with open('output/tfbs_db.pkl', 'wb') as pklFile:
                cPickle.dump(res1,pklFile)
        else:
            print 'Loading precached analysis...'
            with open('output/tfbs_db.pkl', 'rb') as pklFile:
                res1 = cPickle.load(pklFile)

        # Stuff into biclusters
        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, tf(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            b1.addAttribute('tfbs_db',{'tf':r1[1], 'percentTargets':r1[2], 'pValue':r1[3]})
        print 'Done.\n'

        # TFBS_DB has been run on cMonkey run
        c1.tfbs_db = True
        
    # E. 3' UTR PITA
    # If PITA enrichment hasn't been calculated for the biclusters then do so
    if not c1.pita_3pUTR:
        print 'Running PITA on Biclusters:'
        # If this has been run previously just load it up
        if not os.path.exists('output/pita_3pUTR.pkl'):
            # Get ready for multiprocessor goodness
            mgr = Manager()

            # Get a list of all genes in the biclusters
            print 'Get a list of all genes in run...'
            tmpDict = c1.getBiclusters()
            genesInBiclusters = []
            for bicluster in tmpDict:
                genes = tmpDict[bicluster].getGenes()
                for gene in genes:
                    if not gene in genesInBiclusters:
                        genesInBiclusters.append(gene)
            biclusters = mgr.dict(tmpDict)
            del tmpDict

            # Load up PITA miRNA ids
            # If this is the first time then load from the RData file
            if not os.path.exists('miRNA/pita.pkl'):
                print 'Loading PITA predictions...'
                tmpList = []
                tmpDict = {}

                with gzip.open('miRNA/pita_miRNA_sets_geneSymbol.csv.gz','r') as infile:
                    inLines = [i.strip().split(',') for i in infile.readlines() if i.strip()]

                for line in inLines:
                    if line[1] in genesInBiclusters:
                        if not line[1] in tmpList:
                            tmpList.append(line[1])
                        if not line[0] in tmpDict:
                            tmpDict[line[0]] = []
                        tmpDict[line[0]].append(line[1])

                with open('miRNA/pita.pkl','wb') as pklFile:
                    cPickle.dump(tmpDict, pklFile)
                    cPickle.dump(tmpList, pklFile)

            # Otherwise load the dumped pickle file if it exists
            else:
                print 'Loading pickled PITA predictions...'
                with open('miRNA/pita.pkl','rb') as pklFile:
                    tmpDict = cPickle.load(pklFile)
                    tmpList = cPickle.load(pklFile)

            # Setup for analysis
            predDict = mgr.dict(tmpDict)
            pred_totalTargets = mgr.list()
            pred_totalTargets.append(set(tmpList))
            del tmpDict
            del tmpList
            print 'PITA has', len(predDict.keys()),'miRNAs.'

            def clusterHypergeo_pita(biclustId, db = predDict, allGenes=pred_totalTargets[0]):
                # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
                # Take gene list and compute overlap with each miRNA
                genes = allGenes.intersection(biclusters[biclustId].getGenes())
                keys1 = db.keys()
                m1s = []
                q = []
                m = []
                n = []
                k = []
                for m1 in keys1:
                    m1s.append(m1)
                    miRNAGenes = allGenes.intersection(db[m1])
                    q.append(len(set(miRNAGenes).intersection(genes)))
                    m.append(len(miRNAGenes))
                    n.append(len(allGenes)-len(miRNAGenes))
                    k.append(len(genes))
                results = phyper(q,m,n,k)
                min_miRNA = []
                perc_targets = []
                min_pValue = float(1)
                for i in range(1,len(results)):
                    if float(results[i]) <= float(0.05)/float(674) and not q[i]==0 and float(q[i])/float(k[i]) >= 0.1:
                        if min_miRNA==[] or float(results[i]) < min_pValue:
                            min_miRNA = [i]
                            perc_targets = [float(q[i])/float(k[i])]
                            min_pValue = float(results[i])
                        elif float(results[i])==min_pValue:
                            min_miRNA.append(i)
                            perc_targets.append(float(q[i])/float(k[i]))
                print 'Bicluster #', biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA])
                return [biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA]), ' '.join([str(j) for j in perc_targets]), min_pValue]

            # Set allGenes as intersect of pita_totalTargets and genesInBiclusters
            allGenes = pred_totalTargets[0]

            # Run this using all cores available
            print 'Running PITA enrichment analyses...'
            cpus = cpu_count()
            biclustIds = biclusters.keys()
            pool = Pool(processes=cpus)
            res1 = pool.map(clusterHypergeo_pita,biclustIds)
            pool.close()
            pool.join()

            # Dump PITA enrichment results as a pickle file
            with open('output/pita_3pUTR.pkl','wb') as pklFile:
                cPickle.dump(res1, pklFile)
        else:
            print 'Loading precached analysis...'
            with open('output/pita_3pUTR.pkl','rb') as pklFile:
                res1 = cPickle.load(pklFile)

        # Stuff into biclusters
        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += miRNAInDict(m1.lower(),miRNAIDs)
            b1.addAttribute('pita_3pUTR',{'miRNA':r1[1], 'percentTargets':r1[2], 'pValue':r1[3], 'mature_seq_ids':miRNA_mature_seq_ids })
        print 'Done.\n'

        # PITA 3pUTR has been run on cMonkey run
        c1.pita_3pUTR = True

    # F. 3' UTR TargetScan
    # If TargetScan enrichment hasn't been calculated for the biclusters then do so
    if not c1.targetscan_3pUTR:
        print 'Running TargetScan on Biclusters:'
        # If this has been run previously just load it up
        if not os.path.exists('output/targetscan_3pUTR.pkl'):
            # Get ready for multiprocessor goodness
            mgr = Manager()

            # Get a list of all genes in the biclusters
            print 'Get a list of all genes in run...'
            tmpDict = c1.getBiclusters()
            genesInBiclusters = []
            for bicluster in tmpDict:
                genes = tmpDict[bicluster].getGenes()
                for gene in genes:
                    if not gene in genesInBiclusters:
                        genesInBiclusters.append(gene)
            biclusters = mgr.dict(tmpDict)
            del tmpDict

            # Load up TargetScan miRNA ids
            if not os.path.exists('miRNA/targetScan.pkl'):
                print 'Loading TargetScan predictions...'
                tmpList = []
                tmpDict = {}
                with gzip.open('miRNA/targetscan_miRNA_sets_geneSymbol.csv.gz','r') as infile:
                    inLines = [i.strip().split(',') for i in infile.readlines() if i.strip()]

                for line in inLines:
                    if line[1] in genesInBiclusters:
                        if not line[1] in tmpList:
                            tmpList.append(line[1])
                        if not line[0] in tmpDict:
                            tmpDict[line[0]] = []
                        tmpDict[line[0]].append(line[1])

                with open('miRNA/targetScan.pkl','wb') as pklFile:
                    cPickle.dump(tmpDict,pklFile)
                    cPickle.dump(tmpList,pklFile)

            # Otherwise load the dumped pickle file if it exists
            else:
                print 'Loading pickled TargetScan predictions...'
                with open('miRNA/targetScan.pkl','rb') as pklFile:
                    tmpDict = cPickle.load(pklFile)
                    tmpList = cPickle.load(pklFile)

            # Setup for analysis
            predDict = mgr.dict(tmpDict)
            pred_totalTargets = mgr.list()
            pred_totalTargets.append(set(tmpList))
            del tmpDict
            del tmpList
            print 'TargetScan has', len(predDict.keys()),'miRNAs.'

            def clusterHypergeo_ts(biclustId, db = predDict, allGenes = pred_totalTargets[0]):
                # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
                # Take gene list and compute overlap with each miRNA
                genes = allGenes.intersection(biclusters[biclustId].getGenes())
                writeMe = []
                keys1 = db.keys()
                m1s = []
                q = []
                m = []
                n = []
                k = []
                for m1 in keys1:
                    m1s.append(m1)
                    miRNAGenes = allGenes.intersection(db[m1])
                    q.append(len(set(miRNAGenes).intersection(genes)))
                    m.append(len(miRNAGenes))
                    n.append(len(allGenes)-len(miRNAGenes))
                    k.append(len(genes))
                results = phyper(q,m,n,k)
                min_miRNA = []
                perc_targets = []
                min_pValue = float(1)
                for i in range(1,len(results)):
                    if float(results[i]) <= float(0.05)/float(674) and not q[i]==0 and float(q[i])/float(k[i]) >= 0.1:
                        if min_miRNA==[] or float(results[i]) < min_pValue:
                            min_miRNA = [i]
                            perc_targets = [float(q[i])/float(k[i])]
                            min_pValue = float(results[i])
                        elif float(results[i])==min_pValue:
                            min_miRNA.append(i)
                            perc_targets.append(float(q[i])/float(k[i]))
                print 'Bicluster #', biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA])
                return [biclustId, ' '.join([m1s[miRNA] for miRNA in min_miRNA]), ' '.join([str(j) for j in perc_targets]), min_pValue]

            # Set allGenes as intersect of pita_totalTargets and genesInBiclusters
            # allGenes = pred_totalTargets[0].intersection(genesInBiclusters)
            allGenes = pred_totalTargets[0]

            # Run this using all cores available
            print 'Running TargetScan enrichment analyses...'
            cpus = cpu_count()
            biclustIds = biclusters.keys()
            pool = Pool(processes=cpus)
            res1 = pool.map(clusterHypergeo_ts,biclustIds)
            pool.close()
            pool.join()

            # Dump TargetScan enrichment results as a pickle file
            with open('output/targetscan_3pUTR.pkl','wb') as pklFile:
                cPickle.dump(res1,pklFile)
        else:
            print 'Loading precached analysis...'
            with open('output/targetscan_3pUTR.pkl','rb') as pklFile:
                res1 = cPickle.load(pklFile)

        # Stuff into biclusters
        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += miRNAInDict(m1.lower(),miRNAIDs)
            b1.addAttribute('targetscan_3pUTR',{'miRNA':r1[1], 'percentTargets':r1[2], 'pValue':r1[3], 'mature_seq_ids':miRNA_mature_seq_ids })
        print 'Done.\n'

        # TargetScan 3pUTR has been run on cMonkey run
        c1.targetscan_3pUTR = True


    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    if not os.path.exists('output/c1_all.pkl'):
        print 'Dumping Final cMonkey Object:'
        with open('output/c1_all.pkl','wb') as pklFile:
            cPickle.dump(c1, pklFile)
        print 'Done.\n'


#################################################################
## Do postProcessing on cMonkey object                         ##
#################################################################
if not os.path.exists('output/c1_postProc.pkl'):

    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    print 'Loading prechached cMonkey Object (c1_all.pkl):'
    with open('output/c1_all.pkl','rb') as pklFile:
        c1 = cPickle.load(pklFile)
    print 'Done.\n'

    # Load up the expression ratios matrix
    with open(RATIOS_PATH,'r') as ratioFile:
        conditions = [i.strip('"') for i in ratioFile.readline().strip().split('\t')]
        ratios = {}
        for line in ratioFile:
            splitUp = line.strip().split('\t')
            ratios[splitUp[0].strip('"')] = dict(zip(conditions,splitUp[1:]))

    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.genes.txt','w') as cmgFile:
        writeMe = []
        for b1 in c1.getBiclusters():
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getGenes()))
        cmgFile.write('\n'.join(writeMe))

    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.conditions.txt','w') as cmcFile:
        writeMe = []
        for b1 in c1.getBiclusters():
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getConditions()))
        cmcFile.write('\n'.join(writeMe))

    
    # Calculate bicluster eigengene (first principal components)
    if not os.path.exists('output/biclusterEigengenes.csv'):
        ret = subprocess.check_call(['./getEigengene.R',
                                     '-r', RATIOS_PATH,
                                     '-o', 'output'],
                                    stderr=subprocess.STDOUT)
        if ret == 1:
            print "could not create Eigengenes"
            exit(1)

    # Read in bicluster eigengene
    with open('output/biclusterEigengenes.csv','r') as inFile:
        biclustEigengenes = {}
        patients = [i.strip('"') for i in inFile.readline().strip().split(',')]
        patients.pop(0) # Get rid of rowname placeholder
        for line in inFile:
            eigengene = line.strip().split(',')
            bicluster = int(eigengene.pop(0).strip('"'))
            b1 = c1.getBicluster(bicluster)
            b1.addAttribute('pc1', dict(zip(patients, eigengene)))

    # Read in bicluster variance explained
    with open('output/biclusterVarianceExplained.csv','r') as inFile:
        inFile.readline() # Get rid of header
        for line in inFile:
            varExplained = line.strip().split(',')
            bicluster = int(varExplained.pop(0).strip('"'))
            b1 = c1.getBicluster(bicluster)
            b1.addAttribute('pc1.var.exp',varExplained[0])
    
    # Load the phenotype information
    # AGE,chemo_therapy,SURVIVAL,days_to_tumor_progression,SEX.bi,radiation_therapy,DEAD
    phenotypes = {}
    with open('extras/phenotypes.csv','r') as inFile:
        ids = inFile.readline().strip().split(',')[1:]
        for i in ids:
            phenotypes[i] = {}
        for line in inFile:
            splitUp = line.strip().split(',')
            phenotypes[splitUp[0]] = {}
            for i in range(len(ids)):
                phenotypes[ids[i]][splitUp[0]] = splitUp[i+1]

    def postProcess(bicluster):
        def cleanName(name):
            splitUp = name.split('.')
            return splitUp[0]+'.'+splitUp[1]+'.'+splitUp[2]

        attributes = {}
        print ' Postprocessing cluster:', bicluster
        b1 = c1.getBicluster(bicluster)
        attributes['k'] = bicluster
        # Add number of genes and conditions
        attributes['k.rows'] = b1.getNumGenes()
        attributes['k.cols'] = b1.getNumConditions()
        # Get matrix of expression for genes
        genes = b1.getGenes()
        conditions = ratios[genes[0]].keys()
        matrix = [[ratios[gene][condition] for condition in conditions] for gene in genes]
        # Get first principal component variance explained
        fpc = b1.getAttribute('pc1')
        # Corrleation with patient traits
        cleanNames = dict(zip([cleanName(i) for i in conditions],conditions))
        cond2 = set(cleanNames.keys()).intersection(phenotypes['SURVIVAL'].keys())
        pc1_1 = [fpc[cleanNames[i]] for i in cond2]
        for phenotype in ['AGE','SEX.bi','chemo_therapy','radiation_therapy']:
            p1_1 = [phenotypes[phenotype][i] for i in cond2]
            cor1 = correlation(pc1_1, p1_1)
            attributes[phenotype] = dict(zip(['rho','pValue'],cor1))
        # Association of bicluster expression with patient survival
        surv = [phenotypes['SURVIVAL'][i] for i in cond2]
        dead = [phenotypes['DEAD'][i] for i in cond2]
        age = [phenotypes['AGE'][i] for i in cond2]
        s1 = survival(surv, dead, pc1_1, age)
        attributes['Survival'] = dict(zip(['z','pValue'],s1[0]))
        attributes['Survival.AGE'] = dict(zip(['z','pValue'],s1[1]))
        
        return attributes

    if not os.path.exists('output/postProcessed.pkl'):
        # Do post processing
        print 'Do post processing...'
        cpus = cpu_count()
        biclustIds = c1.getBiclusters()
        pool = Pool(processes=cpus)
        res1 = pool.map(postProcess, biclustIds)
        pool.close()
        pool.join()
        print 'Done.\n'

        # Dump res1 into a pkl
        with open('output/postProcessed.pkl','wb') as pklFile:
            cPickle.dump(res1, pklFile)
    else:
        with open('output/postProcessed.pkl','rb') as pklFile:
            res1 = cPickle.load(pklFile)

    # Put results in cMonkey object
    for entry in res1:
        b1 = c1.getBicluster(entry['k'])
        for attribute in entry:
            if not attribute == 'k':
                b1.addAttribute(attribute, entry[attribute])

    #################################################################
    ## TomTom Upstream motifs versus Jaspar and Transfac           ##
    #################################################################
    print 'Running TOMTOM on Upstream Motifs:'
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')

    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    pssms = c1.getPssmsUpstream()
    upstreamMatches = {}

    motif_files = ['motifs/jasparCoreVertebrata_redundant.pkl',
                   'motifs/transfac_2012.1_PSSMs_vertabrate.pkl',
                   'motifs/uniprobePSSMsNonRedundant.pkl',
                   'motifs/selexPSSMsNonRedundant.pkl']
    comparison_pkl_path = 'output/upstreamJasparTransfacComparison.pkl'
    comparison_csv_path = 'output/upstreamComparison_jaspar_transfac.csv'

    if not os.path.exists(comparison_pkl_path):

        target_pssms_in = []
        for motif_file in motif_files:
            with open(motif_file, 'rb') as pklFile:
                pssms = cPickle.load(pklFile)
                for pssm in pssms.values():
                    pssm.setMethod('meme')
                target_pssms_in.append(pssms)

        # Write out results
        with open(comparison_csv_path, 'w') as outFile:
            outFile.write('Motif Name,Original E-Value,Consensus,JASPAR Motif,JASPAR Consensus,TomTom.pValue,TomTom.qValue,Probe In Bicluster,Bicluster Residual')

            # Making MEME formatted files (makeFiles function in utils)
            print 'Making files...'
            for i, target_pssms in enumerate(target_pssms_in):
                makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=pssms.values(),
                          targetPssms=target_pssms.values(), num=i)

            # Run TomTom 
            print 'Comparing Upstream motifs against databases...'
            pool = Pool(processes=cpu_count())
            res1 = pool.map(runTomTom, [i for i in range(len(target_pssms_in))])
            pool.close()
            pool.join()

            print 'Reading in Tomtom run...'
            output_lines = []
            for i in range(len(target_pssms_in)):
                with open('tmp/tomtom_out/tomtom%d.out' % i, 'r') as tomtom_outfile:
                    tomtom_outfile.readline()  # skip header
                    output_lines += [line.strip().split('\t') for line in tomtom_outfile
                                    if float(line.split('\t')[5]) <= 0.05]

            # Now iterate through output and save data
            print 'Now parsing output for %d matches...' % len(output_lines)
            for outputLine in output_lines:
                if len(outputLine) == 10 and float(outputLine[5]) <= 0.05:
                    tfName = outputLine[1].split('_')[0]
                    if not outputLine[0] in upstreamMatches:
                        upstreamMatches[outputLine[0]] = [{'factor':outputLine[1],'confidence':outputLine[5]}]
                    else:
                        upstreamMatches[outputLine[0]].append({'factor':outputLine[1],'confidence':outputLine[5]})

        with open(comparison_pkl_path,'wb') as pklFile:
            cPickle.dump(upstreamMatches, pklFile)
    else:
        print 'Loading precached upstream matches...'
        with open(comparison_pkl_path, 'rb') as pklFile:
            upstreamMatches = cPickle.load(pklFile)

    # Pump into pssms
    matched = []
    for pssmName in upstreamMatches:
        for match in upstreamMatches[pssmName]:
            pssms[pssmName].addMatch(factor=match['factor'], confidence=match['confidence'])
    print 'We matched '+str(len(upstreamMatches))+' upstream motifs.\n'
    

    #################################################################
    ## Prepare to get expanded set of TFs from TFClass             ##
    ## (http://tfclass.bioinf.med.uni-goettingen.de/)              ##
    #################################################################
    print 'Expanding TF factor list with TFClass families...'
    # Read in humanTFs_All.csv with <TF Name>,<Entrez ID>
    tfName2entrezId = {}
    with open('TF/humanTFs_All.csv','r') as inFile:
        inFile.readline() # Get rid of header
        for line in inFile:
            splitUp = line.strip().split(',')
            if splitUp[2] in entrez2id:
                for i in entrez2id[splitUp[2]]:
                    tfName2entrezId[splitUp[0]] = i

    # Read in tfFamilies.csv for expanded list of possible TFs
    tfFamilies = {}
    with open('TF/tfFamilies.csv','r') as inFile:
        inFile.readline() # Get rid of header
        for line in inFile:
            splitUp = line.strip().split(',')
            tmp = []
            for i in splitUp[2].split(' '):
                if i in entrez2id:
                    tmp += entrez2id[i]
            tfFamilies[splitUp[0]] = tmp

    # Add expanded TF regulators
    pssms = c1.getPssmsUpstream()
    for pssm in pssms:
        expandedFactors = {}
        matches = pssms[pssm].getMatches()
        # Collapse matches to a set of entrez IDs and add expanded factors
        if not matches==None:
            for match in matches:
                if match['factor'] in tfName2entrezId:
                    factor = tfName2entrezId[match['factor']]
                    if not factor in expandedFactors:
                        expandedFactors[factor] = [factor]
                        for family in tfFamilies:
                            if factor in tfFamilies[family]:
                                expandedFactors[factor] += tfFamilies[family]
            # Push expanded TF factor list into the PSSM object
            for factor in expandedFactors:
                #print factor,expandedFactors[factor]
                for expandedFactor in list(set(expandedFactors[factor])):
                    pssms[pssm].addExpandedMatch(expandedFactor, factor)
    print 'Finished expanding TF factor list.\n'


    #######################################################################
    ## Filter expanded TFs through correlation with bicluster eigengenes ##
    #######################################################################
    print 'Correlate TFs with eigengenes...'
    
    # Get microarray data
    expData = {}
    with open(ALL_RATIOS_PATH,'r') as inFile:
        names = inFile.readline().strip().split('\t')
        names = [i.strip('"') for i in names]
        for line in inFile:
            splitUp = line.strip().split('\t')
            geneId = splitUp.pop(0).strip('"')
            expData[geneId] = dict(zip(names,splitUp))

    allNames = names

    # [rho, pValue] = correlation(a1,a2)
    biclusters = c1.getBiclusters()
    for b1 in biclusters:
        for pssm in biclusters[b1].getPssmsUpstream():
            factors = pssm.getExpandedMatches()
            compared = {}
            for subset in subsets:
                compared[subset] = []
            # Get maximum amount of correlation positive or negative
            if not factors==None:
                print b1, pssm.getName(), len(factors)
                for factor in factors:
                    for subset in subsets:
                        corMax = []
                        if factor['factor'] in expData.keys():
                            eigengene = biclusters[b1].getAttribute('pc1')
                            if not factor['factor'] in compared[subset]:
                                    compared[subset].append(factor['factor'])
                                    cor1 = correlation([eigengene[i] for i in allNames][subsetsPos[subset][0]:subsetsPos[subset][1]],[expData[factor['factor']][i] for i in allNames][subsetsPos[subset][0]:subsetsPos[subset][1]])
                                    print subset, factor['factor'], cor1
                                    if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                                        corMax = cor1
                        if not corMax==[]:
                            pssm.addCorrelatedMatch(subset,factor['factor'],corMax[0],corMax[1])
    print 'Done.\n'
    
    
    #################################################################
    ## Expand and correlated additional TFs for TFBS_DB            ##
    #################################################################
    print 'Expand and correlate TFBS_DB TFs...'
    # For each bicluster
    for bicluster in c1.getBiclusters():
        b1 = c1.getBicluster(bicluster)
        # Get the tfbs_db attribute and for each TF get the list of expanded factors
        tfs = b1.getAttribute('tfbs_db')
        expandedFactors = {}
        if not tfs==None:
            for tf in tfs['tf'].split(' '):
                if tf[0:2]=='V_':
                    tf = 'V$'+tf[2:]
                # Get the list of expanded factors
                if tf in tfName2entrezId:
                    factor = tfName2entrezId[tf]
                    if not factor in expandedFactors:
                        expandedFactors[factor] = [factor]
                        for family in tfFamilies:
                            if factor in tfFamilies[family]:
                                expandedFactors[factor] += tfFamilies[family]
                        expandedFactors[factor] = list(set(expandedFactors[factor]))
        
        # Push expanded TF factor list into the PSSM object
        if not expandedFactors=={}:
            print factor,expandedFactors
            b1.addAttribute('tfbs_db_expanded',expandedFactors)

    # [rho, pValue] = correlation(a1,a2)
    for bicluster in c1.getBiclusters():
        b1 = c1.getBicluster(bicluster)
        factors = b1.getAttribute('tfbs_db_expanded')
        compared = {}
        for subset in subsets:
            compared[subset] = []
        correlatedFactor = {}
        for subset in subsets:
            correlatedFactor[subset] = []
        # Get maximum amount of correlation positive or negative
        if not factors==None:
            for factor1 in factors:
                for factor2 in factors[factor1]:
                    for subset in subsets:
                        corMax = []
                        if factor2 in expData:
                            eigengene = b1.getAttribute('pc1')
                            if not factor2 in compared:
                                compared[subset].append(factor2)
                                cor1 = correlation([eigengene[i] for i in allNames][subsetsPos[subset][0]:subsetsPos[subset][1]],[expData[factor2][i] for i in allNames][subsetsPos[subset][0]:subsetsPos[subset][1]])
                                print cor1
                                if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                                    corMax = cor1
                            if not corMax==[]:
                                correlatedFactor[subset].append({'factor':factor2,'rho':corMax[0],'pValue':corMax[1]})
        b1.addAttribute('tfbs_db_correlated',correlatedFactor)
    print 'Done.\n'

    #################################################################
    ## Dump first principal components for each bicluster          ##
    #################################################################
    print 'Write biclusterFirstPrincComponents.csv...'
    # Get all first principal components for each bicluster
    fpcWrite = []
    conditions = c1.getBicluster(1).getAttribute('pc1').keys()
    for i in sorted(c1.getBiclusters().keys()):
        pc1 = c1.getBicluster(i).getAttribute('pc1')
        fpcWrite.append(str(i)+','+','.join([str(pc1[j]) for j in conditions]))

    # Dump out file
    with open('output/'+str(fpcFile),'w') as outfile:
        outfile.write('Bicluster,'+','.join([j.strip() for j in conditions])+'\n')
        outfile.write('\n'.join(fpcWrite))

    print 'Done.\n'

    #################################################################
    ## Get permuted p-values for upstream meme motifs              ##
    #################################################################
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')
    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    
    # Compare the random motifs to the original motif in TOMTOM
    permPValues = {}
    matched = 0
    pssms = c1.getPssmsUpstream(de_novo_method='meme')
    if not os.path.exists('output/upstreamMotifPermutedPValues.csv'):
        outFile = open('output/upstreamMotifPermutedPValues.csv','w')
        outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
        pssmsNames = pssms.keys()
        print 'Loading precached random PSSMs...'
        randPssmsDict = {}
        for i in [5,10,15,20,25,30,35,40,45,50,55,60,65]:
            stdout.write(str(i)+' ')
            stdout.flush()
            pklFile = open(str(randPssmsDir)+'/pssms_upstream_'+str(i)+'.pkl','rb')
            randPssmsDict[i] = cPickle.load(pklFile)
            r1 = [randPssmsDict[i][pssm1].setMethod('meme') for pssm1 in randPssmsDict[i]]
            delMes = []
            for randPssm in randPssmsDict[i]:
                if not float(randPssmsDict[i][randPssm].getEValue()) <= float(maxEValue):
                    delMes.append(randPssm)
            for j in delMes:
                del randPssmsDict[i][j]

        print '\nMaking files...'
        for i in range(len(pssms)):
            clustSize = randPSSMClustSize((c1.getBicluster(int(pssmsNames[i].split('_')[0]))).getNumGenes())
            makeFiles(nucFreqs=c1.getNucFreqsUpstream(), queryPssms=[pssms[pssmsNames[i]]],targetPssms=randPssmsDict[clustSize].values(),num=i)

        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.'
        print 'Running TOMTOM to compare PSSMs...'
        pool = Pool(processes=cpus)
        pool.map(runTomTom,range(len(pssms)))
        pool.close()
        pool.join()

        print 'Reading in Tomtom run...'
        for run in range(len(pssms)):
            tomtomPValues = {}
            outputFile = open('tmp/tomtom_out/tomtom'+str(run)+'.out','r')
            output = outputFile.readlines()
            outputFile.close()
            # Now iterate through output and save data
            output.pop(0) # Get rid of header
            while len(output)>0:
                outputLine = output.pop(0).strip().split('\t')
                if len(outputLine)==10:
                    tomtomPValues[outputLine[1]] = float(outputLine[3])
            pValues = tomtomPValues.values()
            similar = 0
            for pValue in pValues:
                if float(pValue) <= float(0.05):
                    similar += 1
            # Write out the results
            mot = outputLine[0].split('_')[1]
            permPValues[outputLine[0]] = { mot+'.consensus':str(pssms[outputLine[0]].getConsensusMotif()), mot+'.permutedEV<=10':str(len(pValues)), mot+'.similar':str(similar), mot+'.permPV':str(float(similar)/float(1000)) }
            pssms[outputLine[0]].setPermutedPValue(str(float(similar)/float(1000)))
            outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].getEValue())+','+str(pssms[outputLine[0]].getConsensusMotif())+','+str(len(pValues))+','+str(similar)+','+str(1000)+','+str(float(similar)/float(1000)))
        outFile.close()
    else:
        print 'Using precalculated upstream p-values...'
        inFile = open('output/upstreamMotifPermutedPValues.csv','r')
        inFile.readline()
        upPValues = [i.strip().split(',') for i in inFile.readlines()]
        inFile.close()
        for line in upPValues:
            pssms[line[0]].setPermutedPValue(str(float(line[7])/float(1000)))
    print 'Done.\n'


    #################################################################
    ## Compare 3' UTR Weeder Motifs to miRBase using miRvestigator ##
    #################################################################
    print 'Running miRvestigator on 3\' UTR Motifs:'
    if not os.path.exists('output/m2m.pkl'):
        print 'Computing miRNA matches...'
        biclusters = c1.getBiclusters()
        pssms = c1.getPssms3pUTR()
        seqs3pUTR = c1.getSeqs3pUTR().values()
        m2m = miRvestigator(pssms.values(),seqs3pUTR,seedModel=[6,7,8],minor=True,p5=True,p3=True,wobble=False,wobbleCut=0.25,baseDir='output',species='mmu')
        with open('output/m2m.pkl','wb') as pklFile:
            cPickle.dump(m2m,pklFile)
    else:
        print 'Loading precached miRNA matches...'
        with open('output/m2m.pkl','rb') as pklFile:
            m2m = cPickle.load(pklFile)
    print 'Done.\n'


    #################################################################
    ## Convert miRNAs and Get permuted p-values for 3' UTR motifs  ##
    #################################################################
    pssms = c1.getPssms3pUTR()
    print 'Loading miRvestigator results...'
    # Convert miRvestigator results
    if not os.path.exists('output/miRvestigatorResults.pkl'):
        inFile = open('output/miRNA/scores.csv','r')
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]
        miRNA_matches = {}
        for line in lines:
            if not line[1]=='NA':
                miRNA_mature_seq_ids = []
                for i in line[1].split('_'):
                    miRNA_mature_seq_ids += miRNAInDict(i.lower(),miRNAIDs)
                miRNA_matches[line[0]] = {'miRNA':line[1],'model':line[2],'mature_seq_ids':miRNA_mature_seq_ids}
                for m1 in miRNA_mature_seq_ids:
                    pssms[line[0]].addMatch(factor=m1, confidence=line[2])
        with open('output/miRvestigatorResults.pkl','wb') as pklFile:
            cPickle.dump(miRNA_matches, pklFile)
    else:
        with open('output/miRvestigatorResults.pkl','rb') as pklFile:
            miRNA_matches = cPickle.load(pklFile)
            for m1 in miRNA_matches:
                for m2 in miRNA_matches[m1]['mature_seq_ids']:
                    #print m1, m2, miRNA_matches
                    pssms[m1].addMatch(factor=m2, confidence=miRNA_matches[m1]['model'])

    # Compile results to put them into the postProcessed
    print 'Get perumted p-values for 3\' UTR motifs...'
    pklFile = open('randPSSMs/weederRand.pkl','rb')
    weederRand = cPickle.load(pklFile)
    pklFile.close()
    pklFile = open('randPSSMs/weederRand_all.pkl','rb')
    weederRand_all = cPickle.load(pklFile)
    pklFile.close()
    clustSizes = sorted(weederRand['8bp'].keys())
    for pssm1 in pssms:
        seqNum = pssms[pssm1].getNumGenes()
        consensus = pssms[pssm1].getConsensusMotif()
        width = len(consensus)
        splitUp = pssm1.split('_')
        clustInd = 5
        for c in clustSizes:
            if seqNum > c:
                clustInd = c
        pValue = float(sum(1 for i in weederRand[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].getEValue())))/float(len(weederRand[str(width)+'bp'][clustInd]))
        pValue_all = float(sum(1 for i in weederRand_all[str(width)+'bp'][clustInd] if float(i) >= float(pssms[pssm1].getEValue())))/float(len(weederRand_all[str(width)+'bp'][clustInd]))
        pssms[pssm1].setPermutedPValue({'pValue':pValue,'pValue_all':pValue_all})
    print 'Done.\n'


    #################################################################
    ## Run replication p-values                                    ##
    #################################################################
    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.genes.txt','w') as cmgFile:
        writeMe = []
        for b1 in c1.getBiclusters():
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getGenes()))
        cmgFile.write('\n'.join(writeMe))

    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.conditions.txt','w') as cmcFile:
        writeMe = []
        for b1 in c1.getBiclusters():
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).getConditions()))
        cmcFile.write('\n'.join(writeMe))

    def runReplication(repScript):
        print '  Replication running for '+repScript+'...'
        repProc = Popen('cd replication_'+repScript+'; R --no-save < replicationDatasetPermutation.R', shell=True, stdout=PIPE, stderr=PIPE)
        output = repProc.communicate()[0].split('\n')

    # Run replication on all datasets
    runEm = []
    for i in ['French','REMBRANDT','GSE7696']:
        if not os.path.exists('output/replicationPvalues_'+i+'.csv'):
            runEm.append(i)
    if len(runEm)>0:
        print 'Run replication..'
        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.'
        pool = Pool(processes=cpus)
        pool.map(runReplication,runEm)
        pool.close()
        pool.join()

    #################################################################
    ## Read in replication p-values                                ##
    #################################################################
    # Read in replication p-values - French Dataset      
    # '','n.rows','overlap.rows','new.resid.norm.gbm','avg.norm.perm.resid.gbm','norm.perm.p.gbm','new.resid.norm.all','avg.norm.perm.resid.all','norm.perm.p.all','pc1.var.exp.gbm','avg.pc1.var.exp.gbm','pc1.perm.p.gbm','pc1.var.exp.all','avg.pc1.var.exp.all','pc1.perm.p.all','survival.gbm','survival.p.gbm','survival.age.gbm','survival.age.p.gbm','survival.all','survival.p.all','survival.age.all','survival.age.p.all'
    print 'Loading replication p-values...'
    inFile = open('output/replicationPvalues_French.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_French',value={'French_new.resid.norm':splitUp[3], 'French_avg.resid.norm':splitUp[4], 'French_norm.perm.p':splitUp[5], 'French_pc1.var.exp':splitUp[9], 'French_avg.pc1.var.exp':splitUp[10], 'French_pc1.perm.p':splitUp[11], 'French_survival':splitUp[15], 'French_survival.p':splitUp[16], 'French_survival.age':splitUp[17], 'French_survival.age.p':splitUp[18]})
        b1.addAttribute(key='replication_French_all',value={'French_all_new.resid.norm':splitUp[6], 'French_all_avg.resid.norm':splitUp[7], 'French_all_norm.perm.p':splitUp[8], 'French_all_pc1.var.exp':splitUp[12], 'French_all_avg.pc1.var.exp':splitUp[13], 'French_all_pc1.perm.p':splitUp[14], 'French_all_survival':splitUp[19], 'French_all_survival.p':splitUp[20], 'French_all_survival.age':splitUp[21], 'French_all_survival.age.p':splitUp[22]})
    inFile.close()
    # Read in replication p-values - REMBRANDT Dataset      
    # "","n.rows","orig.resid","orig.resid.norm","overlap.rows","new.resid","avg.perm.resid","perm.p","new.resid.norm","avg.norm.perm.resid","norm.perm.p","survival","survival.p","survival.age","survival.age.p"
    inFile = open('output/replicationPvalues_REMBRANDT.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_REMBRANDT',value={'REMBRANDT_new.resid.norm':splitUp[3], 'REMBRANDT_avg.resid.norm':splitUp[4], 'REMBRANDT_norm.perm.p':splitUp[5], 'REMBRANDT_pc1.var.exp':splitUp[9], 'REMBRANDT_avg.pc1.var.exp':splitUp[10], 'REMBRANDT_pc1.perm.p':splitUp[11], 'REMBRANDT_survival':splitUp[15], 'REMBRANDT_survival.p':splitUp[16], 'REMBRANDT_survival.age':splitUp[17], 'REMBRANDT_survival.age.p':splitUp[18]})
        b1.addAttribute(key='replication_REMBRANDT_all',value={'REMBRANDT_all_new.resid.norm':splitUp[6], 'REMBRANDT_all_avg.resid.norm':splitUp[7], 'REMBRANDT_all_norm.perm.p':splitUp[8], 'REMBRANDT_all_pc1.var.exp':splitUp[12], 'REMBRANDT_all_avg.pc1.var.exp':splitUp[13], 'REMBRANDT_all_pc1.perm.p':splitUp[14], 'REMBRANDT_all_survival':splitUp[19], 'REMBRANDT_all_survival.p':splitUp[20], 'REMBRANDT_all_survival.age':splitUp[21], 'REMBRANDT_all_survival.age.p':splitUp[22]})
    # Read in replication p-values - GSE7696 Dataset
    # '', 'n.rows','overlap.rows','new.resid.norm.gbm','avg.norm.perm.resid.gbm','norm.perm.p.gbm','pc1.var.exp.gbm','avg.pc1.var.exp.gbm','pc1.perm.p.gbm','survival.gbm','survival.p.gbm','survival.age.gbm','survival.age.p.gbm'
    inFile = open('output/replicationPvalues_GSE7696.csv','r')
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        b1 = c1.getBicluster(int(splitUp[0].replace('"','')))
        b1.addAttribute(key='replication_GSE7696',value={'GSE7696_new.resid.norm':splitUp[3], 'GSE7696_avg.resid.norm':splitUp[4], 'GSE7696_norm.perm.p':splitUp[5], 'GSE7696_pc1.var.exp':splitUp[6], 'GSE7696_avg.pc1.var.exp':splitUp[7], 'GSE7696_pc1.perm.p':splitUp[8], 'GSE7696_survival':splitUp[9], 'GSE7696_survival.p':splitUp[10], 'GSE7696_survival.age':splitUp[11], 'GSE7696_survival.age.p':splitUp[12]})
    inFile.close()
    print 'Done.\n'


    ###########################################################################
    ## Run permuted p-value for variance epxlained first principal component ##
    ###########################################################################
    if not os.path.exists('output/residualPermutedPvalues_permAll.csv'):
        print 'Calculating FPC permuted p-values...'
        enrichProc = Popen("R --no-save < permutedResidualPvalues_permAll_mc.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'

    #################################################################
    ## Read in residual permutations to use for filtering          ##
    #################################################################
    print 'Load residual permuted p-values...'
    with open('output/residualPermutedPvalues_permAll.csv','r') as inFile:
        # "","bicluster","n.rows","n.cols","orig.resid","avg.perm.resid","perm.p","orig.resid.norm","avg.norm.perm.resid","norm.perm.p","pc1.var.exp","avg.pc1.var.exp","pc1.perm.p"
        inFile.readline()
        inLines = inFile.readlines()
        for line in inLines:
            splitUp = line.strip().split(',')
            b1 = c1.getBicluster(int(splitUp[0].strip('"')))
            b1.addAttribute(key='resid.norm.perm.p',value=str(splitUp[9]))
            b1.addAttribute(key='pc1.perm.p',value=str(splitUp[12]))
    print 'Done.\n'


    #################################################################
    ## Run functional enrichment and GO term similarity            ##
    #################################################################
    if not os.path.exists('output/biclusterEnrichment_GOBP.csv'):
        print 'Run functional enrichment...'
        enrichProc = Popen("cd funcEnrichment; R --no-save < enrichment.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'
    if not os.path.exists('output/jiangConrath_hallmarks.csv'):
        print 'Run semantic similarity...'
        enrichProc = Popen("cd funcEnrichment; R --no-save < goSimHallmarksOfCancer.R", shell=True, stdout=PIPE, stderr=PIPE)
        output = enrichProc.communicate()[0]
        print 'Done.\n'

    #################################################################
    ## Read in functional enrichment                               ##
    #################################################################
    print 'Load GO Biological Process functional enrichment...'
    with open('output/biclusterEnrichment_GOBP.csv','r') as inFile:
        inFile.readline() # Get rid of header
        inLines = inFile.readlines()
        lines = [line.strip().split(',') for line in inLines]

    for line in lines:
        b1 = c1.getBicluster(int(line[0].strip('"')))
        b1.addAttribute(key='goTermBP',value=line[2].strip('"').split(';'))
    print 'Done.\n'
    
    #################################################################
    ## Read in hallmarks of cacner                                 ##
    #################################################################
    print 'Load Jiang-Conrath semantic similarity to Hallmarks of Cancer...'
    with open('output/jiangConrath_hallmarks.csv','r') as inFile:
        hallmarks = [i for i in inFile.readline().split(',') if not i.strip('"')=='']
        inLines = inFile.readlines()
        lines = [line.strip().split(',') for line in inLines]
        for line in lines:
            b1 = c1.getBicluster(int(line[0].strip('"')))
            b1.addAttribute(key='hallmarksOfCancer',value=dict(zip(hallmarks,line[1:])))
        print 'Done.\n'

    #################################################################
    ## Save out the final cMonkey object so we don't lose progress ##
    #################################################################
    print 'Dumping Final cMonkey Object:'
    with open('output/c1_postProc.pkl','wb') as pklFile:
        cPickle.dump(c1, pklFile)
    print 'Done.\n'
else:
    print 'Loading from precached cMonkey Object:'
    with open('output/c1_postProc.pkl','rb') as pklFile:
        c1 = cPickle.load(pklFile)
    print 'Done.\n'


#################################################################
## Run NEO and integrate some form of results                  ##
## TODO: Make NEO run on subsets                               ##
#################################################################
if not os.path.exists('output/causality'):
    ## Run the runNEO.R script and do the causality analyses
    print '  Network edge orienting (NEO)...'
    neoProc = Popen('cd NEO; R --no-save < runNEO.R', shell=True, stdout=PIPE, stderr=PIPE)
    output = neoProc.communicate()[0].split('\n')

## Pull together analysis into cohesive output
causalSummary = []
# For each mutation
for dir1 in os.listdir('output/causality'):
    # For each regulator
    if dir1[0:7]=='causal_':
        # For each 
        for file1 in os.listdir('output/causality/'+dir1):
            if file1[0:3]=='sm.':
                with open('output/causality/'+dir1+'/'+file1,'r') as inFile:
                    inLine = inFile.readline() # Get rid of header
                    while 1:
                        inLine = inFile.readline()
                        if not inLine:
                            break
                        splitUp = inLine.strip().split(',')
                        if float(splitUp[6]) >= leo_nb_AtoB and float(splitUp[12]) <= mlogp_M_AtoB:
                            # Somatic Mutation(1), Regulator(3), Biclster(5), leo.nb.AtoB(6), mlogp.M.AtoB(12), PathAB(17), SEPathAB(18), ZPathAB(19), PPathAB(20), BLV.AtoB(25), RMSEA.AtoB(28)
                            causalSummary.append({'Mutation': splitUp[1].strip('"').lstrip('M:'), 'Regulator': splitUp[3].strip('"').lstrip('A:'), 'Bicluster': splitUp[5].strip('"').lstrip('B:bic_'), 'leo.nb.AtoB': splitUp[6], 'mlogp.M.AtoB': splitUp[12], 'PathAB': splitUp[17], 'SEPathAB': splitUp[18], 'ZPathAB': splitUp[19], 'PPathAB': splitUp[20], 'BLV.AtoB': splitUp[25], 'RMSEA.AtoB': splitUp[28]})

## Output:  Somatic Mutation(1), Regulator(3), Biclster(5), leo.nb.AtoB(6), mlogp.M.AtoB(12), PathAB(17), SEPathAB(18), ZPathAB(19), PPathAB(20), BLV.AtoB(25), RMSEA.AtoB(28)
header = ['Mutation', 'Regulator', 'Bicluster', 'leo.nb.AtoB', 'mlogp.M.AtoB', 'PathAB', 'SEPathAB', 'ZPathAB', 'PPathAB', 'BLV.AtoB', 'RMSEA.AtoB']
with open('output/causalitySummary.csv','w') as outFile:
    outFile.write(','.join(header)+'\n')
    outFile.write('\n'.join([','.join([i[j] for j in header]) for i in causalSummary]))

## Dump out correspondent regulators (both mechanistically and causally predicted)
correspondentRegulators = {}
for causalFlow in causalSummary:
    b1 = c1.getBicluster(int(causalFlow['Bicluster']))
    ## Upstream (TFs)
    tfs = []
    # 1. MEME and WEEDER Upstream motifs
    for pssm in b1.getPssmsUpstream():
        for subset in subsets:
            matches = pssm.getCorrelatedMatches(subset)
            if matches:
                for corTf in matches:
                    if corTf['pValue']<=pVCut and abs(corTf['rho'])>=rhoCut:
                        tfs.append(corTf['factor'])
        # 2. TFBS_DB
        for corTf in b1.getAttributes()['tfbs_db_correlated'][subset]:
            if corTf['pValue']<=pVCut and abs(corTf['rho'])>=rhoCut:
                tfs.append(corTf['factor'])
    # 3. Find Correspondent TF regulators
    if causalFlow['Regulator'] in tfs:
        if not int(causalFlow['Bicluster']) in correspondentRegulators:
            correspondentRegulators[int(causalFlow['Bicluster'])] = {'tf':[],'miRNA':[]}
        correspondentRegulators[int(causalFlow['Bicluster'])]['tf'].append(causalFlow['Regulator'])

    ## 3' UTR (miRNA)
    miRNAs = []
    # 1. WEEDER 3'UTR
    for pssm in b1.getPssms3pUTR():
        matches = pssm.getMatches()
        if matches:
            for miR in matches:
                if miR['confidence'] in ['8mer','7mer_a1','7mer_m8']:
                    miRNAs += miRNAIDs_rev[miR['factor']]
    # 2. PITA (not correlated)
    if float(b1.getAttributes()['pita_3pUTR']['pValue'])<=pVCut and float(b1.getAttributes()['pita_3pUTR']['percentTargets'].split(' ')[0])>=percTargets:
        miRNAs += b1.getAttributes()['pita_3pUTR']['miRNA'].split(' ')
    # 3. TargetScan (not correlated)
    if float(b1.getAttributes()['targetscan_3pUTR']['pValue'])<=pVCut and float(b1.getAttributes()['targetscan_3pUTR']['percentTargets'].split(' ')[0])>=percTargets:
        miRNAs += b1.getAttributes()['targetscan_3pUTR']['miRNA'].split(' ')
    # 4. Find Correspondent miRNA regulators
    for miR in miRNAs:
        if  compareMiRNANames(causalFlow['Regulator'].lower(), miR.lower()):
            if not int(causalFlow['Bicluster']) in correspondentRegulators:
                correspondentRegulators[int(causalFlow['Bicluster'])] = {'tf':[],'miRNA':[]}
            correspondentRegulators[int(causalFlow['Bicluster'])]['miRNA'].append(causalFlow['Regulator'])

## Put correspondent regulators into cMonkeyWrapper object
for biclust in correspondentRegulators.keys():
    b1 = c1.getBicluster(biclust)
    b1.addAttribute(key='correspondentRegulators',value=correspondentRegulators[biclust])


#################################################################
## Write out the final post-processed file                     ##
#################################################################
print 'Write postProcessedVFinal.csv...'
postOut = []
hallmarksOfCancer = c1.getBicluster(1).getAttribute('hallmarksOfCancer').keys()
for i in sorted(c1.getBiclusters().keys()):
    writeMe = []
    b1 = c1.getBicluster(i)
    # Write file line by line
    #   a. Bicluster basics:  id, genes, conditions, resid, resid.norm, resid.norm.perm.p
    writeMe += [str(i), # id
                str(b1.getAttribute('k.rows')), # genes
                str(b1.getAttribute('k.cols')), # conditions
                str(b1.getNormResidual()), # normalized residual
                str(b1.getAttribute('resid.norm.perm.p')), # normalized residual permuted p-value
                str(b1.getAttribute('pc1.var.exp')), # Variance explained by first principal component
                str(b1.getAttribute('pc1.perm.p'))] # Variance explained by first principal component
    #   b. Upstream motifs:  meme.motif1.E, meme.motif1.consensus, meme.motif1.matches, meme.motif1.permPV
    motifNames = b1.getPssmsNamesUpstream()
    upstreamMotifs = { 'meme_motif1':None, 'meme_motif2':None, 'weeder_motif1':None, 'weeder_motif2':None }
    for m1 in motifNames:
        splitUp = m1.split('_')
        if splitUp[1]=='motif1' and splitUp[2]=='meme':
            upstreamMotifs['meme_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='meme':
            upstreamMotifs['meme_motif2'] = m1
        if splitUp[1]=='motif1' and splitUp[2]=='weeder':
            upstreamMotifs['weeder_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='weeder':
            upstreamMotifs['weeder_motif2'] = m1
    #   - MEME motifs
    for meme1 in ['meme_motif1','meme_motif2']:
        if not upstreamMotifs[meme1]==None:
            pssm1 = b1.getPssmUpstream(upstreamMotifs[meme1])
            original = []
            matches = 'NA'
            expandedMatches = 'NA'
            correlatedMatches = {}
            originalExpanded = {}
            minCorrelated = {}
            for subset in subsets:
                correlatedMatches[subset] = 'NA'
                originalExpanded[subset] = 'NA'
                minCorrelated[subset] = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([match1['factor'] for match1 in pssm1.getMatches()])
                tmp = pssm1.getExpandedMatches()
                if not tmp==None:
                    expandedMatches = {}
                    for i in tmp:
                        if not i['seedFactor'] in original:
                            original.append(i['seedFactor'])
                        if not i['seedFactor'] in expandedMatches:
                            expandedMatches[i['seedFactor']] = []
                        expandedMatches[i['seedFactor']].append(i['factor'])
                    expandedMatches = ' '.join([seedFactor+':'+';'.join(expandedMatches[seedFactor]) for seedFactor in expandedMatches])
                for subset in subsets:
                    tmp = pssm1.getCorrelatedMatches(subset)
                    if not tmp==None:
                        for match1 in tmp:
                            if match1['pValue']<=pVCut and abs(match1['rho'])>=rhoCut:
                                if correlatedMatches[subset]=='NA':
                                    correlatedMatches[subset] = []
                                correlatedMatches[subset].append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                                if minCorrelated[subset]=='NA' or match1['pValue']<minCorrelated[subset]['pValue']:
                                    minCorrelated[subset] = match1
                        if not correlatedMatches[subset]=='NA':
                            correlatedMatches[subset] = ' '.join(correlatedMatches[subset])
                        if not minCorrelated[subset]=='NA':
                            if minCorrelated[subset]['factor'] in original:
                                originalExpanded[subset] = 'Original'
                            else:
                                originalExpanded[subset] = 'Expanded'
                            minCorrelated[subset] = minCorrelated[subset]['factor']+':'+str(minCorrelated[subset]['rho'])+':'+str(minCorrelated[subset]['pValue'])
            writeMe += ([str(pssm1.getEValue()), # E-value
                        str(pssm1.getPermutedPValue()), # Permuted p-value for motif
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches, # Matches to the motif from TransFac and Jaspar
                        expandedMatches] # Expanded matches usign TFClass
                        + [correlatedMatches[subset]+','+originalExpanded[subset]+','+minCorrelated[subset] for subset in subsets])
        else:
            writeMe += (['NA', # E-value
                        'NA', # Permuted p-value for motif
                        'NA', # Motif consensus sequence
                        'NA', # Matches to the motif from TransFac and Jaspar
                        'NA'] # Expanded matches using TFClass
                        + ['NA,NA,NA' for subset in subsets])
    #   - WEEDER motifs
    for weeder1 in ['weeder_motif1','weeder_motif2']:
        if not upstreamMotifs[weeder1]==None:
            pssm1 = b1.getPssmUpstream(upstreamMotifs[weeder1])
            original = []
            matches = 'NA'
            expandedMatches = 'NA'
            correlatedMatches = {}
            originalExpanded = {}
            minCorrelated = {}
            for subset in subsets:
                correlatedMatches[subset] = 'NA'
                originalExpanded[subset] = 'NA'
                minCorrelated[subset] = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([match1['factor'] for match1 in pssm1.getMatches()])
                tmp = pssm1.getExpandedMatches()
                if not tmp==None:
                    expandedMatches = {}
                    for i in tmp:
                        if not i['seedFactor'] in original:
                            original.append(i['seedFactor'])
                        if not i['seedFactor'] in expandedMatches:
                            expandedMatches[i['seedFactor']] = []
                        expandedMatches[i['seedFactor']].append(i['factor'])
                    expandedMatches = ' '.join([seedFactor+':'+';'.join(expandedMatches[seedFactor]) for seedFactor in expandedMatches])
                for subset in subsets:
                    tmp = pssm1.getCorrelatedMatches(subset)
                    if not tmp==None:
                        for match1 in tmp:
                            if match1['pValue']<=pVCut and abs(match1['rho'])>=rhoCut:
                                if correlatedMatches[subset]=='NA':
                                    correlatedMatches[subset] = []
                                correlatedMatches[subset].append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                                if minCorrelated[subset]=='NA' or match1['pValue']<minCorrelated[subset]['pValue']:
                                    minCorrelated[subset] = match1
                        if not correlatedMatches[subset]=='NA':
                            correlatedMatches[subset] = ' '.join(correlatedMatches[subset])
                        if not minCorrelated[subset]=='NA':
                            if minCorrelated[subset]['factor'] in original:
                                originalExpanded[subset] = 'Original'
                            else:
                                originalExpanded[subset] = 'Expanded'
                            minCorrelated[subset] = minCorrelated[subset]['factor']+':'+str(minCorrelated[subset]['rho'])+':'+str(minCorrelated[subset]['pValue'])
            writeMe += ([str(pssm1.getEValue()), # E-value
                        #str(pssm1.getPermutedPValue()), # Permuted p-value for motif
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches, # Matches to the motif from TransFac and Jaspar
                        expandedMatches] # Expanded matches usign TFClass
                        + [correlatedMatches[subset]+','+originalExpanded[subset]+','+minCorrelated[subset] for subset in subsets])
        else:
            writeMe += (['NA', # E-value
                        #'NA', # Permuted p-value for motif
                        'NA', # Motif consensus sequence
                        'NA', # Matches to the motif from TransFac and Jaspar
                        'NA'] # Expanded matches using TFClass
                        + ['NA,NA,NA' for subset in subsets])
    #   c. Enriched TFs:  TFBS_DB.TFs,TFBS_DB.percTargets,TFBS_DB.pValue
    for association in ['tfbs_db']:
        a1 = b1.getAttribute(association)
        if not a1['tf']=='':
            expandedMatches = 'NA'
            correlatedMatches = {}
            originalExpanded = {}
            minCorrelated = {}
            for subset in subsets:
                correlatedMatches[subset] = 'NA'
                originalExpanded[subset] = 'NA'
                minCorrelated[subset] = 'NA'
            tmp = b1.getAttribute('tfbs_db_expanded')
            if not tmp==None:
                expandedMatches = ' '.join([seedFactor+':'+';'.join(tmp[seedFactor]) for seedFactor in tmp])
            tmp = b1.getAttribute('tfbs_db_correlated')
            for subset in subsets:
                if not tmp[subset]==None:
                    for match1 in tmp[subset]:
                        if match1['pValue']<=pVCut and abs(match1['rho'])>=rhoCut:
                            if correlatedMatches[subset]=='NA':
                                correlatedMatches[subset] = []
                            correlatedMatches[subset].append(match1['factor']+':'+str(match1['rho'])+':'+str(match1['pValue']))
                            if minCorrelated[subset]=='NA' or match1['pValue']<minCorrelated[subset]['pValue']:
                                minCorrelated[subset] = match1
                    if not correlatedMatches[subset]=='NA':
                        correlatedMatches[subset] = ' '.join(correlatedMatches[subset])
                    if not minCorrelated[subset]=='NA':
                        if minCorrelated[subset]['factor'] in original:
                            originalExpanded[subset] = 'Original'
                        else:
                            originalExpanded[subset] = 'Expanded'
                        minCorrelated[subset] = minCorrelated[subset]['factor']+':'+str(minCorrelated[subset]['rho'])+':'+str(minCorrelated[subset]['pValue'])
                if len(correlatedMatches[subset])==0:
                    correlatedMatches[subset] = 'NA'
            writeMe += [str(a1['tf']).replace(';',' '),
                        str(a1['percentTargets']).replace(';',' '),
                        str(a1['pValue']),
                        expandedMatches] + [correlatedMatches[subset] + ',' + originalExpanded[subset] + ',' + minCorrelated[subset] for subset in subsets]
        else:
            writeMe += (['NA','NA','NA','NA'] + ['NA,NA,NA' for subset in subsets])

    #   d. 3' UTR motifs:  weeder.motif1.E, weeder.motif1.permPV, weeder.motif1.permPV_all, weeder.motif1.consensus, weeder.motif1.matches, weeder.motif1.model
    motifNames = b1.getPssmsNames3pUTR()
    p3utrMotifs = { 'weeder_motif1':None, 'weeder_motif2':None }
    for m1 in motifNames:
        splitUp = m1.split('_')
        if splitUp[1]=='motif1' and splitUp[2]=='weeder':
            p3utrMotifs['weeder_motif1'] = m1
        if splitUp[1]=='motif2' and splitUp[2]=='weeder':
            p3utrMotifs['weeder_motif2'] = m1
    #   - WEEDER motifs
    for weeder1 in ['weeder_motif1','weeder_motif2']:
        if not p3utrMotifs[weeder1]==None:
            pssm1 = b1.getPssm3pUTR(p3utrMotifs[weeder1])
            matches = 'NA'
            model = 'NA'
            if not pssm1.getMatches()==None:
                matches = ' '.join([miRNAIDs_rev[match1['factor']] for match1 in pssm1.getMatches()])
                model = pssm1.getMatches()[0]['confidence']
            permutedPValue = 'NA'
            if not pssm1.getPermutedPValue()==None:
                permutedPValue = pssm1.getPermutedPValue()
            writeMe += [str(pssm1.getEValue()), # E-value
                        str(permutedPValue['pValue']), # Permuted p-value for motif
                        str(permutedPValue['pValue_all']), # Permuted p-value for motif (all)
                        pssm1.getConsensusMotif(), # Motif consensus sequence
                        matches, # Matches to the motif to miRBase
                        model] # Model fit by miRvestigator
        else:
            writeMe += ['NA', # E-value
                        'NA', # Permuted p-value for motif
                        'NA', # Permuted p-value for motif (all)
                        'NA', # Motif consensus sequence
                        'NA', # Matches to the motif to miRBase
                        'NA'] # Model fit by miRvestigator

    #   e. Enriched miRNAs:  3pUTR_pita.miRNAs,3pUTR_pita.percTargets,3pUTR_pita.pValue,3pUTR_targetScan.miRNAs,3pUTR_targetScan.percTargets,3pUTR_targetScan.pValue
    for association in ['pita_3pUTR', 'targetscan_3pUTR']:
        a1 = b1.getAttribute(association)
        if not a1['miRNA']=='':
            writeMe += [str(a1['miRNA']).replace(';',' '), str(a1['percentTargets']).replace(';',' '), str(a1['pValue'])]
        else:
            writeMe += ['NA','NA','NA']

    #   f. Correpsondent regulators
    corrRegs = b1.getAttribute('correspondentRegulators')
    if not corrRegs==None:
        # Split up by TF and miRNA
        if not len(corrRegs['tf'])==0:
            writeMe += [' '.join(sorted(set(corrRegs['tf'])))]
        else:
            writeMe += ['NA']
        if not len(corrRegs['miRNA'])==0:
            writeMe += [' '.join(sorted(set(corrRegs['miRNA'])))]
        else:
            writeMe += ['NA']
    else:
        writeMe += ['NA','NA']

    #   g. Associations with traits:  age, sex.bi, chemo_therapy, radiation_therapy
    for association in ['AGE','SEX.bi', 'chemo_therapy','radiation_therapy']:
        ass1 = b1.getAttribute(association)
        writeMe += [str(ass1['rho']), str(ass1['pValue'])]
    surv1 = b1.getAttribute('Survival') 
    survAge1 = b1.getAttribute('Survival.AGE')
    writeMe += [str(surv1['z']), str(surv1['pValue']), str(survAge1['z']), str(survAge1['pValue'])]

    #   h. Replications:  'REMBRANDT_new.resid.norm','REMBRANDT_avg.resid.norm','REMBRANDT_norm.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p'
    replications_French = b1.getAttribute('replication_French')
    replications_REMBRANDT = b1.getAttribute('replication_REMBRANDT')
    replications_French_all = b1.getAttribute('replication_French_all')
    replications_REMBRANDT_all = b1.getAttribute('replication_REMBRANDT_all')
    replications_GSE7696 = b1.getAttribute('replication_GSE7696')
    for replication in ['French_pc1.var.exp','French_avg.pc1.var.exp','French_pc1.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p']:
        writeMe.append(str(replications_French[replication]))
    for replication in ['French_all_pc1.var.exp','French_all_avg.pc1.var.exp','French_all_pc1.perm.p','French_all_survival','French_all_survival.p','French_all_survival.age','French_all_survival.age.p']:
        writeMe.append(str(replications_French_all[replication]))
    for replication in ['REMBRANDT_pc1.var.exp','REMBRANDT_avg.pc1.var.exp','REMBRANDT_pc1.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p']:
        writeMe.append(str(replications_REMBRANDT[replication]))
    for replication in ['REMBRANDT_all_pc1.var.exp','REMBRANDT_all_avg.pc1.var.exp','REMBRANDT_all_pc1.perm.p','REMBRANDT_all_survival','REMBRANDT_all_survival.p','REMBRANDT_all_survival.age','REMBRANDT_all_survival.age.p']:
        writeMe.append(str(replications_REMBRANDT_all[replication]))
    for replication in ['GSE7696_pc1.var.exp','GSE7696_avg.pc1.var.exp','GSE7696_pc1.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p']:
        writeMe.append(str(replications_GSE7696[replication]))

    #   i. Functional enrichment of biclusters using GO term Biological Processes
    bfe1 = b1.getAttribute('goTermBP')
    if bfe1==['']:
        writeMe.append('NA')
    else:
        writeMe.append(';'.join(bfe1))
    
    #   j. Hallmarks of Cancer:  Hanahan and Weinberg, 2011
    bhc1 = b1.getAttribute('hallmarksOfCancer')
    for hallmark in hallmarksOfCancer:
        writeMe.append(str(bhc1[hallmark]))

    #   k. Glioma sub-type enrichment: 'NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM'
    #for overlap in ['NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']:
    #    writeMe.append(str(b1.getAttribute(overlap)))
    
    # Add to the final output file
    postOut.append(deepcopy(writeMe))

with open('output/'+str(postProcessedFile),'w') as postFinal:
    header = ['Bicluster', 'Genes', 'Patients', 'Norm. Residual', 'Norm. Residual Perm. P-Value', 'Var. Exp. First PC', 'Var. Exp. First PC Perm. P-Value'] + \
     ['MEME Motif1 E-Value', 'MEME Motif1 Perm. P-Value', 'Up.MEME Motif1 Consensus', 'Up.MEME Motif1 Matches','Up.MEME Motif1 Expanded Matches', 'Up.MEME Motif1 Correlated Matches', 'Up.MEME Motif1 Original/Expanded', 'Up.MEME Motif1 Minimum Correlated'] + \
     ['Up.MEME Motif2 E-Value', 'Up.MEME Motif2 Perm. P-Value', 'Up.MEME Motif2 Consensus', 'Up.MEME Motif2 Matches', 'Up.MEME Motif2 Expanded Matches', 'Up.MEME Motif2 Correlated Matches', 'Up.MEME Motif2 Original/Expanded', 'Up.MEME Motif2 Minimum Correlated'] + \
     ['Up.WEEDER Motif1 Score', 'Up.WEEDER Motif1 Consensus', 'Up.WEEDER Motif1 Matches', 'Up.WEEDER Motif1 Expanded Matches', 'Up.WEEDER Motif1 Correlated Matches', 'Up.WEEDER Motif1 Original/Expanded', 'Up.WEEDER Motif1 Minimum Correlated'] + \
     ['Up.WEEDER Motif2 Score', 'Up.WEEDER Motif2 Consensus', 'Up.WEEDER Motif2 Matches', 'Up.WEEDER Motif2 Expanded Matches', 'Up.WEEDER Motif2 Correlated Matches', 'Up.WEEDER Motif2 Original/Expanded', 'Up.WEEDER Motif2 Minimum Correlated'] + \
     ['TFBS_DB.TFs', 'TFBS_DB.percTargets', 'TFBS_DB.pValue', 'TFBS_DB.Exapnded Matches', 'TFBS_DB.Correlated Matches', 'TFBS_DB.Original/Expanded','TFBS_DB.Minimum Correlated'] + \
     ['3pUTR.WEEDER Motif1 E-Value', '3pUTR.WEEDER Motif1 Perm. P-Value', '3pUTR.WEEDER Motif1 Perm. P-Value (All)', '3pUTR.WEEDER Motif1 Consensus', '3pUTR.WEEDER Motif1 Matches', '3pUTR.WEEDER Motif1 Model'] + \
     ['3pUTR.WEEDER Motif2 E-Value', '3pUTR.WEEDER Motif2 Perm. P-Value', '3pUTR.WEEDER Motif2 Perm. P-Value (All)', '3pUTR.WEEDER Motif2 Consensus', '3pUTR.WEEDER Motif2 Matches', '3pUTR.WEEDER Motif2 Model'] + \
     ['3pUTR_pita.miRNAs', '3pUTR_pita.percTargets', '3pUTR_pita.pValue'] + \
     ['3pUTR_targetScan.miRNAs', '3pUTR_targetScan.percTargets', '3pUTR_targetScan.pValue'] + \
     ['Correspondent.TFs', 'Correspondent.miRNAs'] + \
     ['Age', 'Age.p', 'Sex', 'Sex.p', 'Chemotherapy', 'Chemotherapy.p', 'RadiationTherapy', 'RadiationTherapy.p', 'Survival', 'Survival.p', 'Survial.covAge', 'Survival.covAge.p'] + \
     ['French_pc1.var.exp','French_avg.pc1.var.exp','French_pc1.perm.p','French_survival','French_survival.p','French_survival.age','French_survival.age.p', 'French_all_pc1.var.exp','French_all_avg.pc1.var.exp','French_all_pc1.perm.p','French_all_survival','French_all_survival.p','French_all_survival.age','French_all_survival.age.p'] + \
     ['REMBRANDT_pc1.var.exp','REMBRANDT_avg.pc1.var.exp','REMBRANDT_pc1.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p', 'REMBRANDT_all_pc1.var.exp','REMBRANDT_all_avg.pc1.var.exp','REMBRANDT_all_pc1.perm.p','REMBRANDT_all_survival','REMBRANDT_all_survival.p','REMBRANDT_all_survival.age','REMBRANDT_all_survival.age.p'] + \
     ['GSE7696_pc1.var.exp','GSE7696_avg.pc1.var.exp','GSE7696_pc1.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p'] + \
     ['GO_Term_BP'] + \
     [i.strip() for i in hallmarksOfCancer]
    postFinal.write(','.join(header)+'\n'+'\n'.join([','.join(i) for i in postOut]))

print 'Done.\n'

