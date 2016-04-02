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
from rpy2.robjects import FloatVector, IntVector, StrVector
from rpy2 import rinterface

# Custom offYerBack libraries
from cMonkeyWrapper import cMonkeyWrapper
from pssm import pssm
from miRvestigator import miRvestigator
import utils
from sys import stdout, exit
import gzip


#################################################################
## rpy2 integration                                            ##
#################################################################

def make_rint_vector(a):
    return IntVector(map(lambda i: rinterface.NA_Integer if i == 'NA' else int(i), a))

def make_rfloat_vector(a):
    return FloatVector(map(lambda f: rinterface.NA_Real if f == 'NA' else float(f), a))


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

# Detenrmine the working directories, we pass this to the
# external R scripts
CURRENT_DIR = os.getcwd()
# note that we need to paste in the '/' at the end because the R scripts
# currently assume it
BASE_DIR = os.path.dirname(os.getcwd()) + '/'
OUTPUT_DIR = os.path.join(CURRENT_DIR, 'output')

SYNONYM_PATH = '../synonymThesaurus.csv.gz'
MIRNA_FASTA_PATH = 'miRNA/hsa.mature.fa'  # Need to update this to the latest database
RATIOS_PATH = '../gbmTCGA_exprMat_medianFiltered.tsv'
ALL_RATIOS_PATH = 'expression/gbmTCGA_exprMat.tsv'
CMONKEY2_RUNDB = '../out/cmonkey_run.db'
REPLICATION_DATASET_NAMES = ['French','REMBRANDT','GSE7696']


g_weeder_args    = None
g_weeder_results = None

#################################################################
## Functions                                                   ##
#################################################################

meme_args = None
cluster_meme_runs = None

def meme(num, seqfile, bgfile, nmotifs, min_motif_width, max_motif_width, revcomp, seed=None):
    """ Run meme and get the output into PSSMs"""
    global cluster_meme_runs

    # Arguments for meme
    args = '%s -bfile %s -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw %d -maxw %d -nmotifs %d' % (seqfile, bgfile, min_motif_width, max_motif_width, nmotifs)

    if revcomp:
        args += ' -revcomp'

    if not seed is None:
        args += ' -cons ' + str(seed)

    print "MEME args: '%s'" % args
    meme_proc = Popen("meme %s" % args, shell=True, stdout=PIPE)
    output = meme_proc.communicate()[0].split('\n')

    PSSMs = []

    for i in range(len(output)):
        desc_comps = output[i].strip().split(' ')

        if desc_comps[0] == 'Motif' and desc_comps[2] == 'position-specific' and desc_comps[3] == 'probability':
            i += 2  # Skip the separator line, go to the summary line
            summary_comps = output[i].strip().split(' ')
            width = int(summary_comps[5])
            sites = int(summary_comps[7])
            evalue = float(summary_comps[9])
            matrix = []

            for j in range(width):
                i += 1
                matrix += [[float(let) for let in output[i].strip().split(' ') if let]]
            PSSMs.append(pssm('%s_motif%s_meme' % ((seqfile.split('_')[1]).split('.')[0], desc_comps[1]),
                              sites, evalue, matrix, [], 'meme'))
    cluster_meme_runs[num] = PSSMs


def run_meme(runarg):
    """Wrapper function to run the meme function using a multiprocessing pool
    """
    global meme_args

    run_num, filepath = runarg
    meme(run_num, filepath, meme_args['bgfile'], meme_args['nmotifs'], meme_args['min_motif_width'],
         meme_args['max_motif_width'], meme_args['revcomp'])


def weeder(bicluster, seqfile, bgfile, size, enriched, revcomp):
    """
    Run weeder and parse its output
    First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
    """
    global g_weeder_results
    print "run weeder on '%s'" % seqfile

    # First run weederTFBS
    weeder_args = "%s %s %s %s" % (seqfile, bgfile, size, enriched)
    if revcomp:
        weeder_args += ' S'
    errout = open('tmp/weeder/stderr.out','w')
    weeder_proc = Popen("weederlauncher %s" % weeder_args, shell=True, stdout=PIPE, stderr=errout)
    output = weeder_proc.communicate()

    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqfile)+'.wee','r')
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
    motif_id = 1
    bicluster_id = int(seqfile.split('_')[1].split('.')[0])
    while 1:
        if len(outLines)<=1:
            break

        if revcomp:
            name = outLines.pop(0).strip() # Get match
            name += '_'+outLines.pop(0).strip()
        else:
            name = outLines.pop(0).strip() # Get match

        if not name.find('(not highest-ranking)') == -1:
            break

        # Get redundant motifs
        outLines.pop(0)
        red_motifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
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
            nums = [float(i.strip()) for i in col.strip().split('\t')[1].split(' ') if i]
            colsum = sum(nums)
            matrix += [[ nums[0] / colsum, nums[1] / colsum, nums[2] / colsum, nums[3] / colsum]]
            col = outLines.pop(0)
        PSSMs += [pssm('%d_motif%d_weeder' % (bicluster_id, motif_id),
                       len(instances), hitBp[len(matrix)][1], matrix, red_motifs, 'weeder')]
        motif_id += 1
    g_weeder_results[bicluster] = PSSMs


def run_weeder(run_arg):
    global g_weeder_args

    run_num, filepath = run_arg
    weeder(run_num, filepath, g_weeder_args['bgfile'], g_weeder_args['size'],
           g_weeder_args['enriched'], g_weeder_args['revcomp'])


def tomtom(num, dist_meth='ed', q_thresh=1, min_overlap=6):
    args = '-dist %s -o tmp/tomtom_out -text -thresh %d -min-overlap %d -verbosity 1 tmp/query%d.tomtom tmp/target%d.tomtom' % (dist_meth, q_thresh, min_overlap, num, num)
    print args

    with open('tmp/stderr_%d.out' % num,'w') as errout:
        tomtom_proc = Popen("tomtom %s" % args, shell=True, stdout=PIPE, stderr=errout)
        with open('tmp/tomtom_out/tomtom%d.out' % num, 'w') as outfile:
            output = tomtom_proc.communicate()[0]
            outfile.write(output)


def run_tomtom(i):
    tomtom(i, 'ed', 1, 6)


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
    cor_test = robj.r['cor.test']
    result = cor_test(make_rfloat_vector(a1), make_rfloat_vector(a2))
    res = dict(zip(result.names, list(result)))
    return res['estimate'][0], res['p.value'][0]


# Compute survival p-value from R
def survival(survival, dead, pc1, age):
    """
    Calculate the survival correlation coefficient and p-value between two variables.
    Input: Four arrays of float or integers.
    Returns:
    """
    surv = robj.r("""
library('survival')
surv <- function(s, dead, pc1, age) {
  scph1 = summary(coxph(Surv(s,dead == 'DEAD') ~ pc1))
  scph2 = summary(coxph(Surv(s,dead == 'DEAD') ~ pc1 + age))
  c(scph1$coef[1,4], scph1$coef[1,5],  scph2$coef[1,4], scph2$coef[1,5])
}
""")
    res = surv(make_rint_vector(survival), StrVector(dead), make_rfloat_vector(pc1), make_rint_vector(age))
    return [[res[0], res[1]], [res[2], res[3]]]


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

def compute_upstream_motifs_meme(c1):
    """
    A. Upstream motifs (MEME) ##
    If MEME hasn't been run on the biclusters upstream sequences then do so
    """
    global meme_args, cluster_meme_runs

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
            run_args = []

            # First make fasta files for all biclusters
            print 'Making Files for MEME Upstream run...'
            for b1 in c1.biclusters:
                seqs = c1.getBiclusterSeqsUpstream(b1)
                if len(seqs) > 0:
                    cluster_filename = 'tmp/meme/fasta/bicluster_%d.fasta' % b1
                    run_args.append((b1, cluster_filename))
                    with open(cluster_filename, 'w') as outfile:
                        outfile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

            # Where all the results will be stored
            cluster_meme_runs = mgr.dict()

            # Parameters to use for running MEME
            meme_args = mgr.dict( { 'bgfile': bgFile, 'nmotifs': nMotifs,
                                    'min_motif_width': motifWidth['upstream'][0],
                                    'max_motif_width': motifWidth['upstream'][1],
                                    'revcomp': revComp['upstream'] } )

            print 'Running MEME on Upstream sequences...'
            print 'There are %d CPUs available.' % cpu_count()
            pool = Pool(processes=cpu_count())
            pool.map(run_meme, run_args)
            pool.close()
            pool.join()

            # Dump weeder results as a pickle file
            with open('output/meme_upstream.pkl','wb') as pklFile:
                cPickle.dump(deepcopy(cluster_meme_runs), pklFile)
        else:
            print 'Loading from precached object...'
            with open('output/meme_upstream.pkl','rb') as pklFile:
                cluster_meme_runs = cPickle.load(pklFile)

        # Add PSSMs to cMonkey object
        print 'Storing output...'
        for i, pssms in cluster_meme_runs.items():
            for pssm1 in pssms:
                b1 = c1.getBicluster(i)
                pssm1.setMethod('meme')
                b1.add_pssm_upstream(pssm1)

        print 'Done with MEMEing.\n'

        # MEME upstream has been run on cMonkey run
        c1.meme_upstream = True


def __compute_motifs_weeder(resultpath, biclusters, add_result, bicluster_seqs, args_dict):
    """Generic motif detection function with weeder
    """
    global g_weeder_args, g_weeder_results

    if not os.path.exists(resultpath):
        if os.path.exists('tmp'):
            rmtree('tmp')
        if not os.path.exists('tmp/weeder/fasta'):
            os.makedirs('tmp/weeder/fasta')

        # Run MEME on all biclusters
        mgr = Manager()
        run_args = []

        # First make fasta files for all biclusters
        for b1 in biclusters:
            seqs = bicluster_seqs(b1)
            if len(seqs) > 0:
                cluster_filename = 'tmp/weeder/fasta/bicluster_%d.fasta' % b1
                run_args.append((b1, cluster_filename))                    
                with open(cluster_filename, 'w') as outfile:
                    outfile.write('\n'.join(['>'+gene+'\n'+seqs[gene] for gene in seqs]))

        # Where all the results will be stored
        g_weeder_results = mgr.dict()

        # Parameters to use for running Weeder
        # Set to run Weeder on 'medium' setting which means 6bp, 8bp and 10bp motifs
        g_weeder_args = mgr.dict(args_dict)

        print 'Running Weeder...'
        print 'There are %d CPUs available.' % cpu_count()
        pool = Pool(processes=cpu_count())
        pool.map(run_weeder, run_args)
        pool.close()
        pool.join()

        # Dump weeder results as a pickle file
        with open(resultpath,'wb') as outfile:
            cPickle.dump(deepcopy(g_weeder_results), outfile)

    else:
        print 'Loading from precached object...'
        with open(resultpath,'rb') as infile:
            g_weeder_results = cPickle.load(infile)

    print 'Storing output...'
    for i, pssms in g_weeder_results.items():
        for p in pssms:
            p.setMethod('weeder')
            add_result(i, p)


def compute_upstream_motifs_weeder(c1):
    if not c1.weeder_upstream:
        __compute_motifs_weeder('output/weeder_upstream.pkl',
                                c1.biclusters,
                                lambda bi, p: c1.getBicluster(bi).add_pssm_upstream(p),
                                lambda bi: c1.getBiclusterSeqsUpstream(bi),
                                { 'bgfile': 'HS', 'size': 'small', 'enriched': 'T50', 'revcomp': True })
        c1.weeder_upstream = True

def compute_3pUTR_weeder(c1):
    if not c1.weeder_3pUTR:
        __compute_motifs_weeder('output/weeder_3pUTR.pkl',
                                c1.biclusters,
                                lambda bi, p: c1.getBicluster(bi).add_pssm_3putr(p),
                                lambda bi: c1.getBiclusterSeqs3pUTR(bi),
                                { 'bgfile': 'HS3P', 'size': 'small', 'enriched': 'T50', 'revcomp': False })
        c1.weeder_3pUTR = True


predDict = None
pred_totalTargets = None
biclusters = None

def cluster_hypergeo(bicluster_id):
    """concurrent computation
    k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
    Take gene list and compute overlap with each miRNA
    """
    global predDict, pred_totalTargets, biclusters
    db = predDict
    all_genes = pred_totalTargets[0]

    genes = all_genes.intersection(biclusters[bicluster_id].genes)
    m1s = []
    q = []
    m = []
    n = []
    k = []
    for m1, m1_genes in db.items():
        m1s.append(m1)
        miRNAGenes = all_genes.intersection(m1_genes)
        q.append(len(set(miRNAGenes).intersection(genes)))
        m.append(len(miRNAGenes))
        n.append(len(all_genes) - len(miRNAGenes))
        k.append(len(genes))

    results = phyper(q,m,n,k)
    min_mirna = []
    perc_targets = []
    min_pvalue = 1.0

    for i in range(1, len(results)):
        if float(results[i]) <= 0.05 / 674.0 and q[i] != 0 and float(q[i]) / float(k[i]) >= 0.1:
            if min_mirna == [] or float(results[i]) < min_pvalue:
                min_mirna = [i]
                perc_targets = [ float(q[i]) / float(k[i]) ]
                min_pvalue = float(results[i])
            elif float(results[i]) == min_pvalue:
                min_mirna.append(i)
                perc_targets.append(float(q[i]) / float(k[i]))

    print 'Bicluster #%d' % bicluster_id, ' '.join([m1s[miRNA] for miRNA in min_mirna])
    return [bicluster_id, ' '.join([m1s[i] for i in min_mirna]),
            ' '.join(map(str, perc_targets)), min_pvalue]


def __bicluster_genes(c1):
    result = set()
    for key, bicluster in c1.biclusters.items():
        for gene in bicluster.genes:
            result.add(gene)
    return result


def __read_predictions(pred_path, pkl_path, genesInBiclusters, manager):
    if not os.path.exists(pkl_path):
        print 'loading predictions...'
        tmp_set = set()
        tmp_dict = {}
        with gzip.open(pred_path, 'r') as infile:
            inLines = [i.strip().split(',') for i in infile.readlines() if i.strip()]

        for line in inLines:
            if line[1] in genesInBiclusters:
                tmp_set.add(line[1])
                if not line[0] in tmp_dict:
                    tmp_dict[line[0]] = []
                tmp_dict[line[0]].append(line[1])

        with open(pkl_path,'wb') as pklFile:
            cPickle.dump(tmp_dict, pklFile)
            cPickle.dump(tmp_set, pklFile)

    # Otherwise load the dumped pickle file if it exists
    else:
        print 'Loading pickled predictions...'
        with open(pkl_path,'rb') as pklFile:
            tmp_dict = cPickle.load(pklFile)
            tmp_set = cPickle.load(pklFile)
    
    predDict = manager.dict(tmp_dict)
    pred_totalTargets = manager.list()
    pred_totalTargets.append(set(tmp_set))
    return predDict, pred_totalTargets


def __compute_enrichment(c1, name,  pkl_path, pred_path, pred_pkl_path):
    """General enrichment analysis function"""
    global predDict, pred_totalTargets, biclusters

    print 'Running %s enrichment on Biclusters:' % name

    if not os.path.exists(pkl_path):
        print 'Get a list of all genes in run...'
        genesInBiclusters = __bicluster_genes(c1)
        mgr = Manager()
        biclusters = mgr.dict(c1.biclusters)
        predDict, pred_totalTargets = __read_predictions(pred_path, pred_pkl_path,
                                                         genesInBiclusters, mgr)
        print '%s prediction has %d TFs.' % (name, len(predDict))
        allGenes = pred_totalTargets[0]

        print 'Running %s enrichment analyses...' % name
        pool = Pool(processes=cpu_count())
        result = pool.map(cluster_hypergeo, c1.biclusters.keys())
        pool.close()
        pool.join()

        with open(pkl_path, 'wb') as pklFile:
            cPickle.dump(result, pklFile)
    else:
        print 'Loading precached analysis...'
        with open(pkl_path, 'rb') as pklFile:
            result = cPickle.load(pklFile)
    return result

def compute_tfbsdb_enrichment(c1):
    """D. Upstream TFBS DB enrichment Analysis"""
    if not c1.tfbs_db:
        res1 = __compute_enrichment(c1, 'TFBS_DB', 'output/tfbs_db.pkl',
                                    'TF/tfbsDb_5000_gs.csv.gz',
                                    'TF/tfbs_db.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, tf(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            b1.add_attribute('tfbs_db', {'tf':r1[1], 'percentTargets': r1[2],
                                         'pValue':r1[3]})
        print 'Done.\n'
        c1.tfbs_db = True


def compute_3pUTR_pita_set_enrichment(c1):
    """E. 3' UTR PITA"""
    if not c1.pita_3pUTR:
        res1 = __compute_enrichment(c1, 'PITA', 'output/pita_3pUTR.pkl',
                                    'miRNA/pita_miRNA_sets_geneSymbol.csv.gz',
                                    'miRNA/pita.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += utils.miRNAInDict(m1.lower(), miRNAIDs)
            b1.add_attribute('pita_3pUTR', {'miRNA': r1[1],
                                            'percentTargets': r1[2],
                                            'pValue': r1[3],
                                            'mature_seq_ids': miRNA_mature_seq_ids})
        print 'Done.\n'
        c1.pita_3pUTR = True


def compute_3pUTR_targetscan_set_enrichment(c1):
    """F. 3' UTR TargetScan"""
    if not c1.targetscan_3pUTR:
        res1 = __compute_enrichment(c1, 'TargetScan', 'output/targetscan_3pUTR.pkl',
                                    'miRNA/targetscan_miRNA_sets_geneSymbol.csv.gz',
                                    'miRNA/targetScan.pkl')

        print 'Storing results...'
        for r1 in res1:
            # r1 = [biclusterId, miRNA(s), Percent Targets, P-Value]
            b1 = c1.getBicluster(r1[0])
            miRNA_mature_seq_ids = []
            for m1 in r1[1]:
                miRNA_mature_seq_ids += utils.miRNAInDict(m1.lower(), miRNAIDs)
            b1.add_attribute('targetscan_3pUTR',
                             {'miRNA': r1[1], 'percentTargets': r1[2],
                              'pValue':r1[3], 'mature_seq_ids': miRNA_mature_seq_ids})
        print 'Done.\n'
        c1.targetscan_3pUTR = True


CM_PKL_PATH = 'output/c1_all.pkl'

def compute_additional_info():
    if not os.path.exists(CM_PKL_PATH):
        c1 = read_cmonkey_run(CMONKEY2_RUNDB)

        ################################################################
        ## Fill in the missing parts                                  ##
        ################################################################
        #  Check to see if all parts are there:                        #
        #   A. Upstream motifs (MEME)                                  #
        #   B. Upstream motif (Weeder)                                 #
        #   C. 3' UTR Weeder-miRvestigator (Weeder)                    #
        #   D. TFBS DB Enrichment                                      #
        #   E. 3' UTR PITA (Set Enrichment)                            #
        #   F. 3' UTR TargetScan (Set Enrichment)                      #
        ################################################################

        compute_upstream_motifs_meme(c1)
        compute_upstream_motifs_weeder(c1)
        compute_3pUTR_weeder(c1)
        compute_tfbsdb_enrichment(c1)
        compute_3pUTR_pita_set_enrichment(c1)
        compute_3pUTR_targetscan_set_enrichment(c1)

        with open(CM_PKL_PATH, 'wb') as outfile:
            cPickle.dump(c1, outfile)
    else:
        print 'Loading prechached cMonkey Object (c1_all.pkl):'
        with open(CM_PKL_PATH, 'rb') as infile:
            c1 = cPickle.load(infile)
    return c1


g_ratios = None
g_phenotypes = None

def post_process(bicluster):
    global g_ratios
    def clean_name(name):
        comps = name.split('.')
        return '%s.%s.%s' % (comps[0], comps[1], comps[2])

    attributes = {}
    print ' Postprocessing cluster:', bicluster
    b1 = c1.getBicluster(bicluster)
    attributes['k'] = bicluster

    # Add number of genes and conditions
    attributes['k.rows'] = b1.num_genes()
    attributes['k.cols'] = b1.num_conditions()

    # Get matrix of expression for genes
    genes = b1.genes
    conditions = g_ratios[genes[0]].keys()
    matrix = [[g_ratios[gene][condition] for condition in conditions] for gene in genes]

    # Get first principal component variance explained
    fpc = b1.attributes['pc1']

    # Corrleation with patient traits
    cleanNames = dict(zip([clean_name(i) for i in conditions], conditions))
    cond2 = set(cleanNames.keys()).intersection(g_phenotypes['SURVIVAL'].keys())
    pc1_1 = [fpc[cleanNames[i]] for i in cond2]

    for phenotype in ['AGE', 'SEX.bi', 'chemo_therapy', 'radiation_therapy']:
        p1_1 = [g_phenotypes[phenotype][i] for i in cond2]
        cor1 = correlation(pc1_1, p1_1)
        attributes[phenotype] = dict(zip(['rho', 'pValue'], cor1))

    # Association of bicluster expression with patient survival
    surv = [g_phenotypes['SURVIVAL'][i] for i in cond2]
    dead = [g_phenotypes['DEAD'][i] for i in cond2]
    age = [g_phenotypes['AGE'][i] for i in cond2]
    s1 = survival(surv, dead, pc1_1, age)
    attributes['Survival'] = dict(zip(['z', 'pValue'], s1[0]))
    attributes['Survival.AGE'] = dict(zip(['z', 'pValue'], s1[1]))

    return attributes

CLUSTER_GENES_PATH = 'output/cluster.members.genes.txt'
CLUSTER_CONDS_PATH = 'output/cluster.members.conditions.txt'
CLUSTER_EIGENGENES_PATH = 'output/biclusterEigengenes.csv'
CLUSTER_VARIANCE_EXPLAINED_PATH = 'output/biclusterVarianceExplained.csv'
PHENOTYPES_PATH = 'extras/phenotypes.csv'
POSTPROCESS_PKL_PATH = 'output/postProcessed.pkl'


def __read_ratios(path, c1):
    print "reading ratios matrix"
    with open(path, 'r') as infile:
        conditions = [i.strip('"') for i in infile.readline().strip().split('\t')]
        ratios = {}
        for line in infile:
            comps = line.strip().split('\t')
            ratios[comps[0].strip('"')] = dict(zip(conditions, comps[1:]))

    print "dump cluster row members"
    with open(CLUSTER_GENES_PATH, 'w') as outfile:
        for b1 in c1.biclusters:
            outfile.write('%s %s\n' % (b1, ' '.join(c1.getBicluster(b1).genes))) 

    print "dump cluster condition members"
    with open(CLUSTER_CONDS_PATH, 'w') as outfile:
        for b1 in c1.biclusters:
            outfile.write('%s %s\n' % (b1, ' '.join(c1.getBicluster(b1).conditions)))
    return ratios


def __get_cluster_eigengenes(ratios_path, c1):
    # Calculate bicluster eigengene (first principal components)
    print "compute bicluster eigengenes"
    if not os.path.exists(CLUSTER_EIGENGENES_PATH):
        ret = subprocess.check_call(['./getEigengene.R',
                                     '-r', ratios_path,
                                     '-o', 'output'],
                                    stderr=subprocess.STDOUT)
        if ret == 1:
            print "could not create Eigengenes"
            exit(1)

    # Read in bicluster eigengene
    print "read bicluster eigengenes"
    with open(CLUSTER_EIGENGENES_PATH, 'r') as infile:
        patients = [i.strip('"') for i in infile.readline().strip().split(',')]
        patients.pop(0) # Get rid of rowname placeholder

        for line in infile:
            eigengene = line.strip().split(',')
            bicluster = int(eigengene.pop(0).strip('"'))
            b1 = c1.getBicluster(bicluster)
            b1.add_attribute('pc1', dict(zip(patients, eigengene)))


def __get_cluster_variance_explained(variance_explained_path, c1):
    print "read bicluster variance explained"
    with open(CLUSTER_VARIANCE_EXPLAINED_PATH, 'r') as infile:
        infile.readline() # Get rid of header
        for line in infile:
            varExplained = line.strip().split(',')
            bicluster = int(varExplained.pop(0).strip('"'))
            b1 = c1.getBicluster(bicluster)
            b1.add_attribute('pc1.var.exp', varExplained[0])


def __get_phenotype_info(phenotypes_path, c1):
    # AGE,chemo_therapy,SURVIVAL,days_to_tumor_progression,SEX.bi,radiation_therapy,DEAD
    print "read phenotype information"
    phenotypes = {}
    with open(phenotypes_path, 'r') as infile:
        ids = infile.readline().strip().split(',')[1:]
        for i in ids:
            phenotypes[i] = {}

        for line in infile:
            splitUp = line.strip().split(',')
            phenotypes[splitUp[0]] = {}
            for i in range(len(ids)):
                phenotypes[ids[i]][splitUp[0]] = splitUp[i+1]
    return phenotypes


def __do_postprocess(postprocess_pkl_path, c1, ratios, phenotypes):
    global g_ratios, g_phenotypes

    if not os.path.exists(postprocess_pkl_path):
        g_ratios = ratios
        g_phenotypes = phenotypes
        print 'Do post processing...'
        pool = Pool(processes=cpu_count())
        res1 = pool.map(post_process, c1.biclusters)
        pool.close()
        pool.join()
        print 'Done.\n'

        # Dump res1 into a pkl
        with open(postprocess_pkl_path, 'wb') as outfile:
            cPickle.dump(res1, outfile)
    else:
        with open(postprocess_pkl_path, 'rb') as infile:
            res1 = cPickle.load(infile)

    # Put results in cMonkey object
    for entry in res1:
        b1 = c1.getBicluster(entry['k'])
        for attribute in entry:
            if not attribute == 'k':
                b1.add_attribute(attribute, entry[attribute])


def __tomtom_upstream_motifs():
    #################################################################
    ## TomTom Upstream motifs versus Jaspar and Transfac           ##
    #################################################################
    print 'Running TOMTOM on Upstream Motifs:'
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')

    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    pssms = c1.pssms_upstream()
    upstreamMatches = {}

    motif_files = ['motifs/jasparCoreVertebrata_redundant.pkl',
                   'motifs/transfac_2012.1_PSSMs_vertabrate.pkl',
                   'motifs/uniprobePSSMsNonRedundant.pkl',
                   'motifs/selexPSSMsNonRedundant.pkl']
    comparison_pkl_path = 'output/upstreamJasparTransfacComparison.pkl'
    comparison_csv_path = 'output/upstreamComparison_jaspar_transfac.csv'

    # WW:
    # we need to read the motif file regardless whether pkl_path exists
    # or we get an exception in the else part
    target_pssms_in = []
    for motif_file in motif_files:
        with open(motif_file, 'rb') as pklFile:
            pssms = cPickle.load(pklFile)
            for pssm in pssms.values():
                pssm.setMethod('meme')
            target_pssms_in.append(pssms)

    if not os.path.exists(comparison_pkl_path):

        # Write out results
        with open(comparison_csv_path, 'w') as outFile:
            outFile.write('Motif Name,Original E-Value,Consensus,JASPAR Motif,JASPAR Consensus,TomTom.pValue,TomTom.qValue,Probe In Bicluster,Bicluster Residual')

            # Making MEME formatted files (makeFiles function in utils)
            print 'Making files...'
            for i, target_pssms in enumerate(target_pssms_in):
                utils.makeFiles(nucFreqs=c1.getNucFreqsUpstream(),
                                queryPssms=pssms.values(),
                                targetPssms=target_pssms.values(),
                                num=i)

            # Run TomTom 
            print 'Comparing Upstream motifs against databases...'
            pool = Pool(processes=cpu_count())
            res1 = pool.map(run_tomtom, [i for i in range(len(target_pssms_in))])
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


def __expand_tf_factor_list():
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
    pssms = c1.pssms_upstream()
    for pssm in pssms:
        expanded_factors = {}
        matches = pssms[pssm].getMatches()
        # Collapse matches to a set of entrez IDs and add expanded factors
        if not matches==None:
            for match in matches:
                if match['factor'] in tfName2entrezId:
                    factor = tfName2entrezId[match['factor']]
                    if not factor in expanded_factors:
                        expanded_factors[factor] = [factor]
                        for family in tfFamilies:
                            if factor in tfFamilies[family]:
                                expanded_factors[factor] += tfFamilies[family]
            # Push expanded TF factor list into the PSSM object
            for factor in expanded_factors:
                for expanded_factor in list(set(expanded_factors[factor])):
                    pssms[pssm].addExpandedMatch(expanded_factor, factor)
    print 'Finished expanding TF factor list.\n'
    return tfName2entrezId, tfFamilies


def __correlate_tfs_with_cluster_eigengenes():
    #######################################################################
    ## Filter expanded TFs through correlation with bicluster eigengenes ##
    #######################################################################
    print 'Correlate TFs with eigengenes...'
    
    # Get microarray data
    exp_data = {}
    with open(ALL_RATIOS_PATH,'r') as inFile:
        names = inFile.readline().strip().split('\t')
        names = [i.strip('"') for i in names]
        for line in inFile:
            splitUp = line.strip().split('\t')
            geneId = splitUp.pop(0).strip('"')
            exp_data[geneId] = dict(zip(names, splitUp))

    all_names = names

    # [rho, pValue] = correlation(a1,a2)
    for b1, bicluster in c1.biclusters.items():
        for pssm in bicluster.pssms_upstream:
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
                        if factor['factor'] in exp_data.keys():
                            eigengene = biclusters[b1].attributes['pc1']
                            if not factor['factor'] in compared[subset]:
                                    compared[subset].append(factor['factor'])
                                    cor1 = correlation([eigengene[i] for i in all_names][subsetsPos[subset][0]:subsetsPos[subset][1]],
                                                       [exp_data[factor['factor']][i] for i in all_names][subsetsPos[subset][0]:subsetsPos[subset][1]])
                                    print subset, factor['factor'], cor1
                                    if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                                        corMax = cor1
                        if not corMax==[]:
                            pssm.addCorrelatedMatch(subset,factor['factor'],corMax[0],corMax[1])
    print 'Done.\n'
    return exp_data, all_names


def __expand_and_correlate_tfbsdb_tfs(tf_name2entrezid, tf_families, exp_data, all_names):
    print 'Expand and correlate TFBS_DB TFs...'

    for bicluster in c1.biclusters:
        b1 = c1.getBicluster(bicluster)
        # Get the tfbs_db attribute and for each TF get the list of expanded factors
        tfs = b1.attributes['tfbs_db']
        expanded_factors = {}
        if not tfs==None:
            for tf in tfs['tf'].split(' '):
                if tf[0:2]=='V_':
                    tf = 'V$'+tf[2:]
                # Get the list of expanded factors
                if tf in tf_name2entrezid:
                    factor = tf_name2entrezid[tf]
                    if not factor in expanded_factors:
                        expanded_factors[factor] = [factor]
                        for family_factors in tf_families.values():
                            if factor in family_factors:
                                expanded_factors[factor] += family_factors
                        expanded_factors[factor] = list(set(expanded_factors[factor]))
        
        # Push expanded TF factor list into the PSSM object
        if len(expanded_factors) > 0:
            print factor, expanded_factors
            b1.add_attribute('tfbs_db_expanded', expanded_factors)

    # [rho, pValue] = correlation(a1,a2)
    for bicluster in c1.biclusters:
        b1 = c1.getBicluster(bicluster)
        factors = b1.attributes['tfbs_db_expanded'] if 'tfbs_db_expanded' in b1.attributes else None
        compared = {}
        for subset in subsets:
            compared[subset] = []
        correlatedFactor = {}
        for subset in subsets:
            correlatedFactor[subset] = []

        # Get maximum amount of correlation positive or negative
        if not factors is None:
            for factor1 in factors:
                for factor2 in factors[factor1]:
                    for subset in subsets:
                        corMax = []
                        if factor2 in exp_data:
                            eigengene = b1.attributes['pc1']
                            if not factor2 in compared:
                                compared[subset].append(factor2)
                                cor1 = correlation([eigengene[i] for i in all_names][subsetsPos[subset][0]:subsetsPos[subset][1]],
                                                   [exp_data[factor2][i] for i in all_names][subsetsPos[subset][0]:subsetsPos[subset][1]])
                                print cor1
                                if corMax==[] or abs(cor1[0])>abs(corMax[0]):
                                    corMax = cor1
                            if not corMax==[]:
                                correlatedFactor[subset].append({'factor':factor2,'rho':corMax[0],'pValue':corMax[1]})
        b1.add_attribute('tfbs_db_correlated', correlatedFactor)
    print 'Done.\n'


def __write_first_principal_components():
    print 'Write biclusterFirstPrincComponents.csv...'
    # Get all first principal components for each bicluster
    fpcWrite = []
    conditions = c1.getBicluster(1).attributes['pc1'].keys()
    for i in sorted(c1.biclusters.keys()):
        pc1 = c1.getBicluster(i).attributes['pc1']
        fpcWrite.append(str(i)+','+','.join([str(pc1[j]) for j in conditions]))

    # Dump out file
    with open('output/'+str(fpcFile),'w') as outfile:
        outfile.write('Bicluster,'+','.join([j.strip() for j in conditions])+'\n')
        outfile.write('\n'.join(fpcWrite))
    print 'Done.\n'


def __get_permuted_pvalues_for_upstream_meme_motifs():
    # Make needed directories
    if os.path.exists('tmp'):
        rmtree('tmp')
    if not os.path.exists('tmp/tomtom_out'):
        os.makedirs('tmp/tomtom_out')
    
    # Compare the random motifs to the original motif in TOMTOM
    permPValues = {}
    matched = 0
    pssms = c1.pssms_upstream(de_novo_method='meme')
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
            clustSize = utils.randPSSMClustSize((c1.getBicluster(int(pssmsNames[i].split('_')[0]))).num_genes())
            utils.makeFiles(nucFreqs=c1.getNucFreqsUpstream(),
                            queryPssms=[pssms[pssmsNames[i]]],
                            targetPssms=randPssmsDict[clustSize].values(),
                            num=i)

        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs available.'
        print 'Running TOMTOM to compare PSSMs...'
        pool = Pool(processes=cpus)
        pool.map(run_tomtom, range(len(pssms)))
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

def __run_mirvestigator_3putr():
    """Compare 3' UTR Weeder Motifs to miRBase using miRvestigator"""
    print 'Running miRvestigator on 3\' UTR Motifs:'
    if not os.path.exists('output/m2m.pkl'):
        print 'Computing miRNA matches...'
        biclusters = c1.biclusters
        pssms = c1.pssms_3putr()
        seqs3pUTR = c1.getSeqs3pUTR().values()
        m2m = miRvestigator(pssms.values(),seqs3pUTR,seedModel=[6,7,8],minor=True,p5=True,p3=True,wobble=False,wobbleCut=0.25,baseDir='output',species='mmu')
        with open('output/m2m.pkl','wb') as pklFile:
            cPickle.dump(m2m,pklFile)
    else:
        print 'Loading precached miRNA matches...'
        with open('output/m2m.pkl','rb') as pklFile:
            m2m = cPickle.load(pklFile)
    print 'Done.\n'

def __convert_mirvestigator_3putr_results():
    """Convert miRNAs and Get permuted p-values for 3' UTR motifs"""

    pssms = c1.pssms_3putr()
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
                    miRNA_mature_seq_ids += utils.miRNAInDict(i.lower(), miRNAIDs)
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

def run_replication(rep_script):
    ## TODO:
    ## these R scripts are external to the project and therefore
    ## we can not generalize them safely. The script has a "loc1" variable which
    ## hardcodes the output path, we need to make both
    ## What actually needs to be done is to move the script into this project and
    ## make it configurable
    print '  Replication running for %s...' % rep_script
    ret = subprocess.check_call('cd replication_%s && Rscript replicationDatasetPermutation.R' % rep_script,
                                stderr=subprocess.STDOUT, shell=True)
    if ret == 1:
        raise Exception('could not run replication - check dependencies')

def __make_replication_pvalues():
    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.genes.txt','w') as cmgFile:
        writeMe = []
        for b1 in c1.biclusters:
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).genes))
        cmgFile.write('\n'.join(writeMe))

    # Dump a file containing all the genes for each cluster
    with open('output/cluster.members.conditions.txt','w') as cmcFile:
        writeMe = []
        for b1 in c1.biclusters:
            writeMe.append(str(b1)+' '+' '.join(c1.getBicluster(b1).conditions))
        cmcFile.write('\n'.join(writeMe))


    # Run replication on all datasets
    run_sets = [name for name in REPLICATION_DATASET_NAMES
                if not os.path.exists('output/replicationPvalues_%s.csv' % name)]

    if len(run_sets) > 0:
        print 'Run replication..'
        # Run this using all cores available
        cpus = cpu_count()
        print 'There are', cpus,'CPUs avialable.'
        pool = Pool(processes=cpus)
        pool.map(run_replication, run_sets)
        pool.close()
        pool.join()


def __read_replication_pvalues():
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
        b1.add_attribute(key='replication_French',value={'French_new.resid.norm':splitUp[3], 'French_avg.resid.norm':splitUp[4], 'French_norm.perm.p':splitUp[5], 'French_pc1.var.exp':splitUp[9], 'French_avg.pc1.var.exp':splitUp[10], 'French_pc1.perm.p':splitUp[11], 'French_survival':splitUp[15], 'French_survival.p':splitUp[16], 'French_survival.age':splitUp[17], 'French_survival.age.p':splitUp[18]})
        b1.add_attribute(key='replication_French_all',value={'French_all_new.resid.norm':splitUp[6], 'French_all_avg.resid.norm':splitUp[7], 'French_all_norm.perm.p':splitUp[8], 'French_all_pc1.var.exp':splitUp[12], 'French_all_avg.pc1.var.exp':splitUp[13], 'French_all_pc1.perm.p':splitUp[14], 'French_all_survival':splitUp[19], 'French_all_survival.p':splitUp[20], 'French_all_survival.age':splitUp[21], 'French_all_survival.age.p':splitUp[22]})
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
        b1.add_attribute(key='replication_REMBRANDT',value={'REMBRANDT_new.resid.norm':splitUp[3], 'REMBRANDT_avg.resid.norm':splitUp[4], 'REMBRANDT_norm.perm.p':splitUp[5], 'REMBRANDT_pc1.var.exp':splitUp[9], 'REMBRANDT_avg.pc1.var.exp':splitUp[10], 'REMBRANDT_pc1.perm.p':splitUp[11], 'REMBRANDT_survival':splitUp[15], 'REMBRANDT_survival.p':splitUp[16], 'REMBRANDT_survival.age':splitUp[17], 'REMBRANDT_survival.age.p':splitUp[18]})
        b1.add_attribute(key='replication_REMBRANDT_all',value={'REMBRANDT_all_new.resid.norm':splitUp[6], 'REMBRANDT_all_avg.resid.norm':splitUp[7], 'REMBRANDT_all_norm.perm.p':splitUp[8], 'REMBRANDT_all_pc1.var.exp':splitUp[12], 'REMBRANDT_all_avg.pc1.var.exp':splitUp[13], 'REMBRANDT_all_pc1.perm.p':splitUp[14], 'REMBRANDT_all_survival':splitUp[19], 'REMBRANDT_all_survival.p':splitUp[20], 'REMBRANDT_all_survival.age':splitUp[21], 'REMBRANDT_all_survival.age.p':splitUp[22]})
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
        b1.add_attribute(key='replication_GSE7696',value={'GSE7696_new.resid.norm':splitUp[3], 'GSE7696_avg.resid.norm':splitUp[4], 'GSE7696_norm.perm.p':splitUp[5], 'GSE7696_pc1.var.exp':splitUp[6], 'GSE7696_avg.pc1.var.exp':splitUp[7], 'GSE7696_pc1.perm.p':splitUp[8], 'GSE7696_survival':splitUp[9], 'GSE7696_survival.p':splitUp[10], 'GSE7696_survival.age':splitUp[11], 'GSE7696_survival.age.p':splitUp[12]})
    inFile.close()
    print 'Done.\n'

def __make_permuted_pvalues():
    ###########################################################################
    ## Run permuted p-value for variance epxlained first principal component ##
    ###########################################################################
    if not os.path.exists('output/residualPermutedPvalues_permAll.csv'):
        print 'Calculating FPC permuted p-values...'
        ret = subprocess.check_call(['./permutedResidualPvalues_permAll_mc.R',
                                     '-b', '..'],
                                    stderr=subprocess.STDOUT)
        if ret == 1:
            print "error in calling R script."
            exit(1)

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
            b1.add_attribute(key='resid.norm.perm.p',value=str(splitUp[9]))
            b1.add_attribute(key='pc1.perm.p',value=str(splitUp[12]))
    print 'Done.\n'


def __make_functional_enrichment_and_go_term_similarity():
    #################################################################
    ## Run functional enrichment and GO term similarity            ##
    #################################################################
    # Note that these are external to the project and have hard-coded paths !!!
    if not os.path.exists('output/biclusterEnrichment_GOBP.csv'):
        print 'Run functional enrichment...'
        ret = subprocess.check_call("cd funcEnrichment && Rscript enrichment.R -o %s" % OUTPUT_DIR,
                                    stderr=subprocess.STDOUT, shell=True)
        if ret == 1:
            raise Exception('could not run functional enrichment')

        print 'Done.\n'

    if not os.path.exists('output/jiangConrath_hallmarks.csv'):
        print 'Run semantic similarity...'
        ret = subprocess.check_call("cd funcEnrichment && Rscript goSimHallmarksOfCancer.R -o %s" % OUTPUT_DIR,
                                    stderr=subprocess.STDOUT, shell=True)
        if ret == 1:
            raise Exception('could not run semantic similarity')

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
        b1.add_attribute(key='goTermBP',value=line[2].strip('"').split(';'))
    print 'Done.\n'


def __add_hallmarks_of_cancer():    
    print 'Load Jiang-Conrath semantic similarity to Hallmarks of Cancer...'
    with open('output/jiangConrath_hallmarks.csv','r') as inFile:
        hallmarks = [i for i in inFile.readline().split(',') if not i.strip('"')=='']
        inLines = inFile.readlines()
        lines = [line.strip().split(',') for line in inLines]
        for line in lines:
            b1 = c1.getBicluster(int(line[0].strip('"')))
            b1.add_attribute(key='hallmarksOfCancer',value=dict(zip(hallmarks,line[1:])))
        print 'Done.\n'

def __write_final_cmonkey_object():
    print 'Dumping Final cMonkey Object:'
    with open('output/c1_postProc.pkl','wb') as pklFile:
        cPickle.dump(c1, pklFile)
    print 'Done.\n'

def perform_postprocessing(c1):
    if not os.path.exists('output/c1_postProc.pkl'):
        ratios = __read_ratios(RATIOS_PATH, c1)
        __get_cluster_eigengenes(RATIOS_PATH, c1)
        __get_cluster_variance_explained(CLUSTER_VARIANCE_EXPLAINED_PATH, c1)
        phenotypes = __get_phenotype_info(PHENOTYPES_PATH, c1)
        __do_postprocess(POSTPROCESS_PKL_PATH, c1, ratios, phenotypes)
        __tomtom_upstream_motifs()
        tf_name2entrezid, tf_families = __expand_tf_factor_list()
        exp_data, all_names = __correlate_tfs_with_cluster_eigengenes()
        __expand_and_correlate_tfbsdb_tfs(tf_name2entrezid, tf_families, exp_data, all_names)
        __write_first_principal_components()
        __get_permuted_pvalues_for_upstream_meme_motifs()
        __run_mirvestigator_3putr()
        __convert_mirvestigator_3putr_results()
        __make_replication_pvalues()
        __read_replication_pvalues()
        __make_permuted_pvalues()
        __make_functional_enrichment_and_go_term_similarity()
        __add_hallmarks_of_cancer()
        __write_final_cmonkey_object()

    else:
        print 'Loading from precached cMonkey Object:'
        with open('output/c1_postProc.pkl','rb') as pklFile:
            c1 = cPickle.load(pklFile)
        print 'Done.\n'
    return c1


def run_neo():
    """Run NEO and integrate some form of results
    TODO: Make NEO run on subsets"""
    if not os.path.exists('output/causality'):
        ## Run the runNEO.R script and do the causality analyses
        print '  Network edge orienting (NEO)...'
        ret = subprocess.check_call("cd NEO && Rscript runNEO.R -b %s" % BASE_DIR,
                                    stderr=subprocess.STDOUT, shell=True)
        if ret == 1:
            raise Exception('could not run causality analyses')

def write_neo_summary():
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
    return causalSummary


def add_correspondent_regulators(c1, causal_summary):
    """Dump out correspondent regulators (both mechanistically and causally predicted)"""
    correspondentRegulators = {}
    for causalFlow in causal_summary:
        b1 = c1.getBicluster(int(causalFlow['Bicluster']))
        ## Upstream (TFs)
        tfs = []
        # 1. MEME and WEEDER Upstream motifs
        for pssm in b1.pssms_upstream:
            for subset in subsets:
                matches = pssm.getCorrelatedMatches(subset)
                if matches:
                    for corTf in matches:
                        if corTf['pValue']<=pVCut and abs(corTf['rho'])>=rhoCut:
                            tfs.append(corTf['factor'])
            # 2. TFBS_DB
            for corTf in b1.attributes['tfbs_db_correlated'][subset]:
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
        for pssm in b1.pssms_3putr:
            matches = pssm.getMatches()
            if matches:
                for miR in matches:
                    if miR['confidence'] in ['8mer','7mer_a1','7mer_m8']:
                        miRNAs += miRNAIDs_rev[miR['factor']]
        # 2. PITA (not correlated)
        if float(b1.attributes['pita_3pUTR']['pValue'])<=pVCut and float(b1.attributes['pita_3pUTR']['percentTargets'].split(' ')[0])>=percTargets:
            miRNAs += b1.attributes['pita_3pUTR']['miRNA'].split(' ')
        # 3. TargetScan (not correlated)
        if float(b1.attributes['targetscan_3pUTR']['pValue'])<=pVCut and float(b1.attributes['targetscan_3pUTR']['percentTargets'].split(' ')[0])>=percTargets:
            miRNAs += b1.attributes['targetscan_3pUTR']['miRNA'].split(' ')
        # 4. Find Correspondent miRNA regulators
        for miR in miRNAs:
            if compareMiRNANames(causalFlow['Regulator'].lower(), miR.lower()):
                if not int(causalFlow['Bicluster']) in correspondentRegulators:
                    correspondentRegulators[int(causalFlow['Bicluster'])] = {'tf':[],'miRNA':[]}
                correspondentRegulators[int(causalFlow['Bicluster'])]['miRNA'].append(causalFlow['Regulator'])

    ## Put correspondent regulators into cMonkeyWrapper object
    for biclust in correspondentRegulators.keys():
        b1 = c1.getBicluster(biclust)
        b1.add_attribute(key='correspondentRegulators',value=correspondentRegulators[biclust])


def write_final_result(c1):
    #################################################################
    ## Write out the final post-processed file                     ##
    #################################################################
    print 'Write postProcessedVFinal.csv...'
    postOut = []
    hallmarksOfCancer = c1.getBicluster(1).attributes['hallmarksOfCancer'].keys()
    for i in sorted(c1.biclusters.keys()):
        writeMe = []
        b1 = c1.getBicluster(i)
        # Write file line by line
        #   a. Bicluster basics:  id, genes, conditions, resid, resid.norm, resid.norm.perm.p
        writeMe += [str(i), # id
                    str(b1.attributes['k.rows']), # genes
                    str(b1.attributes['k.cols']), # conditions
                    str(b1.norm_residual), # normalized residual
                    str(b1.attributes['resid.norm.perm.p']), # normalized residual permuted p-value
                    str(b1.attributes['pc1.var.exp']), # Variance explained by first principal component
                    str(b1.attributes['pc1.perm.p'])] # Variance explained by first principal component
        #   b. Upstream motifs:  meme.motif1.E, meme.motif1.consensus, meme.motif1.matches, meme.motif1.permPV
        motifNames = b1.pssm_names_upstream()
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
                pssm1 = b1.pssm_upstream(upstreamMotifs[meme1])
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
                pssm1 = b1.pssm_upstream(upstreamMotifs[weeder1])
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
            a1 = b1.attributes[association]
            if a1['tf'] != '':
                expandedMatches = 'NA'
                correlatedMatches = {}
                originalExpanded = {}
                minCorrelated = {}
                for subset in subsets:
                    correlatedMatches[subset] = 'NA'
                    originalExpanded[subset] = 'NA'
                    minCorrelated[subset] = 'NA'
                tmp = b1.attributes['tfbs_db_expanded']
                if not tmp is None:
                    expandedMatches = ' '.join([seedFactor+':'+';'.join(tmp[seedFactor]) for seedFactor in tmp])
                tmp = b1.attributes['tfbs_db_correlated']
                for subset in subsets:
                    if not tmp[subset] is None:
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
        motifNames = b1.pssm_names_3putr()
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
                pssm1 = b1.pssm_3putr(p3utrMotifs[weeder1])
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
            a1 = b1.attributes[association]
            if a1['miRNA'] != '':
                writeMe += [str(a1['miRNA']).replace(';',' '), str(a1['percentTargets']).replace(';',' '), str(a1['pValue'])]
            else:
                writeMe += ['NA','NA','NA']

        #   f. Correpsondent regulators
        corrRegs = b1.attributes['correspondentRegulators'] if 'correspondentRegulators' in b1.attributes else None
        if not corrRegs is None:
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
            ass1 = b1.attributes[association]
            writeMe += [str(ass1['rho']), str(ass1['pValue'])]
        surv1 = b1.attributes['Survival']
        survAge1 = b1.attributes['Survival.AGE']
        writeMe += [str(surv1['z']), str(surv1['pValue']), str(survAge1['z']), str(survAge1['pValue'])]

        #   h. Replications:  'REMBRANDT_new.resid.norm','REMBRANDT_avg.resid.norm','REMBRANDT_norm.perm.p','REMBRANDT_survival','REMBRANDT_survival.p','REMBRANDT_survival.age','REMBRANDT_survival.age.p','GSE7696_new.resid.norm','GSE7696_avg.resid.norm','GSE7696_norm.perm.p','GSE7696_survival','GSE7696_survival.p','GSE7696_survival.age','GSE7696_survival.age.p'
        replications_French = b1.attributes['replication_French']
        replications_REMBRANDT = b1.attributes['replication_REMBRANDT']
        replications_French_all = b1.attributes['replication_French_all']
        replications_REMBRANDT_all = b1.attributes['replication_REMBRANDT_all']
        replications_GSE7696 = b1.attributes['replication_GSE7696']
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
        bfe1 = b1.attributes['goTermBP']
        if bfe1 == ['']:
            writeMe.append('NA')
        else:
            writeMe.append(';'.join(bfe1))

        #   j. Hallmarks of Cancer:  Hanahan and Weinberg, 2011
        bhc1 = b1.attributes['hallmarksOfCancer']
        for hallmark in hallmarksOfCancer:
            writeMe.append(str(bhc1[hallmark]))

        #   k. Glioma sub-type enrichment: 'NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM'
        #for overlap in ['NON_TUMOR','ASTROCYTOMA','MIXED','OLIGODENDROGLIOMA','GBM']:
        #    writeMe.append(str(b1.attributes[overlap]))

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

if __name__ == '__main__':
    # Initialize sygnal output directory and conversion dictionaries
    sygnal_init()
    id2entrez, entrez2id = read_synonyms()
    miRNAIDs, miRNAIDs_rev = miRNA_mappings()
    c1 = compute_additional_info()
    c1 = perform_postprocessing(c1)
    run_neo()
    causal_summary = write_neo_summary()
    add_correspondent_regulators(c1, causal_summary)
    write_final_result(c1)

