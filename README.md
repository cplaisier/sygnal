## Systems Genetic Network AnaLysis (SYGNAL) pipeline
We developed the SYstems Genetic Network AnaLysis (SYGNAL) pipeline to integrate correlative, causal and mechanistic inference approaches into a unified framework that systematically infers the causal flow of information from mutations to TFs and miRNAs to perturbed gene expression patterns across patients.

The SYGNAL pipeline is designed to be cloned into a completed [cMonkey<sub>2</sub>](https://github.com/baliga-lab/cmonkey2) run directory. It will then run the SYGNAL pipeline using the cMonkey<sub>2</sub> result database and put all the results into 'sygnal/output' directory, which it will create if not already present. There is fairly through checkpointing as the full SYGNAL pipeline can take some time to run.

### Dependencies
* MEME (http://meme-suite.org/doc/download.html?man_type=web)
* WEEDER (https://github.com/baliga-lab/weeder_patched)
* R (https://cran.r-project.org/)
* rpy2 (http://rpy2.bitbucket.org/)
* AME (http://bioinformatics.org.au/tools/ame/)

Install with:
```
sudo apt-get install r-base r-base-dev
sudo pip install rpy
```
All other programs will have to installed using their installers and following instructions from those software.

### Required data files
As described the SYGNAL pipeline requires user supplied data for patient expression profiles and patient clinical information. It can be run in a modified form without clinical information. Expression data is usually filtered based on differential expression, most variant genes, or expression above a threshold. The reason for this step is that genes with little variance or low expression are not likely to yield much information.

In addition to user supplied expression data and clinical data the pipeline requires additional information in order to run:
* Gene ID thesaurus
* TF-target gene regulatory interactions in JSON format
* TF family information
* miRNA-target gene regulatory interactions in JSON format
* Full expression data including all genes (for TF correlations)
* Expression data and clinical information for replication studies
* Promoter and 3' UTR sequences for organism of study
* Background sequence information for motif callers

### List of parameters
####TF target-gene database construction				
```
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 0.001 -min-overlap 6		
fimo	(4.11.1)	--max-stored-scores 1000000 --verbosity 4 --bgfile bgFile.meme -text --thresh 1e-5 		
```				
####SYGNAL pipeline
```
cmonkey.py	(cMonkey2)	--organism hsa --ratios gbmTCGA_exprMat_medianFiltered.tsv --pipeline pipeline.json --config gbmTCGA.ini --rsat_organism Homo_sapiens --rsat_dir data/Homo_sapiens --string data/string/string.tsv --synonym_file synonymsThesaurus.csv --verbose --case_sensitive
```
#####pipeline.json
```
{
    ""row-scoring"": {
        ""id"": ""combiner"",
        ""function"": { ""module"": ""cmonkey.scoring"", ""class"": ""ScoringFunctionCombiner"" },
        ""args"": {
            ""functions"": [
                { ""id"": ""Rows"",
                  ""function"": { ""module"": ""cmonkey.microarray"", ""class"": ""RowScoringFunction"" }
                },
                { ""id"": ""Networks"",
                  ""function"": { ""module"": ""cmonkey.network"", ""class"": ""ScoringFunction"" }
                },
                { ""id"": ""SetEnrichment"",
                  ""function"": { ""module"": ""cmonkey.set_enrichment"", ""class"": ""ScoringFunction"" }
                }
            ]
        }
    },
    ""column-scoring"": { ""id"": ""Columns"",
                        ""function"": { ""module"": ""cmonkey.scoring"",
                                      ""class"": ""ColumnScoringFunction""} }
}
```

#####gbmTCGA.ini
```
[General]
num_iterations = 2000
num_clusters = 134
normalize_ratios = False
num_cores = 6

[SetEnrichment]
schedule = 1,7
scaling_rvec=seq(1e-5, 1.0, length=num_iterations/2)
map_to_ratio_genes = True
set_types = pita

[SetEnrichment-pita]
set_file = pita_miRNA_sets_geneSymbol.json
weight = 1.0
```
#####Motif searching
```
meme	(4.11.1)	 -bfile bgFile.meme -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw 6 -maxw 12 -nmotifs 2		
weederlauncher	(1.4.2 ISB Patch)	HS small T50 S	For upstream seqeunce searching.	
weederlauncher	(1.4.2 ISB Patch)	HS3P small T50	For 3' UTR seqeunce searching.	
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 1 -min-overlap 6 -verbosity 1		
ame	(Original Version)	--method mhg --scoring totalhits --verbose 1 --bgformat 2 --bgfile bgFile.meme		
```
### Configuration
The SYGNAL pipeline is designed to be configured using the parameters starting at line 58 to 108 of sygnal.py. We will work to move the to a more convenient configuration file in the future.

#### For MEME analysis
MEME_BGFILE       = 'seqs/bgFile.meme'
MEME_NMOTIFS      = 2
MEME_MOTIF_WIDTHS = {'upstream': [6, 12]}
MEME_REVCOMP      = {'upstream': True}

#### Parameters for filtering results
MAX_EVALUE        = 10.0
LEO_NB_ATOB       = 0.5  # Equates to ~3 times better model fit than next best model
MLOGP_M_ATOB      = 0.05
RHO_CUT           = 0.3
PVALUE_CUT        = 0.05
PERC_TARGETS      = 0.1
POSTPROCESSED_FILENAME = 'postProcessed_gbmTCGA_pita.csv'
FPC_FILENAME = 'biclusterFirstPrincComponents.csv'
SUBSETS = ['all'] # Might be nice to include subtypes
SUBSETS_POS = { 'all': [0,422] } # Might be nice to include subtypes
RAND_PSSMS_DIR = 'randPSSMs'

#### Detenrmine the working directories, we pass this to the
#### external R scripts
CURRENT_DIR = os.getcwd()
#### note that we need to paste in the '/' at the end because the R scripts
#### currently assume it
BASE_DIR = os.path.dirname(os.getcwd()) + '/'
OUTPUT_DIR = os.path.join(CURRENT_DIR, 'output')

SYNONYM_PATH = '../synonymThesaurus.csv.gz'
MIRNA_FASTA_PATH = 'miRNA/hsa.mature.fa'  # Need to update this to the latest database
RATIOS_PATH = '../gbmTCGA_exprMat_medianFiltered.tsv'
ALL_RATIOS_PATH = 'expression/gbmTCGA_exprMat.tsv'
CMONKEY2_RUNDB = '../out/cmonkey_run.db'
REPLICATION_DATASET_NAMES = ['French','REMBRANDT','GSE7696']
CM_PKL_PATH = 'output/c1_all.pkl'

MOTIF_FILES = ['motifs/jasparCoreVertebrata_redundant.json',
               'motifs/transfac_2012.1_PSSMs_vertabrate.json',
               'motifs/uniprobePSSMsNonRedundant.json',
               'motifs/selexPSSMsNonRedundant.json']


CLUSTER_GENES_PATH = 'output/cluster.members.genes.txt'
CLUSTER_CONDS_PATH = 'output/cluster.members.conditions.txt'
CLUSTER_EIGENGENES_PATH = 'output/biclusterEigengenes.csv'
CLUSTER_VARIANCE_EXPLAINED_PATH = 'output/biclusterVarianceExplained.csv'
PHENOTYPES_PATH = 'extras/phenotypes.csv'
POSTPROCESS_PKL_PATH = 'output/postProcessed.pkl'

### Output and visualizations
