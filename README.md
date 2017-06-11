## Systems Genetic Network AnaLysis (SYGNAL) pipeline
We developed the SYstems Genetic Network AnaLysis (SYGNAL) pipeline to integrate correlative, causal and mechanistic inference approaches into a unified framework that systematically infers the causal flow of information from mutations to TFs and miRNAs to perturbed gene expression patterns across patients.

The SYGNAL pipeline is designed to be cloned into a completed [cMonkey<sub>2</sub>](https://github.com/baliga-lab/cmonkey2) run directory. It will then run the SYGNAL pipeline using the cMonkey<sub>2</sub> result database and put all the results into 'sygnal/output' directory, which it will create if not already present. There is fairly through checkpointing as the full SYGNAL pipeline can take some time to run.

### Dependencies
#### Applications
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

#### R Dependencies
* WGCNA
* impute
* getopt
* topGO
* org.Hs.eg.db or org.Mm.eg.db depending on species
* GOSim
* multicore

Make sure to set your repositories to include the Bioconductor repositories with setRepositories() and install the apporiate packages through the commands below (you will also have to select the closest mirror as well):

```
setRepositories()
install.packages(c('WGCNA','impute','getopt','topGO','org.Mm.eg.db','org.Hs.eg.db','GOSim','multicore'))
```

In order to run the causality analysis portions of SYGNAL requires the use of the Network Edge Orienting (NEO) package which can be found here:  https://labs.genetics.ucla.edu/horvath/aten/NEO/  The 'neoDecember2015.txt' is the source code needed for running NEO in SYGNAL.

### Required data files
The SYGNAL pipeline requires user supplied data for patient expression profiles and patient clinical information. It can be run in a modified form without clinical information. Expression data is usually filtered based on differential expression, most variant genes, or expression above a threshold. The reason for this step is that genes with little variance or low expression are not likely to yield much information.

* User supplied expression data (microarray or RNAseq) as a tab separated value (TSV) file that has been properly normalized. The expression data should be standardized (mean subtracted and divided by the standard deviation to turn each expression value into a Z score). The header for the expression data is expected to be without the leading tab.

| Cond_1  | ...    | Cond_N   |        |
|:-------:|:------:|:--------:|:------:|
| Gene_1  | Z[1,1] | ...      | Z[1,N] |
| ...     | ...    | ...      | ...    |
| Gene_M  | Z[M,1] | Z[M,N-1] | Z[M,N] |

In addition to user supplied expression data and clinical data the pipeline requires additional information in order to run:
* Gene ID thesaurus - a CSV file with the following format:

| ExprIDs  | Mappings               |
|:--------:|:----------------------:|
| Gene_1   | Gene_1;Alt_1;...;Alt_N |
| ...      | ...;Alt_1;...;Alt_N    |
| Gene_M   | Gene_M;Alt_1;...;Alt_N |

* TF-target gene regulatory interactions (http://tfbsdb.systemsbiology.net/) in the following JSON format:
```
{"TF_ID_1": ["Gene_1", ..., "Gene_M"], ..., "TF_ID_N": ["Gene_1", ..., "Gene_M"]}
```

* TF family information (http://tfclass.bioinf.med.uni-goettingen.de/tfclass)
* miRNA-target gene regulatory interactions (https://genie.weizmann.ac.il/pubs/mir07/mir07_data.html; http://www.targetscan.org/vert_71/) in the following JSON format:
```
{"miRNA_ID_1": ["Gene_1", ..., "Gene_M"], ..., "miRNA_ID_N": ["Gene_1", ..., "Gene_M"]}
```

* Full expression data including all genes (for TF correlations)

* Expression data and clinical information for replication studies
* Promoter and 3' UTR sequences for organism of study ()
* Background sequence information for motif callers

### List of parameters
#### TF target-gene database construction				
```
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 0.001 -min-overlap 6		
fimo	(4.11.1)	--max-stored-scores 1000000 --verbosity 4 --bgfile bgFile.meme -text --thresh 1e-5 		
```				
#### SYGNAL pipeline

```
./sygnal.py [--config <cfg.json>] [--out <output-dir>] [--tmp <tmp-dir>] <cmonkey2-rundb>
```
### Configuration

The SYGNAL pipeline is configured by editing the parameters in the file sygnal_config.json.

##### sygnal_config.json

```
{
    "max-evalue": 10.0,
    "rho-cut": 0.3,
    "pvalue-cut": 0.05,
    "perc-targets": 0.1,
    "leo-nb-atob": 0.5,
    "mlogp-m-atob": 0.05,

    "run-neo": true,

    "synonyms-file": "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/synonymsThesaurus.csv.gz",
    "phenotypes-file": "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/phenotypes_meso.csv",
    "ratios-file": "/proj/omics4tb/cplaisier/cMonkey/meso/exprs/mesothelioma_norm_cv_z.tsv",
    "all-ratios-file": "/proj/omics4tb/cplaisier/cMonkey/meso/exprs_all/mesothelioma_norm.csv",
    "som-muts-file": "/proj/omics4tb/cplaisier/cMonkey/meso/som_muts/somaticMutationsMatrix.csv",
    "mirna-fasta-file": "miRNA/hsa.mature.fa",
    "rand_pssms_dir": "randPSSMs",
    "promoterSeq": "seqs/promoterSeqs_Homo_sapiens.csv.gz",
    "p3utrSeq": "seqs/p3utrSeqs_Homo_sapiens.csv.gz",
    "gene_conv": "extras/entrez2ucsc.csv",
    "replication-dataset-names": ["mesoTCGA"],
    
    "pita": "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/pita_miRNA_sets_entrez_hsa.json",
    "targetscan": "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/targetScan_miRNA_sets_entrez_hsa.json",
    "tfbsdb": "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/tfbsDb_plus_and_minus_5000_entrez.json",

    "meme":  { "upstream": {
                    "nmotifs": 2,
                    "motif-widths": [6, 12],
                    "revcomp": "True",
                    "bgfile": "seqs/bgFile.meme"
                }
            },
    
    "weeder": { "upstream": {
                    "bgfile": "HS",
                    "size": "small",
                    "enriched": "T50",
                    "revcomp": "True"
                },
                "3pUTR": {
                    "bgfile": "HS3P",
                    "size": "small",
                    "enriched": "T50",
                    "revcomp": "False"
                }
            },

    "tomtom": { "upstream":
                {
                    "motif-files": ["/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/jasparCoreVertebrata_redundant.json",
                            "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/transfac_2012.1_PSSMs_vertabrate.json",
                            "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/uniprobePSSMsNonRedundant.json",
                            "/proj/omics4tb/cplaisier/cMonkey/meso/commonFiles/selexPSSMsNonRedundant.json"]
                }
            },
    
    "tfExpansion": { "tfs": "TF/humanTFs_All.csv",
                     "tfFamilies": "TF/tfFamilies.csv"
            },

    "mirvestigator": { "seedModel": [6,7,8],
                       "minor": "True",
                       "minor": "True",
                       "p5": "True",
                       "p3": "True",
                       "wobble": "False",
                       "wobbleCut": 0.25,
                       "species": "hsa"
    },
    
    "first-principal-comps-result-file": "biclusterFirstPrincComponents.csv"
}
```

##### Motif searching
```
meme	(4.11.1)	 -bfile bgFile.meme -nostatus -text -time 600 -dna -maxsize 9999999 -evt 1e9 -mod zoops -minw 6 -maxw 12 -nmotifs 2		
weederlauncher	(1.4.2 ISB Patch)	HS small T50 S	For upstream seqeunce searching.	
weederlauncher	(1.4.2 ISB Patch)	HS3P small T50	For 3' UTR seqeunce searching.	
tomtom	(4.11.1)	-dist ed -o tmp/tomtom_out -text -thresh 1 -min-overlap 6 -verbosity 1		
ame	(Original Version)	--method mhg --scoring totalhits --verbose 1 --bgformat 2 --bgfile bgFile.meme		
```


### Output and visualizations
