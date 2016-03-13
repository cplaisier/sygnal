
# R functions for network edge orienting (NEO)
# This document contains additional network analysis functions that accompany 
# the main function file neo.txt.
# WARNING: many of the functions will not work unless you read neo.txt into the R session.
# The R functions were written by Steve Horvath and Jason Aten.
# To cite this code or the statistical methods please use the following reference:
#Aten JE, Fuller TF, Lusis AJ, Horvath S (2008) Using genetic markers to orient the edges in quantitative trait networks: the NEO software. 
# Tutorials, data and R software code can be downloaded from http://www.genetics.ucla.edu/labs/horvath/aten/NEO/

# Version: March 5, 2008


# The function DirectedNetworkPlot can be used to visualize a directed 
# (causal) network. The rows and columns correspond to traits (nodes).
# The larger the entry CausalScoresAtoB[i,j] the stronger is the 
# evidence that there is a causal edge from trait i to trait j.
# If you want to sort the rows and columns according to an adjacency matrix 
# (undirected network), please either specifyADJ or provide a data frame of 
# traits. ADJ should be a symmetric matrix and its entries should lie betwen 0 and 1.
# The adjacencies encode the pairwise connection strengths between the traits
# The columns of datTraits should correspond to the columns of CausalScoresAtoB.
# Inputs:  causal scores, adjacency matrix
if (exists("DirectedNetworkPlot")) rm(DirectedNetworkPlot);
DirectedNetworkPlot=function(CausalScoresAtoB,  ADJ=NULL, datTraits=NULL , main="", ...) {
if ( dim(as.matrix(CausalScoresAtoB))[[1]] != dim(as.matrix(CausalScoresAtoB))[[2]] ) {stop("Input error. Matrix should be square. dim(as.matrix(CausalScoresAtoB))[[1]] != dim(as.matrix(CausalScoresAtoB))[[2]]")}
if (!is.null(datTraits))  {
if ( dim(as.matrix(datTraits))[[2]] != dim(as.matrix(CausalScoresAtoB))[[2]] ) stop("Input error in function DirectedNetworkPlot. The number of traits does not match the number of columns in the causal score matrix, i.e. dim(as.matrix(datTraits))[[2]] != dim(as.matrix(CausalScoresAtoB))[[2]]. Hint: consider transposing one of the matrices.")
if ( !is.null(ADJ) ) warning("You supplied both an adjacency matrix and a data frame of traits to the function DirectedNetworkPlot. This function will re-define the adjacency matrix as follows: ADJ=abs(cor(datTraits)). Hint: If you want to use the adjacency matrix simply set datTraits=NULL ")
 ADJ=abs(cor(datTraits)) 
}
if ( !is.null(ADJ) ) { if ( dim(as.matrix(ADJ))[[1]] != dim(as.matrix(ADJ))[[2]] ) {stop("Input error. Adjacency matrix ADJ should be square. ")}
if (dim(as.matrix(CausalScoresAtoB))[[1]] != dim(as.matrix(ADJ))[[2]] ) {stop("Input error. Causal Score matrix does not have the same dimensions as the adjacency matrix ADJ")  }
} # end of if ( !is.null(ADJ) )
no.traits=dim(as.matrix(CausalScoresAtoB))[[2]]
if (no.traits != dim(CausalScoresAtoB)[[1]] ) { stop("ERROR: CausalScoresAtoB. no.traits != dim(CausalScoresAtoB)[[1]] ")} 
hier1=NULL
if (!is.null(ADJ))  { 
dissADJ=1-ADJ
hier1=hclust(as.dist(dissADJ), method="average" )
#labeltree=as.character(couleur)
#labelrow  = labeltree
#labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
hier1=as.dendrogram(hier1)
options(expressions = 10000)
heatmap(-as.matrix(CausalScoresAtoB),Rowv=hier1,Colv= hier1, scale="none",revC=T,main=main)
} # end of if (!is.null(ADJ))  
if (is.null(ADJ))  { 
options(expressions = 10000)
heatmap(-as.matrix(CausalScoresAtoB), scale="none",revC=F,main=main, Rowv=NA, Colv=NA)
} # end of if (!is.null(ADJ))  
} #end of function DirectedNetworkPlot









# The following function takes the output of NEO as input and outputs 
# edge orienting scores between the traits (columns) specified in datTraits.
# if the option correlationScores=T then the function outputs the ZEO and BLV scores as well.
if (exists("GetEdgeScores") ) rm(GetEdgeScores) ;   
GetEdgeScores=function(NEOoutput1, datTraits, correlationScores=F){ 
no.traits=dim(as.matrix(datTraits))[[2]]
LEO.NB.OCA.AtoB= data.frame(matrix(NA, nrow=no.traits,ncol=no.traits))
LEO.NB.CPA.AtoB = data.frame(matrix(NA, nrow=no.traits,ncol=no.traits))
ZEO.AtoB= data.frame(matrix(NA, nrow=no.traits,ncol=no.traits))
BLV.AtoB= data.frame(matrix(NA, nrow=no.traits,ncol=no.traits))
dimnames(ZEO.AtoB)[[1]]=dimnames(data.frame(datTraits))[[2]]
dimnames(BLV.AtoB)[[1]]=dimnames(data.frame(datTraits))[[2]]
dimnames(ZEO.AtoB)[[2]]=dimnames(data.frame(datTraits))[[2]]
dimnames(BLV.AtoB)[[2]]=dimnames(data.frame(datTraits))[[2]]
dimnames(LEO.NB.OCA.AtoB)[[1]]=dimnames(data.frame(datTraits))[[2]]
dimnames(LEO.NB.CPA.AtoB)[[1]]=dimnames(data.frame(datTraits))[[2]]
dimnames(LEO.NB.OCA.AtoB)[[2]]=dimnames(data.frame(datTraits))[[2]]
dimnames(LEO.NB.CPA.AtoB)[[2]]=dimnames(data.frame(datTraits))[[2]]
for ( i in c(1:no.traits) ) {
for ( j in c(1:no.traits) ) {
AA=as.matrix(datTraits)[,i]
BB=as.matrix(datTraits)[,j]
abs.cor= abs(cor(AA,BB, use= "p"))
if( abs.cor>.98 ) warning(paste("Traits",  names(data.frame(datTraits))[[i]], "and" , names(data.frame(datTraits))[[j]]  ,  "have absolute correlation higher than .98. No LEO.NB scores will be computed for this edge."  )   )
if ( abs.cor<=0.98){
edge.name=paste(names(data.frame(datTraits))[[i]], "->" , names(data.frame(datTraits))[[j]] )
edgeIndicator=NEOoutput1$stats$edge==edge.name
LEO.NB.OCA.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats[edgeIndicator, "LEO.NB.OCA"][1]))
LEO.NB.CPA.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats[edgeIndicator, "LEO.NB.CPA"][1]))
ZEO.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats[edgeIndicator, "zeo.for"][1]))
BLV.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats[edgeIndicator, "BLV.or.BilayerZscore"][1]))
} # end of if 
} #end of for loop over i
} #end of for loop over j
out1=list(LEO.NB.OCA.AtoB= LEO.NB.OCA.AtoB,  LEO.NB.CPA.AtoB= LEO.NB.CPA.AtoB) 
if (correlationScores ) out1=list(LEO.NB.OCA.AtoB= LEO.NB.OCA.AtoB, LEO.NB.CPA.AtoB= LEO.NB.CPA.AtoB, ZEO.AtoB= ZEO.AtoB, BLV.AtoB= BLV.AtoB )
out1
} # end of function GetEdgeScores







 partialcorOLDversion=function(A,B,C){
corAB=cor(A,B,use="p")
corAC= cor(A,C,use="p")
corBC= cor(B,C,use="p")
if ( is.na(corAC) | is.na(corBC) ) out1=NA
if ( !is.na(corAC) & !is.na(corBC)   & ( abs(corAC)==1 | abs(corBC)==1)   ) out1=0
if (!is.na(corAC) & !is.na(corBC) & abs(corAC) < 1 & abs(corBC) < 1 ) {out1=(corAB-corAC*corBC)/sqrt((1-corAC^2)* (1-corBC^2))}
as.numeric(out1)
} # end of function partialcor




# this function computes the partial correlation coefficient cor(A,B|C)
#, i.e. the correlation between A and B conditional on C.
# One of the entries can be a data frame
if(exists("partialcor") ) rm(partialcor);
partialcor=function(A,B,C){
nA=dim(as.matrix(A))[[1]]
nB=dim(as.matrix(B))[[1]]
nC=dim(as.matrix(C))[[1]]
if (nA != nB | nA != nC ) stop("Input error in function partialcor: nA != nB | nA != nC")
corAB=cor(A,B,use="p")
corAC= cor(A,C,use="p")
corBC= cor(B,C,use="p")
out1=suppressWarnings(
ifelse( is.na(corAC) | is.na(corBC) , NA, 
ifelse( !is.na(corAC) & !is.na(corBC)   & ( abs(corAC)==1 | abs(corBC)==1)  , 0,
(corAB-corAC*corBC)/sqrt((1-corAC^2)* (1-corBC^2)))))
as.numeric(out1)
} # end of function partialcor




# this function caries out the Fisher transform of a correlation coefficient r given
# n are the degrees of freedom, e.g n=no.obs-2 when dealing with a correlation 
# n=no.obs-3 when dealing with a partial correlation 
if (exists("ZFisher") ) rm(ZFisher); 
ZFisher=function(r, n) {
if ( !(n>=0) || is.na(n)  ) stop("Input error in function ZFisher. Missing or not enough degrees of freedom. n<0")
if (sum(abs(as.numeric(r))>1,na.rm=T)>0 ) stop("Input error in function ZFisher. correlation r is larger than 1 or smaller than -1.")
r=ifelse( abs(r)==1, r/1.0000001, r) 
Z=.5*sqrt(n) *log((1+r)/(1-r)) 
Z
} # end of function ZFisher



# The following funcion inputs a data frame of SNPs, and two trait vectors A and B.
# Output is a vector of length=no.SNPs. 
# For each SNP it outputs a score. In the single SNP marker case, an NIS score larger than 1 indicates 
# that the SNP seems to act causally SNP->A->B as opposed to the independence model A<-SNP->B
if(exists("CausalVSIndependentSNPvote") ) rm(CausalVSIndependentSNPvote)
CausalVSIndependentSNPvote=function(MA,A,B,  detailedOutput=F){
A=as.numeric(as.character(A))
B=as.numeric(as.character(B))
no.SNPs=dim(as.matrix(MA))[[2]]
NIS=rep(NA, no.SNPs)
no.obs=length(A) 
if (sum(!is.na(A))<3 | sum(!is.na(B))<3) stop("Input error in the function CausalVSIndependentSNPvote: there are not enough observations in trait A or B: sum(!is.na(A))<3")
if (dim(as.matrix(MA))[[1]] != no.obs) stop("Input error in function CausalVSIndependentSNPvote, the number of observations for trait A does not match the number of rows of the genetic marker data frame MA.")
if (dim(as.matrix(B))[[1]] != no.obs) stop("Input error in function CausalVSIndependentSNPvote, the number of observations for trait A does not match the number of rows of the data frame B.")

dgf1=-3+length(A) 
if (dgf1<0) stop("Warning: in function CausalVSIndependentSNPvote, not enough degrees of freedom.  dgf1=-3+sum(!is.na(M*B*A)) <0 ")
ZABgivenM= ZFisher(r=abs(partialcor(A,B,MA)) ,n=dgf1 )
ZMBgivenA= ZFisher(r=abs(partialcor(MA,B,A)) ,n=dgf1 )
NIS=ZABgivenM-ZMBgivenA
NIS
} # end of function CausalVSIndependentSNPvote






# This function computes the BLV score between a trait A and each trait in data frame B 
# MA denotes a set of candidate pleiotropic markers going into trait A.
# A can be a data frame of traits (columns correspond to traits). 
# B can be a data frame of traits (columns correspond to traits). 
# Output: a matrix of BLV scores where rows correspond to traits and A and columns correspond 
# to traits in B.
# BLV larger than 1.5 suggestion causal flow A to B.
if (exists("BLV") ) rm(BLV);
BLV=function(MA,A,B,testIndependence=T,testMarkerConsistency=T, detailedOutput=F,singleSNPanalysis=F) {
namesA=names(data.frame(A))
namesB= names(data.frame(B))
NumberTraitsA=dim(as.matrix(A))[[2]]
NumberTraitsB=dim(as.matrix(B))[[2]]
BLVoutput=matrix(NA, nrow=NumberTraitsA, ncol=NumberTraitsB) 
for (I in c(1 :NumberTraitsA) ){
AA=as.matrix(A)[,I]
no.obs=length(AA)
if (sum(!is.na(AA))<3 ) stop("Input error in the function BLV: there are not enough observations in trait A: sum(!is.na(A))<3")
if (dim(as.matrix(MA))[[1]] != no.obs) stop("Input error in function BLV, the number of observations for trait A does not match the number of rows of the genetic marker data frame MA.")
if (dim(as.matrix(B))[[1]] != no.obs) stop("Input error in function BLV, the number of observations for trait A does not match the number of rows of the data frame B.")

if ( !singleSNPanalysis ){
for (J in c(1:NumberTraitsB) ) {
BB=as.matrix(B)[,J]
corAB=cor(AA,BB,use="p")
if (testIndependence) {
removeSNPs=CausalVSIndependentSNPvote(MA=MA,A=AA,B=BB)<1;
if (sum(removeSNPs,na.rm=T)>0 ) print(paste("Warning: No. of SNPs removed:", sum(removeSNPs,na.rm=T),"Possible violations: independence model is more plausible than the causal model or the SNP assignment to the trait is not consistent."))
MA=data.frame(MA)[,!removeSNPs | is.na(removeSNPs)] 
} # end of testIndependence

if ( dim(data.frame(MA))[[2]]==0 ) BLVoutput[I,J]=NA

if ( dim(data.frame(MA))[[2]]>0 ) {
BLVvector=rep(NA, dim(as.matrix(MA))[[2]])
for (i in c(1: dim(as.matrix(MA))[[2]]) ) {
M=as.matrix(MA)[,i]
ZMB=ZFisher(r=abs(cor(M,BB,use="p")) , n=-2+sum(!is.na(M*BB))       ) 
ZMBgivenA=ZFisher(r=abs(partialcor(A=M,B=BB, C=AA) ), n=-3+sum(!is.na(M*BB))       ) 
BLVvector[i]=ZMB-ZMBgivenA

if (testMarkerConsistency) {
# now we test whether the marker relates closer to AA than to BB. By assumption, it should
if (  abs(cor(M,AA,use="p")) < abs(cor(M,BB,use="p")) ) 
{ 
if (detailedOutput) {print("Warning: a marker in MA is more closely correlated with a traits in B than with a trait in A.
The corresponding BLV score has been set to -1.")};  BLVvector[i]=-1}
}# end of if testMarkerConsistency

} # end of for loop
BLVoutput[I,J]=mean(BLVvector,na.rm=T)
} # end of if  dim(data.frame(MA))[[2]]>0 ...
} # end of loop over traits in B
} # end of loop over traits in A
BLVoutput=data.frame(BLVoutput)
dimnames(BLVoutput)[[1]]=namesA
dimnames(BLVoutput)[[2]]=namesB
} # end of if ( !singleSNPanalysis )

if ( singleSNPanalysis ){
if (NumberTraitsA>1 | NumberTraitsB>1) stop("Error in the single SNP analysis of the BLV function. 
You have either too many traits (columns) in A or too many columns in B. Only input 1 trait A and 1 trait B.
The single SNP analysis produces a BLV score for every SNP in MA.") 

no.SNPs=dim(as.matrix(MA))[[2]]
BLVoutput=rep(NA, no.SNPs)

for (i in c(1:no.SNPs) ) {
M=as.matrix(MA)[,i]
ZMB=ZFisher(r=abs(cor(M,B,use="p")) , n=-2+sum(!is.na(M*B))       ) 
ZMBgivenA=ZFisher(r=abs(partialcor(A=M,B=B, C=A) ), n=-3+sum(!is.na(M*B))       ) 
BLVoutput[i]=ZMB-ZMBgivenA

if (testMarkerConsistency) {
# now we test whether the marker relates closer to AA than to BB. By assumption, it should
if (  abs(cor(M,A,use="p")) < abs(cor(M,B,use="p")) ) 
{ if (detailedOutput) {print("Warning: a marker in MA is more closely correlated with a traits in B than with a trait in A.
The corresponding BLV score has been set to -1.")};  BLVoutput[i]=-1}
}# end of if testMarkerConsistency
} # end of for loop
} # end of if singleSNPanalysis


BLVoutput
} # end of function BLV




# this function compute the ZEO score where MA denotes 
# candidate pleiotropic markers going into trait A.
# and MB denotes candidate orthogonal causal anchors of B
# it evaluates the direction A->B between the two traits
# threshold is 2.
# A can be a data frame of traits (columns correspond to traits). 
# B can be a data frame of traits (columns correspond to traits). 
# Output: a matrix of ZEO scores where rows correspond to traits and A and columns correspond 
# to traits in B.
# We recommend to set testConfounding=T. In this case,
# the ZEO score will be set to minus .5 if BLV(MA->A->B)<1.5

# We recommend to set testIndependence=T. In this case,
# the ZEO score will be set to minus .5 if BLV(MA->A->B)>1 and BLV(MA->B->A)>1
if (exists("ZEO") ) rm(ZEO);
ZEO=function(MA,MB, A,B, testConfounding=T,testIndependence=T) {
namesA=names(data.frame(A))
namesB= names(data.frame(B))
NumberTraitsA=dim(as.matrix(A))[[2]]
NumberTraitsB=dim(as.matrix(B))[[2]]
ZEO=matrix(NA, nrow=NumberTraitsA, ncol=NumberTraitsB) 
for (I in c(1 :NumberTraitsA) ){
AA=as.matrix(A)[,I]
no.obs=length(AA)
if (dim(as.matrix(MA))[[1]] != no.obs) stop("Input error in function  ZEO, the number of observations for trait A does not match the number of rows of the genetic marker data frame MA.")
if (dim(as.matrix(B))[[1]] != no.obs) stop("Input error in function ZEO, the number of observations for trait A does not match the number of rows of the data frame B.")
noSNPsA=dim(as.matrix(MA))[[2]]
noSNPsB=dim(as.matrix(MB))[[2]]
NumberTraitsB=dim(as.matrix(B))[[2]]
for (J in c(1:NumberTraitsB) ) {
ZEOmatrix=matrix(NA, nrow= noSNPsA, ncol= noSNPsB)
for (i in c(1:noSNPsA) ) {
for (j in c(1:noSNPsB) ) {

BLVforward=BLV(MA=as.matrix(MA)[,i] ,A=AA,B=as.matrix(B)[,J])
BLVreactive=BLV(MA=as.matrix(MA)[,i] ,A=as.matrix(B)[,J], B=AA)
BLVbackward=BLV(MA=as.matrix(MB)[,j] ,A=as.matrix(B)[,J],B=AA)
ZEOmatrix[i,j]=as.numeric(BLVforward- BLVbackward)/2

if( !is.na( BLVforward) & !is.na( BLVbackward) & !is.na( BLVreactive) ){
if (testConfounding & BLVforward<1.5 & BLVbackward< -1.5) {ZEOmatrix[i,j]=-.5 }
if (testIndependence & BLVforward>1 & BLVreactive>1) { ZEOmatrix[i,j]=-.5 } # end of if 
}

}# end for for i 
} # end of for j
ZEO[I,J]=mean(ZEOmatrix, na.rm=T)
}# end of for loop over traits in B
}# end of for loop over traits in B
ZEO=data.frame(ZEO)
dimnames(ZEO)[[1]]=namesA
dimnames(ZEO)[[2]]=namesB
ZEO
} # end of function ZEO









# The following function computes LEO.NB.OCA, LEO.NB.CPA, ZEO and BLV scores for input data frames A and B
# It outputs edge scores for all edges between A and B that pass the trait.cor.threshold.
# It outputs edge scores for both possible directions: AtoB and BtoA.
# Output is a list 
# We recommend to set testConfounding=T, which sets the ZEO and the LEO.NB.OCA to -0.5 if 
# if the BLV and the CPA score indicate confounding.
# if the option correlationScores=T then the function outputs the ZEO and BLV scores as well.
if (exists("NEO.scores")) rm("NEO.scores")
NEO.scores=function(datSNP ,A ,B=NULL, detailed.output=F, trait.cor.threshold=0.1, top.N.snps.per.trait=4,testConfounding=T, correlationScores=F) {
if (is.null(B) ) B=A
datSNP=data.frame(datSNP)

names(datSNP)=paste("SNP.", names(datSNP), sep="")
no.samples=dim(as.matrix(A))[[1]]
no.traitsA=dim(as.matrix(A))[[2]]
no.traitsB=dim(as.matrix(B))[[2]]

if ( no.samples<3) stop("Input error in function LEO.scores: no.samples (no of entries) of trait A is less than 3")
if (dim(as.matrix(B))[[1]]!=no.samples) stop("Input error in function LEO.scores dim(as.matrix(B))[[1]]!=no.samples")
if (dim(as.matrix(datSNP))[[1]] !=no.samples) stop("Input error in function LEO.scores 
dim(as.matrix(datSNP))[[1]] !=no.samples")

LEO.NB.OCA.AtoB= data.frame(matrix(NA, nrow=no.traitsA,ncol=no.traitsB))
LEO.NB.CPA.AtoB = data.frame(matrix(NA, nrow=no.traitsA,ncol=no.traitsB))
ZEO.AtoB= data.frame(matrix(NA, nrow=no.traitsA,ncol=no.traitsB))
BLV.AtoB= data.frame(matrix(NA, nrow=no.traitsA,ncol=no.traitsB))

dimnames(ZEO.AtoB)[[1]]=dimnames(data.frame(A))[[2]]
dimnames(BLV.AtoB)[[1]]=dimnames(data.frame(A))[[2]]
dimnames(ZEO.AtoB)[[2]]=dimnames(data.frame(B))[[2]]
dimnames(BLV.AtoB)[[2]]=dimnames(data.frame(B))[[2]]
dimnames(LEO.NB.OCA.AtoB)[[1]]=dimnames(data.frame(A))[[2]]
dimnames(LEO.NB.CPA.AtoB)[[1]]=dimnames(data.frame(A))[[2]]
dimnames(LEO.NB.OCA.AtoB)[[2]]=dimnames(data.frame(B))[[2]]
dimnames(LEO.NB.CPA.AtoB)[[2]]=dimnames(data.frame(B))[[2]]

for ( i in c(1:no.traitsA) ) {
for ( j in c(1:no.traitsB) ) {
AA=as.matrix(A)[,i]
BB=as.matrix(B)[,j]
abs.cor= abs(cor(AA,BB, use= "p"))
if( abs.cor>.95 ) warning(paste("Traits",  names(A)[[i]], "and" , names(B)[[j]]  ,  "have absolute correlation higher than .95. No LEO.NB scores will be computed for this edge."  )   )
if (abs.cor>= trait.cor.threshold & abs.cor<=0.95)
{datCombined=data.frame(datSNP,A=AA,B=BB)
       pm=neo.get.param()
       pm$no.log=T
       pm$quiet=T
       pm$cor.dep.th=  trait.cor.threshold
       pm$top.N.snps.per.trait= top.N.snps.per.trait
       pm$do.m1m2.average = F
  NEOoutput1=neo(datCombined,pm=pm)
if (detailed.output) {print(NEOoutput1$stats[, c("LEO.NB.OCA", "LEO.NB.CPA", "ZEO", "zeo.for", "BLV.or.BilayerZscore" , "PearsonCor", "abs.pcor", "EdgeSurvives.Final.filter")])}
LEO.NB.OCA.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats["A -> B", "LEO.NB.OCA"][1]))
LEO.NB.CPA.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats["A -> B", "LEO.NB.CPA"][1]))
ZEO.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats["A -> B", "zeo.for"][1]))
BLV.AtoB[i,j]= as.numeric(as.character(NEOoutput1$stats["A -> B", "BLV.or.BilayerZscore"][1]))
} # end of if 
} #end of for loop over i
} #end of for loop over j
if (testConfounding ) {
LEO.NB.OCA.AtoB[ !is.na(LEO.NB.OCA.AtoB) &  (LEO.NB.OCA.AtoB>0)  & (LEO.NB.CPA.AtoB< -.5 | is.na(LEO.NB.CPA.AtoB))  ]=-.5
ZEO.AtoB[ !is.na(ZEO.AtoB)  & ZEO.AtoB>1 &  (BLV.AtoB<0 | is.na(BLV.AtoB))]=-.5
} # end of if testConfounding
out1= list(LEO.NB.OCA.AtoB= LEO.NB.OCA.AtoB,
LEO.NB.CPA.AtoB= LEO.NB.CPA.AtoB)
if (correlationScores) {out1= list(LEO.NB.OCA.AtoB= LEO.NB.OCA.AtoB,
LEO.NB.CPA.AtoB= LEO.NB.CPA.AtoB,
ZEO.AtoB= ZEO.AtoB,
BLV.AtoB= BLV.AtoB
)}
out1
} # end of function NEO.scores


