library(GOSim)

l1 = list()
l1$SelfSufficiencyInGrowthSignals = c('GO:0009967','GO:0030307','GO:0008284','GO:0045787','GO:0007165')
l1$InsensitivityToAntigrowthSignals = c('GO:0009968','GO:0030308','GO:0008285','GO:0045786','GO:0007165')
l1$EvadingApoptosis = c('GO:0043069','GO:0043066','GO:0045768')
l1$LimitlessReplicativePotential = c('GO:0001302','GO:0032206','GO:0090398')
l1$SustainedAngiogenesis = c('GO:0045765','GO:0045766','GO:0030949','GO:0001570')
l1$TissueInvasionAndMetastasis = c('GO:0042060','GO:0007162','GO:0033631','GO:0044331','GO:0001837','GO:0016477','GO:0048870','GO:0007155')
l1$GenomeInstabilityAndMutation = c('GO:0051276','GO:0045005','GO:0006281')
l1$TumorPromotingInflammation = c('GO:0002419','GO:0002420','GO:0002857','GO:0002842','GO:0002367','GO:0050776')
l1$ReprogrammingEnergyMetabolism = c('GO:0006096','GO:0071456')
l1$EvadingImmuneDetection = c('GO:0002837','GO:0002418','GO:0002367','GO:0050776')

d0 = c('GO:0009887','GO:0009605','GO:0048870','GO:0050896','GO:0007202','GO:0009611','GO:0006928','GO:0006950','GO:0048513','GO:0007166','GO:0014070')
d0 = c('GO:0007155','GO:0001568','GO:0001525','GO:0048514')

d1 = read.csv('../output/biclusterEnrichment_GOBP.csv',header=T,row.names=1)
l2 = list()
for(cluster in rownames(d1)) {
    l2[[cluster]] = strsplit(as.character(d1[cluster,2]),';')[[1]]
}

hallmarks = matrix(ncol=length(names(l1)), nrow=length(names(l2)), dimnames=list(names(l2), names(l1)))
for(cluster in names(l2)) {
    if (!(length(l2[[cluster]])==0)) {
        for(hallmark in names(l1)) {
            d2 = getTermSim(c(l2[[cluster]],l1[[hallmark]]),method='JiangConrath')
            hallmarks[cluster,hallmark] = max(d2[1:length(l2[[cluster]]),-(1:length(l2[[cluster]]))])
        }
    }
}

write.csv(hallmarks,'../output/jiangConrath_hallmarks.csv')

