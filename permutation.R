#!/usr/bin/env Rscript

#spike only
in_vitro_sites=c('68','69','95','211','484','501','569','613','624','655',
                 '681','708','709','723','728','729','767','985','1128','1219')
B117=c('69','70','144','501','570','681','716','982','1118')
B1351=c('18','80','215','242','243','245','246','417','484','701')
P1=c('18','20','26','138','190','417','484','501','655','1027','1176')

 comb_muts=unique(c(B117,P1,B1351))


comb_muts=sort(as.numeric(comb_muts))
Spikesites=2:1272 #ignore start & stop codons- makes it minisculely more conservative

sim_num=100000000
sim_hits=rep(0,sim_num)
for(i in 1:sim_num){
  sim_hits[i]=length(intersect(comb_muts,sample(Spikesites,20,replace=F)  ))
}
#gives P value
length(which(sim_hits>=length(intersect(in_vitro_sites,comb_muts))))/sim_num
#plots null distribution base R
hist(sim_hits[1:sim_num],main=paste(sim_num,'simulated mutational histories'),xlab='Number of VOC overlapping mutations in simulation run',ylab='Number of times x overlap observed',
     ylim=c(0,sim_num))

