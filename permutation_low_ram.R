#!/usr/bin/env Rscript

#spike only
in_vitro_sites=c('68','69','95','211','484','501','569','613','624','655',
                 '681','708','709','723','728','729','767','985','1128','1219')
B117=c('69','70','144','501','570','681','716','982','1118')
B1351=c('18','80','215','242','243','245','246','417','484','701')
P1=c('18','20','26','138','190','417','484','501','655','1027','1176')

comb_muts=unique(c(B117,P1,B1351))


comb_muts=sort(as.numeric(comb_muts))
Spikesites=2:1273 #ignore start & stop codons- makes it minisculely more conservative

sim_num=100000000
sim_hits=rep(0,length(Spikesites))
for(i in 1:sim_num){
  overlap=length(intersect(comb_muts,sample(Spikesites,length(in_vitro_sites),replace=F)))
  sim_hits[overlap+1]=sim_hits[overlap+1]+1
  if(!as.logical(i%%(sim_num/10))){
    print(paste(100*i/sim_num,'% complete'))
  }
}

#gives P value
sum(sim_hits[(length(intersect(in_vitro_sites,comb_muts))+1):length(sim_hits)])/sim_num
#plots null distribution base R
plot(1:length(sim_hits)-1,sim_hits,xlim=c(0,head(which(sim_hits==0),1)))

#to save for later
#write.csv(cbind(1:length(sim_hits)-1,sim_hits),file='~/sim_hits.table.csv',row.names = F)

#plots null distribution ggplot2
library(ggplot2);library(RColorBrewer)

sim_hits=data.frame(overlap=1:length(sim_hits)-1,sim_hits=sim_hits)

print(head(sim_hits,which(sim_hits[,2]==0)[1]))

ggplot(sim_hits,aes(x=overlap,y=sim_hits*10)) + 
  geom_col(fill=brewer.pal(3,'Dark2')[1]) +
  scale_y_log10(breaks=c(10^(0:round(log10(sim_num*10)))),
                labels=as.character(c(0,10^(0:round(log10(sim_num))))),##slightly complicated work around because ggplot hates x^0 being >0
                limits=c(1,sim_num*10))+
  theme_minimal()+xlab('Number of VOC overlapping mutations in simulation run')+
  ylab('Number of times x overlap observed (log10+1)')+
  scale_x_continuous(breaks=0:(head(which(sim_hits[,2]==0),1)),
                     limits=c(-0.5,head(which(sim_hits[,2]==0),1)-1.5))+
  ggtitle(paste(sim_num,'simulated mutational histories'))


