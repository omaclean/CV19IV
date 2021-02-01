#!/usr/bin/env Rscript

in_vitro_sites=c('68','69','95','211','484','501','569','613','624','655',
                 '681','708','709','723','728','729','767','985','1128','1219')
B117=c('69','70','144','501','570','681','716','982','1118')
P1=c('18','20','26','138','190','417','484','501','655','1027')
B1351=c('18','80','215','242','243','245','246','417','484','701')

# download table from http://cov-glue.cvr.gla.ac.uk/#/replacement 
# in this case filter was set to 'Virus genome region' matches 'S' 
#later data pulls may produce different results as GLUE pulls new gisaid data
dat=read.csv('~/Downloads/replacements.csv')

#########set thresholds etc.
minimum_seq_num=5
maxlength=1273 #change for use on other genes
dat=dat[dat$numSeqs>=minimum_seq_num,]
interval=20

################ calculate diversity
hits=table(dat$codonNumber)
div=rep(NA,(maxlength-interval+1))
for(i in 1:(maxlength-interval+1)){
  div[i]=sum(hits[names(hits)%in%(i:(i+interval-1))])/interval
}

############################### plot diversity

#libraries for colour scale label & palette
library(DescTools);library(RColorBrewer) 

plot(in_vitro_sites,rep(-0.05,length(in_vitro_sites)),col=alpha('#000000',.5),lwd=2,pch='|',
     ylim=c(0,max(div)+0.1),xlim=c(0,maxlength),ylab='',xlab='',cex=2.2,axes=F)
par(new=T)
#generate continuous colour palette for bars
pal=colorRampPalette(c("blue","#02CAE1","#FA7D14","red"))

pal2=pal(round(max(div)*1000))
plot(1:(maxlength-interval+1)+interval/2,div,
     ylim=c(0,max(div)+0.1),type="h",xlim=c(0,maxlength),xlab="central spike amino acid number in sliding window",
     ylab="average amino acid replacements observed per site",
     col=pal2[round(div*1000)],
     main=paste('Diversity in ',interval,' amino acid sliding window \n(calculated from ',nrow(dat),' mutations seen in minimum of ',minimum_seq_num,' sequences)',sep=''),axes=F)

ColorLegend(x=1273+5,y=max(div)+0.1,cols=pal2,labels=c(0:ceiling(max(div))),width = 10,cex=1)


par(new=T)
plot(c(B117,P1,B1351),rep(-0.07,length(c(B117,P1,B1351))),col=alpha(brewer.pal(3,'Dark2')[3],.4),lwd=2,pch=2,
     ylim=c(0,max(div)+0.1),xlim=c(0,maxlength),ylab='',xlab='',cex=1.2,axes=F)

axis(1,labels=c(0,1273),at=c(0,1273))
axis(2,labels=0:ceiling(max(div)),at=0:ceiling(max(div)))
legend(legend='VOC defining mutations',x=900,y=2.5,bty='n',col=alpha(brewer.pal(3,'Dark2')[3],.9),pch=2,cex=1.2)
legend(legend='in vitro mutations',x=900,y=2.2,bty='n',col=alpha('#000000',.9),pch='|',cex=1.2)

