# This script is for visualizing ibs dendrograms, calling clonal groups, PCoA, ADMIXTURE and Fst plotting
#Used in tandem with ANGSD files generated from DR_ANGSD.txt
library(WGCNA)
library(sparcl)
# reading list of bam files = order of samples in IBS matrix
bams=read.table("bams",header=F)[,1]
bams=sub("\\.bam","",bams,perl=T)

ma=as.matrix(read.table('Spp_out_whole.covMat'))
dimnames(ma)=list(bams,bams)

##Below is a manual creation of metadata and grouping of samples. Can also use the DR_metadata file in repo to color samples by site. 
sample=sapply(strsplit(split='_',x=colnames(ma)),FUN='[',1)
site_number=sapply(strsplit(split='_',x=colnames(ma)),FUN='[',4)
meta=as.data.frame(cbind(sample,site_number))
rownames(meta)=colnames(ma)
#View bams and assign a number 1-5 based on the site which it originated from
meta$site_number<-c(3,2,4,2,3,2,3,4,5,5,1,3,5,5,1,1,4,1,3,1,1,4,3,5) #This is only an example, order is dependent on your bams file
#Changing site numbers to names
library(dplyr)
value_map <- c("1" = "EP", "2" = "FM", "3" = "PL", "4" = "PA", "5" = "CR","6"="NA")
meta <- meta %>%
  mutate(site = value_map[as.character(site_number)])
meta$site<-as.character(meta$site)
site=meta$site
##

# setting up colors for plotting
palette(rainbow(length(unique(meta$site))))
colors=as.numeric(as.factor(meta$site))
colpops=as.numeric(as.factor(sort(unique(meta$site))))

##### IBS Dendrogram and Genotyping #####
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
hc=hclust(as.dist(ma),"ave")
quartz()
plot(hc,cex=0.5)  # this shows how similar unique MLGs are
quartz()
ColorDendrogram(hclust(as.dist(ma),"ave"), y = colors, labels = bams, branchlength=0.05)
#h= horizontal line used for calling clones.
h=0.18 # change this to be appropriate
abline(h=h,col="black",lty=3) 
quartz.save(file='Spp_dendrogramIBS.pdf',type='pdf')

#Generating genotype file
cc=cutree(hc,h=h) #change threshold
cn=cc
# save sample to genotype info
cn_df=as.data.frame(cn)
colnames(cn_df)[1]='genotype'
length(unique(cn_df$genotype)) 
write.csv(cn_df,'Spp_sample_to_genotype.csv')

#Create a .csv file with sample names (column 1), genotype (column 2), qrank (column 3), and site (column 4)
PastGeno<-read.csv("Past_genotype_meta.csv")
View(PastGeno)
geno_final=PastGeno
unique(PastGeno$genotype)
##Choosing which samples to keep from each genotype
library(dplyr)
samples2keep = PastGeno %>% group_by(genotype) %>% slice(which.max(qrank))
# write.csv based on samples2keep 
# Final bams list should have one sample per genotype before moving onto ADMIXTURE and MDS plots.

##### STOP #####
# Before moving to PCoA plots
  # For P. astreoides: Be sure that all clones have been removed, then go to step 4 in DR_ANGSD.txt (ngsLD) 
  # For O. faveolata: Go to step 4 in DR_ANGSD.txt (ngsLD)
  # For A. agaricites: Be sure that all clones have been removed.
# Need to complete steps 4 and 5 in DR_ANGSD.txt first







##### PCoA plots #####
library(vegan)
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site,conds)
# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
library(adegenet) 
cmd=pp0  # change to cmd=pp to see constrained ordination (data projection to maximize by-site separation)
plot(cmd,choices=axes2plot,display="sites",type="n") 
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

quartz.save(file='MDSplot_scaled.pdf',type='pdf')
dev.off()

quartz()
# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
quartz.save(file='MDSplot_unscaled.pdf',type='pdf')
dev.off()

# see what pca looks like
pca=prcomp(ma,scale. = TRUE)
pca_df=as.data.frame(pca$x)
pca_df$site=site
library(ggplot2)
quartz()
ggplot(pca_df,aes(x=PC1,y=PC2,color=site))+geom_point(size=2)+stat_ellipse()


##### ADMIXTURE #####
#Plotting for nsgADMIX results
setwd('/path/to/admix/folder')
# population structure with NGSadmix
# first choose a K value -- ran NGSadmix 10x for each K 
# adapted from: https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html
#read in the data
data<-list.files('/Users/sibelle/Desktop/CEE/DR_POP_GEN/Past_postngsLD/PastLD_ngsADMIX/',pattern = ".log", full.names = T)

#use lapply to read in all our log files at once
#have to change the bounds depending on number of files (data) - ex 20 runs x 10 K values = 200 bound limit
bigData<-lapply(1:200, FUN = function(i) readLines(data[i]))

# find the line that starts with "best like=" or just "b"
library(stringr)
#this will pull out the line that starts with "b" from each file and return it as a list
#have to change the bounds depending on number of files (bigData)
foundset<-sapply(1:200, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])
#now we need to pull out the first number in the string (max likelihood value), we'll do this with the function sub
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

#make a dataframe with an index corresponding to K values
logs<-data.frame(K = rep(1:10, each=20))
#add to it likelihood values
logs$like<-sapply(strsplit(split='=',as.character(foundset)), FUN='[',2)
logs$like<-as.numeric(sapply(strsplit(split=' after ',logs$like), FUN='[',1))

#### choosing a K for ngsADMIX results
## following Evanno 2005 https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2005.02553.x

library(ggplot2)
library(dplyr)
likes=as.data.frame(logs %>% group_by(K) %>% summarize(mean=mean(like),sd=sd(like)))


plot1=ggplot(likes,aes(x=K,y=mean))+geom_point()+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=1)+ylab('L(k)')+
  geom_line()

# rate of change (Lprime) is calculated as L(k)-L(k-1) aka the difference in likelihoods between Ks
Lprime=c()
K=c(1:10)
for (i in K) {
  delL=mean(logs[logs$K==i,]$like) - mean(logs[logs$K==(i-1),]$like)
  Lprime=c(Lprime,delL)
}
Lprime=as.data.frame(cbind(K,Lprime))

plot2=ggplot(Lprime,aes(x=K,y=Lprime))+geom_point()+geom_line()+ylab('L\'(k)')

# now second order rate of change Lpp) which is Lprime(L+1)-Lprime(k)
Lpp=c()
K=c(1:10)
for (i in K) {
  delLprime=abs(Lprime[Lprime$K==(i+1),]$Lprime - Lprime[Lprime$K==(i),]$Lprime)
  Lpp=c(Lpp,delLprime)
}
Lpp=as.data.frame(cbind(K,Lpp))
plot3=ggplot(Lpp,aes(x=K,y=Lpp))+geom_point()+geom_line()+ylab('L\'\'(k)')

# finally delta k which is Lpp / sd(like)
sd=tapply(logs$like, logs$K, FUN= function(x) sd(abs(x)))
delK=as.data.frame(cbind(Lpp,sd))
delK$delK=abs(delK$Lpp)/abs(delK$sd)
plot4=ggplot(delK,aes(x=K,y=delK))+geom_point()+geom_line()+ylab('delK')

quartz()
library(gridExtra)
library(grid)
title='NGSadmix Spp' # change title based on dataset
grid.arrange(plot1,plot2,plot3,plot4+ylim(0,60000),ncol=2,top=textGrob(title))
quartz.save(file='NGSadmix_Spp.pdf',type='pdf')
dev.off()

##### Plotting best K #####
setwd('/Users/sibelle/Desktop/CEE/DR_POP_GEN/Admixture_Spp/')

# assembling the input table
dir="/Users/sibelle/Desktop/CEE/DR_POP_GEN/Admixture_Spp/" # path to input files

inName="Spp_final_k4run20.qopt"

npops=4
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

# can use meta dataframe made earlier
head(meta)
#Make a data.frame with one column as sample an one column as site
inds=as.character(meta$sample)
inds==bams
#rename column 1 header to inds and column 2 to pops
inds2pop=meta[,-2]
colnames(inds2pop)=c('inds','pops')
write.table(inds2pop,quote=F,row.names=F,col.names=F,file=paste(dir,'inds2pops',sep=''))

##Can jump to here if just changing inName and npops for plot
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)
names(i2p)=c("ind","pop")
#i2p=i2p[!i2p$ind %in% c('CC7_host.bam','TX_host.bam',"SA_508_host.bam"),] # no labs or SAs
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look
tbl$pop=as.factor(tbl$pop)

#Can graph 3 Ks on same plot with par(mfrow=c(3,1))
source("/Users/sibelle/Desktop/CEE/DR_POP_GEN/plot_admixture.R")
quartz()
par(mfrow=c(3,1))
ords=plotAdmixture(data=tbl,npops=npops,angle=0,vshift=0,hshift=0)
quartz.save('Spp_final_NGSadmix_K2to4.pdf',type='pdf')

##### Plotting Fst #####
# Remove all same-site fst values form Fst.txt file before reading into R (ex remove PA_PA row)
# Be sure that Fst file is formatted correctly - header should look like this
# PA_PL	0.017051	0.027778
# PA_FM	0.002063	0.0027765
setwd('/Path/to/Fst/file/and/bams/file')
# read in Fst
Fst_df=read.delim('Ofav_FstFinal_LD.txt', header = F, sep = "\t",col.names = c('PopPair','meanFst','weightFst'))
# subset out pop against itself i.e. negative values
Fst_df=subset(Fst_df, weightFst >= 0) #Don't run this for Ofav - contains some negative Fst values we want to keep

# make pops into two columns
Fst_df$PopPair=as.character(Fst_df$PopPair)
Fst_df$pop1=sapply(strsplit(split='_',x=Fst_df$PopPair),FUN='[',1)
Fst_df$pop2=sapply(strsplit(split='_',x=Fst_df$PopPair),FUN='[',2)
Fst_df$pop1=factor(Fst_df$pop1, levels = c('PA','PL','FM','CR','EP'))
Fst_df$pop2=factor(Fst_df$pop2, levels = c('EP','CR','FM','PL','PA'))
quartz()
library(corrplot)
library(reshape2)
cor_flat=subset(Fst_df,select=c(pop1,pop2,weightFst))
corr = dcast(cor_flat,pop1 ~ pop2, value.var ='weightFst')
rownames(corr)=corr$pop1
corr$pop1=NULL
corr=as.matrix(corr)
corr=corr[,c(5:1)]
quartz()
corrplot(corr, type = 'lower', is.corr = FALSE, diag = FALSE, addCoef.col = 'black', tl.col = 'black', tl.cex = 2,col.lim = c(-0.09, 0.12))
quartz.save('Ofav_Fst.pdf',type='pdf')


#Done!