#############################################################################
############################# Figure 2 and S2 ###############################
#############################################################################

############################### Figure 2a ###################################
lof_com<-data.frame(catss=c("LoFs","LoFs","LoFs","LoFs","LoFs"),
                   cats=c("stop_gained","stop_lost","start_lost","splice_acceptor","splice_donor"),
                   nums=c(8434,19,2310,5776,10628))

lof_com$cats<-factor(lof_com$cats,levels=c("start_lost","stop_lost","stop_gained",
                                           "splice_donor","splice_acceptor"))

p<-ggplot(lof_com, aes(y =catss,x = nums,fill =cats ))+
  geom_bar(stat ="identity",width=0.5,position ="stack")+  
  labs(x = "Derived allel frequency(DAF)",y = "Propotion")+                        
  ##guides(fill = guide_legend(reverse = F))+ 
  ##scale_x_continuous(breaks=frq_long$freq, labels = frq_long$freq)+
  scale_fill_manual(values = c("#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF"))+
  theme_classic()+
  theme(legend.title = element_blank(),             
        legend.text = element_text(size = 8,color = "black"),       
        legend.position = "top",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

ggsave("los_composition.pdf", plot=p,width =12, height = 4, units = "cm")
ggsave("los_composition.png", plot=p,width =12, height = 4, units = "cm",dpi = 680)
######################## Figure S2a ##########################################
setwd("D:/vcf_for_gte/losfunction/16_coding_snp")
library(ggplot2)
library(ggsci)
library(viridis)
library(MASS)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

coding_cats<-c("stop_gained","stop_lost","start_lost","missense","synonymous","stop_retained",
               "initiator_codon","splice_acceptor","splice_donor")


pcadd_st<-data.frame(ID=NA,p_mean=NA,p_se=NA,nums=NA,sums=NA)

for(i in coding_cats){
  tmp<-read.csv(paste(i,".pcadd.txt",sep=""),header = F,sep="\t")
  tmp$density <- get_density(tmp$V3, tmp$V2, h = c(1, 1), n = 100)
  tmp_sum<-sum(tmp$V2)
  tmp_mean<-mean(tmp$V2)
  tmp_se<-sd(tmp$V2)/sqrt(nrow(tmp))
  tmp_all<-c(i,tmp_mean,tmp_se,nrow(tmp),tmp_sum)
  pcadd_st<-rbind(pcadd_st,tmp_all)
}

pcadd_st<-na.omit(pcadd_st)
pcadd_st$p_mean<-as.numeric(pcadd_st$p_mean)
pcadd_st$p_se<-as.numeric(pcadd_st$p_se)
pcadd_st$nums<-as.numeric(pcadd_st$nums)
pcadd_st$sums<-as.numeric(pcadd_st$sums)

pcadd_st$ID<-factor(pcadd_st$ID,levels=c("synonymous","missense","start_lost","stop_lost",
                                         "stop_gained","splice_donor","splice_acceptor",
                                         "stop_retained","initiator_codon"))
pcadd_st<-pcadd_st[order(pcadd_st$ID),]
lof_pcadd<-sum(pcadd_st$sums[5:9])/sum(pcadd_st$nums[5:9])

p<-ggplot(pcadd_st,aes(ID, log10(nums),fill=ID))+
  geom_bar(stat = "identity",width=0.4)+
  scale_fill_manual(values = c("#00468BFF","#0099B4FF","#925E9FFF","#FDAF91FF",
                               "#AD002AFF","#ADB6B6FF","#1B1919FF","#42B540FF",
                               "#ED0000FF"))+
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()+labs(x="",y="log10(Numbers)")+
  ###scale_y_reverse(expand=c(0,0))+
  theme(legend.position = "", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        
        ##axis.text.x=element_text(size = 12,color = "black",angle = 30),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.title.y=element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x= element_blank()
  )
ggsave("coding_all.bar.pcadd.png",plot=p,width = 6,height =2,dpi=680)
ggsave("coding_all.bar.pcadd.pdf",plot=p,width = 6,height =2) 

p1<-ggplot(pcadd_st, aes(ID, p_mean,color=ID))+ geom_point(size=2)+
  scale_color_manual(values = c("#00468BFF","#0099B4FF","#925E9FFF","#FDAF91FF",
                                "#AD002AFF","#ADB6B6FF","#1B1919FF","#42B540FF",
                                "#ED0000FF"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), 
                width=0.4,position = position_dodge(.9)) + 
  scale_x_discrete(labels=c("Synonymous","Missense","Start_lost","Stop_lost","Stop_gained",
                            "Splice_donor","Splice_acceptor",
                            "Stop_retained","Initiator_codon"))+
  labs(y="pCADD scores",x="") + 
  theme_classic()+ylim(0,18)+
  theme(legend.position = "",
        axis.text.x = element_text(size = 14,color = "black",angle=30,vjust=0.7,hjust=0.5),
        ##axis.ticks.x = element_blank(),
        ##axis.line.y  = element_blank(),
        axis.text.y = element_text(size = 14,color = "black"),
        axis.title.y= element_text(size = 14,colour="black")
  )
ggsave("coding_all.pcadd.png",width = 6,height =3,dpi=680)
ggsave("coding_all.pcadd.pdf",width = 6,height =3)

########################### Figure S2b #####################################
library(RIdeogram)
library(circlize)
library(rsvg)
setwd("D:/vcf_for_gte/losfunction/RIdeogram")
pig_karyotype <- read.table("pig_karyotype", sep = "\t", header = T, stringsAsFactors = F)
snp_info <-read.table("vep_snpeff_rm_coovar_lof_plot.txt", sep = "\t", header = T, stringsAsFactors = F)
snp_info<-subset(snp_info,,select=c("Chr","Start","End"))
snp_density <- genomicDensity(snp_info,window.size = 1e6,count_by = "number",overlap = F)

colnames(snp_density)<-c("Chr","Start","End","Value")
ideogram(karyotype = pig_karyotype,overlaid = snp_density,width=170)
convertSVG("chromosome.svg",file="final_lof.pdf",device = "pdf",width = 16,
           height = 12)
rsvg_pdf("chromosome.svg", "final_lof.pdf",width = 12,height = 16)
snp_density <-snp_density[order(-snp_density$value),]
snp_density_top<-snp_density[1:round(2252*0.01),]
write.csv(snp_density_top,"snp_density_top1.csv",quote = F)
########################### Figure S2c #####################################
setwd("D:/vcf_for_gte/losfunction/01_pca/lof")
nums1<-length(group.info$GROUP[group.info$GROUP=="ASD"])
nums2<-length(group.info$GROUP[group.info$GROUP=="ASW"])
nums3<-length(group.info$GROUP[group.info$GROUP=="EUD"])
nums4<-length(group.info$GROUP[group.info$GROUP=="EUW"])
nums5<-length(group.info$GROUP[group.info$GROUP=="SUI"])
res<-read.csv(paste("tsne_umap_res/GTE.1602.pcs","6",".","tsne",".csv",sep=""),header=T,sep = ",")
colnames(res)<-c("ID","PC1","PC2")
res_all<-merge(group.info,res,by="ID")
group_new<-read.csv("C:/vcf_for_gte/01_pca/tsne_umap_res/new_pop.list",header = F,sep = "\t")
res_all$GROUP1<-"N"
for (i in group_new$V1){
  res_all[res_all$ID == i,]$GROUP1<-res_all[res_all$ID == i,]$GROUP
}

ggplot(res_all, aes(x=PC1, y=PC2,color=GROUP1)) +
  geom_point(size=1.5)+ theme_bw()+xlab("t-SNE-1") +ylab("t-SNE-2")+
  scale_color_manual(values = c("#FE7D7D","#FEC655","#75C4FE","#56B956","#ADB6B6FF"))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(0.5,1.0),legend.direction = "horizontal",
        legend.key.height=unit(0.3,"line"),
        legend.key.width=unit(0.3,"line"),
        axis.text.x = element_text(size = 18,color = "black"),
        axis.text.y = element_text(size = 18,color = "black"),
        axis.title = element_text(size = 18,color = "black"),
        legend.text = element_text(size = 13,color = "black",hjust=0))
ggsave("C:/vcf_for_gte/losfunction/01_pca/new_old_tsne_6.png",width = 16, units = "cm", height = 12, dpi = 680)
ggsave("C:/vcf_for_gte/losfunction/01_pca/new_old_tsne_6.pdf",width =16, units = "cm", height = 12)

######################### Figure 2b #########################################
library(tidyr)
setwd("D:/vcf_for_gte/losfunction/16_coding_snp/cat_frq")
lofs<-c("stop_gained","stop_lost","start_lost","splice_acceptor","splice_donor")
frqs<-c(0.001,0.005,0.01,0.05,0.1,0.3,0.5,0.7,0.9,1.0)
lofs_st<-data.frame(popss=NA,cats=NA,nums=NA)
for(lof in lofs){
  tmp<-read.csv(paste(lof,".final_lof.frq",sep=""),header = F,sep="\t")
  colnames(tmp)<-c("CHR","POS","N","N_CHR","A1","A2")
  tmp<-tmp[2:nrow(tmp),]
  tmp1<-separate(tmp,A2,into=c("DA","DAF"),sep=":")
  tmp1$DAF<-as.numeric(tmp1$DAF)
  tmp1<-na.omit(tmp1)
  tmp2<-tmp1[tmp1$DAF>0.0,]
  frqtmp<-""
  for(frq in frqs){
    if(frqtmp == ""){
      nums<-nrow(tmp2[tmp2$DAF<=frq,])/nrow(tmp2)
      frq1<-c(lof,paste("(0.000",",",frq,"]",sep=""),nums)
      lofs_st<-rbind(lofs_st,frq1)
      frqtmp<-frq
    }else{
      nums<-nrow(tmp2[tmp2$DAF>frqtmp&tmp2$DAF<=frq,])/nrow(tmp2)
      frq1<-c(lof,paste("(",frqtmp,",",frq,"]",sep=""),nums)
      lofs_st<-rbind(lofs_st,frq1)
      frqtmp<-frq
    }
  }
}

for(i in c("synonymous","missense")){
  tmp<-read.csv(paste("D:/vcf_for_gte/losfunction/16_coding_snp/",i,".pcadd.txt",sep=""),header = F,sep="\t")
  colnames(tmp)<-c("ID","pcadd","frq","cats")
  tmp<-tmp[tmp$frq>0.0,]
  frqtmp<-""
  for(frq in frqs){
    if(frqtmp == ""){
      nums<-nrow(tmp[tmp$frq<=frq,])/nrow(tmp)
      frq1<-c(i,paste("(0.000",",",frq,"]",sep=""),nums)
      lofs_st<-rbind(lofs_st,frq1)
      frqtmp<-frq
    }else{
      nums<-nrow(tmp[tmp$frq>frqtmp&tmp$frq<=frq,])/nrow(tmp)
      frq1<-c(i,paste("(",frqtmp,",",frq,"]",sep=""),nums)
      lofs_st<-rbind(lofs_st,frq1)
      frqtmp<-frq
    }
  }
}

lofs_st<-na.omit(lofs_st)
lofs_st$nums<-as.numeric(lofs_st$nums)
lofs_st$popss<-factor(lofs_st$popss,levels = c("synonymous","missense",
                                               "start_lost","stop_lost","stop_gained","splice_donor","splice_acceptor"))
lofs_st$cats<-factor(lofs_st$cats,levels = c("(0.000,0.001]","(0.001,0.005]","(0.005,0.01]","(0.01,0.05]",
                                             "(0.05,0.1]","(0.1,0.3]","(0.3,0.5]","(0.5,0.7]","(0.7,0.9]","(0.9,1]"))
p<-ggplot(lofs_st, aes(y =nums,x = cats,fill =popss ))+
  geom_bar(stat ="identity",width=0.5,position = "dodge")+  
  labs(x = "Derived Allele Frequency(DAF)",y = "Propotion")+                        
  ##guides(fill = guide_legend(reverse = F))+ 
  scale_fill_manual(values =c("#00468BFF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF",
                              "#ADB6B6FF","#1B1919FF"))+
  theme_classic()+
  theme(legend.title = element_blank(),             
        legend.text = element_text(size = 8,color = "black"),       
        legend.position = "right",
        axis.text.x = element_text(size=14,color = "black",angle = 30,vjust=0.65),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size = 14,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"))
ggsave("freq_los_in_5.pdf", plot=p,width =14, height = 8, units = "cm")
ggsave("freq_los_in_5.png", plot=p,width =14, height = 8, units = "cm",dpi = 680)
############################### Figure S2d ##################################
library(ggplot2)
library(tidyverse)
setwd("D:/vcf_for_gte/losfunction/06_pop_frq")
pops<-c("ASIA_D","ASIA_W","EURO_D","EURO_W")
frqs_st<-data.frame(popss=NA,cats=NA,nums=NA)
for(pop in pops){
  tmp<-read.csv(paste("final_lof",pop,"frq",sep = "."),header = F,sep = "\t")
  colnames(tmp)<-c("CHR","POS","N","N_CHR","A1","A2")
  tmp<-tmp[2:nrow(tmp),]
  tmp1<-separate(tmp,A2,into=c("DA","DAF"),sep=":")
  tmp1$DAF<-as.numeric(tmp1$DAF)
  tmp1<-na.omit(tmp1)
  tmp2<-tmp1[tmp1$DAF>0.0,]
  frqs<-c(0.001,0.005,0.01,0.05,0.1,0.3,0.5,0.7,0.9,1.0)
  frqtmp<-""
  for(frq in frqs){
    if(frqtmp == ""){
      nums<-nrow(tmp2[tmp2$DAF<=frq,])/nrow(tmp2)
      frq1<-c(pop,paste("(0.000",",",frq,"]",sep=""),nums)
      frqs_st<-rbind(frqs_st,frq1)
      frqtmp<-frq
    }else{
      nums<-nrow(tmp2[tmp2$DAF>frqtmp&tmp2$DAF<=frq,])/nrow(tmp2)
      frq1<-c(pop,paste("(",frqtmp,",",frq,"]",sep=""),nums)
      frqs_st<-rbind(frqs_st,frq1)
      frqtmp<-frq
    }
  }
}

frqs_st<-na.omit(frqs_st)
frqs_st$nums<-as.numeric(frqs_st$nums)
frqs_st$cats<-factor(frqs_st$cats,levels = c("(0.000,0.001]","(0.001,0.005]","(0.005,0.01]","(0.01,0.05]",
                                             "(0.05,0.1]","(0.1,0.3]","(0.3,0.5]","(0.5,0.7]","(0.7,0.9]","(0.9,1]"))
p<-ggplot(frqs_st, aes(y = nums,x = cats,fill =popss ))+
  geom_bar(stat ="identity",position = "dodge",width=0.5)+  
  labs(x = "Derived allel frequency(DAF)",y = "Propotion")+                        
  scale_fill_manual(values = c("#FE7D7D","#FEC655","#75C4FE","#56B956"))+
  theme_classic()+
  theme(legend.title = element_blank(),             
        legend.text = element_text(size = 8,color = "black"),       
        legend.position = "top",
        axis.text.x = element_text(size=14,color = "black",angle = 30,vjust=0.65),
        axis.text.y = element_text(size=14,color = "black"),
        axis.title.x = element_text(size = 14,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"))

p
ggsave("freq_los_in_pops.pdf", plot=p,width =12, height = 8, units = "cm")
ggsave("freq_los_in_pops.png", plot=p,width =12, height = 8, units = "cm",dpi = 680)
############################### Figure 2c ###################################
setwd("D:/vcf_for_gte/losfunction/04_roh/pcadd/")
library(ggpubr)
i="plink"
los_pcadd<-read.csv(paste("GTE.pops.roh.los_all.",i,".static",sep=""),header = F,sep=" ")
los_pcadd<-los_pcadd[los_pcadd$V8 !=0,]

pops_4<-c("ASIA_W","ASIA_D","EURO_D","EURO_W")
los_pcadd_4<-los_pcadd[los_pcadd$V1 %in% pops_4,]

los_pcadd_4$inroh<-los_pcadd_4$V4/los_pcadd_4$V8
los_pcadd_4$notinroh<-(los_pcadd_4$V3-los_pcadd_4$V4)/(los_pcadd_4$V7-los_pcadd_4$V8)
los_pcadd_4_l<-melt(los_pcadd_4,id.vars = "V1",measure.vars = c("inroh","notinroh"),
                    variable.name = "los_pcadd",value.name = "pcadd_val")

ylab<-bquote(Nums[LoFs]/Nums[SNPs])
ylim1 = boxplot.stats(los_pcadd_4_l$pcadd_val)$stats[c(1, 5)]
p<-ggboxplot(los_pcadd_4_l,x="V1",y="pcadd_val",color = "los_pcadd",outlier.shape = NA)+
  stat_compare_means(aes(group=los_pcadd), label = "p.signif",
                     method = "t.test",label.y = ylim1[2]*1.9)+
  labs(x="",y=ylab)+
  scale_fill_lancet(guide="none")+
  coord_cartesian(ylim = c(0,ylim1[2]*2))+
  scale_x_discrete(labels=c("ASD","ASW","EUD","EUW"))+
  theme_classic()+
  scale_color_discrete("Regions",breaks = c("inroh","notinroh"),
                       labels=c("ROH","None ROH"))+
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 8,color = "black"),
  )
ggsave("roh_los_in_pop4_plink_box_all.png",plot=p,width =6, height =8, units = "cm",dpi = 680)
ggsave("roh_los_in_pop4_plink_box_all.pdf",plot=p,width =6, height =8, units = "cm")

los_pcadd_4$homoinroh<-los_pcadd_4$V6/los_pcadd_4$V4
los_pcadd_4$homonotinroh<-(los_pcadd_4$V5-los_pcadd_4$V6)/(los_pcadd_4$V3-los_pcadd_4$V4)
los_pcadd_4_homo<-melt(los_pcadd_4,id.vars = "V1",measure.vars = c("homoinroh","homonotinroh"),
                       variable.name = "los_pcadd",value.name = "pcadd_val")

ylim1 = boxplot.stats(los_pcadd_4_homo$pcadd_val)$stats[c(1, 5)]

ylab<-bquote(Nums[LoFs(DAF==1)]/Nums[LoFs(DAF>0)])

p<-ggboxplot(los_pcadd_4_homo,x="V1",y="pcadd_val",color = "los_pcadd",outlier.shape = NA)+
  stat_compare_means(aes(group=los_pcadd),method = "t.test",label = "p.signif",label.y = 0.98)+
  labs(x="",y=ylab)+
  ##ylim(0,50)+
  scale_fill_lancet(guide="none")+
  coord_cartesian(ylim = c(ylim1[1]*0.25,1))+
  ##scale_fill_manual(values = c("#FE7D7D","#FEC655","#75C4FE","#56B956"),
  ##labels = c("ASD","ASW","EUD","EUW"),guide="none")+
  scale_x_discrete(labels=c("ASD","ASW","EUD","EUW"))+
  theme_classic()+
  scale_color_discrete("Regions",breaks = c("homoinroh","homonotinroh"),
                       labels=c("ROH","None ROH"))+
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 8,color = "black"),
  )
ggsave("roh_homolos_in_pop4_plink_box_all.png", plot=p,width =6, height =8, units = "cm",dpi = 680)
ggsave("roh_homolos_in_pop4_plink_box_all.pdf", plot=p,width =6, height =8, units = "cm")

############################### Figure 2d ###################################
eqtl_egene_0.05<-c("Blastomere","Cartilage","Duodenum","Kidney","Lymph_node","Morula","Oocyte","Pituitary","Placenta")
sqtl_egene_0.05<-c("Adipose","Blood","Brain","Embryo","Hypothalamus","Ileum","Jejunum","Large_intestine",
                   "Liver","Lung","Muscle","Ovary","Small_intestine","Spleen","Testis","Uterus")
setwd("D:/vcf_for_gte/losfunction/08_los")
library(ggplot2)
library(reshape2) 
library(ggsci)
tmp<-read.csv("QTL.all.sig.txt",header = F,sep = "\t")
tmp$all<-(tmp$V6/tmp$V5)/(tmp$V4/tmp$V3)
tmp$splice<-(tmp$V8/tmp$V7)/(tmp$V4/tmp$V3)
tmp$codes<-((tmp$V6-tmp$V8)/(tmp$V5-tmp$V7))/(tmp$V4/tmp$V3)
tmp_e<-subset(tmp,V1=="eQTL",select = c("V1","V2","all","splice","codes"))
tmp_s<-subset(tmp,V1=="sQTL",select = c("V1","V2","all","splice","codes"))
tmp_e<-tmp_e[!(tmp_e$V2 %in% eqtl_egene_0.05),]
tmp_s<-tmp_s[(tmp_s$V2 %in% sqtl_egene_0.05),]
tmp_e_l<-melt(tmp_e,id.vars = c("V1","V2"),
              measure.vars = c("all","splice","codes"),
              variable.name = "catagory",
              value.name = "enrichment")
tmp_s_l<-melt(tmp_s,id.vars = c("V1","V2"),
              measure.vars = c("all","splice","codes"),
              variable.name = "catagory",
              value.name = "enrichment")
tmp_l<-rbind(tmp_e_l,tmp_s_l)

ggplot(tmp_l, aes(x = V1, y = enrichment)) +
  geom_violin(aes(fill = catagory), trim = FALSE)+
  scale_fill_manual(values = c("grey","79cc3dff","#00468BFF"))+
  theme_classic()+
  labs(x="",y="")+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 17),
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("LoFs.enrichment.esQTL.png",width = 5.5,height =3.5,dpi=680)
ggsave("LoFs.enrichment.esQTL.pdf",width = 5.5,height =3.5) 
############################### Figure S2e ##################################
setwd("D:/vcf_for_gte/losfunction/04_roh/pcadd/")
library(ggpubr)
i="plink"
los_pcadd<-read.csv(paste("GTE.pops.roh.los_all.",i,".static",sep=""),header = F,sep=" ")
los_pcadd<-los_pcadd[los_pcadd$V8 !=0,]
pops_4<-c("ASIA_W","ASIA_D","EURO_D","EURO_W")
los_pcadd_4<-los_pcadd[los_pcadd$V1 %in% pops_4,]
los_pcadd_4$froh<-los_pcadd_4$V9/2265774640
ggplot(los_pcadd_4,aes(x=V1,y=froh,fill=V1))+geom_boxplot()+
  scale_fill_manual(values = c("#FE7D7D","#FEC655","#75C4FE","#56B956"))+
  scale_x_discrete(labels = c("ASD","ASW","EUD","EUW"))+ 
  labs(x="",y=expression(bolditalic("F"[ROH])))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 17,color = "black"),
        axis.text.y = element_text(size = 17,color = "black"),
        axis.title.y = element_text(size = 17,color = "black"),
        legend.position = "none")
ggsave("froh_in_pop4.png", width =8, height =8, units = "cm",dpi = 680)
ggsave("froh_in_pop4.pdf", width =8, height =8, units = "cm")

########################### Figure S2f ######################################
library(ggplot2)
setwd("D:\\vcf_for_gte\\losfunction\\08_los\\qtl")
tmp<-read.csv("LoFs_nums_egene_in_QTL.csv",header = T)

Add_R <- function(dataframe,x,y,factor){
  cor <- data.frame()
  dataframe[,factor] <- as.factor(dataframe[,factor])
  lev <- levels(dataframe[,factor])
  for (i in c(1:length(lev))) {
    name <- lev[i]
    data <- dataframe[which(dataframe[,factor] == name),]
    lm <- summary(lm(data,formula = data[,y]~data[,x]))
    r_squared <- round(lm$r.squared,2)
    inter <- round(lm$coefficients[1,1],2)
    coefficients <- round(lm$coefficients[1,2],2)
    max_x <- max(data[,x])
    max_y <- max(data[,y])
    if(inter>0){
      eq <- substitute(""~R^2~"="~a~","~hat(y)~" = "~b%.%x+c~ "",list(a = r_squared,b = coefficients,c = inter))
    }else{
      inter <- abs(inter)
      eq <- substitute(""~R^2~"="~a~","~hat(y)~" = "~b%.%x-c~"",list(a = r_squared,b = coefficients,c = inter))
    }
    cor <- rbind(cor,cbind(rsqua = r_squared,coef = coefficients,intercept = inter,max_x = max_x,max_y = max_y,exp = ""))
    exp <- as.character(as.expression(eq))
    cor$exp[i] <- exp
    row.names(cor)[i] <- name
  }
  for (i in  c(1:5)){
    cor[,i] <- as.numeric(cor[,i])
  }
  return(cor)
}

df <- Add_R(tmp,"percentage_lofs","percentage_egen","catagory")
tmp$egene_sample<-tmp$percentage_egene/tmp$sample_nums
tmp<-tmp[order(tmp$egene_sample,decreasing = T),]
tmp1<-tmp[tmp$catagory=="eQTL",]

ggplot(tmp,aes(percentage_lofs,percentage_egene))+
  geom_point(size = 3,aes(color = catagory,shape = catagory,fill = catagory))+
  geom_smooth(aes(color = catagory,fill =catagory),method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+
  scale_color_manual(values = c("slateblue2","tomato2"))+
  xlab("Sig LoFs pecentage")+ylab("Sig egene percentage")+
  ##geom_text(data = df,aes(max_x,max_y,label = exp),vjust = 3,hjust =0.8,size = 3.5,parse = T,color=c("slateblue2","tomato2"))+
  ##coord_cartesian(xlim = c(4,9),expand = F,ylim = c(2,4.7))+
  theme_classic(base_size = 28)+
  theme(panel.background = element_blank(),
        legend.position = "top",
        axis.line = element_line(colour = "black",size = 0.5),
        panel.grid.major.y = element_line(colour = "grey",linetype = 2,size=0.5))
ggsave("cor_LoFs_egene_nums.png",width = 12,height =8,dpi=680)
ggsave("cor_LoFs_egene_nums.pdf",width = 12,height =8) 

############################ Figure 2e ######################################
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggsci)
library(patchwork)
setwd("D:\\vcf_for_gte\\losfunction\\08_los\\qtl")
egene_file<-read.csv("egene.enrich.gene.34_1.txt",header = F,sep=" ")
egene_file$foldchang<-(egene_file$V7/egene_file$V6)/(egene_file$V3/egene_file$V2)
egene_file$foldchang1<-(egene_file$V8/egene_file$V7)/(egene_file$V4/egene_file$V3)
egene_file$foldchang2<-(egene_file$V9/egene_file$V7)/(egene_file$V5/egene_file$V3)
egene_file$pval<-""
egene_file$pval1<-""
egene_file$pval2<-""
egene_file$ave_len<-egene_file$V10/egene_file$V6
for(i in 1:nrow(egene_file)){
  egene_file[i,]$pval<-fisher.test(matrix(c(egene_file[i,]$V7,egene_file[i,]$V3,egene_file[i,]$V6,egene_file[i,]$V2),ncol = 2,nrow = 2),
                                   alternative = "two.sided")$p.value
  egene_file[i,]$pval1<-fisher.test(matrix(c(egene_file[i,]$V8,egene_file[i,]$V4,egene_file[i,]$V7,egene_file[i,]$V3),ncol = 2,nrow = 2),
                                    alternative = "two.sided")$p.value
  egene_file[i,]$pval2<-fisher.test(matrix(c(egene_file[i,]$V9,egene_file[i,]$V5,egene_file[i,]$V7,egene_file[i,]$V3),ncol = 2,nrow = 2),
                                    alternative = "two.sided")$p.value
}
egene_file$fdr<-p.adjust(egene_file$pval,method = "fdr")
egene_file$fdr1<-p.adjust(egene_file$pval1,method = "fdr")
egene_file$fdr2<-p.adjust(egene_file$pval2,method = "fdr")

egene_file<-subset(egene_file,V1 != 0,)
a<-cor.test(egene_file$foldchang,egene_file$V1)
p_txt<-paste("R=",round(a$estimate[["cor"]],2),",p=",round(a$p.value,8),sep = "")
egene_file$foldchang_p<-egene_file$foldchang/egene_file$ave_len
b<-cor.test(egene_file$foldchang_p,egene_file$V1)
p_txt<-paste("R=",round(b$estimate[["cor"]],2),",p=",round(b$p.value,5),sep = "")

p<-ggplot(egene_file, aes(V1,foldchang))+
  geom_point(size=6,aes(alpha=0.9))+
  geom_smooth(method = lm)+
  ##geom_text(aes(4,2.093188e-05),label=p_txt,color="black",size=6)+
  labs(y="Enrichment",x="Nums of tissues shared by egene") +
  theme_classic()+ scale_x_continuous(breaks=c(seq(1,34,3)),labels = c(seq(1,34,3)))+
  geom_hline(yintercept = 1.0,linetype=2,colour="red")+
  theme_classic()+
  theme(##panel.background = element_blank(),
    legend.position = "",
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black")) 
ggsave("enrich_by_tissues_egenes_all.png",width = 6,height =5,dpi=680)
ggsave("enrich_by_tissues_egenes_all.pdf",width = 6,height =5) 

######################## Figure 2g ##########################################
tmp2<-subset(egene_file,,select = c(V1,foldchang1))
tmp2$cats<-"egene overlapped with one LoF"
colnames(tmp2)[2]<-"foldchang2"
tmp3<-subset(egene_file,,select = c(V1,foldchang2))
tmp3$cats<-"egene overlapped with more than one LoF"
tmp4<-rbind(tmp2,tmp3)
cor.test(tmp3$V1,tmp3$foldchang2)
cor.test(tmp2$V1,tmp2$foldchang2)
p1<-ggplot(tmp4, aes(V1,foldchang2))+
  geom_point(size = 3,aes(color = cats,shape = cats,fill = cats))+
  geom_smooth(aes(color = cats,fill =cats),method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+
  scale_color_manual(values = c("slateblue2","tomato2"))+
  ##geom_text(aes(4,2.093188e-05),label=p_txt,color="black",size=6)+
  labs(y="Enrichment",x="Nums of tissues shared by egene") + 
  theme_classic()+ scale_x_continuous(breaks=c(seq(1,34,3)),labels = c(seq(1,34,3)))+
  geom_hline(yintercept = 1.0,linetype=2,colour="red")+
  theme(##panel.background = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black"),
  )
ggsave("enrich_by_tissues_egenes_all_1_2.png",width = 6,height =5,dpi=680)
ggsave("enrich_by_tissues_egenes_all_1_2.pdf",width = 6,height =5) 
############################# Figure 2f #####################################
setwd("D:/vcf_for_gte/losfunction/08_los/gwas")
library(ggplot2)
tmp<-read.csv("LoFs.13.enrichment",header = F,sep = "\t")
tmp$V1<-gsub("M_","",tmp$V1)
tmp$V2<-as.numeric(gsub("<","",tmp$V2))
tmp$V5<-(2^tmp$V2)
tmp$V6<-(2^tmp$V3)
tmp$V7<-(2^tmp$V4)

ggplot(tmp,aes(x=V6, y=V1))+
  geom_point(size=1.5)+
  geom_segment(aes(x=V5,xend=V7,y=V1,yend=V1),cex=0.2)+
  geom_vline(xintercept=2,color="red",linetype="dashed")+
  scale_y_discrete(position = 'right')+
  labs(y=NULL)+theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("LoFs.enrichment.13.png",width = 3,height =2.5,dpi=680)
ggsave("LoFs.enrichment.13.pdf",width = 3,height =2.5) 
############################## Figure 2g and S2h ############################
setwd("D:\\vcf_for_gte\\losfunction\\08_los\\gwas")
library(stringr)
library(patchwork)
library(ggplot2)
tmp<-read.csv("all_sig_snp_regions_eqtl.txt",header = F,sep="\t")
for(i in unique(tmp$V2)){
  i<-"17_50768323"
  j<-unlist(strsplit(i,"_"))
  Chr<-j[1]
  chr<-gsub("chr","",Chr)
  
  pos<-as.numeric(j[2])
  
  snp_id<-i
  tmp1<-subset(tmp,V2 == i,)
  
  tmp1$colrs<-""
  tmp1$sizs<-""
  for(j in 1:nrow(tmp1)){
    if(tmp1$V6[j]< 5e-08 & tmp1$V4[j] == pos){
      tmp1$colrs[j]<-1
      tmp1$sizs[j]<-3
    } else if(tmp1$V6[j]<5e-08 & tmp1$V4[j] != pos){
      tmp1$colrs[j]<-2
      tmp1$sizs[j]<-1
    }else{
      tmp1$colrs[j]<-3
      tmp1$sizs[j]<-1
    }
  }
  num1<-length(unique(tmp1$V1)) 
  p<-ggplot(tmp1, aes(V4/1000000,-log10(V6)))+
    geom_point(aes(color=colrs,size=sizs,shape=sizs),alpha=0.8)+
    scale_color_manual(values = c("#ED0000FF","#0099B4FF","#ADB6B6FF"))+
    scale_size_manual(values = c(1,6))+
    ###scale_shape_manual(values = c(16,25))+
    labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
    theme_classic()+
    geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
    facet_wrap(.~V1,scales = "free_y",ncol = 1)+
    theme(##panel.background = element_blank(),
      legend.position = "",
      ##panel.grid.major.y = element_line(colour = "grey",linetype = 2),
      axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))),
      axis.title.y = element_text(size = rel(2),hjust = 0.5),
      axis.title.x = element_text(size = rel(2),hjust = 0.5),
      axis.text.x = element_text(size = rel(2),hjust = 0.5,vjust=0.9),
      axis.text.y = element_text(hjust = 1,size = rel(2)),
      axis.ticks = element_line(size = rel(1.3)),
      ##plot.title = element_text(size = rel(1.8)),
      plot.margin = margin(15,9,9,30),
      strip.text = element_text(size = 15)
    ) 
  
  if(file.exists(snp_id)){
    
    tmp2<-read.table(snp_id,header = F,sep="\t")
    split_b<-str_split(tmp2$V3,"_")
    tmp2$pos<-as.numeric(sapply(split_b,"[",2))
    tmp2$chr<-sapply(split_b,"[",1)
    
    tmp2$GROUP<-paste(tmp2$V1,tmp2$V2,sep=":")
    tmp2$colrs<-""
    tmp2$sizs<-""
    
    for(j in 1:nrow(tmp2)){
      if(tmp2$pos[j] == pos){
        tmp2$colrs[j]<-1
        tmp2$sizs[j]<-3
      }else{
        tmp2$colrs[j]<-3
        tmp2$sizs[j]<-1
      }
    }
    
    p1<-ggplot(tmp2, aes(pos/1000000,-log10(V8)))+
      geom_point(aes(color=colrs,size=sizs,shape=sizs),alpha=0.8)+
      scale_color_manual(values = c("#ED0000FF","#ADB6B6FF"))+
      scale_size_manual(values = c(1,6))+
      ###scale_shape_manual(values = c(16,25))+
      labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
      theme_classic()+facet_wrap(.~GROUP,scales = "free_y",ncol = 1)+
      ##geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
      theme(##panel.background = element_blank(),
        legend.position = "",
        ##panel.grid.major.y = element_line(colour = "grey",linetype = 2),
        axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))),
        axis.title.y = element_text(size = rel(2),hjust = 0.5),
        axis.title.x = element_text(size = rel(2),hjust = 0.5),
        axis.text.x = element_text(size = rel(2),hjust = 0.5,vjust=0.9),
        axis.text.y = element_text(hjust = 1,size = rel(2)),
        axis.ticks = element_line(size = rel(1.3)),
        ##plot.title = element_text(size = rel(1.8)),
        plot.margin = margin(15,9,9,30),
        strip.position = element_text(size = 15),
        
      )  
    
    nums<-length(unique(tmp2$GROUP))
    
    p3<-p/p1+plot_layout(heights= c(num1, nums/1.5))
    
    ggsave(paste(i,"_gwas_1.png",sep=""),plot=p3,width = 6,height =10,dpi=680,limitsize = FALSE)
    ggsave(paste(i,"_gwas_1.pdf",sep=""),plot=p3,width = 6,height =10,limitsize = FALSE) 
    
  }else{
    
    ggsave(paste(i,"_gwas_1.png",sep=""),plot=p,width = 12,height =3*num1,dpi=680,limitsize = FALSE)
    ggsave(paste(i,"_gwas_1.pdf",sep=""),plot=p,width = 12,height =3*num1,limitsize = FALSE) 
  }  
  
}
######################### Figure 2h #########################################
library(ggplot2)
setwd("D:\\vcf_for_gte\\losfunction\\08_los\\gwas")
tmp<-data.frame(pops=c("ASD","ASD","ASW","ASW","EUD","EUD","EUW","EUW"),
                allele=c("T","C","T","C","T","C","T","C"),
                propo=c(0.944527,0.0554734,1,0,0.746491,0.253509,0.764706,0.235294))
ggplot(data=tmp,aes(propo,pops,,fill=allele))+geom_bar(stat = "identity")+
  theme_classic()+scale_x_continuous(expand = c(0,0))
ggsave("four_pops.allele_chr17_5076832.pdf",width = 6,height = 3)
########################## Figure S2i #####################################
library(tidyverse)
library(ggsci)
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/08_los/ase")
ase<-read.table("all.sample.select.site3.txt",header=T,sep="\t")
ase$ids<-paste(ase$Chr,ase$pos,sep="_")
for(i in unique(ase$ids)){
  
  tmp<-subset(ase,ids==i,)
  tmp <- tmp %>%
    mutate(sizes1 = case_when(  
      FDR > 0.05 ~ 1, 
      FDR<= 0.05  ~ 2))
  tmp$colr<-""
  for(j in 1:nrow(tmp)){
    if(tmp[j,]$FDR>0.05){
      tmp[j,]$colr<-"AAA"
    }else{
      tmp[j,]$colr<-tmp[j,]$Tissue
    }
  }
  
  tmp$sizes1<-as.factor(tmp$sizes1)
  a<-max(tmp$acount)
  b<-max(tmp$bcount)
  if(a>=b){
    ggplot(tmp)+geom_point(aes(x=acount,y=bcount,color=colr,shape=sizes1,size=sizes1))+
      scale_color_manual(values = c("#66666699","#FF000099","#FF990099","#00FF0099","#6699FF99",
                                    "#CC33FF99","#99991E99","#FF00CC99","#CC000099","#FFCCCC99",
                                    "#CCFF0099","#35800099","#0000CC99","#99CCFF99","#00FFFF99",
                                    "#CC99FF99","#79CC3D99","#1F77B4FF","#FF7F0EFF","#2CA02CFF",
                                    "#D62728FF","#9467BDFF","#E377C2ff","#BCBD22FF","#17BECFFF",
                                    "#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF","#C5B0D5FF",
                                    "#F7B6D2FF","#DBDB8DFF","#9EDAE5FF","#D62728FF","#9467BDFF"))+
      geom_abline( aes(intercept = 0,slope = 1),color="red",linetype="dashed",size=0.5)+
      ylim(0,a)+xlim(0,a)+xlab("reference allel count") +ylab("alternative allel count")+
      scale_shape_manual(values = c(20,42))+
      scale_size_manual(values = c(6,20))+
      theme_classic(base_size = 10)+theme_classic(base_size = 10)+
      guides(color=guide_legend(override.aes = list(size=5),
                                keywidth=0.1,
                                keyheight=0.1,
                                default.unit="inch"))
    
  }else{
    ggplot(tmp)+geom_point(aes(x=acount,y=bcount,color=colr,shape=sizes1,size=sizes1))+
      scale_color_manual(values = c("#66666699","#FF000099","#FF990099","#00FF0099","#6699FF99",
                                    "#CC33FF99","#99991E99","#FF00CC99","#CC000099","#FFCCCC99",
                                    "#CCFF0099","#35800099","#0000CC99","#99CCFF99","#00FFFF99",
                                    "#CC99FF99","#79CC3D99","#1F77B4FF","#FF7F0EFF","#2CA02CFF",
                                    "#D62728FF","#9467BDFF","#E377C2ff","#BCBD22FF","#17BECFFF",
                                    "#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF","#C5B0D5FF",
                                    "#F7B6D2FF","#DBDB8DFF","#9EDAE5FF","#D62728FF","#9467BDFF"))+
      geom_abline( aes(intercept = 0,slope = 1),color="red",linetype="dashed",size=0.5)+
      ylim(0,b)+xlim(0,b)+xlab("reference allel count") +ylab("alternative allel count")+
      scale_shape_manual(values = c(20,42))+
      scale_size_manual(values = c(6,20))+
      theme_classic(base_size = 10)+
      guides(color=guide_legend(override.aes = list(size=5),
                                keywidth=0.1,
                                keyheight=0.1,
                                default.unit="inch"))
    
  }
  ggsave(paste("ase",i,".png",sep=""),width = 12,height =8,dpi=680)
  ggsave(paste("ase",i,".pdf",sep=""),width = 12,height =8)     
}
