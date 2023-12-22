#############################################################################
######################## Figure 3 and S3 ####################################
#############################################################################

######################## Figure 3a ##########################################
library(ggplot2)
library(ggsci)
setwd("D:/vcf_for_gte/losfunction/07_utr")
all_pcadd<-read.csv("utr_pcadd",header = F,sep = "\t")
all_pcadd$V1[all_pcadd$V1=="5_prime_UTR_variant"]<-"5\'UTR"
all_pcadd$V1[all_pcadd$V1=="missense_variant"]<-"Missense"
all_pcadd$V1[all_pcadd$V1=="synonymous_variant"]<-"Synonymous"
lei<-c("LoFs","uAUG_gained","uAUG_lost","uSTOP_gained","uSTOP_lost","Missense","Synonymous","5\'UTR")
pcadd_st<-data.frame(ID=NA,p_mean=NA,p_se=NA)
for(i in lei){
  tmp<-subset(all_pcadd,V1==i,)
  tmp_mean<-mean(tmp$V3)
  tmp_se<-sd(tmp$V3)/sqrt(nrow(tmp))
  tmp_all<-c(i,tmp_mean,tmp_se)
  pcadd_st<-rbind(pcadd_st,tmp_all)
}
pcadd_st<-na.omit(pcadd_st)
pcadd_st$p_mean<-as.numeric(pcadd_st$p_mean)
pcadd_st$p_se<-as.numeric(pcadd_st$p_se)
pcadd_st$ID<-factor(pcadd_st$ID,
                    c("LoFs","Missense","Synonymous","uAUG_gained","uAUG_lost","uSTOP_gained","uSTOP_lost","5\'UTR"))
ggplot(pcadd_st, aes(ID, p_mean,color=ID))+ geom_point(size=5)+
  scale_color_manual(values = c("#7E6148FF","#0099B4FF","#00468BFF","#E64B35FF","#00A087FF",
                                "#F39B7FFF","#8491B4FF","#FF990099"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), 
                width=0.4,position = position_dodge(.9)) + labs(y="",x="") + 
  theme_classic()+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.6),
        axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("all.pcadd.png",width = 9,height =7,dpi=680)
ggsave("all.pcadd.pdf",width =9,height =7)
######################### Figure 3b #########################################
setwd("D:/vcf_for_gte/losfunction/07_utr")
strength_st<-data.frame(GROUP=NA,ID=NA,p_mean=NA,p_se=NA)
for(k in utrs ){
  utr_aug<-read.csv(paste(k,".pcadd",sep=""),header=F,sep="\t")
  strengths<-c("Weak","Moderate","Strong")
  for(i in strengths){
    tmp<-subset(utr_aug,V6==i,)
    tmp_mean<-mean(tmp$V4)
    tmp_se<-sd(tmp$V4)/sqrt(nrow(tmp))
    tmp_all<-c(k,i,tmp_mean,tmp_se)
    strength_st<-rbind(strength_st,tmp_all)
  }
}
strength_st<-na.omit(strength_st)
strength_st$p_mean<-as.numeric(strength_st$p_mean)
strength_st$p_se<-as.numeric(strength_st$p_se)
strength_st$ID<-factor(strength_st$ID,levels=c("Weak","Moderate","Strong"))
ggplot(strength_st, aes(ID, p_mean,,color=GROUP))+ geom_point(size=2.5)+
  ##scale_color_lancet()+
  scale_color_manual(values = c("#E64B35FF","#00A087FF","#F39B7FFF","#8491B4FF"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
  labs(y="",x="") + 
  theme_classic()+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="LoFs"],linetype=2,colour="#00468BFF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Missense"],linetype=2,colour="#ED0000FF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Synonymous"],linetype=2,colour="#42B540FF")+
  theme(legend.position = "",
        strip.text = element_blank(),
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black")
  )+ facet_wrap( ~ GROUP,nrow = 1,ncol = 4)
scales = "free_y"
ggsave("utr.four.strength.png",width = 9,height =7,dpi=680)
ggsave("utr.four.strength.pdf",width = 9,height =7)
############################## Figure 3c ####################################
setwd("D:/vcf_for_gte/losfunction/07_utr")
catagory_st<-data.frame(GROUP=NA,ID=NA,p_mean=NA,p_se=NA)
for(k in utrs ){
  utr_aug<-read.csv(paste(k,".pcadd",sep=""),header=F,sep="\t")
  utr_aug$V7[utr_aug$V7=="inFrame_oORF"]<-"InFrame_oORF"
  cats<-unique(utr_aug$V7)
  for(i in cats){
    tmp<-subset(utr_aug,V7==i,)
    tmp_mean<-mean(tmp$V4)
    tmp_se<-sd(tmp$V4)/sqrt(nrow(tmp))
    tmp_all<-c(k,i,tmp_mean,tmp_se)
    catagory_st<-rbind(catagory_st,tmp_all)
  }
}
catagory_st<-na.omit(catagory_st)
catagory_st$p_mean<-as.numeric(catagory_st$p_mean)
catagory_st$p_se<-as.numeric(catagory_st$p_se)
ggplot(catagory_st, aes(ID, p_mean,,color=GROUP))+ geom_point(size=2.5)+
  ##scale_color_lancet()+
  scale_color_manual(values = c("#E64B35FF","#00A087FF","#F39B7FFF","#8491B4FF"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
  labs(y="",x="") + 
  theme_classic()+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="LoFs"],linetype=2,colour="#00468BFF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Missense"],linetype=2,colour="#ED0000FF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Synonymous"],linetype=2,colour="#42B540FF")+
  theme(legend.position = "",
        strip.background =element_blank(),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 90),
        axis.text.y = element_text(size = 17,color = "black")
  )+ facet_wrap( ~ GROUP,scales="free",nrow = 1,ncol = 4)  
##facet_grid(~GROUP,scales="free",space="free_x")
ggsave("utr.four.catagory.png",width = 9,height =7,dpi=680)
ggsave("utr.four.catagory.pdf",width = 9,height =7)
######################## Figure S3a ########################################
setwd("D:/vcf_for_gte/losfunction/07_utr")
distance_st<-data.frame(GROUP=NA,ID=NA,p_mean=NA,p_se=NA)
for(k in utrs){
  utr_aug<-read.csv(paste(k,".pcadd",sep=""),header=F,sep="\t")
  tmp<-subset(utr_aug,V5<50,)
  tmp_mean<-mean(tmp$V4)
  tmp_se<-sd(tmp$V4)/sqrt(nrow(tmp))
  tmp_all<-c(k,"<50bp",tmp_mean,tmp_se)
  distance_st<-rbind(distance_st,tmp_all)
  tmp<-subset(utr_aug,V5>=50,)
  tmp_mean<-mean(tmp$V4)
  tmp_se<-sd(tmp$V4)/sqrt(nrow(tmp))
  tmp_all<-c(k,">=50bp",tmp_mean,tmp_se)
  distance_st<-rbind(distance_st,tmp_all)
}

distance_st<-na.omit(distance_st)
distance_st$p_mean<-as.numeric(distance_st$p_mean)
distance_st$p_se<-as.numeric(distance_st$p_se)
ggplot(distance_st, aes(ID, p_mean,,color=GROUP))+ geom_point(size=2.5)+
  ##scale_color_lancet()+
  scale_color_manual(values = c("#E64B35FF","#00A087FF","#F39B7FFF","#8491B4FF"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
  labs(y="",x="") + 
  theme_classic()+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="LoFs"],linetype=2,colour="#00468BFF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Missense"],linetype=2,colour="#ED0000FF")+
  ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Synonymous"],linetype=2,colour="#42B540FF")+
  theme(legend.position = "",
        strip.background =element_blank(),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black")
  )+ facet_wrap( ~ GROUP,nrow = 1,ncol = 4)

ggsave("utr.four.distance.png",width = 9,height =7,dpi=680)
ggsave("utr.four.distance.pdf",width = 9,height =7)

########################## Figure 3d and S3b ############################### 
setwd("D:/vcf_for_gte/losfunction/11_noncoding/all_cats")
tmp<-fread("noncoding_cats_pcadd_enrich.txt")
tmp$V2<-factor(tmp$V2,levels = c("0.01","0.025","0.05","0.075","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5",
                                 "0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1"))
tmp$foldchang<-(tmp$V8/tmp$V7)/(tmp$V6/tmp$V5)
for(i in unique(tmp$V1)){
  tmp1<-subset(tmp,V1==i,)
  fc_st<-data.frame(percs=NA,regs=NA,p_mean=NA,p_se=NA)
  for(j in unique(tmp1$V2)){
    for(k in unique(tmp1$V4)){
      tmp2<-subset(tmp1,V2==j & V4== k,)
      tmp_mean<-mean(tmp2$foldchang)
      tmp_se<-sd(tmp2$foldchang)/sqrt(nrow(tmp2))
      tmp_all<-c(j,k,tmp_mean,tmp_se)
      fc_st<-rbind(fc_st,tmp_all)
    }
  }
  
  fc_st<-na.omit(fc_st)
  fc_st$p_mean<-as.numeric(fc_st$p_mean)
  fc_st$p_se<-as.numeric(fc_st$p_se)
  fc_st$regs<-factor(fc_st$regs,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                           "E12","E13","E14","E15"))
  ggplot(fc_st, aes(percs, p_mean,group=regs,color=regs))+ geom_line(size=0.5) +
    geom_point(size=1.5)+
    ##scale_color_lancet()+
    scale_color_manual(values = c("Red","DeepPink", "ForestGreen","Green3", "Green1","Yellow", "Goldenrod3","Gold2",
                                  "Gold1", "LightGoldenrod1","SkyBlue", "Chocolate1","grey41", "grey61", "black"))+
    geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
    labs(y="",x="") + 
    theme_classic()+
    theme(legend.position = "right",
          axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
          axis.text.y = element_text(size = 17,color = "black")
    )
  ggsave(paste(i,".noncoding_enrichment.png",sep=""),width = 8,height =3.5,dpi=680)
  ggsave(paste(i,".noncoding_enrichment.pdf",sep=""),width = 8,height =3.5)
  
  fc_st1<-subset(fc_st,regs=="E1"|regs=="E12",)
  ggplot(fc_st1, aes(percs, p_mean,group=regs,color=regs))+ geom_line(size=0.5) +
    geom_point(size=1.5)+
    ##scale_color_lancet()+
    scale_color_manual(values = c("Red","Chocolate1"))+
    geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
    geom_hline(yintercept = 1.0,linetype="dashed")+
    labs(y="",x="") + 
    theme_classic()+
    ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="LoFs"],linetype=2,colour="#00468BFF")+
    ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Missense"],linetype=2,colour="#ED0000FF")+
    ##geom_hline(yintercept = pcadd_st$p_mean[pcadd_st$ID=="Synonymous"],linetype=2,colour="#42B540FF")+
    theme(legend.position = "right",
          axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
          axis.text.y = element_text(size = 17,color = "black")
    )
  ggsave("all.noncoding_enrichment1.112.png",width = 9,height =3.5,dpi=680)
  ggsave("all.noncoding_enrichment1.112.pdf",width = 9,height =3.5) 
}

####################### Figure 3e  ########################################
library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd("D:/vcf_for_gte/losfunction/07_utr/enrich")
utr<-read.csv("utr.5.enrichment.txt",header=F,sep=" ")
utr$V3<-factor(utr$V3,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                 "E12","E13","E14","E15"))
utr$foldchange<-0
utr$pval<-0
for(i in unique(utr$V1)){
  for(j in unique(utr$V2)){
    for(k in unique(utr$V3)){
      l1<-utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$V5
      l2<-utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$V7
      l3<-utr[utr$V1 == i & utr$V2== j & utr$V3=="E15",]$V5
      l4<-utr[utr$V1 == i & utr$V2== j & utr$V3=="E15",]$V7
      utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$foldchange<-(l2/l1)/(l4/l3)
      utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$pval<-fisher.test(matrix(c(l2,l4,l1,l3),ncol = 2,nrow = 2),
                                                                    alternative = "two.sided")$p.value
    }
  }
}

utr3<-subset(utr,V3 != "E15",)
utr3$fdr<-p.adjust(utr3$pval,method = "fdr")
utr3 <- utr3 %>%
  mutate(text = case_when(  # 一定要 get ??? case_when() 函数奥秘
    fdr< 0.05 ~ "*",
    fdr >=0.05 ~" "))
utr4<-subset(utr3,V1!="UTR",)
utr5<-subset(utr3,V1=="UTR",)
ggplot(utr5, aes(V3, foldchange,color=V3)) + 
  geom_boxplot() +
  scale_x_discrete(labels=c("TssA","TssAHet","TxFlnk","TxFlnkWk","TxFlnkHet","EnhA","EnhAMe","EnhAWk",
                            "EnhAHet","EnhPois","ATAC_ls","TssBiv","Repr","ReprWk"))+
  scale_color_manual(values = c("Red","DeepPink", "ForestGreen","Green3", "Green1","Yellow", "Goldenrod3","Gold2",
                                "Gold1", "LightGoldenrod1","SkyBlue", "Chocolate1","grey41", "grey61", "black"))+
  theme_classic() + 
  theme(legend.position = "",
        axis.title.x=element_blank(), # 去掉 title
        axis.text.x = element_text(angle = 30, hjust = 1,size = 14), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14))
ggsave("HIUTR.enrich.box.png", width = 6, height = 4, dpi = 680)
ggsave("HIUTR.enrich.box.pdf", width = 6, height = 4)
########################## Figure S3c #######################################
library(tidyverse)
library(ggplot2)
library(ggpubr)
setwd("D:/vcf_for_gte/losfunction/07_utr/enrich")
utr<-read.csv("utr.5.enrichment.txt",header=F,sep=" ")
utr$V3<-factor(utr$V3,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                 "E12","E13","E14","E15"))
utr$foldchange<-0
utr$pval<-0
for(i in unique(utr$V1)){
  for(j in unique(utr$V2)){
    for(k in unique(utr$V3)){
      l1<-utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$V5
      l2<-utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$V7
      l3<-utr[utr$V1 == i & utr$V2== j & utr$V3=="E15",]$V5
      l4<-utr[utr$V1 == i & utr$V2== j & utr$V3=="E15",]$V7
      utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$foldchange<-(l2/l1)/(l4/l3)
      utr[utr$V1 == i & utr$V2== j & utr$V3== k,]$pval<-fisher.test(matrix(c(l2,l4,l1,l3),ncol = 2,nrow = 2),
                                                                    alternative = "two.sided")$p.value
    }
  }
}
utr3<-subset(utr,V3 != "E15",)
ggplot(utr3, aes(V2, V3)) + 
  geom_tile(aes(fill = foldchange ), colour = "grey", size = 0.1)+
  scale_fill_gradientn(colors = c("#90EE90","white","#EA2E2D")) +
  scale_y_discrete(labels=c("TssA","TssAHet","TxFlnk","TxFlnkWk","TxFlnkHet","EnhA","EnhAMe","EnhAWk",
                            "EnhAHet","EnhPois","ATAC_ls","TssBiv","Repr","ReprWk"))+
  scale_x_discrete(labels=c("Adi.","Cec.","Cer.","Col.","Cor.","Duo.","Hyp.","Ile.",
                            "Jej.","Liv.","Lug.","Mus.","Spl.","Sto."),position = "top")+
  geom_text(aes(label=text),col ="black",size = 5,vjust = 0.5) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x ???
        axis.title.y=element_blank(), # 去掉 y ???
        axis.text.x = element_text(angle = 90, hjust = 0,size = 10, color="black"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, color="black"),
        strip.text.x  = element_text(size=14)) +
  labs(fill =paste0("* FDR< 0.05","\n\n","Foldchange")) +   
  facet_wrap(.~V1,nrow=1,strip.position = "bottom")
ggsave("UTR.enrich4_E15.png", width = 16, height = 4, dpi = 680)
ggsave("UTR.enrich4_E15.pdf", width = 16, height = 4)

###################### Figure 3f ############################################
eqtl_egene_0.05<-c("Blastomere","Cartilage","Duodenum","Kidney","Lymph_node","Morula","Oocyte","Pituitary","Placenta")
sqtl_egene_0.05<-c("Adipose","Blood","Brain","Embryo","Hypothalamus","Ileum","Jejunum","Large_intestine",
                   "Liver","Lung","Muscle","Ovary","Small_intestine","Spleen","Testis","Uterus")
setwd("D:/vcf_for_gte/losfunction/07_utr/qtl")
library(ggplot2)
tmp<-read.csv("utr.esqtl.enrichment",header = F,sep = " ")
tmp$all<-(tmp$V7/tmp$V6)/(tmp$V5/tmp$V4)

tmp_e<-subset(tmp,V1=="eQTL",select = c("V1","V2","V3","all"))
tmp_s<-subset(tmp,V1=="sQTL",select = c("V1","V2","V3","all"))

tmp_e<-tmp_e[!(tmp_e$V2 %in% eqtl_egene_0.05),]
tmp_s<-tmp_s[(tmp_s$V2 %in% sqtl_egene_0.05),]

tmp_l<-rbind(tmp_e,tmp_s)

ggplot(tmp_l,aes(V1,all,fill=V3)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#FF990099","#E64B35FF","#00A087FF",
                               "#F39B7FFF","#8491B4FF"))+
  labs(y="",x="")+
  theme_classic()+
  labs(x="",y="")+guides(alpha="none",)+
  theme(legend.position = "top",
        axis.text.y = element_text(size = 17,color = "black"),
        axis.text.x= element_text(size = 17,color = "black")
  )
ggsave("UTR.enrichment.esQTL.boxplot.png",width = 6,height =5,dpi=680)
ggsave("UTR.enrichment.esQTL.boxplot.pdf",width = 6,height =5) 

ggplot(tmp_l,aes(V3,all,color=V3)) +
  geom_point(size=1.0,position = position_dodge2(0.8),fill="black")+
  geom_line(group=1)+
  scale_color_manual(values = c("#FF990099","#E64B35FF","#00A087FF",
                                "#F39B7FFF","#8491B4FF"))+
  labs(y="",x="")+
  theme_classic()+
  labs(x="",y="")+guides(alpha="none",)+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 17,color = "black")
  )+facet_grid(vars(V1),vars(V2),scales = "fixed")

ggsave("UTR.enrichment.esQTL.png",width = 17,height =3.5,dpi=680)
ggsave("UTR.enrichment.esQTL.pdf",width = 17,height =3.5) 

########################### Figure 3g #####################################
setwd("D:/vcf_for_gte/losfunction/07_utr")
library(ggplot2)
tmp<-read.csv("utr.enrichment.gwas.txt",header = F,sep = "\t")
tmp$V1<-gsub("M_","",tmp$V1)
tmp$V2<-as.numeric(gsub("<","",tmp$V2))
tmp$V5<-2^(tmp$V2)
tmp$V6<-2^(tmp$V3)
tmp$V7<-2^(tmp$V4)
ggplot(tmp,aes(x=V6, y=V1))+
  geom_point(size=1.5)+
  geom_segment(aes(x=V5,xend=V7,y=V1,yend=V1),cex=0.2)+
  geom_vline(xintercept=1,color="red",linetype="dashed")+
  scale_y_discrete(position = 'right')+
  labs(y=NULL)+theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("utrs.enrichment.13.png",width = 2.5,height =2,dpi=680)
ggsave("utrs.enrichment.13.pdf",width = 2.5,height =2) 

######################## Figure 3h #########################################
library(ggplot2)
library(stringr)
library(patchwork)
setwd("D:/vcf_for_gte/losfunction/07_utr/gwas")
snp1<-read.csv("M_TLWT_BA.hiutr.chr2_87000000-89050000",header=F,sep="\t")
snp1$colrs<-1
snp1[snp1$V2==87778783,]$colrs<-2
snp1$colrs<-as.factor(snp1$colrs)
"#0099B4FF"
p<-ggplot(snp1, aes(V2/1000000,-log10(V4)))+
  geom_point(aes(color=colrs,size=colrs,shape=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ADB6B6FF","#ED0000FF"))+
  scale_size_manual(values = c(1,6))+
  ###scale_shape_manual(values = c(16,25))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+
  geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
  coord_cartesian(xlim = c(87000000/1000000, 89050000/1000000))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
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

tmp2<-read.table("2_87778783",header = F,sep="\t")
split_b<-str_split(tmp2$V3,"_")
tmp2$pos<-as.numeric(sapply(split_b,"[",2))
tmp2$chr<-sapply(split_b,"[",1)

tmp2$GROUP<-paste(tmp2$V1,tmp2$V2,sep=":")
tmp2$colrs<-""
tmp2$sizs<-""

for(j in 1:nrow(tmp2)){
  if(tmp2$pos[j] == 87778783){
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
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  coord_cartesian(xlim = c(87000000/1000000, 89050000/1000000))+
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
    strip.text = element_text(size = 15),
  )  
p3<-p/p1+plot_layout(heights= c(1, 1.5))

ggsave("M_TLWT_BA.hiutr.chr2_87000000-89050000.png",plot=p3,width = 12,height =6,dpi=680,limitsize = FALSE)
ggsave("M_TLWT_BA.hiutr.chr2_87000000-89050000.pdf",plot=p3,width = 12,height =6,limitsize = FALSE) 
##################################### Figure S3d #############################
ase<-read.csv("utr_2_87778783_ase.csv",header=T,sep=",")
ase$size2<-"1"
ase[ase$FDR<0.05,]$size2<-"2"
ggplot(ase, aes(x=Tissue, y=ref_ratio)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill=size2,color=size2),
               binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c("blue","red"))+
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.y = element_text(size = 17,hjust = 0.5),
    axis.title.x = element_text(size = 17,hjust = 0.5),
    axis.text.x = element_text(size = 17,hjust = 0.5,vjust=0.5,angle = 30),
    axis.text.y = element_text(hjust = 1,size = 17),
  ) 
ggsave("utr_2_87778783_ase.png",width = 6,height =5,dpi=680)
ggsave("utr_2_87778783_ase.pdf",width = 6,height =5) 