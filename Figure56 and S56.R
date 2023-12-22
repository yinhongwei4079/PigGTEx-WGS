#############################################################################
###################### Figure 5,6 S5 and S6 #################################
#############################################################################

all_col<-c("Adipose"="#CC66FF","Artery"="#AAAAFF","Blastocyst"="#99FF00", "Blastomere"="#4DAF4A",
           "Blood"="#FF0000","Brain"="#EEEE00","Cartilage"="#7777FF","Cecum"="#66ad05","Cerebellum"="#a7d1fa",
           "Colon"="#8EABD2","Cortex"="#2bb1ff","Duodenum"="#00B0F0",
           "Embryo"="#1B9E77","Fetal_thymus"="#FF5555","Frontal_cortex"="#FFD700","Heart"="#33CCCC",
           "Hypothalamus"="#FDFDBF","Ileum"="#8EA9DB","Jejunum"="#B4C6E7","Kidney"="#8B0F55",
           "Large_intestine"="#386CB0","Liver"="#E2EFDA","Lung"="#7570B3","Macrophage"="#FFCCCC",
           "Milk"="#727F3F","Lymph_node"="#FFAA99","Morula"="#A6D854","Muscle"="#AAEEFF",
           "Oocyte"="#33DD33","Ovary"="#99BB88","Pituitary"="#FFDD99","Placenta"="#22FFDD",
           "Small_intestine"="#0000FF","Spleen"="#FF6600","Stomach"="#8bfdc7","Synovial_membrane"="black",
           "Testis"="#A6CEE3","Uterus"="#A6761D")
all_col<-as.data.frame(all_col)
#################### Figure 5a ##############################################
library(ggplot2)
library(ggsci)
setwd("D:/vcf_for_gte/losfunction/05_baseji")
roc<-read.csv("all_tissue.roc.txt",header = F,sep = "\t")
roc$sample<-paste(roc$V1,roc$V2,sep="_")

roc<-subset(roc,V1!="Cerebellum",)

ggplot(roc,aes(V3,V4,group=sample,color=V1))+
  geom_line(size=0.3)+
  scale_color_manual(values = c("maroon","darkgreen","red","blue","green","orange","purple","pink","brown",
                                "navy","#4B0082","#808000","darkblue"))+
  labs(y="True positive rate",x="False positive rate") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black"),
    axis.title = element_text(size = 17,color = "black")
  )
ggsave("all_tissue.roc.png",width = 6,height =5,dpi=680)
ggsave("all_tissue.roc.pdf",width = 6,height =5)    

#################### Figure S5b ###########################################
setwd("D:/vcf_for_gte/losfunction/05_baseji")
pr<-read.csv("all_tissue.pr.txt",header = F,sep = "\t")
pr<-subset(pr,V1!="Cerebellum",)
pr$sample<-paste(pr$V1,pr$V2,sep="_")
ggplot(pr,aes(V4,V3,group=sample,color=V1))+
  geom_line(size=0.3)+
  scale_color_manual(values = c("maroon","darkgreen","red","blue","green","orange","purple","pink","brown",
                                "navy","#4B0082","#808000","darkblue"))+
  labs(x="Recall",y="Precision") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black"),
    axis.title = element_text(size = 17,color = "black")
  )
ggsave("all_tissue.pr.png",width = 6,height =5,dpi=680)
ggsave("all_tissue.pr.pdf",width = 6,height =5)  
############################### Figure 5b ##################################
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji")
pos<-read.csv("all.targets.preds.anno.bed",header = F,sep = "\t")
pos<-subset(pos,V1 != "Cerebellum", )

region_cor<-data.frame(regions=NA,samples=NA,cors=NA)
for(region in unique(pos$V9)){
  tmp<-subset(pos,V9==region,)
  for(tissue in unique(tmp$V1)){
    tmp1<-subset(tmp,V1==tissue)
    t0<-cor.test(tmp1$V5,tmp1$V6)
    t1<-cor.test(tmp1$V7,tmp1$V8)
    l<-c(region,paste(tissue,"t0",sep=""),t0$estimate[['cor']])
    l1<-c(region,paste(tissue,"t1",sep=""),t1$estimate[['cor']])
    region_cor<-rbind(region_cor,l)
    region_cor<-rbind(region_cor,l1)
  }
}
region_cor<-na.omit(region_cor)
region_cor$cors<-as.numeric(region_cor$cors)
ebtop<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
ebbottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
ggplot(region_cor,aes(regions,cors)) +
  geom_point(size=0.5)+
  geom_jitter(width = 0.2)+
  stat_summary(geom = "bar",fun = "mean",alpha=0.8,
               position = position_dodge(0.9),width=0.6,fill="#00468BFF")+
  stat_summary(geom = "errorbar",
               fun.min = ebbottom,
               fun.max = ebtop,
               position = position_dodge(0.9),
               width=0.4,size=0.8)+  
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+labs(y="Correlation(R2)",x="")+
  theme(
    axis.text.x = element_text(size = 17,color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 17,color = "black"),
    axis.title = element_text(size = 17,color = "black"))
ggsave("all_tissue.region_cor.png",width = 6,height =5,dpi=680)
ggsave("all_tissue.region_cor.pdf",width = 6,height =5) 
####################### Figure S5c ##########################################
setwd("D:/vcf_for_gte/losfunction/05_baseji/test")
library(ggplot2)
library(ggsci)
roc<-read.csv("all_adipose.roc.txt",header = F,sep = "\t")
roc$sample<-paste(roc$V1,roc$V2,sep="_")
ggplot(roc,aes(V3,V4,group=sample,color=V1))+
  geom_line(size=0.5)+
  scale_color_lancet()+
  labs(y="True positive rate",x="False positive rate") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 17,color = "black"),
    axis.text.y = element_text(size = 17,color = "black"),
    axis.title = element_text(size = 17,color = "black")
  )
ggsave("adipose.zhao.roc.png",width = 6,height =5,dpi=680)
ggsave("adipose.zhao.roc.pdf",width = 6,height =5)

################### Figure 5c ################################################
setwd("D:/vcf_for_gte/losfunction/05_baseji/tf")
library(ggplot2)
library(data.table)
library(ggunchained)
library(ggpubr)
library(Hmisc)
tmp<-fread("all_tissue.chr1-18.aver.sad.tf10.txt",header=F,sep="\t")
tmp<-subset(tmp,V1 != "Cerebellum")
tmp$V7<-log10(abs(tmp$V6))
##pnorm(1.645,mean=0,sd=1)-pnorm(-1.645,mean=0,sd=1)确定几倍标准差的置信区间
###top1%区间
sum_stat<-data.frame(tissues=NA,max_0.95=NA)
tissues<-unique(tmp$V1)
for(i in tissues){
  tmp1<-subset(tmp,V1==i & V2=="backgroup",)
  set = density(tmp1$V7)
  max_midu = set$x[which.max(set$y)]
  right = tmp1[tmp1$V7 > max_midu,]$V7
  max = max_midu + sd(c(right,max_midu*2-right))*2.37
  sum_stat<-rbind(sum_stat,c(i,max))
}
sum_stat<-na.omit(sum_stat)
sum_stat$max_0.95<-as.numeric(sum_stat$max_0.95)
colnames(sum_stat)[1]<-"V1"
ggplot(tmp,aes(V1,V7,fill=V2))+geom_split_violin()+
  stat_summary(fun = mean,size=0.8,geom = "point",
               position = position_dodge(width = 0.4))+
  stat_compare_means(data = tmp, aes(x = V1,y = V7),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "-")),
                     label = "p.signif",label.y = max(tmp$V7),hide.ns = F)+
  scale_fill_manual(values = c('#FB5554','#868B31'))+
  labs(y="Predicted effect",x="Tissues")+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
        legend.position = 'top')
ggsave("all_tisuue.snp.effect.violin.png",width = 9,height = 5,dpi=680)
ggsave("all_tisuue.snp.effect.violin.pdf",width = 9,height = 5)
####################### Figure 5d and S5d ###################################
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji/upsteam")
tmp<-read.csv("eqtl.up2kb.enrich1.txt",header = F,sep = " ")
tmp$foldchange<-(tmp$V7/tmp$V6)/(tmp$V5/tmp$V4)
tmp<-subset(tmp,V2 != "Cerebellum")

tmp1<-read.csv("sqtl.up2kb.enrich1.txt",,header = F,sep = " ")
tmp1$foldchange<-(tmp1$V7/tmp1$V6)/(tmp1$V5/tmp1$V4)
tmp1<-subset(tmp1,V2 != "Cerebellum")

eqtl_egene_0.05<-c("Blastomere","Cartilage","Duodenum","Kidney","Lymph_node","Morula","Oocyte","Pituitary","Placenta")
tmp_e<-tmp[!(tmp$V1 %in% eqtl_egene_0.05),]
sqtl_egene_0.05<-c("Adipose","Blood","Brain","Embryo","Hypothalamus","Ileum","Jejunum","Large_intestine",
                   "Liver","Lung","Muscle","Ovary","Small_intestine","Spleen","Testis","Uterus")
tmp_s<-tmp1[(tmp1$V1 %in% sqtl_egene_0.05),]
tmp_e$cata<-"eQTL"
tmp_s$cata<-"sQTL"
tmp3<-rbind(tmp_e,tmp_s)
tmp4<-subset(tmp3,V3==0.1 |V3==0.2 |V3==0.3,)
tmp4$V3<-as.factor(tmp4$V3)
ggplot(tmp4,aes(V3,foldchange,color=cata)) + geom_boxplot()+
  geom_hline(yintercept =1.0,linetype=2,colour="red")+
  ##scale_x_continuous(expand = c(0.01, 0))+
  ##3scale_y_discrete(expand = c(0,0))+
  theme_classic()+
  labs(y="Enrichment",x="")+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black"))
ggsave("upstream_enrichment3_boxplot.png",width =6,height =5,dpi=680)
ggsave("upstream_enrichment3_boxplot.pdf",width = 6,height =5)

tmp3$V3<-as.factor(tmp3$V3)
ggplot(tmp3,aes(V3,foldchange,color=cata)) + geom_boxplot()+
  geom_hline(yintercept =1.0,linetype=2,colour="red")+
  ##scale_x_continuous(expand = c(0.01, 0))+
  ##3scale_y_discrete(expand = c(0,0))+
  theme_classic()+
  labs(y="Enrichment",x="")+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black"))
ggsave("upstream_enrichment3_boxplot_all.png",width =11,height =4,dpi=680)
ggsave("upstream_enrichment3_boxplot_all.pdf",width = 11,height =5)
########################### Figure 5f and S5e ###############################
library(data.table)
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji/all_sad")
tmp<-read.csv("all.sad.daf_1.txt",header = F,sep = " ")
tmp$foldchang<-(tmp$V7/tmp$V6)/(tmp$V5/tmp$V4)
tmp1<-subset(tmp,V3==0.1,)
tmp1$V1<-as.factor(tmp1$V1)
mean(tmp1[tmp1$V3==0.1 & tmp1$V1 ==1.00,]$foldchang)
ggplot(tmp1, aes(V1,foldchang)) + 
  geom_boxplot()+
  geom_hline(yintercept = 1.0,linetype="dashed",color="red")+
  ##scale_fill_gradient2(low="#90EE90",mid="white",high="#F7AA97",midpoint = 0) +
  ##geom_text(aes(label=round(tmp1$perces,2)),col ="black",size = 5,vjust = 0.5) +
  theme_classic() + 
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.title.y=element_blank(), # 去掉 y ???
        axis.text.x = element_text(angle = 90, vjust = 0.5,size = 14, color="black"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, color="black")) 
ggsave('all.tissue.all.cata.sad.enrich.png',width = 12,  height = 4, dpi = 680)
ggsave('all.tissue.all.cata.sad.enrich.pdf',width = 12,  height = 4)

tmp$V1<-as.factor(tmp$V1)
ggplot(tmp, aes(V1,foldchang)) + 
  geom_boxplot()+
  geom_hline(yintercept = 1.0,linetype="dashed",color="red")+
  ##scale_fill_gradient2(low="#90EE90",mid="white",high="#F7AA97",midpoint = 0) +
  ##geom_text(aes(label=round(tmp1$perces,2)),col ="black",size = 5,vjust = 0.5) +
  theme_classic() + 
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.title.y=element_blank(), # 去掉 y ???
        axis.text.x = element_text(angle = 90, vjust = 0.5,size = 14, color="black"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, color="black"))+
  facet_wrap(vars(V3),ncol=1,strip.position = "right",scales = "free_y")
ggsave('all.tissue.all.cata.sad.enrich.all.png',width = 12,  height = 9, dpi = 680)
ggsave('all.tissue.all.cata.sad.enrich.all.pdf',width = 12,  height = 9)
########################## Figure 5e ########################################
library(data.table)
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji/all_sad")
a<-list.files('D:/vcf_for_gte/losfunction/05_baseji/all_sad',pattern="*_sad.txt")
sad_st<-data.frame(tissues=NA,catagory=NA,p_mean=NA,p_se=NA)
for (i in a){
  tissue<-unlist(strsplit(i,"_"))[1]
  tmp<-fread(i)
  for(j in unique(tmp$V5)){
    tmp1<-subset(tmp,V5 == j,)
    tmp_mean<-mean(tmp1$V4)
    tmp_se<-sd(tmp1$V4)/sqrt(nrow(tmp1))
    l<-c(tissue,j,tmp_mean,tmp_se)
    sad_st<-rbind(sad_st,l)
  }
}
sad_st<-na.omit(sad_st)
sad_st$p_mean<-as.numeric(sad_st$p_mean)
sad_st$p_se<-as.numeric(sad_st$p_se)
sad_st$catagory<-factor(sad_st$catagory,levels = c("synonymous_variant","missense_variant","LoFs",
                                                   "non-CDS","CDS","uAUG_gained","uAUG_lost",
                                                   "uSTOP_gained","uSTOP_lost","5_utr"))
ggplot(sad_st, aes(tissues, p_mean,fill=catagory,group=catagory)) + 
  geom_bar(stat = "identity",position="dodge")+
  scale_fill_manual(values = c("#00468BFF","#0099B4FF","#7E6148FF","79cc3dff","#9900ccff","#E64B35FF","#00A087FF",
                               "#F39B7FFF","#8491B4FF","#FF990099"))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black"))
ggsave('all.tissue.all.cata.sad.png',width = 6,  height = 4, dpi = 680)
ggsave('all.tissue.all.cata.sad.pdf',width = 6,  height = 4) 

##################### Figure 6a and S6a ######################################
library(ggplot2)
library(data.table)
library(ggsci)
setwd("D:/vcf_for_gte/losfunction/05_baseji/qtl_snp")
tmp1<-fread("eqtl.snp.enrich1.txt")
tmp2<-fread("sqtl.snp.enrich1.txt")

tmp1$foldchang<-(tmp1$V7/tmp1$V6)/(tmp1$V5/tmp1$V4)
tmp2$foldchang<-(tmp2$V7/tmp2$V6)/(tmp2$V5/tmp2$V4)

eqtl_egene_0.05<-c("Blastomere","Cartilage","Duodenum","Kidney","Lymph_node","Morula","Oocyte","Pituitary","Placenta")
tmp_e<-tmp1[!(tmp1$V1 %in% eqtl_egene_0.05),]

sqtl_egene_0.05<-c("Adipose","Blood","Brain","Embryo","Hypothalamus","Ileum", "Jejunum","Large_intestine",
                   "Liver","Lung","Muscle","Ovary","Small_intestine","Spleen","Testis","Uterus")
tmp_s<-tmp2[(tmp2$V1 %in% sqtl_egene_0.05),]

tmp_e$cats<-"eQTL"
tmp_s$cats<-"sQTL"
tmp_c<-rbind(tmp_e,tmp_s)
tmp_c<-subset(tmp_c,V2!="Cerebellum",)

fd_st_c<-data.frame(cats=NA,tissues=NA,sad=NA,mean_fd=NA,sd_fd=NA)
for(i in unique(tmp_c$cats)){
  for(j in unique(tmp_c$V2)){
    for(k in unique(tmp_c$V3)){
      tmp_m<-subset(tmp_c,cats==i & V2 == j & V3 ==k ,)
      tmp_mean<-mean(tmp_m$foldchang)
      tmp_sd<-sd(tmp_m$foldchang)/sqrt(nrow(tmp_m))
      l<-c(i,j,k,tmp_mean,tmp_sd)
      fd_st_c<-rbind(fd_st_c,l)
    }
  }
}

fd_st_c<-na.omit(fd_st_c)
fd_st_c$sad<-as.factor(fd_st_c$sad)

fd_st_c$mean_fd<-as.numeric(fd_st_c$mean_fd)
fd_st_c$sd_fd<-as.numeric(fd_st_c$sd_fd)

ggplot(fd_st_c,aes(tissues,mean_fd,fill=sad,color=sad)) +
  geom_point(size=.2,position = position_dodge2(0.8))+
  geom_errorbar(aes(ymin=mean_fd - sd_fd, ymax=mean_fd + sd_fd),size=0.1, 
                width=0.2,position = position_dodge(.8)) + labs(y="",x="")+
  scale_colour_ucscgb()+
  theme_classic()+
  labs(x="Enrichment",y="")+guides(alpha="none",)+
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black")
  )+facet_wrap(.~cats,nrow=1,strip.position = "top",scales = "fixed")

ggsave("esQTL.reg.sad.noncoding_enrichment1.egene0.05.png",width = 8,height =3.5,dpi=680)
ggsave("esQTL.reg.sad.noncoding_enrichment1.egene0.05.pdf",width = 8,height =3.5) 

tmp_c$V3<-as.factor(tmp_c$V3)
ggplot(tmp_c,aes(V3,foldchang)) +
  geom_boxplot()+
  theme_classic()+
  geom_hline(yintercept =1.0,linetype=2,colour="red")+
  labs(x="Enrichment",y="")+guides(alpha="none",)+
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black")
  )+facet_wrap(.~cats,nrow=1,strip.position = "top",scales = "fixed")
ggsave("esQTL.reg.sad.noncoding_enrichment1.egene0.05.boxplot.png",width = 8,height =3.5,dpi=680)
ggsave("esQTL.reg.sad.noncoding_enrichment1.egene0.05.boxplot.pdf",width = 8,height =3.5)

################################ Figure 6b #################################
library(dplyr)
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji/gwas")
fgwas<-read.csv("traits_13.in_top0.01.res.txt",header = F,sep="\t")
f_max<-fgwas %>% group_by(V1) %>% slice_max(V4)
f_max$V1<-gsub("M_","",f_max$V1)
f_max$V3<-as.numeric(gsub("<","",f_max$V3))
f_max$V5<-as.numeric(f_max$V5)
f_max$V9<-paste(f_max$V1,"(",f_max$V2,")",sep = "")
f_max$V6<-2^(f_max$V3)
f_max$V7<-2^(f_max$V4)
f_max$V8<-2^(f_max$V5)

ggplot(f_max,aes(x=V7, y=V9))+
  geom_point(size=1.5)+
  geom_segment(aes(x=V6,xend=V8,y=V9,yend=V9),cex=0.2)+
  geom_vline(xintercept=1,color="red",linetype="dashed")+
  scale_y_discrete(position = 'right')+
  labs(y=NULL)+theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(0, 30))

ggsave("sad.enrichment.13.png",width = 2.5,height =2,dpi=680)
ggsave("sad.enrichment.13.pdf",width = 2.5,height =2)
###################### Figure 6c and S6b ####################################
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/05_baseji/gwas")
snp1<-read.csv("BFT.chr1_164800000-166050000",header=F,sep="\t")
snp3<-read.csv("BFT.chr7_29750000-31200000",header=F,sep="\t")
snp1$colrs<-1
snp1[snp1$V5==165048629,]$colrs<-2
snp1[snp1$V5==165083511,]$colrs<-2
snp1[snp1$V5==165119144,]$colrs<-2
snp1[snp1$V5==165775340,]$colrs<-2
snp1$colrs<-as.factor(snp1$colrs)

p<-ggplot(snp1, aes(V5/1000000,-log10(V7)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ADB6B6FF","#ED0000FF"))+
  scale_size_manual(values = c(1,4))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~V1,scales = "free_y",ncol = 1)+
  coord_cartesian(xlim = c(164800000/1000000, 166050000/1000000))+
  geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
  theme(
    legend.position = "",
    axis.title.y = element_text(size = 17,hjust = 0.5),
    axis.title.x = element_text(size = 17,hjust = 0.5),
    axis.text.x = element_text(size = 17,hjust = 0.5,vjust=0.9),
    axis.text.y = element_text(hjust = 1,size = 17),
    strip.text = element_text(size = 15)
  )


tmp2<-read.table("1_164800000",header = F,sep="\t")
split_b<-str_split(tmp2$V3,"_")
tmp2$pos<-as.numeric(sapply(split_b,"[",2))
tmp2$chr<-sapply(split_b,"[",1)

tmp2$GROUP<-paste(tmp2$V1,tmp2$V2,sep=":")
tmp2$colrs<-3
tmp2[tmp2$V1=="Muscle" & tmp2$pos==165083511,]$colrs<-1
tmp2[tmp2$V1=="Muscle" & tmp2$pos==165119144,]$colrs<-1
tmp2[tmp2$V1=="Muscle" & tmp2$pos==165775340,]$colrs<-1
tmp2[tmp2$V1=="Hypothalamus" & tmp2$pos==165048629,]$colrs<-1
tmp2[tmp2$V1=="Hypothalamus" & tmp2$pos==165083511,]$colrs<-1
tmp2[tmp2$V1=="Hypothalamus" & tmp2$pos==165119144,]$colrs<-1
tmp2[tmp2$V1=="Hypothalamus" & tmp2$pos==165775340,]$colrs<-1
tmp2$colrs<-as.factor(tmp2$colrs)

p1<-ggplot(tmp2, aes(pos/1000000,-log10(V8)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ED0000FF","#ADB6B6FF"))+
  scale_size_manual(values = c(6,1))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~GROUP,scales = "free_y",ncol = 1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  coord_cartesian(xlim = c(164800000/1000000, 166050000/1000000))+
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

ggsave("BFT.chr1_164800000-166050000.sad.png",plot=p3,width = 12,height =6,dpi=680)
ggsave("BFT.chr1_164800000-166050000.sad.pdf",plot=p3,width = 12,height =6)

snp3<-read.csv("BFT.chr7_29750000-31200000",header=F,sep="\t")
snp3$colrs<-1
snp3[snp3$V5==29825689,]$colrs<-2
snp3[snp3$V5==30205575,]$colrs<-2
snp3[snp3$V5==30363233,]$colrs<-2
snp3[snp3$V5==30389304,]$colrs<-2
snp3$colrs<-as.factor(snp3$colrs)

p<-ggplot(snp3, aes(V5/1000000,-log10(V7)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ADB6B6FF","#ED0000FF"))+
  scale_size_manual(values = c(1,4))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~V1,scales = "free_y",ncol = 1)+
  coord_cartesian(xlim = c(29750000/1000000, 31200000/1000000))+
  geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
  theme(
    legend.position = "",
    axis.title.y = element_text(size = 17,hjust = 0.5),
    axis.title.x = element_text(size = 17,hjust = 0.5),
    axis.text.x = element_text(size = 17,hjust = 0.5,vjust=0.9),
    axis.text.y = element_text(hjust = 1,size = 17),
  ) 

tmp3<-read.table("7_29750000",header = F,sep="\t")

split_b<-str_split(tmp3$V3,"_")
tmp3$pos<-as.numeric(sapply(split_b,"[",2))
tmp3$chr<-sapply(split_b,"[",1)

tmp3$GROUP<-paste(tmp3$V1,tmp3$V2,sep=":")
tmp3$colrs<-3
tmp3[tmp3$V1=="Blastocyst" & tmp3$pos==29825689 & tmp3$V2 =="ENSSSCG00000001533",]$colrs<-1
tmp3[tmp3$V1=="Heart" & tmp3$pos==30389304 & tmp3$V2 =="ENSSSCG00000001533",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==29825689 & tmp3$V2 =="ENSSSCG00000001533",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==30363233 & tmp3$V2 =="ENSSSCG00000001533",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==30389304 & tmp3$V2 =="ENSSSCG00000001533",]$colrs<-1

tmp3[tmp3$V1=="Testis" & tmp3$pos==29825689 & tmp3$V2 =="ENSSSCG00000027053",]$colrs<-1
tmp3[tmp3$V1=="Testis" & tmp3$pos==30363233 & tmp3$V2 =="ENSSSCG00000027053",]$colrs<-1
tmp3[tmp3$V1=="Testis" & tmp3$pos==30389304 & tmp3$V2 =="ENSSSCG00000027053",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==30389304 & tmp3$V2 =="ENSSSCG00000027053",]$colrs<-1
tmp3<-subset(tmp3,GROUP != "Blastocyst:ENSSSCG00000027053", )
tmp3$colrs<-as.factor(tmp3$colrs)

p1<-ggplot(tmp3, aes(pos/1000000,-log10(V8)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ED0000FF","#ADB6B6FF"))+
  scale_size_manual(values = c(4,1))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~GROUP,scales = "free_y",ncol = 1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  coord_cartesian(xlim = c(29750000/1000000, 31200000/1000000))+
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
p3<-p/p1+plot_layout(heights= c(1, 3))

ggsave("BFT.chr7_29750000-31200000.sad.png",width = 6,height=9,dpi=680)
ggsave("BFT.chr7_29750000-31200000.sad.pdf",width = 6,height=9,dpi=680)
###################### Figure S6c and 6d ###################################
ase<-read.csv("BFT.chr1_164800000-166050000.ase.csv",header=T,sep=",")
ase1<-read.csv("BFT.chr7_29750000-31200000.ase.csv",header=T,sep=",")

ase$pos<-as.factor(ase$pos)
ase$size2<-"1"
ase[ase$FDR<0.05,]$size2<-"2"

ase1$size2<-"1"
ase1[ase1$FDR<0.05,]$size2<-"2"

ggplot(ase, aes(x=pos, y=ref_ratio)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill=size2,color=size2),
               binaxis = 'y', stackdir = 'center', dotsize = 1)+
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

ggsave("BFT.chr1_164800000-166050000.ase.sad.png",width = 6,height =5,dpi=680)
ggsave("BFT.chr1_164800000-166050000.ase.sad.pdf",width = 6,height =5) 

ase1$pos<-as.factor(ase1$pos)
ggplot(ase1, aes(x=pos, y=ref_ratio)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill=size2,color=size2),
               binaxis = 'y', stackdir = 'center', dotsize = 1)+
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

ggsave("BFT.chr7_29750000-31200000.ase.sad.png",width = 6,height =5,dpi=680)
ggsave("BFT.chr7_29750000-31200000.ase.sad.pdf",width = 6,height =5) 
