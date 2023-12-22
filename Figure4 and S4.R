#############################################################################
######################## Figure 4 and S4 ####################################
#############################################################################

######################## Figure S4a ########################################
library(ggplot2)
library(ggsci)
library(scales)
setwd("D:/vcf_for_gte/losfunction/11_noncoding")

snp_comp<-data.frame(catss=c("SNPs","SNPs"),catogory=c("In Coding Region","In Non-coding Region"),nums=c(2253251,193918480))
snp_comp$var1<-round(snp_comp$nums/sum(snp_comp$nums),3)
snp_comp$var2<-paste(snp_comp$catogory,"(",snp_comp$var1*100,"%)",sep="")
ggplot( snp_comp, aes( x= "",y = var1, fill = catogory))+
  geom_bar( stat ="identity")+
  coord_polar(theta="y")+
  scale_fill_manual(values = c("79cc3dff","#9900ccff"))+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        axis.title = element_blank())
ggsave('all_composit_vep.png',width = 5,  height = 4, dpi = 680)
ggsave('all_composit_vep.pdf',width = 5,  height = 4)  

###################### Figure 4a ###########################################
library(data.table)
library(ggplot2)
library(ggpubr)
setwd("D:/vcf_for_gte/losfunction/11_noncoding")
tmp<-fread("all.no.non.pcadd")
tmp1<-fread("synonymous.missense.pcadd.bed")
tmp2<-rbind(tmp,tmp1)

pcadd_st<-data.frame(ID=NA,p_mean=NA,p_se=NA,nums=NA,sums=NA)
for(i in unique(tmp2$V5)){
  tmp3<-subset(tmp2,V5==i,)
  tmp_sum<-sum(tmp3$V4)
  tmp_mean<-mean(tmp3$V4)
  tmp_se<-sd(tmp3$V4)/sqrt(nrow(tmp3))
  tmp_all<-c(i,tmp_mean,tmp_se,nrow(tmp3),tmp_sum)
  pcadd_st<-rbind(pcadd_st,tmp_all)
}

pcadd_st<-na.omit(pcadd_st)
pcadd_st$p_mean<-as.numeric(pcadd_st$p_mean)
pcadd_st$p_se<-as.numeric(pcadd_st$p_se)
pcadd_st$ID<-factor(pcadd_st$ID,levels = c("Missense","Synonymous","non-CDS","CDS"))

ggplot(pcadd_st, aes(ID, p_mean,fill=ID))+ geom_bar(stat ="identity",width = 0.6,position = "dodge")+
  scale_fill_manual(values=c("#0099B4FF","#00468BFF","#9900ccff","79cc3dff"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), 
                width=0.5,position = "dodge",size=0.2) + labs(y="pCADD scores",x="") + 
  theme_classic()+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle=0,hjust = 1.0),
        ##axis.ticks.x = element_blank(),
        ##axis.line.y  = element_blank(),
        axis.text.y = element_text(size = 17,color = "black"),
        axis.title.y= element_text(size = 17,colour="black")
  )
ggsave('noncoding_CDS_vep.png',width = 5,  height = 4, dpi = 680)
ggsave('noncoding_CDS_vep.pdf',width = 5,  height = 4) 
########################### Figure S4b #####################################
setwd("D:/vcf_for_gte/losfunction/11_noncoding")
tmp<-read.csv("noncoding_pcadd",header=F,sep="\t")
tmp[tmp$V1=="non_coding_transcript_exon",]$V1<-"non_coding_transcript"

pcadd_st<-data.frame(ID=NA,p_mean=NA,p_se=NA,nums=NA,sums=NA)
for(i in unique(tmp$V1)){
  tmp1<-subset(tmp,V1==i,)
  tmp_sum<-sum(tmp1$V2)
  tmp_mean<-mean(tmp1$V2)
  tmp_se<-sd(tmp1$V2)/sqrt(nrow(tmp1))
  tmp_all<-c(i,tmp_mean,tmp_se,nrow(tmp1),tmp_sum)
  pcadd_st<-rbind(pcadd_st,tmp_all)
}
l<-c("Missense",11.7450614336233,0.0101033475124342,1,1)
k<-c("Synonymous",2.6501815556071,0.00333649618503868,1,1)
pcadd_st<-rbind(pcadd_st,l,k)
pcadd_st<-na.omit(pcadd_st)
pcadd_st$p_mean<-as.numeric(pcadd_st$p_mean)
pcadd_st$p_se<-as.numeric(pcadd_st$p_se)
pcadd_st$nums<-as.numeric(pcadd_st$nums)
pcadd_st$ID<-factor(pcadd_st$ID,levels = c("3_prime_UTR","5_prime_UTR","upstream_gene","downstream_gene",
                                             "intergenic","intron","non_coding_transcript",
                                             "mature_miRNA","splice_acceptor","splice_donor","stop_lost",
                                             "stop_retained","Synonymous","Missense"))
ggplot(pcadd_st, aes(ID, p_mean,color=ID))+ geom_point(size=2.5)+
  scale_color_manual(values=c("#FF0000FF","#FF9900FF","#FFCC00FF","#00FF00FF",
                                "#6699FFFF","#CC33FFFF","#999999FF",
                                "#FF00CCFF","#1B1919FF","#ADB6B6FF","#FDAF91FF",
                                "#42B540FF","#00468BFF","#0099B4FF"))+
    ### stat_summary(fun.data=mean_cl_normal,geom="pointrange", 
    ##size=1.0,position = position_dodge2(0.8))+
    geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), 
                  width=0.5,position = position_dodge(.8)) + labs(y="pCADD scores",x="") + 
    theme_classic()+
    theme(legend.position = "",
          axis.text.x = element_text(size = 17,color = "black",angle=270,hjust = 1.0),
          ##axis.ticks.x = element_blank(),
          ##axis.line.y  = element_blank(),
          axis.text.y = element_text(size = 17,color = "black"),
          axis.title.y= element_text(size = 17,colour="black"))
ggsave("noncoding_all.pcadd.png",width = 5,height =4,dpi=680)
ggsave("noncoding_all.pcadd.pdf",width = 5,height =4)
############################## Figure 4b and S4c ##############################
library(ggplot2)
library(data.table)
setwd("D:/vcf_for_gte/losfunction/11_noncoding/")
tmp<-read.csv("D:/vcf_for_gte/losfunction/11_noncoding/intron.all.enrichment2.txt",header = F,sep=" ")
tmp$V1<-factor(tmp$V1,levels=c("top0.10","top0.20","top0.30","top0.40","top0.50","top0.60"
                               ,"top0.70","top0.80","top0.90","top1.00"))

tmp$foldchange<-(tmp$V7/tmp$V6)/(tmp$V5/tmp$V4)

fc_st<-data.frame(percs=NA,regs=NA,p_mean=NA,p_se=NA)
for(j in unique(tmp$V1)){
  for(k in unique(tmp$V3)){
    tmp2<-subset(tmp,V1==j & V3== k,)
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
ggsave("all.noncoding_enrichment1.png",width = 9,height =3.5,dpi=680)
ggsave("all.noncoding_enrichment1.pdf",width = 9,height =3.5) 

fc_st1<-subset(fc_st,regs=="E1"|regs=="E12",)
ggplot(fc_st1, aes(percs, p_mean,group=regs,color=regs))+ geom_line(size=0.5) +
  geom_point(size=1.5)+
  ##scale_color_lancet()+
  scale_color_manual(values = c("Red","Chocolate1"))+
  geom_errorbar(aes(ymin=p_mean - p_se, ymax=p_mean + p_se), width=0.2) + 
  geom_hline(yintercept = 1.0,linetype="dashed")+
  labs(y="",x="") + 
  theme_classic()+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("all.noncoding_enrichment1.112.png",width = 9,height =3.5,dpi=680)
ggsave("all.noncoding_enrichment1.112.pdf",width = 9,height =3.5) 

###################### Figure 4c and S4d #####################################
setwd("D:/vcf_for_gte/losfunction/07_utr")
library(data.table)
library(ggplot2)
tmp<-fread("noncoding_cats_pcadd_enrich.1.txt")
tmp$V2<-factor(tmp$V2,levels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))
tmp$foldchang<-(tmp$V8/tmp$V7)/(tmp$V6/tmp$V5)
fc_st2<-subset(tmp,V2==0.1,)
fc_st3<-subset(fc_st2,V1 != "non_coding_transcript_exon_variant" & V1 != "non_coding_transcript_variant",)
fc_st3$V4<-factor(fc_st3$V4,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                       "E12","E13","E14","E15"))
ggplot(fc_st3, aes(V4, foldchang,group=V1,color=V1))+ 
  stat_summary(aes(V4, foldchang,group=V1),fun=mean,
               geom = "pointrange",position = position_dodge(0.5),size=0.3,
               fun.min = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.max = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  geom_hline(yintercept = 1.0,linetype="dashed")+
  scale_color_manual(values=c("#FF0000FF","#FF9900FF","#00FF00FF","#6699FFFF","#CC33FFFF","#FFCC00FF"))+
  labs(y="",x="") + 
  theme_classic()+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("all.noncoding_enrichment.png",width = 9,height =3.5,dpi=680)
ggsave("all.noncoding_enrichment.pdf",width = 9,height =3.5) 

fc_st4<-subset(fc_st3,V4=="E1" |V4=="E12",)
ggplot(fc_st4, aes(V4, foldchang,group=V1,color=V1))+ 
  stat_summary(aes(V4, foldchang,group=V1),size=0.1,fun=mean,
               geom = "pointrange",position = position_dodge(0.3),
               fun.min = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.max = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  scale_color_manual(values=c("#FF0000FF","#FF9900FF","#00FF00FF","#6699FFFF","#CC33FFFF","#FFCC00FF"))+
  geom_hline(yintercept = 1.0,linetype="dashed")+
  labs(y="",x="") + 
  theme_classic()+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("all.noncoding_enrichment.112.png",width = 4.5,height =3.5,dpi=680)
ggsave("all.noncoding_enrichment.112.pdf",width = 4.5,height =3.5) 

########################### Figure 4d and 4e ################################
library(ggplot2)
library(dplyr)
library(ggsci)
qtl<-c("sqtl","eqtl")
i<-"eqtl"
setwd(paste("D:/vcf_for_gte/losfunction/11_noncoding/",i,sep = ""))
file_dir<-paste("all.",i,"_noncoding_pcadd_reg1.txt",sep = "")
tmp<-read.csv(file_dir,header = F,sep = " ")
tmp$V7<-factor(tmp$V7,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                 "E12","E13","E14","E15"))
tmp$foldchang<-(tmp$V10/tmp$V9)/(tmp$V3/tmp$V4)
eqtl_egene_0.05<-c("Blastomere","Cartilage","Duodenum","Kidney","Lymph_node","Morula","Oocyte","Pituitary","Placenta")
tmp_e<-tmp[!(tmp$V1 %in% eqtl_egene_0.05),]

i<-"sqtl"
setwd(paste("D:/vcf_for_gte/losfunction/11_noncoding/",i,sep = ""))
file_dir<-paste("all.",i,"_noncoding_pcadd_reg1.txt",sep = "")
tmp<-read.csv(file_dir,header = F,sep = " ")
tmp$V7<-factor(tmp$V7,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                 "E12","E13","E14","E15"))
tmp$foldchang<-(tmp$V10/tmp$V9)/(tmp$V3/tmp$V4)
sqtl_egene_0.05<-c("Adipose","Blood","Brain","Embryo","Hypothalamus","Ileum","Jejunum","Large_intestine",
                   "Liver","Lung","Muscle","Ovary","Small_intestine","Spleen","Testis","Uterus")
tmp_s<-tmp[(tmp$V1 %in% sqtl_egene_0.05),]

tmp_e$cats<-"eQTL"
tmp_s$cats<-"sQTL"
tmp_c<-rbind(tmp_e,tmp_s)
fd_st_c<-data.frame(cats=NA,pcadd=NA,reg=NA,mean_fd=NA,sd_fd=NA)
for(i in unique(tmp_c$cats)){
  for(j in unique(tmp_c$V6)){
    for(k in unique(tmp_c$V7)){
      tmp_m<-subset(tmp_c,cats==i & V6 == j & V7 ==k ,)
      tmp_mean<-mean(tmp_m$foldchang)
      tmp_sd<-sd(tmp_m$foldchang)/sqrt(nrow(tmp_m))
      l<-c(i,j,k,tmp_mean,tmp_sd)
      fd_st_c<-rbind(fd_st_c,l)
    }
  }
}
fd_st_c$reg<-factor(fd_st_c$reg,levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11",
                                           "E12","E13","E14","E15"))
fd_st_c$pcadd<-as.factor(fd_st_c$pcadd)
fd_st_c<-na.omit(fd_st_c)
fd_st_c$mean_fd<-as.numeric(fd_st_c$mean_fd)
fd_st_c$sd_fd<-as.numeric(fd_st_c$sd_fd)
fd_st_c$pcadd<-as.factor(fd_st_c$pcadd)
ggplot(fd_st_c,aes(pcadd,mean_fd,fill=reg,color=reg)) +
  geom_point(size=1.0)+
  ##geom_errorbar(aes(ymin=mean_fd - sd_fd, ymax=mean_fd + sd_fd),width=0.4) + 
  labs(y="",x="")+
  geom_line(aes(pcadd,mean_fd,group=reg),size=0.25)+
  scale_color_manual(values = c("Red","DeepPink", "ForestGreen","Green3", "Green1","Yellow", "Goldenrod3","Gold2",
                                "Gold1", "LightGoldenrod1","SkyBlue", "Chocolate1","grey41", "grey61", "black"))+
  theme_classic()+
  labs(x="Enrichment",y="")+guides(alpha="none",)+
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black")
  )+facet_wrap(.~cats,nrow=1,strip.position = "top",scales = "fixed")

ggsave("all.reg.pcadd.noncoding_enrichment2.egene0.05.png",width = 8,height =3.5,dpi=680)
ggsave("all.reg.pcadd.noncoding_enrichment2.egene0.05.pdf",width = 8,height =3.5) 
fd_st_c1<-subset(fd_st_c,reg=="E1"|reg=="E12",)

ggplot(fd_st_c1,aes(pcadd,mean_fd,fill=cats,color=cats)) +
  geom_point(size=1.0)+
  geom_errorbar(aes(ymin=mean_fd - sd_fd, ymax=mean_fd + sd_fd),
                width=0.4) + 
  labs(y="",x="")+
  geom_line(aes(pcadd,mean_fd,group=cats),size=0.25)+
  theme_classic()+
  labs(x="Enrichment",y="")+guides(alpha="none",)+
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7,angle = 30),
        axis.text.y = element_text(size = 17,color = "black")
  )+facet_wrap(.~reg,nrow=1,strip.position = "top",scales = "fixed")

ggsave("E112.reg.pcadd.noncoding_enrichment2.egene0.05.png",width = 9,height =3.5,dpi=680)
ggsave("E112.reg.pcadd.noncoding_enrichment2.egene0.05.pdf",width = 9,height =3.5)

######################## Figure S4f #########################################
library(ggridges)
tmp_muscle<-subset(tmp_c,V1=="Muscle")
tmp_muscle$V6<-as.factor(tmp_muscle$V6)
ggplot(tmp_muscle,aes(foldchang,V7,fill=cats,color=cats,alpha=0.8)) + geom_density_ridges()+
  scale_x_continuous(expand = c(0.01, 0))+
  scale_y_discrete(expand = c(0,0),labels=c("TssA","TssAHet","TxFlnk","TxFlnkWk","TxFlnkHet","EnhA","EnhAMe","EnhAWk",
                                            "EnhAHet","EnhPois","ATAC_ls","TssBiv","Repr","ReprWk","Qui"))+
  theme_classic()+
  labs(x="Enrichment",y="")+guides(alpha="none",)+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(size = 17,color = "black",vjust=0.7),
        axis.text.y = element_text(size = 17,color = "black")
  )
ggsave("Muscle.reg.pcadd.noncoding_enrichment1.egene0.05.png",width = 5.5,height =8,dpi=680)
ggsave("Muscle.reg.pcadd.noncoding_enrichment1.egene0.05.pdf",width = 5.5,height =8) 
######################### Figure 4e ##########################################
setwd("D:/vcf_for_gte/losfunction/11_noncoding/gwas")
library(dplyr)
library(ggplot2)
fgwas<-read.csv("noncoding.fgwas.result.txt",header = F,sep = "\t")
fgwas<-fgwas[fgwas$V4!="<-20",]
f_max<-fgwas %>% group_by(V1) %>% slice_max(V5)
f_max$V4<-as.numeric(f_max$V4)
f_max$V6<-as.numeric(f_max$V6)
f_max$V9<-paste(f_max$V1,"(",f_max$V2,"-",f_max$V3,")",sep = "")
f_max$V10<-2^(f_max$V4)
f_max$V11<-2^(f_max$V5)
f_max$V12<-2^(f_max$V6)
ggplot(f_max,aes(x=V11, y=V9))+
  geom_point(size=1.5)+
  geom_segment(aes(x=V10,xend=V12,y=V9,yend=V9),cex=0.2)+
  geom_vline(xintercept=2,color="red",linetype="dashed")+
  scale_y_discrete(position = 'right')+
  labs(y=NULL)+theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(0, 35))

ggsave("noncoding.enrichment.13.gwas.png",width = 2.2,height =2,dpi=680)
ggsave("noncoding.enrichment.13.gwas.pdf",width = 2.2,height =2) 
############################ Figure 4f ######################################
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/11_noncoding/gwas")
snp1<-read.csv("ADG.noncoding.chr1_52100000-53450000.region",header=F,sep="\t")
snp1$colrs<-1
snp1[snp1$V5==52433942,]$colrs<-2
snp1[snp1$V5==52991009,]$colrs<-2
snp1[snp1$V5==52991100,]$colrs<-2
snp1[snp1$V5==52991128,]$colrs<-2
snp1[snp1$V5==52991154,]$colrs<-2
snp1[snp1$V5==52991210,]$colrs<-2
snp1[snp1$V5==52996873,]$colrs<-2
snp1[snp1$V5==52997058,]$colrs<-2
snp1[snp1$V5==52997128,]$colrs<-2
snp1[snp1$V5==53255886,]$colrs<-2
snp1[snp1$V5==53255888,]$colrs<-2
snp1[snp1$V5==53256516,]$colrs<-2
snp1[snp1$V5==53256959,]$colrs<-2
snp1[snp1$V5==53257017,]$colrs<-2
snp1[snp1$V5==53257274,]$colrs<-2
snp1[snp1$V5==53257383,]$colrs<-2

snp1$colrs<-as.factor(snp1$colrs)

p<-ggplot(snp1, aes(V5/1000000,-log10(V7)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ADB6B6FF","#ED0000FF"))+
  scale_size_manual(values = c(1,4))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~V1,scales = "free_y",ncol = 1)+
  coord_cartesian(xlim = c(52100000/1000000, 53450000/1000000))+
  geom_hline(yintercept = -log10(5e-08),linetype=2,colour="red")+
  theme(
    legend.position = "",
    axis.title.y = element_text(size = 17,hjust = 0.5),
    axis.title.x = element_text(size = 17,hjust = 0.5),
    axis.text.x = element_text(size = 17,hjust = 0.5,vjust=0.9),
    axis.text.y = element_text(hjust = 1,size = 17),
  )  



tmp3<-read.table("1_52100000",header = F,sep="\t")

split_b<-str_split(tmp3$V3,"_")
tmp3$pos<-as.numeric(sapply(split_b,"[",2))
tmp3$chr<-sapply(split_b,"[",1)

tmp3$GROUP<-paste(tmp3$V1,tmp3$V2,sep=":")
tmp3$colrs<-3
tmp3[tmp3$V1=="Liver" & tmp3$pos==52991128 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Liver" & tmp3$pos==52991154 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Liver" & tmp3$pos==52997058 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Liver" & tmp3$pos==52997128 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Liver" & tmp3$pos==52997178 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Liver" & tmp3$pos==53255888 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1
tmp3[tmp3$V1=="Brain" & tmp3$pos==53255888 & tmp3$V2 =="ENSSSCG00000004289",]$colrs<-1

tmp3[tmp3$V1=="Muscle" & tmp3$pos==52991128 & tmp3$V2 =="ENSSSCG00000004281",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==52991154 & tmp3$V2 =="ENSSSCG00000004281",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==52997058 & tmp3$V2 =="ENSSSCG00000004281",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==52997128 & tmp3$V2 =="ENSSSCG00000004281",]$colrs<-1
tmp3[tmp3$V1=="Muscle" & tmp3$pos==52997178 & tmp3$V2 =="ENSSSCG00000004281",]$colrs<-1
tmp3$colrs<-as.factor(tmp3$colrs)

p1<-ggplot(tmp3, aes(pos/1000000,-log10(V8)))+
  geom_point(aes(color=colrs,size=colrs),alpha=0.8)+
  scale_color_manual(values = c("#ED0000FF","#ADB6B6FF"))+
  scale_size_manual(values = c(4,1))+
  labs(y=expression(-log10(italic(P))),x=expression(Position(Mb))) + 
  theme_classic()+facet_wrap(.~GROUP,scales = "free_y",ncol = 1)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  coord_cartesian(xlim = c(52100000/1000000, 53450000/1000000))+
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



ggsave("ADG.noncoding.chr1_52100000-53450000.sad.png",plot=p3,width = 8,height =6,dpi=680)
ggsave("ADG.noncoding.chr1_52100000-53450000.sad.pdf",plot=p3,width = 8,height =6)
################################ Figure S4g ##################################
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/11_noncoding/gwas")
ase<-read.csv("ase.csv",header = T,sep=",")
ase$pos<-as.factor(ase$pos)
ase$size2<-"1"
ase[ase$FDR<0.05,]$size2<-"2"
ggplot(ase, aes(x=pos, y=ref_ratio)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill=size2,color=size2),
               binaxis = 'y', stackdir = 'center', dotsize = 1.0)+
  scale_fill_manual(values = c("blue","red"))+
  scale_color_manual(values = c("blue","red"))+
  theme_classic()+
  theme(
    legend.position = "",
    axis.title.y = element_text(size = 17,hjust = 0.5),
    axis.title.x = element_text(size = 17,hjust = 0.5),
    axis.text.x = element_text(size = 17,hjust = 0.5,vjust=0.5,angle = 30),
    axis.text.y = element_text(hjust = 1,size = 17),
  ) 

ggsave("ASE.noncoding.chr1_52100000-53450000.mai.png",width = 6,height =5,dpi=680)
ggsave("ASE.noncoding.chr1_52100000-53450000.mai.pdf",width = 6,height =5)