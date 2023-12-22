#############################################################################
##################### Figure 1 and S1 #######################################
#############################################################################

####### Figure 1b and S1b (sex coverage and grouo)###########################
library(ggplot2)
library(patchwork)
library(ggsci)
library(scales)
setwd("D:/vcf_for_gte/losfunction/01_pca")
sex_data<-data.frame(group=c("new","new","database","database"),
                     cats=c("Male","Female","Male","Female"),nums=c(226,284,586,721))
cov_data<-data.frame(group=c(rep("new",5),rep("database",5)),
                     cats=c("<5X","5~10X","10~20X","20~30X",">=30X","<5X","5~10X","10~20X","20~30X",">=30X"),
                     nums=c(0,159,293,4,54,213,236,493,260,105))
pop_data<-data.frame(group=c(rep("new",5),rep("database",5)),
                     cats=c("ASD","ASW","EUD","EUW","SUI","ASD","ASW","EUD","EUW","SUI"),
                     nums=c(370,33,89,9,0,413,57,766,46,45))
cov_data$cats<-factor(cov_data$cats,levels=c("<5X","5~10X","10~20X","20~30X",">=30X"))
sex_data$group<-factor(sex_data$group,levels = c("new","database"))
cov_data$group<-factor(cov_data$group,levels = c("new","database"))
pop_data$group<-factor(pop_data$group,levels = c("new","database"))
p<-ggplot(sex_data,aes(cats,nums,fill=group))+geom_bar(stat = "identity",width = 0.3,color="black")+
  theme_classic()+ylab("")+xlab("Sex")+labs(title ="Sex Distribution")+
  scale_fill_manual(values=c("white","grey"),
                    label=c("Our study","Public"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1200))+
  ##geom_text(aes(label=nums,vjust=3,hjust=0.5,size=17,colour="#1B1919FF",fontface="bold"))+
  scale_colour_manual(values=c("white", "white"))+
  guides(colour="none",size="none",fill=guide_legend(title = NULL))+
  theme(legend.position = "top",
        axis.text.x = element_text(size = 17,color = "black"),
        axis.text.y = element_text(size = 17,color = "black"),
        axis.title.x = element_text(size = 17,color = "black"),
        plot.title = element_text(size = 17,hjust=0.5)
  )
p1<-ggplot(cov_data,aes(cats,nums,fill=group))+geom_bar(stat = "identity",width = 0.5,color="black")+
  theme_classic()+ylab("")+xlab("Depth")+labs(title ="Depth Distribution")+
  scale_fill_manual(values=c("white","grey"),
                    label=c("Existing","New generated"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,840))+
  ##geom_text(aes(label=nums,vjust=3,hjust=0.5,size=17,colour="#1B1919FF",fontface="bold"))+
  scale_colour_manual(values=c("white", "white"))+
  guides(colour="none",size="none",fill=guide_legend(title = NULL))+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle=30,vjust=0.6),
        axis.text.y = element_text(size = 17,color = "black"),
        axis.title.x = element_text(size = 17,color = "black"),
        plot.title = element_text(size = 17,hjust=0.5)
  )
p2<-ggplot(pop_data,aes(cats,nums,fill=group))+geom_bar(stat = "identity",width = 0.5,color="black")+
  theme_classic()+ylab("")+xlab("Population")+labs(title ="Population Composition")+
  scale_fill_manual(values=c("white","grey"),
                    label=c("Existing","New generated"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,920))+
  ##geom_text(aes(label=nums,vjust=2,hjust=0.5,size=17,colour="#1B1919FF",fontface="bold"))+
  scale_colour_manual(values=c("white", "white"))+
  guides(colour="none",size="none",fill=guide_legend(title = NULL))+
  theme(legend.position = "",
        axis.text.x = element_text(size = 17,color = "black",angle = 30,vjust=0.6),
        axis.text.y = element_text(size = 17,color = "black"),
        axis.title.x = element_text(size = 17,color = "black"),
        plot.title = element_text(size = 17,hjust=0.5)
  )
p2+p1+p
ggsave("sex_cov_pop.png",width = 18,height =6,dpi=680)
ggsave("sex_cov_pop.pdf",width = 18,height =6)

############### Figure S1a( indel length distribution) #######################
etwd("D:/vcf_for_gte/losfunction/indel_sta/GTE.1602.chr1_18.indel.out")
indel_len<-read.csv("indels.0.dat",header = T,sep="\t")
colnames(indel_len)<-c("length","counts")

ggplot(indel_len,aes(length,counts))+geom_bar(stat = "identity",width = 0.5)+
  theme_classic()+ylab("Count")+xlab("InDel Length")+
  scale_fill_manual(values="#42B540FF")+
  scale_y_continuous(expand = c(0,0))+xlim(-50,50)+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 16,color = "black"),
        axis.text.y = element_text(size = 16,color = "black"),
        axis.title = element_text(size = 16,color = "black"))
ggsave("indel_len.png",width = 8,height =6,dpi=680)
ggsave("indel_len.png.pdf",width = 8,height =6)

############## Figure 1c and S1c (TSNE) #####################################
setwd("D:/vcf_for_gte/01_pca")
library(ggplot2)
group.info<-read.csv("D:/vcf_for_gte/group_info/GTE_pop_clean_0715_five_population",header = T,sep = "\t")
nums1<-length(group.info$GROUP[group.info$GROUP=="ASD"])
nums2<-length(group.info$GROUP[group.info$GROUP=="ASW"])
nums3<-length(group.info$GROUP[group.info$GROUP=="EUD"])
nums4<-length(group.info$GROUP[group.info$GROUP=="EUW"])
nums5<-length(group.info$GROUP[group.info$GROUP=="SUI"])
for(i in 2:14){
    res<-read.csv(paste("tsne_umap_res/GTE.1602.pcs",i,".",meth,".csv",sep=""),header=T,sep = ",")
    colnames(res)<-c("ID","PC1","PC2")
    res_all<-merge(group.info,res,by="ID")
      ggplot(res_all, aes(x=PC1, y=PC2,color=GROUP)) +
        geom_point(size=1.5)+ theme_bw()+xlab("t-SNE-1") +ylab("t-SNE-2")+
        scale_color_manual(values = c("#FE7D7D","#FEC655","#75C4FE","#56B956","#DA70D6"),
                           labels=c(paste("ASD(",nums1,")",sep = ""),paste("ASW(",nums2,")",sep = ""),
                                    paste("EUD(",nums3,")",sep = ""),paste("EUW(",nums4,")",sep = ""),
                                    paste("SUI(",nums5,")",sep = "")))+
        theme_classic()+
        theme(legend.title=element_blank(),
              legend.position=c(0.5,1.0),legend.direction = "horizontal",
              legend.key.height=unit(0.3,"line"),
              legend.key.width=unit(0.3,"line"),
              axis.text.x = element_text(size = 18,color = "black"),
              axis.text.y = element_text(size = 18,color = "black"),
              axis.title = element_text(size = 18,color = "black"),
              legend.text = element_text(size = 13,color = "black",hjust=0))

    ggsave(paste(meth,"/GTE.1602.pcs",i,".tsne.",'png',sep =""),width = 14, units = "cm", height = 12, dpi = 680)
    ggsave(paste(meth,"/GTE.1602.pcs",i,".tsne.",'pdf',sep =""),width = 14, units = "cm", height = 12)
}
setwd("D:/vcf_for_gte/losfunction/01_pca/lof")
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

################# Figure S1d and 1d (pie for snp and indel)#################
library(eulerr)
library(patchwork)
library(ggplot2)
setwd("D:/vcf_for_gte/losfunction/02_TFBS")
snp_vd<-euler(c(Our_study=32432306,dbSNP=21535320,"Our_study&dbSNP"=35348566))
p<-plot(snp_vd,fills=list(fill=c("white","grey","grey")),
        labels=list(col="white",font=2,size=8),edeges=FALSE,quantities = TRUE)

indel_vd<-euler(c(Our_study=7812667,dbSNP=4265533,"Our_study&dbSNP"=931954))
p1<-plot(indel_vd,fills=list(fill=c("white","grey","grey")),
         labels=list(col="white",font=2,size=8),edeges=FALSE,quantities = TRUE)
p
p1
p2<-p+p1
ggsave("newsample.ven.snp.png",plot=p,width = 18, units = "cm", height = 12, dpi = 680)
ggsave("newsample.ven.snp.pdf",plot=p,width = 18, units = "cm", height =12)
ggsave("newsample.ven.indels.png",plot=p1,width = 18, units = "cm", height = 12, dpi = 680)
ggsave("newsample.ven.indels.pdf",plot=p1,width = 18, units = "cm", height =12)

snp_vd<-euler(c(Our_study=155648239,dbSNP=16360394,"Our_study&dbSNP"=40523492))
p<-plot(snp_vd,fills=list(fill=c("white","grey","grey")),
        labels=list(col="white",font=2,size=8),edeges=FALSE,quantities = TRUE)

indel_vd<-euler(c(Our_study=39332667,dbSNP=2956790,"Our_study&dbSNP"=2240697))
p1<-plot(indel_vd,fills=list(fill=c("white","grey","grey")),
         labels=list(col="white",font=2,size=8),edeges=FALSE,quantities = TRUE)
p+p1
ggsave("ven.snp.png",plot=p,width = 18, units = "cm", height = 12, dpi = 680)
ggsave("ven.snp.pdf",plot=p,width = 18, units = "cm", height =12)
ggsave("ven.indels.png",plot=p1,width = 18, units = "cm", height = 12, dpi = 680)
ggsave("ven.indels.pdf",plot=p1,width = 18, units = "cm", height =12)

########################Figure 1e and S1e ###################################
setwd("C:/vcf_for_gte/losfunction/binary_snp_anno")
library(ggplot2)
library(ggsci)
library(reshape2)
library(plotrix)
options(scipen=200)
#####vep_all #################################
snp_pos<-read.csv("gte_snp_anno_vep.csv",header = T)
snp_pos$GROUP<-factor(snp_pos$GROUP,levels = c("our","dbsnp"))
snp_pos$category<-factor(snp_pos$category,
                         levels=c("intergenic","intron","prime_5_UTR","prime_3_UTR","synonymous","missense","splice","nonsense"))

p<-ggplot(data=snp_pos,mapping=aes(x=category,y=nums1,fill=GROUP))+
  geom_bar(stat="identity",position = position_dodge(0.9))+
  labs(y=expression(paste("SNP Numbers(",10^3,")",sep="")),x="") + theme_classic()+
  ##scale_y_continuous(position = "right")+
  scale_fill_manual(values=c("#00468BFF","#42B540FF"))+
  scale_x_discrete(labels= c("Intergenic","Intron","5'UTR","3'UTR","Synonymous","Missense","Essential splice","Nonsense"))+
  theme(axis.text.x = element_text(size = 15,color = "black",angle = 30,hjust = 1),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 15,color = "black")
  )

p2<- gg.gap(plot = p,
            segments = list(c(50,130),c(2900,15000)),
            tick_width = c(30, 700, 30000),
            ylim = c(1, 99000),
            rel_heights=c(0.20,0,0.2,0,0.2)
)
ggsave("gte_snp_pos_vep_all.png", plot=p2,width =24, height =18, units = "cm",dpi = 680)
ggsave("gte_snp_pos_vep_all.pdf", plot=p2,width =24, height =18, units = "cm")

############vep for sus group################################################
snp_pos<-read.csv("gte_rmout_snp_anno_vep.csv",header = T)
snp_pos$GROUP<-factor(snp_pos$GROUP,levels = c("our","dbsnp"))
snp_pos$category<-factor(snp_pos$category,
                         levels=c("intergenic","intron","prime_5_UTR","prime_3_UTR","synonymous","missense","splice","nonsense"))
p<-ggplot(data=snp_pos,mapping=aes(x=category,y=nums1,fill=GROUP))+
  geom_bar(stat="identity",position = position_dodge(0.9))+ 
  labs(y=expression(paste("SNP numbers(",10^3,")",sep="")),x="") + theme_classic()+
  ##scale_y_continuous(position = "right")+
  scale_fill_manual(values=c("#00468BFF","#42B540FF"))+
  scale_x_discrete(labels= c("Intergenic","Intron","5'UTR","3'UTR","Synonymous","Missense","Essential splice","Nonsense"))+
  theme(axis.text.x = element_text(size = 15,color = "black",angle = 30,hjust = 1),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title.y = element_text(size = 15,color = "black")
  )

p2<- gg.gap(plot = p,
            segments = list(c(30,130),c(1300,14000)),
            tick_width = c(20, 400, 8000),
            ylim = c(1, 42000),
            rel_heights=c(0.20,0,0.2,0,0.2)
)
ggsave("gte_rmout_snp_pos_vep_all.png", plot=p2,width =24, height =18, units = "cm",dpi = 680)
ggsave("gte_rmout_snp_pos_vep_all.pdf", plot=p2,width =24, height =18, units = "cm")

###########################Figure S1f and S1g ################################
library(reshape2)
library(ggplot2)
library(patchwork)
setwd("C:/vcf_for_gte/losfunction/freq")
i<-"gte"
new_pi<-read.csv(paste(i,".new_pop.windowed.pi",sep=""),header=T,sep="\t")
old_pi<-read.csv(paste(i,".old_pop.windowed.pi",sep=""),header=T,sep="\t")
new_frq<-read.csv(paste("freq_stat.",i,".new_pop.frq1",sep=""),header=T,sep="\t")
old_frq<-read.csv(paste("freq_stat.",i,".old_pop.frq1",sep=""),header=T,sep="\t")

p<-ggplot(new_frq,aes(x=MAF))+geom_histogram(aes(y=..density..),binwidth = 0.01,
                                             color="#ADB6B6FF",fill="#00468BFF")+
  geom_density(alpha=0.8,fill="#00468BFF")+theme_classic()+
  geom_vline(xintercept = median(new_frq$MAF),color="#00468BFF",linetype="dashed")+
  theme_classic(base_size = 15)

p1<-ggplot(old_frq,aes(x=MAF))+geom_histogram(aes(y=..density..),binwidth = 0.01,
                                              color="#ADB6B6FF",fill="#42B540FF")+
  geom_density(alpha=0.8,fill="#42B540FF")+theme_classic()+
  geom_vline(xintercept = median(old_frq$MAF),color="#42B540FF",linetype="dashed")+
  theme_classic(base_size = 15)

p2<-ggplot(new_pi,aes(x=PI))+geom_histogram(aes(y=..density..),binwidth = 0.00001,
                                            color="#ADB6B6FF",fill="#00468BFF")+
  geom_density(alpha=0.8,fill="#00468BFF")+theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  geom_vline(xintercept = median(new_pi$PI),color="#00468BFF",linetype="dashed")+
  theme_classic(base_size = 15)

p3<-ggplot(old_pi,aes(x=PI))+geom_histogram(aes(y=..density..),binwidth = 0.00001,
                                            color="#ADB6B6FF",fill="#42B540FF")+
  geom_density(alpha=0.8,fill="#42B540FF")+theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  geom_vline(xintercept = median(old_pi$PI),color="#42B540FF",linetype="dashed")+
  theme_classic(base_size = 15)

q<-(p|p2)/(p1|p3)
median(new_frq$MAF)
median(old_frq$MAF)
median(new_pi$PI)
median(old_pi$PI)
ggsave("gte_pi_new_old.pdf",plot=q,width =16, height = 12, units = "cm")
ggsave("gte_pi_new_old.png",plot=q,width =16, height = 12, units = "cm",dpi = 680)


