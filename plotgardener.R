library(plotgardener)
library(RColorBrewer)### display.brewer.all()
library(GenomicFeatures)
library(RMariaDB)
library(AnnotationForge)
library(biomaRt)
library(dplyr)
library(ggplot2)

## Create a plotgardener page
txdb <- makeTxDbFromEnsembl(organism="Sus scrofa", release=100)
mart <- useEnsembl("genes", "sscrofa_gene_ensembl", version=100)

genes <- 
  getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name")) |>
  as_tibble() |>
  dplyr::rename(
    GID="ensembl_gene_id",
    GENENAME="external_gene_name") |>
  dplyr::mutate(
    SYMBOL=case_when(
      GENENAME == "" ~ GID,
      GENENAME %in% unique(GENENAME[duplicated(GENENAME)]) ~ GID,
      TRUE ~ GENENAME),
    GENEID=GID)

# Forge the OrgDb.
makeOrgPackage(
  gene_info=genes,
  version="0.0.1",
  maintainer="############",
  author="############",
  outputDir=".",
  tax_id=9823,
  genus="Sus",
  species="scrofa"
)


install.packages("D:/ad_gte/org.Sscrofa.eg.db", repos=NULL, type="source")
library("org.Sscrofa.eg.db")

ss_assembly <- assembly(
  Genome="Sscrofa11.1",
  TxDb=txdb,
  gene.id.column="GENEID",
  OrgDb=org.Sscrofa.eg.db,
  display.column="SYMBOL"
)

setwd("E:/Hi-c")

region<-"7_29750000_31200000"
chr<-unlist(strsplit(region,"_"))[1]
chrstart<-as.numeric(unlist(strsplit(region,"_"))[2])
chrend<-as.numeric(unlist(strsplit(region,"_"))[3])

pdf("ADG.chr1_52100000-53450000.liver.sad.pdf",width = 12,height =4)
pageCreate(width = 7, height = 3, default.units = "inches")

hic_file<-readHic(file="EUD_LargeWhite_PEFs_mapq30_8106.hic",chrom=chr,assembly = ss_assembly,
                  chromstart = chrstart, chromend = chrend,resolution = 50000)
## Text section labelhttp://127.0.0.1:8189/graphics/plot_zoom_png?width=1200&height=900
plotText(label = "C", fontsize = 24,
         x = 0.1, y = 0.1, just = c("left", "top"), default.units = "inches")
## Set genomic and dimension parameters in a `params` object
params_c <- pgParams(chrom = chr, chromstart = chrstart, chromend = chrend, 
                     assembly = ss_assembly,
                     x = 0, width = 7.0, default.units = "inches")
params_d <- pgParams(chrom = chr, chromstart = chrstart, chromend = chrend, 
                     assembly = ss_assembly,
                     x = 0, width = 7.0, default.units = "inches")
## Plot Hi-C triangle
hic_gm <- plotHicTriangle(data = hic_file, params = params_c,
                           resolution = 50000,palette = colorRampPalette(brewer.pal(n = 9, "Reds")),
                          y = 2.5, height = 2.5, just = c("left", "bottom"))
## Annotate Hi-C heatmap legend
annoHeatmapLegend(plot = hic_gm, fontsize = 7, 
                  x = 2.5, y = 0.5, width = 0.5, height = 1.5, 
                  just = c("right", "top"), default.units = "inches")

pdf("BFT.chr7_29750000_31200000.muscle.gene.pdf",width = 12,height =4)
pageCreate(width = 9, height = 1, default.units = "inches")
## Plot genes
params_d <- pgParams(chrom = chr, chromstart = chrstart, chromend = chrend, 
                     assembly = ss_assembly,
                     x = 0, width = 9, default.units = "inches")
genes_gm <- plotGenes(params = params_d, stroke = 1,fontsize = 14,
                      fontcolor = c("#A9A9A9", "#000000"),
                      fill = c("#A9A9A9", "#000000"),
                      strandLabels = T,
                      y = 0, height = 1.0)
## Annotate genome label
annoGenomeLabel(plot = genes_gm, params = params_d, 
                scale = "Mb", fontsize = 14,y = 1.0)
pageGuideHide()

dev.off()


pdf("BFT.chr7_29750000_31200000.muscle.gene.pdf",width = 12,height =4)
pageCreate(width = 4, height = 1, default.units = "inches")
## Plot genes
params_d <- pgParams(chrom = chr, chromstart = chrstart, chromend = chrend, 
                     assembly = ss_assembly,
                     x = 0, width = 4, default.units = "inches")
genes_gm <- plotGenes(params = params_d, stroke = 1,fontsize = 14,
                      fontcolor = c("#A9A9A9", "#000000"),
                      fill = c("#A9A9A9", "#000000"),
                      strandLabels = T,
                      y = 0, height = 1.0)
## Annotate genome label
annoGenomeLabel(plot = genes_gm, params = params_d, 
                scale = "Mb", fontsize = 14,y = 1.0)
pageGuideHide()
dev.off()


