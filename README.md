# PigGTEx-WGS
Figure 1, 2, 3, 4, 5, 6, S1, S2, S3, S4, S5 and S6 were the origin R codes of plot for our papers.

WGS_pipeline.sh is the pipeline we deal with the WGS sample.

itools.sh and bam_stat.sh is the codes we statics the bam files by iTools Xamtools and samtools.

umap_data.py is the python codes to deal with the bi-allelic SNPs to run PCA, TSNE and UMAP.

select_ancestral_allel.py is the python codes that we identified the ancestral allele based on method of Bianco, E(A deep catalog of autosomal single nucleotide variation in the pig).

add_ancestral.vcf.sh is codes that we added the anceatral allele to VCF file with bcftools.

annotate_binary_snp was the codes to annotate the bi-allelic SNPs with CooVar, SNPEFF, VEP and pCADD.

generate_tf_select.sh was the codes to generate the TF based on ATAC-seq by Hommer2.
