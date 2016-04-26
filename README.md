JGM_QC
======

##**Pipeline for QC analysis of NGS libraries**

Current version 0.3, is only compatible for JAX clusters, has many dependencies including:

	fastqc/0.11.2
	bowtie2/2.0.6
	samtools/0.1.19
	tophat/2.0.7
	R/3.1.1
	rsem/1.2.12
	cufflinks/2.2.1 
	bwa/0.7.9a java/1.8.0_73 
	bismark/0.13.0
	trimmomatic/
	xenome/

Species supported as of current version:

	Human - hg19 (GRCh37), transcriptome version: ensembl v70
	Mouse - mm10 (GRCm38), transcriptome version: ensembl v74
	Rat   - rn5, transcriptome version: custom built

Pipelines available:

	mRNA - SE and PE
	mRNA_Unstranded
	mRNA_Tophat
	tRNA
	ATAC
	ChIP
	RRBS - NuGEN and custom
	WGS
	PDX_mRNA
	PDX_WGS
	
Metrics reported: