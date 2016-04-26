# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 19:04:14 2015

@author: mby

version 1 of JGM-QC freeze on November 17th 2015. This version will be tested for full functionality and
any changes necessary will be reflected in version 2.

"""

import os
import sys
import subprocess
import shlex
import argparse
import pysam
import time

#R1=''
#R2=''
#SPECIES=''
#LIBRARY_TYPE=''
samplename=''

'''ChiP-Seq'''

def chip_seq(R1,R2,SPECIES,output_folder):
    if R2 != None:
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
        status_align=run_bowtie_PE(trimR1,trimR2,SPECIES,output_folder)
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        stats=alignment_stats(status_sort),''

        remove_files(trimR1)
        remove_files(trimR2)
        remove_files(alignfile)
        remove_files(status_dups[2])        
    
    else:
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        status_align=run_bowtie_SE(trimR1,SPECIES,output_folder)
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert='',''
        stats=alignment_stats(status_sort),''
       
        remove_files(trimR1)
        remove_files(alignfile)
        remove_files(status_dups[2])
    
    return status_trim,status_align,status_dups,status_insert,stats
'''RNA-Seq Stranded'''

def mrna_stranded(R1,R2,SPECIES,folder):
    global samplename
    output_folder=folder
    status_trim=trimmer(output_folder,R1,R2)
    trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
    trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
    status_align=run_tophat(trimR1,trimR2,output_folder,SPECIES,'fr-firststrand')
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_rnametrics=rnaseq_metrics(output_folder,'fr-firststrand', SPECIES,status_sort)
    stats=alignment_stats(status_sort),''
    status_cufflinks='',''
    #cufflinks(output_folder,SPECIES,'fr-firststrand',status_sort)


    remove_files(trimR1)
    remove_files(trimR2)
    remove_files(alignfile)
    remove_files(status_dups[2])
    remove_files(output_folder+'unmapped.bam')
    remove_files(status_dups[2].replace('.bam','.bai'))
    
    return status_trim,status_align,status_dups,status_insert,status_rnametrics,status_cufflinks,stats

def mrna_stranded_RSEM(R1,R2,SPECIES,folder):
    global samplename
    
    if R2 != None:
        output_folder=folder
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
        status_align=rsem(trimR1,trimR2,output_folder,SPECIES,'fr-firststrand')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-firststrand', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''
    
        remove_files(trimR1)
        remove_files(trimR2)
        remove_files(alignfile)
        remove_files(status_dups[2])
        remove_files(alignfile.replace('genome','transcript'))
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.marked.bai'))        
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam.bai'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.bai'))
            

    else:
        output_folder=folder
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        status_align=rsem_se(trimR1,output_folder,SPECIES,'fr-firststrand')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-firststrand', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''
        status_insert='',''
        
        remove_files(trimR1)
        remove_files(alignfile)
        remove_files(alignfile.replace('genome','transcript'))
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.marked.bai'))        
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam.bai'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.bam.bai'))
        
    return status_trim,status_align,status_dups,status_insert,status_rnametrics,stats
###WORKING###

def pdx_mrna_stranded(R1,R2,SPECIES,folder):
    global samplename
    output_folder=folder
    status_trim=trimmer(output_folder,R1,R2)
    trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
    trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
    
    cmd='/data/mby/QC_PIPELINE/bin/xenome/1.0.0/xenome classify -T 12 -P /data/shared/research_pipelines_reference_data/human/RNA/UpdatedTransIndex/Index_with_Human_unplaced/trans_hg19_NOD_based_on_mm10_k25 --pairs --host-name mouse --graft-name human --output-filename-prefix %s%s -i %s -i %s' %(output_folder,samplename,trimR1,trimR2)
    cmd=shlex.split(cmd)    
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    classify_stats=open(output_folder+'classification_stats','w')

    for line in stdout:
        classify_stats.writelines(line)
    classify_stats.close()
    
    classify_R1=output_folder+samplename+'_human_1.fastq'
    classify_R2=output_folder+samplename+'_human_2.fastq'
    status_align=run_tophat(classify_R1,classify_R2,output_folder,SPECIES,'fr-firststrand')
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_rnametrics=rnaseq_metrics(output_folder,'fr-firststrand', SPECIES,status_sort)
    stats=alignment_stats(status_sort),''
    
    remove_files(classify_R1)
    remove_files(classify_R2)
    remove_files(alignfile)
    remove_files(status_dups[2])
    remove_files(status_dups[2].replace('.bam','.bai'))
    remove_files(classify_R1.replace('human','mouse'))
    remove_files(classify_R2.replace('human','mouse'))  
    remove_files(classify_R1.replace('human','ambiguous'))
    remove_files(classify_R2.replace('human','ambiguous'))
    remove_files(classify_R1.replace('human','both'))
    remove_files(classify_R2.replace('human','both'))
    remove_files(classify_R1.replace('human','neither'))
    remove_files(classify_R2.replace('human','neither'))
    remove_files(output_folder+'unmapped.bam')

    return status_trim,status_align,status_dups,status_insert,status_rnametrics,stats
        

'''RNA-seq Unstranded'''

def mrna_unstranded(R1,R2,SPECIES,folder):
    if R2 != None:
        output_folder=folder
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
        status_align=run_tophat_SE(trimR1,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        status_cufflinks=''
        #cufflinks(output_folder,SPECIES,strand,status_sort)
        stats=alignment_stats(status_sort),''
    
        remove_files(trimR1)
        remove_files(trimR2)
        remove_files(alignfile)
        remove_files(status_dups[2])
    
        return status_trim,status_align,status_dups,status_insert,status_rnametrics,status_cufflinks,stats
    else:
        output_folder=folder
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        status_align=run_tophat_SE(trimR1,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        status_cufflinks=''
        #cufflinks(output_folder,SPECIES,'fr-unstranded',status_sort)
        stats=alignment_stats(status_sort),''
    
        remove_files(trimR1)
        remove_files(alignfile)
        remove_files(status_dups[2])
        remove_files(status_dups[2].replace('.bam','.bam.bai'))        
    
        return status_trim,status_align,status_dups,status_insert,status_rnametrics,status_cufflinks,stats

def mrna_unstranded_RSEM(R1,R2,SPECIES,folder):
    global samplename
    output_folder=folder
    if R2 != None:
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        trimR2=output_folder+R2.split('/')[-1]+'_filtered_trimmed'
        status_align=rsem(trimR1,trimR2,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''

        remove_files(trimR1)
        remove_files(trimR2)
        remove_files(alignfile)
        remove_files(status_dups[2])
        remove_files(alignfile.replace('genome','transcript'))
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.marked.bai'))        
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam.bai'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.bam.bai'))
    else:
        status_trim=trimmer(output_folder,R1,R2)
        trimR1=output_folder+R1.split('/')[-1]+'_filtered_trimmed'
        status_align=rsem_se(trimR1,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''
        status_insert='',''
        
        remove_files(trimR1)
        remove_files(alignfile)
        remove_files(alignfile.replace('genome','transcript'))
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.marked.bai'))        
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam.bai'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.bam.bai'))
        
    return status_trim,status_align,status_dups,status_insert,status_rnametrics,stats
        
    
def mrna_unstranded_RSEM_trim(R1,R2,SPECIES,folder):
    global samplename
    output_folder=folder
    if R2 != None:
        paired1=output_folder+samplename+'_R1.fastq_filtered'
        paired2=output_folder+samplename+'_R2.fastq_filtered'
        unpair1=output_folder+samplename+'_R1.fastq_unpaired'
        unpair2=output_folder+samplename+'_R2.fastq_unpaired'
        
        cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
        PE %s %s \
        %s %s %s %s \
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' %(R1,R2,paired1,unpair1,paired2,unpair2)
        print "Trimming:", cmd
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout, stderr=proc.communicate()
        
        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
        status_trim=error_logs,out_logs
    
        status_align=rsem(paired1,paired2,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''
        
        remove_files(paired1)
        remove_files(paired2)
        remove_files(unpair1)
        remove_files(unpair2)
        remove_files(alignfile)
        remove_files(status_dups[2])
    else:
        paired1=output_folder+samplename+'_R1.fastq_filtered'
        
        cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
        SE %s\
        %s\
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
        ILLUMINACLIP:/data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/adapters/SmartSeq_Adapter.fa:2:30:10' %(R1,paired1)
        print "Trimming:", cmd
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout, stderr=proc.communicate()
        
        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
        status_trim=error_logs,out_logs    
        
        status_align=rsem_se(paired1,output_folder,SPECIES,'fr-unstranded')
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
        stats=alignment_stats(status_sort),''
        status_insert='',''
        
        remove_files(paired1)
        remove_files(alignfile)
        remove_files(alignfile.replace('genome','transcript'))
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.marked.bai'))        
        remove_files(alignfile.replace('genome.bam','transcript.sorted.bam.bai'))
        remove_files(alignfile.replace('genome.bam','genome.sorted.bam.bai'))
 
    
    return status_trim,status_align,status_dups,status_insert,status_rnametrics,stats
    
def mrna_unstranded_trim(R1,R2,SPECIES,folder):
    global samplename
    output_folder=folder
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
    ILLUMINACLIP:/data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print "Trimming:", cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    
    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs

    status_align=run_tophat(paired1,paired2,output_folder,SPECIES,'fr-unstranded')
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_rnametrics=rnaseq_metrics(output_folder,'fr-unstranded', SPECIES,status_sort)
    stats=alignment_stats(status_sort),''
    
    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])    
    
    return status_trim,status_align,status_dups,status_insert,status_rnametrics,stats

'''ATAC-Seq'''

def atac(R1,R2,SPECIES,output_folder):
    global samplename
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    
    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs
    
    #trimR1=paired1.replace('.fastq_filtered','.trim.fastq')
    #trimR2=paired2.replace('.fastq_filtered','.trim.fastq')
    
    #cmd='python /data/mby/QC_PIPELINE/bin/pyadapter_trim.py %s %s' %(paired1,paired2)
    #os.system(cmd)
    
    status_align=run_bowtie_PE(paired1,paired2,SPECIES,output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    stats=alignment_stats(status_sort),''

    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])

    return status_trim,status_align,status_dups,status_insert,stats

def atac_bwa(R1,R2,SPECIES,output_folder):
    global samplename
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    
    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs
    
    #trimR1=paired1.replace('.fastq_filtered','.trim.fastq')
    #trimR2=paired2.replace('.fastq_filtered','.trim.fastq')
    
    #cmd='python /data/mby/QC_PIPELINE/bin/pyadapter_trim.py %s %s' %(paired1,paired2)
    #os.system(cmd)
    
    status_align=run_bwa_PE(paired1,paired2,SPECIES,output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    stats=alignment_stats(status_sort),''

    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])

    return status_trim,status_align,status_dups,status_insert,stats
    
def rrbs(R1,R2,SPECIES,output_folder):
    if R2 is None:
        global samplename
        cmd='/data/mby/QC_PIPELINE/bin/trim_galore_zip/trim_galore --rrbs --path_to_cutadapt /data/mby/QC_PIPELINE/bin/cutadapt %s -o %s' %(R1,output_folder)
        print 'Trimmming:', cmd    
        
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()

        for files in os.listdir(output_folder):
            if files.endswith('trimmed.fq') or files.endswith('trimmed.fq.gz'):
                trimR1=os.path.join(output_folder,files)

        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
        status_trim=error_logs,out_logs
        status_align=run_bismark_SE(trimR1,SPECIES,"",output_folder)
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        stats=alignment_stats(alignfile),''

        cmd='java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar \
        INPUT=%s \
        OUTPUT=%s \
        SO=queryname' %(alignfile,alignfile.replace('.sam','.sorted.sam'))
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
        fileout=alignfile.replace('.sam','.querysorted.sam')

        methylcalls=bismark_methylation_extractor(R1,R2,output_folder,fileout,SPECIES)
        remove_files(trimR1)
        remove_files(trimR1+'_unmapped_reads.fq')       
        remove_files(trimR1+'_ambiguous_reads.fq')
        remove_files(alignfile)
        remove_files(status_dups[2])
        remove_files(fileout)
        
        cpgreport=fileout.replace('.sam','.CpG_report.txt')
        methylkit_report=cpgreport.replace('.txt','_MethylKit.txt')
        
        cmd = '''awk  '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print "" $1,$2,$3,$4/($4+$5),$4+$5;}' %s > %s''' %(cpgreport, methylkit_report)
        print cmd
        os.system(cmd)
        
        return status_trim,status_align,status_dups,methylcalls,stats
    
    else:
        global samplename
               
        cmd='/data/mby/QC_PIPELINE/bin/trim_galore_zip/trim_galore --paired --rrbs --path_to_cutadapt /data/mby/QC_PIPELINE/bin/cutadapt %s %s -o %s' %(R1,R2,output_folder)
        print 'Trimmming:', cmd            
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()

        for files in os.listdir(output_folder):
            if files.endswith('val_1.fq') or files.endswith('val_1.fq.gz'):
                trimR1=os.path.join(output_folder,files)
            elif files.endswith('val_2.fq') or files.endswith('val_2.fq.gz'):
                trimR2=os.path.join(output_folder,files)
        print trimR1, trimR2
        
        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
        status_trim=error_logs,out_logs
        status_align=run_bismark(trimR1,trimR2,SPECIES,"",output_folder)
        alignfile=status_align[2]
        status_sort=sort(alignfile)
        status_dups=duplicates(output_folder,status_sort)
        status_insert=insert_size(output_folder,status_sort)
        stats=alignment_stats(status_sort),''
        
        '''Resort aligned file by queryname required for methylation extractor input'''        
        
        cmd='java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar \
        INPUT=%s \
        OUTPUT=%s \
        SO=queryname' %(alignfile,alignfile.replace('.sam','.querysorted.sam'))
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
        fileout=alignfile.replace('.sam','.querysorted.sam')

        methylcalls=bismark_methylation_extractor(R1,R2,output_folder,fileout,SPECEIS)
        remove_files(trimR1)
        remove_files(trimR2)
        remove_files(alignfile)
        remove_files(status_dups[2])
        remove_files(fileout)

        cpgreport=fileout.replace('.sam','.CpG_report.txt')
        methylkit_report=cpgreport.replace('.txt','_MethylKit.txt')

        cmd = '''awk  '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print "" $1,$2,$3,$4/($4+$5),$4+$5;}' %s > %s''' %(cpgreport, methylkit_report)
        print cmd
        os.system(cmd)
        
        return status_trim,status_align,status_dups,status_insert,methylcalls,stats


def nugen_rrbs(R1,R2,INDEX,SPECIES,output_folder):
    global samplename

    cmd='/data/mby/QC_PIPELINE/bin/trim_galore_zip/trim_galore --paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC --path_to_cutadapt /data/mby/QC_PIPELINE/bin/cutadapt -o %s %s %s' %(output_folder,R1,R2)
    print 'Trimmming:', cmd   
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs
    
    for files in os.listdir(output_folder):
        if files.endswith('val_1.fq') or files.endswith('val_1.fq.gz'):
            trimR1=os.path.join(output_folder,files)
        elif files.endswith('val_2.fq') or files.endswith('val_2.fq.gz'):
            trimR2=os.path.join(output_folder,files)
    print trimR1, trimR2

    cmd='python /data/mby/QC_PIPELINE/bin/trimRRBSdiversityAdaptCustomers.py -1 %s -2 %s' %(trimR1,trimR2)
    print cmd
    os.system(cmd)
    
    divtrimR1=trimR1.replace('.fq','_trimmed.fq')
    divtrimR2=trimR2.replace('.fq','_trimmed.fq')
    
    status_align=run_bismark(divtrimR1,divtrimR2,SPECIES,"",output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    stats=alignment_stats(status_sort),''
    
    cmd='/data/mby/QC_PIPELINE/bin/strip_bismark_sam.sh %s' %(alignfile)
    print cmd
    os.system(cmd)
    
    stripped_samfile=alignfile.replace('.sam','.sam_stripped.sam')
    cmd='python /data/mby/QC_PIPELINE/bin/nudup.py -2 -f %s -o %s%s.nudup %s' %(INDEX,output_folder,samplename,stripped_samfile)
    print cmd
    os.system(cmd)
    
    sorted_dedup_bam=output_folder+samplename+'.nudup.sorted.dedup.bam'
    cmd='samtools sort -n  %s %s' %(sorted_dedup_bam,sorted_dedup_bam.replace('dedup.bam','query.sorted'))
    print cmd
    os.system(cmd)
    
    fileout=sorted_dedup_bam.replace('dedup.bam','query.sorted.bam')
    methylcalls=bismark_methylation_extractor(R1,R2,output_folder,fileout,SPECIES)
    remove_files(trimR1)
    remove_files(trimR2)
    remove_files(divtrimR1)
    remove_files(divtrimR2)
    remove_files(divtrimR1+'_unmapped_reads_1.fq')
    remove_files(divtrimR2+'_unmapped_reads_2.fq')
    remove_files(divtrimR1+'_ambiguous_reads_1.fq')
    remove_files(divtrimR2+'_ambiguous_reads_2.fq')
    remove_files(alignfile)
    remove_files(status_sort)
    remove_files(stripped_samfile)
    remove_files(status_dups[2])
    remove_files(sorted_dedup_bam)
    remove_files(sorted_dedup_bam.replace('dedup','markdup'))
    
    cpgreport=fileout.replace('.bam','.CpG_report.txt')
    methylkit_report=cpgreport.replace('.txt','_MethylKit.txt')

    cmd = '''awk  '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print "" $1,$2,$3,$4/($4+$5),$4+$5;}' %s > %s''' %(cpgreport, methylkit_report)
    print cmd
    os.system(cmd)

    return status_trim,status_align,status_dups,status_insert,methylcalls,stats

    
def bismark_methylation_extractor(R1,R2,outfolder,samin,SPECIES):
    
    genomefolder=get_bismark_index(SPECIES)
    
    if R2 is None:
        cmd='bismark_methylation_extractor --single-end --multicore 8 --ignore 0 --ignore_3prime 0 --bedgraph --comprehensive --cytosine_report --genome_folder %s --output %s %s' %(genomefolder,outfolder,samin)
        print 'Methylation Extractor:' ,cmd
        
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
            
        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line        
    else:
        cmd='bismark_methylation_extractor --paired-end --multicore 8 --no_overlap --ignore 0 --ignore_r2 0 --ignore_3prime 0 --ignore_3prime_r2 0 --bedgraph --comprehensive --cytosine_report --genome_folder %s --output %s %s' %(genomefolder,outfolder,samin)
        print 'Methylation Extractor:' ,cmd
        
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
            
        error_logs,out_logs='',''
        for line in stderr:
           error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
    
    return error_logs,out_logs


def wgs(R1,R2,SPECIES,output_folder):
    global samplename
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
    ILLUMINACLIP:/data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    
    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs
    
    status_align=run_bwa_PE(paired1,paired2,SPECIES,output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_wgs=wgs_metrics(output_folder,SPECIES,status_dups[2])
    stats=alignment_stats(status_sort),''

    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])
    
    return status_trim,status_align,status_dups,status_insert,status_wgs,stats


def wes(R1,R2,SPECIES,output_folder):
    global samplename
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
    ILLUMINACLIP:/data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()

    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs
    
    status_align=run_bwa_PE(paired1,paired2,SPECIES,output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_wgs=wgs_metrics(output_folder,SPECIES,status_dups[2])
    stats=alignment_stats(status_sort),''

    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])
    
    return status_trim,status_align,status_dups,status_insert,status_wgs,stats

def pdx_wgs(R1,R2,SPECIES,output_folder):

    global samplename
    paired1=output_folder+samplename+'_R1.fastq_filtered'
    paired2=output_folder+samplename+'_R2.fastq_filtered'
    unpair1=output_folder+samplename+'_R1.fastq_unpaired'
    unpair2=output_folder+samplename+'_R2.fastq_unpaired'
    
    cmd ='java -jar /data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/trimmomatic-0.32.jar\
    PE %s %s \
    %s %s %s %s \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
    ILLUMINACLIP:/data/mby/QC_PIPELINE/bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10' %(R1,R2,paired1,unpair1,paired2,unpair2)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    
    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    status_trim=error_logs,out_logs

    cmd='/data/mby/QC_PIPELINE/bin/xenome/1.0.0/xenome classify -T 12 -P /data/shared/research_pipelines_reference_data/human/RNA/UpdatedTransIndex/Index_with_Human_unplaced/trans_hg19_NOD_based_on_mm10_k25 --pairs --host-name mouse --graft-name human --output-filename-prefix %s%s -i %s -i %s' %(output_folder,samplename,paired1,paired2)
    cmd=shlex.split(cmd)    
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr=proc.communicate()
    classify_stats=open(output_folder+'classification_stats','w')

    for line in stdout:
        classify_stats.writelines(line)
    classify_stats.close()
    
    classify_R1=output_folder+samplename+'_human_1.fastq'
    classify_R2=output_folder+samplename+'_human_2.fastq'

    status_align=run_bwa_PE(classify_R1,classify_R2,SPECIES,output_folder)
    alignfile=status_align[2]
    status_sort=sort(alignfile)
    status_dups=duplicates(output_folder,status_sort)
    status_insert=insert_size(output_folder,status_sort)
    status_wgs=wgs_metrics(output_folder,SPECIES,status_dups[2])
    stats=alignment_stats(status_sort),''

    remove_files(paired1)
    remove_files(paired2)
    remove_files(unpair1)
    remove_files(unpair2)
    remove_files(alignfile)
    remove_files(status_dups[2])
    remove_files(classify_R1.replace('human','mouse'))
    remove_files(classify_R2.replace('human','mouse'))  
    remove_files(classify_R1.replace('human','ambiguous'))
    remove_files(classify_R2.replace('human','ambiguous'))
    remove_files(classify_R1.replace('human','both'))
    remove_files(classify_R2.replace('human','both'))
    remove_files(classify_R1.replace('human','neither'))
    remove_files(classify_R2.replace('human','neither'))

    
    return status_trim,status_align,status_dups,status_insert,status_wgs,stats

    
    
'''Fastqc'''

def fastqc(R1,R2,folder):
    if R2 != None:
        fqcfolder=folder+'fastqc/'
        os.mkdir(fqcfolder)
        cmd='fastqc -o %s %s %s' %(fqcfolder, R1, R2)
        print cmd
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
#        print 'fastqc:', cmd
    else:
        fqcfolder=folder+'/'+'fastqc/'
        os.mkdir(fqcfolder)
        cmd='fastqc -o %s %s' %(fqcfolder, R1)
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
        print 'fastqc:', cmd


'''Trimming and Filtering Reads'''
    
def trimmer(odir,R1,R2):
    if R2 != None:
        cmd='python /data/mby/QC_PIPELINE/bin/filter_trim.py -d %s %s %s' %(odir,R1,R2)
    #    cmd=shlex.split(cmd)
    #    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE, shell=True)
        print 'Trimmed:', cmd
        os.system(cmd)
    #    stdout, stderr = proc.communicate()
        error_logs,out_logs='',''    
    #    for line in stderr:
    #        print line
    #        error_logs=error_logs+line
    #    for line in stdout:
    #        print line
    #        out_logs=out_logs+line
        
    else:
        cmd='python /data/mby/QC_PIPELINE/bin/filter_trim_can_handle_SE.py -d %s %s' %(odir,R1)
        print 'Trimmed:', cmd
        error_logs,out_logs='',''  
        os.system(cmd)
        
    return error_logs,out_logs

'''Alignment Functions'''
    
def run_tophat(R1,R2,output_folder,SPECIES,strand):
    index=get_index_path(SPECIES)
    cmd='tophat -p 16 -z 0 -o %s -a 6 -m 2 --min-intron-length 70 -g 20 --library-type %s --transcriptome-index %s -n 2 %s %s %s' %(output_folder,strand,index[1],index[0],R1,R2)
    cmd=shlex.split(cmd)
    print 'Tophat:', cmd
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    bamfile=output_folder+'/'+'accepted_hits.bam'
    global samplename
    renamebam=output_folder+samplename+'.bam'
    cmd='mv %s %s' %(bamfile,renamebam) 
    os.system(cmd)
    return error_logs,out_logs,renamebam
    
def run_tophat_SE(R1,output_folder,SPECIES,strand):
    global samplename
    index=get_index_path(SPECIES)
    cmd='tophat -p 16 -z 0 -o %s -a 6 -m 2 --min-intron-length 70 -g 20 --library-type %s --transcriptome-index %s -n 2 %s %s' %(output_folder,strand,index[1],index[0],R1)
    cmd=shlex.split(cmd)
    print 'Tophat:', cmd
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    bamfile=output_folder+'/'+'accepted_hits.bam'
    renamebam=output_folder+samplename+'.bam'
    cmd='mv %s %s' %(bamfile,renamebam) 
    os.system(cmd)
    return error_logs,out_logs,renamebam
    
def rsem(R1,R2,output_folder,SPECIES,strand):
    error_logs,out_logs='',''
    index=get_rsem_index(SPECIES)
    
    if strand == 'fr-firststrand':
        strand='0'
    elif strand == 'fr-unstranded':
        strand='0.5'
    
    global samplename
    renamebam=output_folder+samplename
    cmd='rsem-calculate-expression --paired-end --output-genome-bam -p 16 --forward-prob %s --bowtie2 %s %s %s %s' %(strand,R1,R2,index[0],renamebam)
    print 'RSEM:', cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    return error_logs,out_logs,renamebam+'.genome.bam'
    
def rsem_se(R1,output_folder,SPECIES,strand):
    error_logs,out_logs='',''
    index=get_rsem_index(SPECIES)
    
    if strand == 'fr-firststrand':
        strand='0'
    elif strand == 'fr-unstranded':
        strand='0.5'
    
    global samplename
    renamebam=output_folder+samplename
    cmd='rsem-calculate-expression --fragment-length-mean 280 --fragment-length-sd 50 --seed-length 25 --output-genome-bam -p 16 --forward-prob %s --bowtie2 %s %s %s' %(strand,R1,index[0],renamebam)
    print 'RSEM:', cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    return error_logs,out_logs,renamebam+'.genome.bam'
    

def run_bowtie_PE(R1,R2,SPECIES,output_folder):
    index=get_bowtie_index(SPECIES)
    global samplename
    alignfile=output_folder+samplename+'.sam'
    
    cmd='bowtie2 -p 8 -x %s -1 %s -2 %s --no-unal -S %s' %(index[0],R1,R2,alignfile)
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    print 'Bowtie2:', cmd
    return error_logs,out_logs,alignfile
    
def run_bwa_PE(R1,R2,SPECIES,output_folder):
    index=get_bwa_index(SPECIES)
    global samplename
    alignfile=output_folder+samplename+'.sam'
    
    cmd='bwa mem -M -t 8 %s %s %s > %s' %(index[0],R1,R2,alignfile)
    os.system(cmd)
    error_logs,out_logs='',''    
    print 'BWA:', cmd
    return error_logs,out_logs,alignfile

def run_bowtie_SE(R1,SPECIES,output_folder):
    index=get_bowtie_index(SPECIES)
    global samplename
    alignfile=output_folder+samplename+'.sam'
    
    cmd='bowtie2 -p 8 -x %s -U %s --no-unal -S %s' %(index[0],R1,alignfile)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    print 'Bowtie2:', cmd
    return error_logs,out_logs,alignfile

def run_bismark(R1,R2,SPECIES,directionality,output_folder):
    index=get_bismark_index(SPECIES)
    samfile=R1+'_bismark_bt2_pe.sam'

    cmd='bismark --bowtie2 -p 16 %s --ambiguous --unmapped %s -1 %s -2 %s -o %s' %(directionality,index,R1,R2,output_folder)
    print 'Align:', cmd    
    
    cmd=shlex.split(cmd)    
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    
    return error_logs,out_logs, samfile

def run_bismark_SE(R1,SPECIES,directionality,output_folder):
    index=get_bismark_index(SPECIES)
    samfile=R1+'_bismark_bt2.sam'

    cmd='bismark --bowtie2 -p 16 %s --ambiguous --unmapped %s %s -o %s' %(directionality,index,R1,output_folder)
    print 'Align:', cmd    
    
    cmd=shlex.split(cmd)    
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    error_logs,out_logs='',''
    for line in stderr:
       error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    
    return error_logs,out_logs, samfile    

'''Run Metrics'''

def insert_size(output_folder,marked_file):
    if marked_file.endswith('.sam'):
        pdf=marked_file.replace('.sam','.insert-sizes.pdf')
        output=marked_file.replace('.sam','.insert-size.metrics')
    elif marked_file.endswith('.bam'):
        pdf=marked_file.replace('.bam','.insert-sizes.pdf')
        output=marked_file.replace('.bam','.insert-size.metrics')
    
    cmd= 'java -jar /opt/compsci/picard/1.95/CollectInsertSizeMetrics.jar \
    INPUT=%s \
    HISTOGRAM_FILE=%s\
    OUTPUT=%s ' %(marked_file,pdf,output)
    print 'InsertSize:', cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    return error_logs,out_logs

def duplicates(output_folder,alignfile):
    global LIBRARY_TYPE
#    if LIBRARY_TYPE == 'ATAC':
#        if alignfile.endswith('.sam'):
#            marked_file=alignfile.replace('.sam','.marked.sam')
#            output=alignfile.replace('.sam','.duplicates.metrics')
#        elif alignfile.endswith('.bam'):
#            marked_file=alignfile.replace('.bam','.marked.bam')
#            output=alignfile.replace('.bam','.duplicates.metrics')
#        cmd='java -jar /opt/compsci/picard/1.95/MarkDuplicates.jar \
#        INPUT=%s \
#        OUTPUT=%s \
#        CREATE_INDEX=TRUE \
#        REMOVE_DUPLICATES=true \
#        METRICS_FILE=%s' %(alignfile,marked_file,output)
#        cmd=shlex.split(cmd)
#        print cmd
#        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
#        stdout, stderr = proc.communicate()
#        error_logs,out_logs='',''    
#        for line in stderr:
#            error_logs=error_logs+line
#        for line in stdout:
#            out_logs=out_logs+line
#        return error_logs,out_logs,marked_file    
#    else:
        
    if alignfile.endswith('.sam'):
        marked_file=alignfile.replace('.sam','.marked.sam')
        output=alignfile.replace('.sam','.duplicates.metrics')
    elif alignfile.endswith('.bam'):
        marked_file=alignfile.replace('.bam','.marked.bam')
        output=alignfile.replace('.bam','.duplicates.metrics')
    cmd='java -jar /opt/compsci/picard/1.95/MarkDuplicates.jar \
    INPUT=%s \
    OUTPUT=%s \
    CREATE_INDEX=TRUE \
    REMOVE_DUPLICATES=TRUE \
    METRICS_FILE=%s' %(alignfile,marked_file,output)
    print 'MarkDuplicates:',cmd
    cmd=shlex.split(cmd)

    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    return error_logs,out_logs,marked_file

def wgs_metrics(output_folder,SPECIES,alignfile):
    refFasta=get_reference(SPECIES)[0]

    if alignfile.endswith('.sam'):
        output=alignfile.replace('marked.sam','.wgs.metrics')
    elif alignfile.endswith('.bam'):
        output=alignfile.replace('marked.bam','.wgs.metrics')
    cmd='java -jar /data/mby/QC_PIPELINE/bin/picard-tools-2.1.0/picard.jar CollectWgsMetrics \
    INPUT=%s \
    OUTPUT=%s \
    REFERENCE_SEQUENCE=%s' %(alignfile,output,refFasta)
    print 'WGS Metrics:',cmd
    cmd=shlex.split(cmd)

    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line
    return error_logs,out_logs

def rnaseq_metrics(output_folder, STRAND,SPECIES,marked_file):
    
    index=get_reference(SPECIES)
    rrna_interval_file=get_rrna_interval(SPECIES)
    ref_flat=get_ref_flat(SPECIES)

    if marked_file.endswith('.sam'):
        pdf=marked_file.replace('.sam','.rnaseq.metrics.pdf')
        output=marked_file.replace('.sam','.rnaseq.metrics')
    elif marked_file.endswith('.bam'):
        pdf=marked_file.replace('.bam','.rnaseq.metrics.pdf')
        output=marked_file.replace('.bam','.rnaseq.metrics')

    
    if STRAND == 'fr-firststrand':
        strand='SECOND_READ_TRANSCRIPTION_STRAND'
    elif STRAND == 'fr-unstranded':
        strand='NONE'
    
    if SPECIES == 'Rat' or SPECIES == 'Mouse' or SPECIES == 'HumanArray':
        cmd='java -jar /opt/compsci/picard/1.95/CollectRnaSeqMetrics.jar \
        REF_FLAT=%s \
        STRAND_SPECIFICITY=%s \
        CHART_OUTPUT=%s \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
        INPUT=%s \
        OUTPUT=%s' %(ref_flat,strand,pdf,marked_file,output)
        print 'RNASeq-Metrics:', cmd
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        error_logs,out_logs='',''    
        for line in stderr:
            error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line

    else:
        
        cmd='java -jar /opt/compsci/picard/1.95/CollectRnaSeqMetrics.jar \
        REF_FLAT=%s \
        STRAND_SPECIFICITY=%s \
        RIBOSOMAL_INTERVALS=%s\
        CHART_OUTPUT=%s \
        METRIC_ACCUMULATION_LEVEL=ALL_READS \
        INPUT=%s \
        OUTPUT=%s \
        REFERENCE_SEQUENCE=%s' %(ref_flat,strand,rrna_interval_file,pdf,marked_file,output,index[-1])
        print 'RNASeq-Metrics:', cmd        
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        error_logs,out_logs='',''    
        for line in stderr:
            error_logs=error_logs+line
        for line in stdout:
            out_logs=out_logs+line
    
    return error_logs, out_logs

def alignment_stats(alignfile):
    obj=ParseBAM(alignfile)
    stats=obj.stat()
    return stats

class ParseBAM:
    multi_hit_tags=['H0','H1','H2','IH','NH']
    qcut=30
    def __init__(self,inputFile):
        try:
            self.samfile = pysam.Samfile(inputFile,'rb')
            if len(self.samfile.header) ==0:
                print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                sys.exit(1)
            self.bam_format = True
            self.stats=[]
        except:
            self.samfile = pysam.Samfile(inputFile,'r')
            if len(self.samfile.header) ==0:
                print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                sys.exit(1)
            self.bam_format = False
            self.stats=[]
    def stat (self,q_cut=30):
        R_total=0
        R_qc_fail=0
        R_duplicate=0
        R_nonprimary=0
        R_unmap =0
        R_multipleHit=0
        R_uniqHit=0	#all the following count should be sum to uniqHit
        R_read1=0
        R_read2=0
        R_reverse =0
        R_forward=0
        R_nonSplice=0
        R_splice=0
        R_properPair =0 
        R_pair_diff_chrom = 0
        R_mitochondria = 0
        R_X=0
        R_Y=0
        try:
            while(1):
                flag=0
                aligned_read = self.samfile.next()
                R_total +=1
                if aligned_read.is_qcfail:			#skip QC fail read
                    R_qc_fail +=1
                    continue
                if aligned_read.is_duplicate:		#skip duplicate read
                    R_duplicate +=1
                    continue
                if aligned_read.is_secondary:		#skip non primary hit
                    R_nonprimary +=1
                    continue
                if aligned_read.is_unmapped:		#skip unmap read
                    R_unmap +=1
                    continue		
                if aligned_read.mapq < q_cut:
                    R_multipleHit +=1
                    continue						#skip multiple map read				
                if aligned_read.mapq >= q_cut:
                    R_uniqHit +=1
                if aligned_read.is_read1:
                    R_read1 +=1
                if aligned_read.is_read2:
                    R_read2 +=1
                if aligned_read.is_reverse:
                    R_reverse +=1
                else:
                    R_forward +=1
                if aligned_read.is_proper_pair:
                    R_properPair +=1
                    R_read1_ref = self.samfile.getrname(aligned_read.tid)
                    R_read2_ref = self.samfile.getrname(aligned_read.rnext)
                    if R_read1_ref != R_read2_ref:
                        R_pair_diff_chrom +=1
                else:
                    pass
                
#                if aligned_read.is_paired:
#                    if aligned_read.is_secondary:
#                        pass
#                    else:
#                        if aligned_read.mapq>=q_cut:
#                            if self.samfile.getrname(aligned_read.tid) == 'chrM' or self.samfile.getrname(aligned_read.tid) == 'MT':
#                                R_mitochondria+=1
#                            elif self.samfile.getrname(aligned_read.tid) == 'chrX' or self.samfile.getrname(aligned_read.tid) == 'X':
#                                R_X+=1
#                            elif self.samfile.getrname(aligned_read.tid) == 'chrY' or self.samfile.getrname(aligned_read.tid) == 'Y':
#                                R_Y+=1
#
                if self.samfile.getrname(aligned_read.tid) == 'chrM' or self.samfile.getrname(aligned_read.tid) == 'MT':
                    R_mitochondria+=1
                elif self.samfile.getrname(aligned_read.tid) == 'chrX' or self.samfile.getrname(aligned_read.tid) == 'X':
                    R_X+=1
                elif self.samfile.getrname(aligned_read.tid) == 'chrY' or self.samfile.getrname(aligned_read.tid) == 'Y':
                    R_Y+=1

        except StopIteration:
            pass
        line='Total_Alignments:\t%d\nQC_Fail:\t%d\nDuplicates:\t%d\nSecondaryAlignments:\t%d\nMultipleAlignments:\t%d\nUniqueAlignments:\t%d\nProperPair:\t%d\nMitochondria:\t%d\nchrX:\t%d\nchrY:\t%d\n' \
        %(R_total,R_qc_fail,R_duplicate,R_nonprimary,R_multipleHit,R_uniqHit,R_properPair,R_mitochondria,R_X,R_Y)
        return line

'''Necessary Files'''

'''[species]=(<BOWTIE_GENOME_INDEX>,<BOWTIE_TRANSCRIPTOME_INDEX>,<GFF>,<GENOME.fa>)'''

def get_index_path(SPECIES):
    indexes={'Human':('/data/mby/INDEXES/HUMAN/male.spikeins.hg19','/data/mby/INDEXES/HUMAN/gencode.v19.annotation.tRNA.male.spikeins','/data/mby/INDEXES/HUMAN/gencode.v19.annotation.tRNA.male.spikeins.gff','/data/mby/INDEXES/HUMAN/RSEM_INDEX/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ERCC92.fa'),
             'Mouse':('/data/mby/INDEXES/MOUSE/C57BL6J','<T-INDEX>','<GFF>','/data/mby/INDEXES/MOUSE/C57BL6J.fa'),
             'Rat':('/data/mby/INDEXES/RAT/rn5.ERCC92','<T-INDEX>','<GFF>','/data/mby/INDEXES/RAT/rn5.ERCC92.fasta')
             }
    return indexes[SPECIES]

def get_bowtie_index(SPECIES):
    indexes={'Mouse':('/data/mby/INDEXES/MOUSE/BOWTIE_INDEX/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed',),
             'Human':('/data/mby/INDEXES/HUMAN/male.spikeins.hg19',),
            }
    return indexes[SPECIES]

def get_reference(SPECIES):
    global LIBRARY_TYPE
    if 'mRNA_Tophat' in LIBRARY_TYPE or 'PDX-mRNA' in LIBRARY_TYPE:
        reference={'Human':('/data/mby/INDEXES/HUMAN/male.spikeins.hg19.fa',),
                   'Mouse':('/data/mby/INDEXES/MOUSE/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92.fasta',)
                   }
        return reference[SPECIES]
    elif 'WGS' in LIBRARY_TYPE:
        reference={'Human':('/data/mby/INDEXES/HUMAN/BWA/hg19.fa',),
                   'Mouse':('/data/mby/INDEXES/MOUSE/BWA_INDEX/mm9.fa',),
                   'Rat':('/data/mby/INDEXES/RAT/rn5.ERCC92.fasta',),
                   'Yeast':('/data/mby/INDEXES/YEAST/genome.fa',),
                   }
        return reference[SPECIES]

    else:
        reference={'Human':('/data/mby/INDEXES/HUMAN/RSEM_INDEX/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ERCC92.fa',),
                   'Mouse':('/data/mby/INDEXES/MOUSE/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92.fasta',),
                   'Rat':('/data/mby/INDEXES/RAT/rn5.ERCC92.fasta',),
                   'Yeast':('/data/mby/INDEXES/YEAST/genome.fa',),
                   'HumanArray':('/data/mby/INDEXES//HUMAN/RSEM_INDEX_WITH_ARRAY_CONTROL/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ArrayControl.fa')
                   }
        return reference[SPECIES]

def get_bismark_index(SPECIES):
    indexes={'Human':'/data/mby/INDEXES/HUMAN/BISMARK_INDEX/',
             'Mouse':'/data/shared/research_pipelines_reference_data/mouse/Bisulphite/Bismark_Index/WholeGenomeFasta/',
             'Yeast':'/data/mby/INDEXES/YEAST/BISMARK_INDEX/',
             'HumanYeast':'/data/mby/INDEXES/OTHER/BISMARK_INDEX_sacCer3_hg19/'
             }
    return indexes[SPECIES]

def get_ref_flat(SPECIES):
    ref_flat={'Human':'/data/mby/INDEXES/HUMAN/gencode_v19_refFlat_UCSC.RefFlat',
              'Mouse':'/data/mby/INDEXES/MOUSE/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92.refFlat',
              'Rat':'/data/mby/INDEXES/RAT/rn5.ERCC92.added3genes.refFlat',
              'HumanArray':'/data/mby/INDEXES/HUMAN/gencode_v19_refFlat_UCSC.RefFlat'
              }
    return ref_flat[SPECIES]

def get_rsem_index(SPECIES):
    indexes={'Human':('/data/mby/INDEXES/HUMAN/RSEM_INDEX/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ERCC92',),
             'Mouse':('/data/mby/INDEXES/MOUSE/RSEM_INDEX/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92',),
             'Rat':('/data/mby/INDEXES/RAT/RSEM_INDEX/rn5.ERCC92.added3genes',),
             'HumanArray':('/data/mby/INDEXES//HUMAN/RSEM_INDEX_WITH_ARRAY_CONTROL/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added..ArrayControl',),
             }
    return indexes[SPECIES]
    
def get_bwa_index(SPECIES):
    indexes={'Human':('/data/mby/INDEXES/HUMAN/BWA/hg19.fa',),
             'Mouse':('/data/mby/INDEXES/MOUSE/BWA_INDEX/mm9.fa',)}
    return indexes[SPECIES]
    

def get_rrna_interval(SPECIES):
    global LIBRARY_TYPE
    if 'Tophat' in LIBRARY_TYPE or 'PDX' in LIBRARY_TYPE:
        intervals={'Human':'/data/mby/INDEXES/HUMAN/Ribosomal_Interval_List_HEADER',
                   'Mouse':'/data/mby/INDEXES/MOUSE/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92.rRNA.intervals',
                   'Rat':'/data/mby/INDEXES/RAT/rn5.ERCC92.rRNA.intervals',
                   }
        return intervals[SPECIES]        
    else:
        intervals={'Human':'/data/mby/INDEXES/HUMAN/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ERCC92.ribosomal.intervals',
                   'Mouse':'/data/mby/INDEXES/MOUSE/Mus_musculus_GRCm38_74_dna_primary_assembly_patch_removed.ERCC92.rRNA.intervals',
                   'Rat':'/data/mby/INDEXES/RAT/rn5.ERCC92.added3genes.rRNA.intervals',
                   'HumanArray':'/data/mby/INDEXES/HUMAN/RSEM_INDEX_WITH_ARRAY_CONTROL/Homo_sapiens_GRCh37_70_dna_primary_assembly_chr_added.ArrayControl.rRNA.intervals',
                   }
        return intervals[SPECIES]

'''mRNA Seq Counts'''
    
#==============================================================================
# def htseq_count(output_folder, SPECIES):
#     inputbam=output_folder+'/'+output_folder.split('/')[-1]+'.marked.bam's
#     index=get_index_path(SPECIES)
#     cmd='htseq-count -f bam %s %s > %s' %(inputbam, index[2],inputbam.replace('.bam','.gene_counts'))
#     cmd=shlex.split()
#     proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
#     proc.communicate()
#     
#==============================================================================

def cufflinks(output_folder,SPECIES,strand,status_sort):
    index=get_index_path(SPECIES)
    refGTF=index[-2]
    refFasta=get_reference(SPECIES)
    cmd='cufflinks -I 30000 -p 12 -g %s --max-bundle-frags 10000000000000 --library-type %s -o %s %s' %(refGTF,strand,output_folder,status_sort)
    print cmd
    cmd=shlex.split(cmd)
    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    proc.communicate()
    error_logs,out_logs='',''    
    for line in stderr:
        error_logs=error_logs+line
    for line in stdout:
        out_logs=out_logs+line

    return error_logs, out_logs

'''Generate Logs'''
    
def logs(output_folder,logs):
    logout=output_folder+'.pipeline.olog'
    logerr=output_folder+'.pipeline.elog'
    logstats=output_folder+'alignment.stats'
    print logout, logerr
    lo=open(logout,'w')
    le=open(logerr,'w')
    ls=open(logstats,'w')
    for i in logs:
        lo.writelines('################BEGIN QC PIPELINE STDOUT#####################\n')
        for line in i[1]:
            lo.writelines(line)
        lo.writelines('\n\n\n\n\n')
    lo.close()
    
    for i in logs:
        le.writelines('################BEGIN QC PIPELINE STDERR#####################\n')
        for line in i[0]:
            le.writelines(line)
        le.writelines('\n\n\n\n\n')
    le.close()
    
    ls.writelines(logs[-1][0])
    ls.close()
    
    
def logs(output_folder,logs):
    logout=output_folder+'.pipeline.olog'
    logerr=output_folder+'.pipeline.elog'
    logstats=output_folder+'alignment.stats'
    print logout, logerr
    lo=open(logout,'w')
    le=open(logerr,'w')
    ls=open(logstats,'w')
    for i in logs:
        lo.writelines('################BEGIN QC PIPELINE STDOUT#####################\n')
        for line in i[1]:
            lo.writelines(line)
        lo.writelines('\n\n\n\n\n')
    lo.close()
    
    for i in logs:
        le.writelines('################BEGIN QC PIPELINE STDERR#####################\n')
        for line in i[0]:
            le.writelines(line)
        le.writelines('\n\n\n\n\n')
    le.close()
    
    ls.writelines(logs[-1][0])
    ls.close()
    
    
def load_modules():
    #os.system('module load python/2.7.8')
    os.system('module load bowtie2/2.0.6')
    os.system('module load samtools/1.1')
    os.system('module load tophat/2.0.7')
    os.system('module load cufflinks/2.2.1')
    os.system('module load java/1.7.0')
    os.system('module load BEDtools/2.22.0')
    os.system('module load bwa/0.7.9a')
    os.system('module load fastqc/0.11.2')
    os.system('module load bismark/0.13.0')
    #os.system('module list')

def get_sample_name(R1):
    if R1.endswith('_R1_001.fastq.gz'):
        samplename=R1.split('/')[-1].replace('_R1_001.fastq.gz','')
    elif R1.endswith('_R1_001.fastq'):
        samplename=R1.split('/')[-1].replace('_R1_001.fastq','')
    elif R1.endswith('_R1.fastq.gz'):
        samplename=R1.split('/')[-1].replace('_R1.fastq.gz','')
    elif R1.endswith('_R1.fastq'):
        samplename=R1.split('/')[-1].replace('_R1.fastq','')
    elif R1.endswith('.fastq'):
        samplename=R1.split('/')[-1].replace('.fastq','')
    else:
        samplename=R1.split('/')[-1]
    print samplename
    return samplename

def sort(alignfile):
    if alignfile.endswith('.sam'):
        cmd='java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar \
        INPUT=%s \
        OUTPUT=%s \
        SO=coordinate' %(alignfile,alignfile.replace('.sam','.sorted.sam'))
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
        fileout=alignfile.replace('.sam','.sorted.sam')
    elif alignfile.endswith('.bam'):
        cmd='samtools sort %s %s' %(alignfile, alignfile.replace('.bam','.sorted'))
        cmd=shlex.split(cmd)
        proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.communicate()
        fileout=alignfile.replace('.bam','.sorted.bam')
    return fileout

''' READ METRICS and STATS FILE FOR REPORT GENERATION

ALL Library Types: Trim_stats (from filter_trim.py or Trimmomatic)
                   Mapping Stats (from Duplicates Metrics)
                   Duplicate Rates (from Duplicates Mertics)
                   ChrM Mapping (ParseBam Generated Alignment Stats)

For Paired End:    Insert Size (from Insert Size Metrics)

For mRNA-Seq:      RNA-Seq Stranded

'''
def get_read_filter_stats(output_folder,R1,R2):
    total_reads=None
    if R2 is None:
        for filein in os.listdir(output_folder): 
            if filein.endswith('_stat'):
                statsfile=output_folder+filein
                with open(statsfile,'r') as fin:
                    for line in fin:
                        if line.startswith("Total number of reads"):
                            total_reads=line.strip().split('\t')[1]
                        if line.startswith("Total number of HQ filtered reads"):
                            trim_filtered=line.strip().split('\t')[1]

            elif filein.endswith('_trimming_report.txt'):
                statsfile=output_folder+filein
                with open(statsfile,'r') as fin:
                    for line in fin:
                        if line.endswith("sequences processed in total"):
                            total_reads=line.split(' ')[0]
                        if line.startswith("Number of sequence pairs removed"):
                            trim_filtered=line.split(':')[1].split('(')[0]

        if total_reads is None:
            statsfile=output_folder+'.pipeline.elog'
            with open(statsfile,'r') as fin:
                for line in fin:
                    if line.startswith('Input Reads'):
                        total_reads=str(int(line[line.index('Reads:')+len('Reads:'):line.index('Surviving:')]))
                        trim_filtered=str(int(line[line.index('Surviving:')+len('Surviving:'):line.index('(')])) 
                        break
                    
    else:
        for filein in os.listdir(output_folder): 
            if filein.endswith('_stat'):
                statsfile=output_folder+filein
                with open(statsfile,'r') as fin:
                    for line in fin:
                        if line.startswith("Total number of reads"):
                            total_reads=line.strip().split('\t')[1]
                        if line.startswith("Total number of HQ filtered reads"):
                            trim_filtered=line.strip().split('\t')[1] 
            
            elif filein.endswith('_trimming_report.txt'):
                if 'R3' in filein:
                    statsfile=output_folder+filein
                    with open(statsfile,'r') as fin:
                        for line in fin:
                            if line.startswith("Total number of sequences"):
                                total_reads=line.split(':')[1]
                            if line.startswith("Number of sequence pairs removed"):
                                trim_filtered=line.split(':')[1].split('(')[0]
                        
        if total_reads is None:
            statsfile=output_folder+'.pipeline.elog'
            with open(statsfile,'r') as fin:
                for line in fin:
                    if line.startswith('Input Read Pairs'):
                        total_reads=str(int(line[line.index('Pairs:')+len('Pairs:'):line.index('Both Surviving:')]))
                        trim_filtered=str(int(line[line.index('Both Surviving:')+len('Both Surviving:'):line.index('(')])) 
                        break

    return total_reads,trim_filtered
    
def get_duplicate_stats(output_folder,R1,R2):
    singletons=None
    if R2 is None:
         for filein in os.listdir(output_folder):
            if filein.endswith('duplicates.metrics'):
                dupsfile=output_folder+filein
                dups=[]
                with open(dupsfile) as fin:
                    for line in fin:
                        line=line.strip()
                        if line and not (line.startswith('#')):
                            dups.append(line.split('\t'))
                for i,header in enumerate(dups[0]):
                    if header in ['PERCENT_DUPLICATION']:
                        duprate=dups[1][i]
                    elif header in['UNPAIRED_READS_EXAMINED']:
                        mapped_reads=dups[1][i]
    else:
         for filein in os.listdir(output_folder):
            if filein.endswith('duplicates.metrics'):
                dupsfile=output_folder+filein
                dups=[]
                with open(dupsfile) as fin:
                    for line in fin:
                        line=line.strip()
                        if line and not (line.startswith('#')):
                            dups.append(line.split('\t'))
                for i,header in enumerate(dups[0]):
                    if header in ['PERCENT_DUPLICATION']:
                        duprate=dups[1][i]
                    elif header in['UNPAIRED_READS_EXAMINED']:
                        singletons=dups[1][i]
                    elif header in['READ_PAIRS_EXAMINED']:
                        mapped_reads=dups[1][i]
    return duprate,mapped_reads,singletons
    
def get_insert_size_stats(output_folder,R1,R2):
    for filein in os.listdir(output_folder):
        if filein.endswith('insert-size.metrics'):
            insertfile=output_folder+filein
            insert=[]
            with open(insertfile) as fin:
                for line in fin:
                    line=line.strip()
                    if line and not (line.startswith('#')):
                        insert.append(line.split('\t'))
            for i,header in enumerate(insert[0]):
                if header in ['MEDIAN_INSERT_SIZE']:
                    insert_size=insert[1][i]
                elif header in['STANDARD_DEVIATION']:
                    insert_size_SD=insert[1][i]
    return insert_size,insert_size_SD

def get_rnaseq_stats(output_folder,R1,R2):
    for filein in os.listdir(output_folder):

        if filein.endswith('rnaseq.metrics'):
            rnafile=output_folder+filein
            rna=[]
            with open(rnafile) as fin:
                for line in fin:
                    line=line.strip()
                    if line and not (line.startswith('#')):
                        rna.append(line.split('\t'))
            for i,header in enumerate(rna[0]):
                if header in ['PCT_RIBOSOMAL_BASES']:
                    pct_rrna=rna[1][i]
                elif header in ['PCT_CODING_BASES']:
                    pct_coding=rna[1][i]
                elif header in ['PCT_UTR_BASES']:
                    pct_utr=rna[1][i]
                elif header in ['PCT_INTRONIC_BASES']:
                    pct_intron=rna[1][i]
                elif header in ['PCT_INTERGENIC_BASES']:
                    pct_intergenic=rna[1][i]
                elif header in ['PCT_CORRECT_STRAND_READS']:
                    pct_strand=rna[1][i]
                elif header in ['MEDIAN_5PRIME_BIAS']:
                    med_5bias=rna[1][i]
                elif header in ['MEDIAN_3PRIME_BIAS']:
                    med_3bias=rna[1][i]
                elif header in['MEDIAN_5PRIME_TO_3PRIME_BIAS']:
                    med_53bias=rna[1][i]
                elif header in['PCT_CORRECT_STRAND_READS']:
                    pct_strand=rna[1][i]

    return pct_coding,pct_utr,pct_intron,pct_intergenic,pct_rrna,med_5bias,med_3bias,med_53bias,pct_strand

def get_wgs_stats(output_folder,R1,R2):
    for filein in os.listdir(output_folder):

        if filein.endswith('wgs.metrics'):
            wgsfile=output_folder+filein
            wgs=[]
            with open(wgsfile) as fin:
                for line in fin:
                    line=line.strip()
                    if line and not (line.startswith('#')):
                        wgs.append(line.split('\t'))
            for i,header in enumerate(wgs[0]):
                if header in ['GENOME_TERRITORY']:
                    genome_territory=wgs[1][i]
                elif header in ['MEAN_COVERAGE']:
                    mean_coverage=wgs[1][i]
                elif header in ['SD_COVERAGE']:
                    sd_coverage=wgs[1][i]
                elif header in ['PCT_EXC_TOTAL']:
                    pct_exc=wgs[1][i]
                elif header in ['PCT_10X']:
                    pct_10x=wgs[1][i]
                elif header in ['PCT_20X']:
                    pct_20x=wgs[1][i]
                elif header in ['PCT_30X']:
                    pct_30x=wgs[1][i]
                elif header in ['PCT_60X']:
                    pct_60x=wgs[1][i]
                elif header in['HET_SNP_SENSITIVITY']:
                    het_snp_sensitivty=wgs[1][i]

    return genome_territory,mean_coverage,sd_coverage,pct_exc,pct_10x,pct_20x,pct_30x,pct_60x,het_snp_sensitivty


def chrM(output_folder):
    for filein in os.listdir(output_folder):

        if filein.endswith('.stats'):
            chrmfile=output_folder+filein
            with open(chrmfile) as fin:
                for line in fin:
                    line=line.strip()
                    if line.startswith('Mito'):
                        chrm_reads=line.split(':')[-1].replace('\t','')
                    elif line.startswith('chrX'):
                        chrx_reads=line.split(':')[-1].replace('\t','')
                    elif line.startswith('chrY'):
                        chry_reads=line.split(':')[-1].replace('\t','')
    return chrm_reads,chrx_reads,chry_reads

def generate_reports_new(output_folder, R1,R2, ltype):
    global samplename
    pct_coding=None
    mean_coverage=None
    
    if R2 is None:
        total_reads,trim_filtered=get_read_filter_stats(output_folder,R1,R2)
        duprate,mapped_reads,singletons=get_duplicate_stats(output_folder,R1,R2)
        chrM_reads,chrx_reads,chry_reads=chrM(output_folder)
        if 'mRNA' in ltype or 'tRNA' in ltype:
            pct_coding,pct_utr,pct_intron,pct_intergenic,pct_rrna,med_5bias,med_3bias,med_53bias,pct_strand=get_rnaseq_stats(output_folder,R1,R2)
        elif 'WGS' in ltype:
            genome_territory,mean_coverage,sd_coverage,pct_exc,pct_10x,pct_20x,pct_30x,pct_60x,het_snp_sensitivty=get_wgs_stats(output_folder,R1,R2)

    else:
        total_reads,trim_filtered=get_read_filter_stats(output_folder,R1,R2)
        duprate,mapped_reads,singletons=get_duplicate_stats(output_folder,R1,R2)
        insert_size,insert_size_SD=get_insert_size_stats(output_folder,R1,R2)
        chrM_reads,chrx_reads,chry_reads=chrM(output_folder)
        if 'mRNA' in ltype or 'tRNA' in ltype:
            pct_coding,pct_utr,pct_intron,pct_intergenic,pct_rrna,med_5bias,med_3bias,med_53bias,pct_strand=get_rnaseq_stats(output_folder,R1,R2)
        elif 'WGS' in ltype:
            genome_territory,mean_coverage,sd_coverage,pct_exc,pct_10x,pct_20x,pct_30x,pct_60x,het_snp_sensitivty=get_wgs_stats(output_folder,R1,R2)

    if pct_coding is None and mean_coverage is None:
        if R2 is None:
            headerline='Total_reads\t%s\nTrimmed_filtered_reads\t%s\nMapped_reads\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\n' %(total_reads,trim_filtered,mapped_reads,duprate,chrM_reads,chrx_reads,chry_reads)
        else:
            headerline='Total_pairs\t%s\nTrimmed_filtered_pairs\t%s\nMapped_pairs\t%s\nSingletons\t%s\nMedian_Insert_Size\t%s\nInset_Size_SD\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\n' %(total_reads,trim_filtered,mapped_reads,singletons,insert_size,insert_size_SD,duprate,chrM_reads,chrx_reads,chry_reads)      
    elif pct_coding is not None:
        if R2 is None:
            headerline='Total_reads\t%s\nTrimmed_filtered_reads\t%s\nMapped_reads\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\nPercent_Coding\t%s\nPercent_UTR\t%s\nPercent_Intron\t%s\nPercent_Intergenic\t%s\nPercent_RRNA\t%s\nMedian_5`_Bias\t%s\nMedian_3`_Bias\t%s\nMedian_5`to3`_Bias\t%s\nPercent_Correct_strand\t%s\n' %(total_reads,trim_filtered,mapped_reads,duprate,chrM_reads,chrx_reads,chry_reads,pct_coding,pct_utr,pct_intron,pct_intergenic,pct_rrna,med_5bias,med_3bias,med_53bias,pct_strand)
        else:
            headerline='Total_pairs\t%s\nTrimmed_filtered_pairs\t%s\nMapped_pairs\t%s\nSingletons\t%s\nMedian_Insert_Size\t%s\nInset_Size_SD\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\nPercent_Coding\t%s\nPercent_UTR\t%s\nPercent_Intron\t%s\nPercent_Intergenic\t%s\nPercent_RRNA\t%s\nMedian_5`_Bias\t%s\nMedian_3`_Bias\t%s\nMedian_5`to3`_Bias\t%s\nPercent_Correct_strand\t%s\n' %(total_reads,trim_filtered,mapped_reads,singletons,insert_size,insert_size_SD,duprate,chrM_reads,chrx_reads,chry_reads,pct_coding,pct_utr,pct_intron,pct_intergenic,pct_rrna,med_5bias,med_3bias,med_53bias,pct_strand)
    elif mean_coverage is not None:
        if R2 is None:
            headerline='Total_reads\t%s\nTrimmed_filtered_reads\t%s\nMapped_reads\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\nGenome_Territory\t%s\nMean_Coverage\t%s\nSD_Coverage\t%s\nPercent_Excluded\t%s\nPercent_10X_Coverage\t%s\nPercent_20X_Coverage\t%s\nPercent_30X_Coverage\t%s\nPercent_60X_Coverage\t%s\nHet_SNP_sensitivity\t%s\n' %(total_reads,trim_filtered,mapped_reads,duprate,chrM_reads,chrx_reads,chry_reads,genome_territory,mean_coverage,sd_coverage,pct_exc,pct_10x,pct_20x,pct_30x,pct_60x,het_snp_sensitivty)
        else:
            headerline='Total_pairs\t%s\nTrimmed_filtered_pairs\t%s\nMapped_pairs\t%s\nSingletons\t%s\nMedian_Insert_Size\t%s\nInset_Size_SD\t%s\nPercent_Duplication\t%s\nchrM_Mapped\t%s\nchrX_Mapped\t%s\nchrY_Mapped\t%s\nGenome_Territory\t%s\nMean_Coverage\t%s\nSD_Coverage\t%s\nPercent_Excluded\t%s\nPercent_10X_Coverage\t%s\nPercent_20X_Coverage\t%s\nPercent_30X_Coverage\t%s\nPercent_60X_Coverage\t%s\nHet_SNP_sensitivity\t%s\n' %(total_reads,trim_filtered,mapped_reads,singletons,insert_size,insert_size_SD,duprate,chrM_reads,chrx_reads,chry_reads,genome_territory,mean_coverage,sd_coverage,pct_exc,pct_10x,pct_20x,pct_30x,pct_60x,het_snp_sensitivty)
    
    logstats=output_folder+samplename+'.QC.metrics'
    fout=open(logstats,'w')
    fout.writelines(headerline)
    fout.close()


def populate_database(output_folder,ltype,sampleid):
    global samplename
    if sampleid==None:
        pass
    else:
        import MySQLdb
        db=MySQLdb.connect(host='seqdma01',port=3306, user='root', passwd='seqdma321!',db='JAXSEQDB')
        cur=cur = db.cursor()
        
        metrics={}
        
        for filein in os.listdir(output_folder):
    
            if filein.endswith('.QC.metrics'):
                metricsfile=output_folder+filein
                with open(metricsfile) as fin:
                    for line in fin:
                        line=line.strip()
                        metrics[line.split('\t')[0]]=line.split('\t')[1]
        
        string="INSERT INTO SequencedMetadata (SD_Name_ID, \
        SD_Total_Read_Count,\
        SD_Reads_PF,\
        SD_Percentage_Mapping_Rate,\
        SD_Avg_Insert_Size, \
        SD_Percentage_Duplication_Rate,\
        SD_Target_Detection_Count, \
        SD_rRNA_Read_Count, \
        SD_mtDNA_Read_Count, \
        SD_Intonic_Read_Count,\
        SD_Exonic_Read_Count,\
        SD_Flaggs, \
        SD_RunFolder, \
        SD_DateTime, \
        SD_FASTQ_Name, \
        SampleMetadata_SM_ID,\
        SD_Archived) \
        VALUES (%s,%s,%s,%s,%s,%s,%s,)" \
        %('','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',metrics['Total_Reads'],)

def remove_files(filein):
    cmd = 'rm %s' %(filein)
    print 'Remove File:', cmd
    os.system(cmd)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--read1',
                        dest='read1',
                        required=True,
                        help='Read 1 Fastq.gz file',
                        default=None)
    parser.add_argument('-2', '--read2',
                        dest='read2',
                        required=False,
                        help='Read 2 Fastq.gz file',
                        default=None)
    parser.add_argument('-b', '--index_read',
                        dest='index_read',
                        required=False,
                        help='Index Read.gz file for Nugen RRBS',
                        default=None)
    parser.add_argument('-s','--species',
                       dest='species',
                       help='Species',
                       default=None,
                       required=True,
                       choices=['Human','Mouse','Rat','Yeast','HumanYeast','HumanArray'])
    parser.add_argument('-l','--library_type',
                       dest='libtype',
                       help='Library Type',
                       default=None,
                       required=True,
                       choices=['mRNA_Tophat',
                       'mRNA',
                       'mRNA_Unstranded',
                       'ChiP',
                       'ATAC',
                       'RRBS',
                       'WGS',
                       'WES',
                       'PDX-mRNA',
                       'PDX-WGS',
                       'PDX-tRNA',
                       'tRNA_Tophat',
                       'tRNA',
                       'mRNA_trim',
                       'mRNA_Unstranded_Tophat',
                       'mRNA_Unstranded_Tophat_trim'])
    parser.add_argument('-o','--output-folder',
                       dest='outfolder',
                       help='OutputFolder',
                       default=None,
                       required=True)
    parser.add_argument('-i','--sample-metadata-id',
                        dest='sampleId',
                        help='Sample ID for seqdma',
                        default=None,
                        required=False)
                        
    options = parser.parse_args()
    global LIBRARY_TYPE    
    R1=options.read1
    R2=options.read2
    SPECIES=options.species
    LIBRARY_TYPE=options.libtype
    INDEX_FASTQ=options.index_read
    
    if R2 is None:
        if LIBRARY_TYPE=='ChiP' or LIBRARY_TYPE=='mRNA_Unstranded' or LIBRARY_TYPE=='RRBS' or LIBRARY_TYPE=='mRNA_Unstranded_Tophat'  or LIBRARY_TYPE=='mRNA' or LIBRARY_TYPE=='mRNA_trim':
            pass
        else:
            parser.print_help()
            print 'SingleEnd is currently not supported with this Library Type'
            sys.exit(0)
    
    if SPECIES == 'HumanYeast' or SPECIES == 'YEAST':
        if LIBRARY_TYPE == 'RRBS':
            pass
        else:
            print 'Species Analysis is currently not supported with this library'
            sys.exit(0)
    
    global samplename
    samplename=get_sample_name(R1)
    
#    load_modules()
    
    if options.outfolder.endswith('/'):
        alignfolder=options.outfolder
    else:
        alignfolder=options.outfolder+'/'
    #os.mkdir(falignfolder)
#    alignfolder='/scratch/%s/' %(samplename)
#    os.mkdir(alignfolder)
    
    statusa=fastqc(R1,R2,alignfolder)
        
    if LIBRARY_TYPE=='mRNA_Tophat' or LIBRARY_TYPE=='tRNA_Tophat':
        allstatus=mrna_stranded(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='ChiP':
        allstatus=chip_seq(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='mRNA' or LIBRARY_TYPE=='tRNA':
        allstatus=mrna_stranded_RSEM(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)    
    elif LIBRARY_TYPE=='ATAC':
        allstatus=atac_bwa(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)  
    elif LIBRARY_TYPE=='PDX-mRNA' or LIBRARY_TYPE=='PDX-tRNA':
        allstatus=pdx_mrna_stranded(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='mRNA_Unstranded' or LIBRARY_TYPE=='tRNA_Unstranded':
        allstatus=mrna_unstranded_RSEM(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='mRNA_trim':
        allstatus=mrna_unstranded_RSEM_trim(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='mRNA_Unstranded_Tophat' or LIBRARY_TYPE=='tRNA_Unstranded_Tophat':
        allstatus=mrna_unstranded(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='mRNA_Unstranded_Tophat_trim':
        allstatus=mrna_unstranded_trim(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='RRBS':
        if INDEX_FASTQ is None:
            allstatus=rrbs(R1,R2,SPECIES,alignfolder)
            logs(alignfolder,allstatus)
            generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
        else:
            allstatus=nugen_rrbs(R1,R2,INDEX_FASTQ,SPECIES,alignfolder)
            logs(alignfolder,allstatus)
            generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)
    elif LIBRARY_TYPE=='WGS':
        allstatus=wgs(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)            
    elif LIBRARY_TYPE=='PDX-WGS':
        allstatus=pdx_wgs(R1,R2,SPECIES,alignfolder)
        logs(alignfolder,allstatus)
        generate_reports_new(alignfolder,R1,R2,LIBRARY_TYPE)            

    else:
        print 'Sorry - Not Implemented'
        sys.exit(0)

if __name__ == '__main__':
    main()
        
