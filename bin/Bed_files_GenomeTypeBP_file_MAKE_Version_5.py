#!/usr/bin/env python

import subprocess
import pandas as pd
import argparse
from subprocess import run, PIPE, Popen
import os
from Bio import SeqIO

# Argument inputs
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--refFlat", help="attach *.refFlat.txt file")
parser.add_argument("-n", "--fna_file", help="attach FNA file of the genome in fasta")
parser.add_argument("-o", "--output_dir", help="specify the output directory, default = ./genomebp", default="./genomebp")
parser.add_argument("-g", "--gene_data", help="attach *.gene_data.txt file")

args=parser.parse_args()
refFlat_txt = args.refFlat
fna_file = args.fna_file
output_dir = args.output_dir
gene_data_txt = args.gene_data

# Argument inputs

# bedtool functions
def bedtools_sort_n_merge(gene_bed_file_path):
    run([f"sort -k1,1 -k2,2n {gene_bed_file_path} > {gene_bed_file_path.replace('.bed', '.sorted.bed.tmp')}"], shell=True)
    run([
            f"bedtools merge -i {gene_bed_file_path.replace('.bed', '.sorted.bed.tmp')} > {gene_bed_file_path.replace('.bed', '.sorted.merged.bed.tmp')}"],
        shell=True)
    return f"{gene_bed_file_path.replace('.bed', '.sorted.merged.bed.tmp')}"

def bedtools_subtract(file_A, file_B, file_Out):
    run([f"bedtools subtract -a {file_A} -b {file_B} > {file_Out}"], shell=True)


## Make the output folder
if os.path.isdir(output_dir) == False:
    os.makedirs(output_dir)


## Generate the bed file for all the genes by processing refFlat file
gene_bed_file_path=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_genes.bed.tmp")
df_gene_data_txt=pd.read_csv(gene_data_txt,sep="\t")
df_all_genes=pd.DataFrame()
df_all_genes[0]=df_gene_data_txt["Chrom"]
df_all_genes[1]=df_gene_data_txt["Start"]
df_all_genes[2]=df_gene_data_txt["End"]
df_all_genes.to_csv(gene_bed_file_path,index=False,header=False,sep="\t") # Write the dataframe to a tsv bed file
merged_allgene_bedfile=bedtools_sort_n_merge(gene_bed_file_path) # Sort the bed file and run merge on it


## Generate the bed file for all CDS regions from the refFlat file. (NOTE: This is not the true CDS! The true CDS is generated later in this script):
df_refflat_txt=pd.read_csv(refFlat_txt,sep="\t",header=None)
CDS_bed_file_path=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_CDS_VerRefFlat.bed.tmp")
df_all_CDS=pd.DataFrame()
df_all_CDS[0]=df_refflat_txt[2]
df_all_CDS[1]=df_refflat_txt[6]
df_all_CDS[2]=df_refflat_txt[7]

df_all_CDS_filtered = df_all_CDS[(df_all_CDS[1]!=0) | (df_all_CDS[2]!=0)] # This is to filter out the Chr 0 0 lines, so that bedtools will be happy later on
df_all_CDS_filtered.to_csv(CDS_bed_file_path,index=False,header=False,sep="\t")
# Sort the bed file and run merge on it
merged_CDS_bedfile_verRefFlat = bedtools_sort_n_merge(CDS_bed_file_path)


## Generate the bed file for all Exon regions:
df_refFlat_txt=pd.read_csv(refFlat_txt,header=None,sep="\t")
Exon_bed_file_path=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_Exon.bed.tmp")
with open(Exon_bed_file_path,"w") as newfile:
    for index, line in df_refFlat_txt.iterrows():
        exon_star_list = line[9].rstrip(",").split(",")
        exon_end_list = line[10].rstrip(",").split(",")
        chrom_name = line[2]
        strand_name = line[3]
        for num, pos in enumerate(exon_star_list):
            if int(pos) != 0 | int(exon_end_list[num]) != 0: # This is to filter out the Chr 0 0 lines, so that bedtools will be happy later on
                newfile.write(f"{chrom_name}\t{pos}\t{exon_end_list[num]}\t{line[0]}|{line[1]}\t1\t{strand_name}\n") # Columns: Chr  start   end geneName|transcriptName
# Sort the bed file and run merge on it
merged_Exon_bedfile = bedtools_sort_n_merge(Exon_bed_file_path)


## Generate the bed files for all Introns:
intron_bed_file=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_Intron.bed.tmp")
bedtools_subtract(merged_allgene_bedfile,merged_Exon_bedfile,intron_bed_file)
# Sort the bed file and run merge on it
merged_Inton_bedfile = bedtools_sort_n_merge(intron_bed_file)


## Generate the bed files for all UTRs:
bedtools_subtract(merged_Exon_bedfile, merged_CDS_bedfile_verRefFlat, os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_UTRs.bed.tmp"))
# Sort the bed file and run merge on it
merged_UTRs_bedfile = bedtools_sort_n_merge(os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_UTRs.bed.tmp"))

## Generate the real CDS regions!! Load all beds, gtf and fasta to IGV to see the reason behind this process
real_CDS_bedfile=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_actual_CDS.bed.tmp") # This is becuause I found that the real CDS = EXONS (from RefFlat) - UTRs (Calculated) above
bedtools_subtract(merged_Exon_bedfile, merged_UTRs_bedfile, real_CDS_bedfile)

## Generate a bed file for transcripts. This is purly for downstream analysis
df_transcript_bed=df_refflat_txt.loc[:,[2,4,5,1,3]]
transcript_bed_file=os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".all_transcripts.bed.tmp")
df_transcript_bed.to_csv(transcript_bed_file,sep="\t",header=None, index=False)


######################## ############ ############  Now that we have built all the required bed files, Time to calcualte the base distribution of different regions  ############ ############ ############ ############ ############

def base_count_sum(bed_file,base_type): # This function calculate the sum of base pairs for the genome and for individual contigs, and return a dataframe
    if os.path.getsize(bed_file) > 0:
        df_bed=pd.read_csv(bed_file,header=None,sep="\t")
        #make a distance column
        df_bed[3]=df_bed[2].astype(int)-df_bed[1].astype(int)+1
        whole_genome_count=df_bed[3].sum() # this is to count the base of this type within the whole genome
        chroms=df_bed[0].unique()
        chrom_list=[]
        base_sum_list=[]
        for chrom in chroms:
            base_sum_list.append(df_bed[df_bed[0]==chrom][3].sum())
            chrom_list.append(chrom)
        index=["Genome"]+chrom_list
        base_counts=[whole_genome_count]+base_sum_list
    else: 
        base_counts = [0]
        index=["Genome"]
    df_out=pd.DataFrame(index=index)
    df_out[base_type]=base_counts
    # print(df_out)
    return df_out

def inter_geneic_dis_cal(df_contig):
    end_loc = df_contig.loc[:, 2].astype(int).tolist()
    end_loc_col = [df_contig[1].astype(int).tolist()[0]] + end_loc[:-1]  # This is to add a column to the gene bed dataframe so i can calculate the inter-gene distances faster
    df_contig = df_contig.assign(pre_end = end_loc_col)
    df_contig.loc[:, 1].astype(int)
    df_contig = df_contig.assign( inter_dist = df_contig[1] - df_contig["pre_end"])
    return df_contig["inter_dist"].sum()

df_CDS=base_count_sum(real_CDS_bedfile,"PC-CDS")
df_Intron=base_count_sum(merged_Inton_bedfile,"INTRON")
df_UTR=base_count_sum(merged_UTRs_bedfile,"PC-UTR")

#calculate the intergeneic distance:( this takes the most time)
df_genes_bed=pd.read_csv(merged_allgene_bedfile,header=None,sep="\t")
chroms = list(df_genes_bed[0].unique())
inter_gene_base_count=[]
for chrom in chroms:
    df_chrom=df_genes_bed[df_genes_bed[0] == chrom]
    inter_gene_base_count.append(inter_geneic_dis_cal(df_chrom))
inter_gene_base_count_genome=sum(inter_gene_base_count)
inter_gene_base_count_final=[inter_gene_base_count_genome] + inter_gene_base_count
df_interGene=pd.DataFrame(index=["Genome"] + chroms)
df_interGene["INTERGENE"] = inter_gene_base_count_final

# Calculate the chromosome length
contig_list=[]
contig_length_list=[]
for record in SeqIO.parse(fna_file,"fasta"):
    contig_list.append(record.id)
    contig_length_list.append(len(record.seq))
genome_length=sum(contig_length_list)
df_total_length=pd.DataFrame(index=["Genome"]+contig_list)
df_total_length["TOTAL"]=[genome_length] + contig_length_list

df_final = pd.concat([df_total_length, df_interGene, df_Intron, df_UTR, df_CDS], axis=1).fillna(0).astype(int)

# Calculate ribosome length
# This part depends on rRNA annotation, will parse all the rRNA information from the gene_data.txt file
df_rRNA=df_gene_data_txt[df_gene_data_txt["Biotype"] == "rRNA"]
if len(df_rRNA) > 0:
    df_rRNA = df_rRNA.assign(rRNA_len = df_rRNA["End"] - df_rRNA["Start"])
    rRNA_sum = df_rRNA["rRNA_len"].sum()
else:
    rRNA_sum = 0
df_final["RIBOSOME"] = [rRNA_sum] + [0] * (len(df_final)-1)

# Make the NC-EXON column
# This is mostly for the R script that generate the secundo report, the NC-EXON adds all columns together to the genome/contig length
df_CSV_ready = df_final.assign(NCEXON = df_final["TOTAL"] - df_final["INTERGENE"] - df_final["INTRON"] - df_final["PC-UTR"] - df_final["PC-CDS"] - df_final["RIBOSOME"])
df_CSV_ready.rename(columns={"NCEXON": "NC-EXON"}, inplace=True)
df_CSV_ready = df_CSV_ready.iloc[[0]]
df_CSV_ready.to_csv(os.path.join(output_dir,os.path.basename(refFlat_txt).split(".refFlat.")[0]+".GenomeBpTypes.txt"), sep="\t")

#Remove all the tmp files
subprocess.run([f"rm {output_dir}/*.tmp"],shell=True)