#!/usr/bin/env python

import pandas as pd
import os
import argparse
pd.set_option('display.max_columns', None)

parser = argparse.ArgumentParser(description='This script generate *.gene_data.txt file from a gtf file')
parser.add_argument('-g', '--gtf_file', help= "Input gtf file associated to genome")
parser.add_argument('-d', '--gene_data_file', help= "Input .gene_data.txt file generated from genome")
parser.add_argument('-r', '--refFlat_file', help= "Input .refFlat.txt file generated from genome gtf")
parser.add_argument('-o', '--output_dir', default="./", help="Output directory, default=cwd, default='./'")
args=parser.parse_args()
gtf_file = args.gtf_file
gene_data_file=args.gene_data_file
refFlat_file=args.refFlat_file
out_dir = args.output_dir
# gtf_file="mm39.Ens_106.gtf"
# gene_data_file="../gene_data_txt_make/mm39.Ens_106.gene_data.txt"
# refFlat_file="../GenomeBPcount_fileMake/mm39.Ens_106.refFlat.txt"
# out_dir="transOut_test"

if os.path.isdir(out_dir) == False: os.mkdir(out_dir)

df_gtf = pd.read_csv(gtf_file, sep="\t", header=None, comment='#')
df_gene_data = pd.read_csv(gene_data_file, sep="\t")
df_refFlat = pd.read_csv(refFlat_file, header=None, sep="\t")
df_transcript=df_gtf[df_gtf[2] == "transcript"]

if len(df_transcript) > 0:
    df_transcript_copy=df_transcript.copy()
    df_transcript_copy["Transcript_ID"] = df_transcript_copy[8].str.extract('transcript_id ([^;]+)')
    df_transcript_copy["Transcript_ID"] = df_transcript_copy["Transcript_ID"].str.replace('"','') # Get rid of the quotation marks
    df_transcript_copy["Gene_ID"] = df_transcript_copy[8].str.extract('gene_id ([^;]+)')
    df_transcript_copy["Gene_ID"] = df_transcript_copy["Gene_ID"].str.replace('"', '')  # Get rid of the quotation marks
    df_transcript_copy = df_transcript_copy.assign(cDNA_len = df_transcript_copy[4] - df_transcript_copy[3] + 1)
    df_transcript_copy["transcript_Biotype"] = df_transcript_copy[8].str.extract('transcript_biotype ([^;]+)')
    df_transcript_copy["transcript_Biotype"] = df_transcript_copy["transcript_Biotype"].str.replace('"','')  # Get rid of the quotation marks
    df_refFlat = df_refFlat.loc[:,[1,6,7,8]] #Take what is needed from the refFlat file into mem, this is to make the script faster
    df_refFlat = df_refFlat.rename(columns={1:"Transcript_ID", 6:"CDS_Start", 7:"CDS_End", 8:"N_Exons"})
    df_refFlat = df_refFlat.assign(CDS_Len=df_refFlat["CDS_End"] - df_refFlat["CDS_Start"])
    df_allInfo = pd.merge(df_transcript_copy, df_refFlat, on="Transcript_ID", how='left') # merge the two dataframes based on the df_transcript_copy
    df_gene_data = df_gene_data.loc[:,["Gene_ID","Genomic_Len","Name"]] #Take what is needed from the gene_data.txt file into mem
    df_allInfo = pd.merge(df_allInfo, df_gene_data, on="Gene_ID", how='left')
    df_allInfo = df_allInfo.rename(columns={0:"Chrom", 3:"Start", 4:"End", 6:"Strand"})
    df_allInfo = df_allInfo.loc[:,["Transcript_ID", "Name", "Gene_ID", "Chrom", "Start", "End", "Strand", "Genomic_Len", "cDNA_len", "N_Exons", "transcript_Biotype", "CDS_Len", "CDS_Start", "CDS_End"]]
    transcript_data_txt_fileName = os.path.join(out_dir,os.path.basename(gene_data_file).split(".gene_data")[0]+".transcript_data.txt")
    df_allInfo.to_csv(transcript_data_txt_fileName, sep="\t", index=None)

else: 
    print(f"Transcript information absent/incomplete in gtf file: {gtf_file}, cannot generate transcript_data.txt file")
    transcript_data_txt_fileName = os.path.join(out_dir,os.path.basename(gene_data_file).split(".gene_data")[0]+".transcript_data.txt")
    with open(transcript_data_txt_fileName, "w") as new:
        new.write(f"Transcript information absent/incomplete in gtf file: {gtf_file}, cannot generate transcript_data.txt file")
