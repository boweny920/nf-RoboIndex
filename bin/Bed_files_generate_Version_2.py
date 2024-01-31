#!/usr/bin/env python

import subprocess
import pandas as pd
import argparse
from subprocess import run, PIPE, Popen
import os
pd.set_option('display.max_columns', None)


parser = argparse.ArgumentParser(description='This script generate *.gene_data.txt file from a gtf file')
parser.add_argument('-g', '--gtf_file', help="Input gtf file associated to genome")
parser.add_argument('-d', '--id', help='ID of the genome index being processed, should be the same as the one specified in the sample sheet')
parser.add_argument('-v' , '--annotation_version', help='Version of genome annotation (usually ensembl')
parser.add_argument('-o', '--output_dir', help="Output directory ")
args=parser.parse_args()
gtf_file = args.gtf_file
out_dir = args.output_dir
bedfile_base_name = f"{args.id}.{args.annotation_version}" # e.g. hg38.Ens_106

if os.path.isdir(out_dir) == False: os.mkdir(out_dir)

df_gtf = pd.read_csv(gtf_file, header=None, sep="\t", comment='#')
df_gtf[0] = df_gtf[0].astype(str)
df_gtf[3] = df_gtf[3] - 1 #Make sure the start of the locations are subtracted by 1 to fit the bed file format (Zhang.Ning input)

##Make the exon bed file & uxon bed files##
df_exon=df_gtf[df_gtf[2] == "exon"]
if len(df_exon)>0:
    df_exon_copy=df_exon.copy() # This is to get pandas to shutup with the SettingWithCopyWarning
    df_exon_copy["exonID"] = df_exon_copy[8].str.extract('exon_id ([^;]+)')
    df_exon_copy["exonID"] = df_exon_copy["exonID"].str.replace('"','') # Get rid of the quotation marks
    df_exon_copy["geneID"] = df_exon_copy[8].str.extract('gene_id ([^;]+)') # This is for Uxon bed file making
    df_exon_copy["geneID"] = df_exon_copy["geneID"].str.replace('"','') # Get rid of the quotation marks
    df_exon_copy["chr_gene_ID"] = df_exon_copy[0] + "|" + df_exon_copy["geneID"] # This is for getting uxon regions using bedtools merge
    df_exon_copy=df_exon_copy.assign(exonGeneID=df_exon_copy["geneID"] + "|" + df_exon_copy["exonID"])
    df_exon_copy[9] = [1] * len(df_exon)
    df_exon_bed = df_exon_copy.loc[:,[0,3,4,"exonGeneID",9,6]]
    exon_bed_fileName = os.path.join(out_dir, bedfile_base_name + ".exons.bed")
    df_exon_bed.to_csv(exon_bed_fileName, index=None, header=None, sep="\t")
    #Making Uxons from here
    df_uxon_bed_tmp = df_exon_copy.loc[:,["chr_gene_ID",3,4]]
    df_uxon_bed_tmp_fileName = os.path.join(out_dir, bedfile_base_name + ".uxon.bed.tmp")
    df_uxon_bed_tmp.to_csv(df_uxon_bed_tmp_fileName, index=None, header=None, sep="\t")
    df_uxon_bed_tmp_merged_fileName = os.path.join(out_dir, bedfile_base_name + ".uxon.bed.merged.tmp")
    with open(df_uxon_bed_tmp_merged_fileName,"w") as newfile: # Generate the merged exons based on geneID and chrom name
        subprocess.run([f"sort -k1,1 -k2,2n {df_uxon_bed_tmp_fileName}| bedtools merge"] , shell=True, stdout=newfile) # each feature is listed only once!
    df_uxons=pd.read_csv(df_uxon_bed_tmp_merged_fileName,header=None,sep="\t")
    df_uxons_copy=df_uxons.copy()
    df_uxons_copy["chr"] = df_uxons_copy[0].str.split("|").str[0]
    df_uxons_copy["geneID"] =  df_uxons_copy[0].str.split("|").str[1]
    df_uxons_copy = df_uxons_copy.loc[:,["chr",1,2,"geneID"]]
    uxon_bed_fileName_tmp = os.path.join(out_dir, bedfile_base_name + ".uxons.bed.tmp")
    df_uxons_copy.to_csv(uxon_bed_fileName_tmp, sep="\t", index=None, header=None)
    uxon_bed_fileName = os.path.join(out_dir, bedfile_base_name + ".uxons.bed")
    with open(uxon_bed_fileName,"w") as newfile:
        subprocess.run([f"sort -k1,1 -V -k2,2n {uxon_bed_fileName_tmp} "], shell=True,
                       stdout=newfile)  # each feature is listed only once!

    subprocess,run([f"rm {out_dir}/*.tmp"],shell=True)
else: print("GTF file does not contain exon information, thus cannot generate associated bed file")

##Make the gene bed file & intergenic bed file##
df_gene=df_gtf[df_gtf[2] == "gene"]
if len(df_gene) > 0:
    df_gene_copy=df_gene.copy()
    df_gene_copy[8]=df_gene_copy[8].str.extract('gene_id ([^;]+)')
    df_gene_copy[8] = df_gene_copy[8].str.replace('"','') # Get rid of the quotation marks
    df_gene_copy[9] = [1] * len(df_gene_copy)
    df_gene_copy=df_gene_copy.loc[:,[0,3,4,8,9,6]]
    gens_bed_fileName=os.path.join(out_dir, bedfile_base_name+".genes.bed")
    df_gene_copy.to_csv(gens_bed_fileName, index=None, header=None, sep="\t")
    # Make the intergenic bed from here
    end_loc = df_gene_copy.loc[:, 4].tolist()
    end_loc_col = [df_gene_copy[3].tolist()[0]] + end_loc[:-1]  # Making columns of the end locations of the previous genes
    end_chrom=df_gene_copy.loc[:,0].tolist()
    end_chrom_col=[df_gene_copy[0].tolist()[0]] + end_chrom[:-1] # Making columns of the chrom names of the previous genes
    df_gene_copy = df_gene_copy.assign(pre_end=end_loc_col)
    df_gene_copy = df_gene_copy.assign(pre_end_chrom=end_chrom_col)
    df_gene_copy = df_gene_copy[df_gene_copy[0] == df_gene_copy["pre_end_chrom"]] # This is to make sure the intergenic distances stay within the associated chrom
    df_interGene=df_gene_copy.loc[:,[0,"pre_end",3]]
    df_interGene=df_interGene.iloc[1:,:]
    intergens_bed_fileName = os.path.join(out_dir, bedfile_base_name + ".intergenic.bed")
    df_interGene.to_csv(intergens_bed_fileName, index=None, header=None, sep="\t")
else: print("GTF file does not contain gene information, thus cannot generate associated bed file")

##Make the transcript bed file & intron bed file ##
df_transcript=df_gtf[df_gtf[2] == "transcript"]
if len(df_transcript) > 0:
    df_transcript_copy=df_transcript.copy()
    df_transcript_copy["transID"]=df_transcript_copy[8].str.extract('transcript_id ([^;]+)')
    df_transcript_copy["transID"] = df_transcript_copy["transID"].str.replace('"','') # Get rid of the quotation marks
    df_transcript_copy["geneID"] = df_transcript_copy[8].str.extract('gene_id ([^;]+)')
    df_transcript_copy["geneID"] = df_transcript_copy["geneID"].str.replace('"', '')  # Get rid of the quotation marks
    df_transcript_copy[9] = [1] * len(df_transcript_copy)
    df_transcript_copy = df_transcript_copy.assign(geneTrans_ID=df_transcript_copy["geneID"] + "|" + df_transcript_copy["transID"])
    df_transcript_bed = df_transcript_copy.loc[:,[0,3,4,"geneTrans_ID",9,6]]
    transcript_bed_fileName=os.path.join(out_dir, bedfile_base_name+".transcripts.bed")
    df_transcript_bed.to_csv(transcript_bed_fileName, index=None, header=None, sep="\t")
    ##Make the intron bed files from here
    if len(df_exon) > 0:
        df_transcript_copy["chr_transID_strand"] =  df_transcript[0] + "|" + df_transcript_copy["transID"] + "|" + df_transcript_copy[6]
        df_transcript_tmp=df_transcript_copy.loc[:,["chr_transID_strand",3,4]]
        df_exon_copy["transID"] = df_exon[8].str.extract('transcript_id ([^;]+)')
        df_exon_copy["transID"]=df_exon_copy["transID"].str.replace('"','')
        df_exon_copy["chr_transID_strand"] = df_exon[0] + "|" + df_exon_copy["transID"] + "|" + df_exon_copy[6] # Working with the previous copy df to save mem
        df_exon_tmp = df_exon_copy.loc[:, ["chr_transID_strand", 3, 4]]
        df_exon_tmp_filename = os.path.join(out_dir, bedfile_base_name+".exon.bed.tmp")
        df_transcript_tmp_filename = os.path.join(out_dir, bedfile_base_name +".transcripts.bed.tmp")
        df_exon_tmp.to_csv(df_exon_tmp_filename, index=None, header=None, sep="\t")
        df_transcript_tmp.to_csv(df_transcript_tmp_filename, index=None, header=None, sep="\t")
        with open(os.path.join(out_dir,"intron.tmp"),"w") as newfile:
            subprocess.run([f"bedtools2 subtract -a {df_transcript_tmp_filename} -b {df_exon_tmp_filename}"], shell=True, stdout=newfile) # Subtract transcript with exons to get introns # Using the the new version of bedtools
        if os.path.getsize(os.path.join(out_dir,"intron.tmp")) > 0:
            df_introns=pd.read_csv(os.path.join(out_dir,"intron.tmp"), sep="\t", header=None)
            df_introns_copy=df_introns.copy()
            df_introns_copy["chr"] = df_introns_copy[0].str.split("|").str[0]
            df_introns_copy["transID"] = df_introns_copy[0].str.split("|").str[1]
            df_introns_copy["strand"] = df_introns_copy[0].str.split("|").str[2]
            df_introns_copy["score"] = [1] * len(df_introns_copy)
            df_introns_copy = df_introns_copy.loc[:, ["chr", 1, 2, "transID","score","strand"]]
            intron_bed_fileName = os.path.join(out_dir, bedfile_base_name + ".introns.bed")
            df_introns_copy.to_csv(intron_bed_fileName, sep="\t", index=None, header=None)
        subprocess, run([f"rm {out_dir}/*.tmp"], shell=True)
else: print("GTF file does not contain transcript information, thus cannot generate associated bed file")

##Make the UTR bed file
df_utr=df_gtf[(df_gtf[2] == "three_prime_utr") | (df_gtf[2] == "five_prime_utr") ]
if len(df_utr) > 0:
    df_utr_copy=df_utr.copy()
    df_utr_copy["transID"] = df_utr_copy[8].str.extract('transcript_id ([^;]+)')
    df_utr_copy["geneID"] = df_utr_copy[8].str.extract('gene_id ([^;]+)')
    df_utr_copy = df_utr_copy.assign(UTR_ID=df_utr_copy["geneID"] + "|" + df_utr_copy["transID"] + "|" + df_utr_copy[2])
    df_utr_copy["UTR_ID"] = df_utr_copy["UTR_ID"].str.replace('"','') # Get rid of the quotation marks
    df_utr_copy[9] = [1] * len(df_utr_copy)
    df_utr_bed = df_utr_copy.loc[:,[0,3,4,"UTR_ID",9,6]]
    df_utr_bed_fileName=os.path.join(out_dir, bedfile_base_name + ".UTRs.bed")
    df_utr_bed.to_csv(df_utr_bed_fileName, sep="\t", index=None, header=None)
else: print("GTF file does not contain UTR information, thus cannot generate associated bed file")
