#!/usr/bin/env python

import pandas as pd
import os
import argparse
import numpy as np
import subprocess

#Input arguments
parser = argparse.ArgumentParser(description='This script generate *.gene_data.txt file from a gtf file')
parser.add_argument('-g', '--gtf_file', help= "Input gtf file associated to genome")
parser.add_argument('-o', '--output_dir', default="./" ,help="Output directory, default=cwd, default='./'")
parser.add_argument('-e', '--ensembl_database', action='store',default=False, type=str,help="Provide -e/--ensembl_databases [an ensemble_database_name in ensembl_datasets.csv] if input data is of ensembl annotation")
args=parser.parse_args()

gtf_file = args.gtf_file
out_dir = args.output_dir
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
df_gtf = pd.read_csv(gtf_file, sep="\t", header=None, comment='#')

#Collect gene related information from input gtf file
gene_dic={}
for index, row in df_gtf.iterrows():
    info_list=row[8].split(";")
    gene_id=[v for v in info_list if "gene_id" in v][0].split()[-1].replace('"','')
    gene_strand=row[6]
    chrom_name=row[0]
    #Get gene name
    if len([v for v in info_list if "gene_name" in v]) > 0:
        gene_name=[v for v in info_list if "gene_name" in v][0].split()[-1].replace('"','')
    else: gene_name=""
    #Get gene source
    if len([v for v in info_list if "gene_source" in v]) > 0:
        gene_source=[v for v in info_list if "gene_source" in v][0].split()[-1].replace('"','')
    else: gene_source = "unknown"
    #Get gene biotype
    if len([v for v in info_list if "gene_biotype" in v]) > 0:
        gene_biotype=[v for v in info_list if "gene_biotype" in v][0].split()[-1].replace('"','')
    else: gene_biotype = "Unknown"
    #Get the positions for exons, CDS, and transcripts
    if "exon" in row[2]:
        exon_pos = (row[3],row[4])
    else: exon_pos=None  ## REMEMBER TO REMOVE THE NONE LATER!!!
    if "CDS" in row[2]:
        CDS_pos = (row[3],row[4])
    else: CDS_pos = None
    if "transcript" in row[2]:
        transcript_pos = (row[3],row[4])
    else: transcript_pos = None

    if gene_id not in gene_dic:
        gene_dic.setdefault(gene_id,{"gene_name":gene_name, "Chrom":chrom_name ,"gene_source":gene_source, "gene_biotype":gene_biotype,
                                     "strand":gene_strand, "exon_pos":{exon_pos}, "CDS_pos": {CDS_pos}, "transcript_pos": {transcript_pos}})
    else:
        gene_dic[gene_id]["exon_pos"].add(exon_pos)
        gene_dic[gene_id]["CDS_pos"].add(CDS_pos)
        gene_dic[gene_id]["transcript_pos"].add(transcript_pos)

#Making the columns of the fina table
Gene_ID_list=[]
Gene_Name_list=[]
Gene_ID_Name_list=[]
Chrom_name_list=[]
gene_start_list=[]
gene_end_list=[]
gene_strand_list=[]
gene_length_list=[]
exon_len_list=[]
CDS_len_list=[]
N_Trans_list=[]
N_Exons_list=[]
BioType_list=[]
CDS_pos_start_list=[]
CDS_pos_end_list=[]

def merge_positions_and_length_Cal(pos_list):
    #Merge the overlap positions of CDS, EXONS, Transcriptions and etc
    positions = [v for v in pos_list if v != None]
    if len(positions) > 0:
        # sort positions by start time
        sorted_positions = sorted(positions, key=lambda pos: pos[0])
        # merge overlapping positions
        merged_positions = [sorted_positions[0]]
        for pos in sorted_positions[1:]:
            last_pos = merged_positions[-1]
            if pos[0] <= last_pos[1]:
                merged_positions[-1] = (last_pos[0], max(pos[1], last_pos[1]))
            else:
                merged_positions.append(pos)
                # Calculate the total length after merging the locations
        list_of_length = [v[1] - v[0] + 1 for v in merged_positions]
        return sum(list_of_length), len(sorted_positions), sorted_positions
    else:
        return np.nan, np.nan, [None]

def gene_start_end_calculate(exons_pos,key,gtf_file):
    # From the gtf of mm39 from ensumble, the last exon location is the end of the gene, the start of the first exon is the start of the gene
    positions = [v for v in exons_pos if v != None]
    if len(positions) > 0:
        sorted_positions = sorted(positions, key=lambda pos: pos[0])
        return [sorted_positions[0][0], sorted_positions[-1][-1]]
    else: # This is for genes in gtf that do not have exon information, they may be psedo genes
        start_end=subprocess.run([f"grep {key} {gtf_file} | sort -u | cut -f 4,5"],shell=True,stdout=subprocess.PIPE).stdout.decode('utf-8').split()
        if len(start_end) <0 :
            return [np.nan, np.nan]
        else:
            return [int(v) for v in start_end]

for key in gene_dic:
    Gene_ID_list.append(key)
    Gene_Name_list.append(gene_dic[key]["gene_name"])
    if gene_dic[key]["gene_name"] != "": Gene_ID_Name_list.append(key + "|" + gene_dic[key]["gene_name"])
    else: Gene_ID_Name_list.append(key)
    Chrom_name_list.append((gene_dic[key]["Chrom"]))
    gene_start_end=gene_start_end_calculate(gene_dic[key]["exon_pos"],key,gtf_file)
    gene_start_list.append(gene_start_end[0])
    gene_end_list.append(gene_start_end[1])
    gene_strand_list.append(gene_dic[key]["strand"])
    gene_length_list.append(gene_start_end[1]-gene_start_end[0] + 1)
    exon_len, exon_Num, exon_pos_sorted = merge_positions_and_length_Cal(gene_dic[key]["exon_pos"])
    CDS_len, CDS_Num, CDS_pos_sorted = merge_positions_and_length_Cal(gene_dic[key]["CDS_pos"])
    trans_len, trans_Num, trans_pos_sorted = merge_positions_and_length_Cal(gene_dic[key]["transcript_pos"])
    exon_len_list.append(exon_len)
    CDS_len_list.append(CDS_len)
    N_Trans_list.append(trans_Num)
    N_Exons_list.append(exon_Num)
    BioType_list.append(gene_dic[key]["gene_biotype"])
    CDS_pos_start_tmp = [v[0] for v in CDS_pos_sorted if v != None]
    CDS_pos_end_tmp = [v[1] for v in CDS_pos_sorted if v != None]
    CDS_pos_start_list.append(",".join(str(v) for v in CDS_pos_start_tmp))
    CDS_pos_end_list.append(",".join(str(v) for v in CDS_pos_end_tmp))

#Making the columns into a dataframe
df_gene_data=pd.DataFrame()
df_gene_data["Gene_ID"] = Gene_ID_list
df_gene_data["Name"] = Gene_Name_list
df_gene_data["Gene_ID|Name"] = Gene_ID_Name_list
df_gene_data["Chrom"] = Chrom_name_list
df_gene_data["Start"] = gene_start_list
df_gene_data["End"] = gene_end_list
df_gene_data["Strand"] = gene_strand_list

df_gene_data["Genomic_Len"] = gene_length_list
df_gene_data["Genomic_Len"] = df_gene_data["Genomic_Len"].fillna("0").astype(int) # This is to get rid of the decimals (unnecessary, but pretty)
df_gene_data["Exonic_Len"] = exon_len_list
df_gene_data["Exonic_Len"] = df_gene_data["Exonic_Len"].fillna("0").astype(int) # This is to get rid of the decimals
df_gene_data["CDS_Len"] = CDS_len_list
df_gene_data["CDS_Len"] = df_gene_data["CDS_Len"].fillna("0").astype(int) # This is to get rid of the decimals

df_gene_data["N_Trans"] = N_Trans_list
df_gene_data["N_Exons"] = N_Exons_list
df_gene_data["Biotype"] = BioType_list
df_gene_data["CDS_Genomic_Start"] = CDS_pos_start_list
df_gene_data["CDS_Genomic_End"] = CDS_pos_end_list


## Get the annotation information for ensembl organisms, this will also generate a table of geneIDs and gene descriptions (.biomart.annotation.tsv)
if args.ensembl_database != False:
    biomart_annotation_out_file = os.path.join(out_dir, os.path.basename(gtf_file).rstrip(".gtf")+".biomart.annotation.tsv")
    ensembl_database = args.ensembl_database
    cmd = [
        f"wget -O {biomart_annotation_out_file} 'http://www.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"{ensembl_database}\" interface = \"default\" ><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"description\" /></Dataset></Query>'"]
    # The annotation fetched from biomart is a tsv table with the following columns: geneID description ensembl_transcript_id   ensembl_peptide_id
    try:
        subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
    except subprocess.CalledProcessError as e:
        print(f"Error in Biomart annotation information fetching:\n{e}")
    else:
        print(f"Biomart annotation information saved in {biomart_annotation_out_file}")

    df_biomart_annotation = pd.read_csv(biomart_annotation_out_file, sep="\t", names=["ensembl_gene_id", "description"])

    description_list = []
    for index, row in df_gene_data.iterrows():
        description_list.append("|".join(str(v).replace('nan', "") for v in df_biomart_annotation[
            df_biomart_annotation["ensembl_gene_id"] == row["Gene_ID"]]["description"].tolist()))
    df_gene_data["Gene_Description"] = description_list
else:
    df_gene_data["Gene_Description"] = [""] * len(df_gene_data)

df_gene_data.to_csv(os.path.join(out_dir, os.path.basename(gtf_file).replace(".gtf",".gene_data.txt")), sep="\t", index=False)