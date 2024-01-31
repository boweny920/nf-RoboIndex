#!/usr/bin/env python

### IMPORTANT ###
# Before running this script,
# run: java -Xmx10g -jar /n/apps/CentOS7/bin/picard.jar CreateSequenceDictionary R=reference.fa O=reference.dict
# to generate picard dictionary file
### ###

import pandas as pd
import argparse
import os
import subprocess
pd.set_option('display.max_columns', None)

parser = argparse.ArgumentParser(description='This script generate *.riboList.default.txt file from associated *.gene_data.txt file')
parser.add_argument('-g', '--gene_data', help= "Input .gene_data.txt file associated to genome")
parser.add_argument('-o', '--output_dir', default= "./", help="Output directory, default=cwd, default='./'")
parser.add_argument('-p', '--picard_dict', help= "Input dict file associated to genome made by picard CreateSequenceDictionary")
args=parser.parse_args()
geneData_file = args.gene_data
out_dir = args.output_dir
fa_dict = args.picard_dict
# geneData_file = "../gene_data_txt_make/mm39.Ens_106.gene_data.txt"
# out_dir = "rRNA_test_ver1"
# fa_dict = "./mm39.dict"

if os.path.isdir(out_dir) == False: os.mkdir(out_dir)

df_geneData = pd.read_csv(geneData_file, sep="\t")
df_geneData_rRNA = df_geneData[df_geneData["Biotype"] == "rRNA"]
if len(df_geneData_rRNA) > 0:
    df_geneData_rRNA = df_geneData_rRNA.loc[:,["Chrom","Start","End","Strand","Gene_ID"]]
    rRNA_fileNAME=os.path.join(out_dir,os.path.basename(geneData_file).split(".gene_data")[0]+".riboList.default.txt")
    df_geneData_rRNA.to_csv(rRNA_fileNAME,header=None,sep="\t", index=None)

    headers=subprocess.run(["awk -F \"\t\" '{print $1\"\t\"$2\"\t\"$3}' %s"%fa_dict],stdout=subprocess.PIPE, shell=True).stdout
    with open(rRNA_fileNAME,"r+") as newfile:
        content=newfile.read()
        newfile.seek(0,0)
        newfile.write(headers.decode('utf-8') + content)
else:
    print("No rRNA information present in input .gene_data.txt file")
    rRNA_fileNAME=os.path.join(out_dir,os.path.basename(geneData_file).split(".gene_data")[0]+".riboList.default.txt")
    headers=subprocess.run(["awk -F \"\t\" '{print $1\"\t\"$2\"\t\"$3}' %s"%fa_dict],stdout=subprocess.PIPE, shell=True).stdout
    with open(rRNA_fileNAME, "w") as newfile:
        newfile.write(headers.decode('utf-8'))