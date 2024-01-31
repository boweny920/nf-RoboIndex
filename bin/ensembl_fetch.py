#!/usr/bin/env python

import argparse
import requests
import pandas as pd
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_sheet", help="attach sample-Sheet")
parser.add_argument("-d", "--output_index_directory", help="directory of where all the index folders will be stored")
args=parser.parse_args()
df_sampleSheet = pd.read_csv(args.sample_sheet)
output_index_directory = args.output_index_directory

def download_file(url, filename):
    response = requests.get(url)
    if response.status_code == 200:
        with open(filename, 'wb') as file:
            file.write(response.content)
        print(f'Download complete: {filename}')
        subprocess.run(['gunzip',filename])
        # subprocess.run(['rm',filename])
        if ".gtf" in filename:
            with open('tmp.gtf',"w") as newfile:
                newfile.write(subprocess.run(['grep','-v','#',filename.replace('.gz','')],stdout=subprocess.PIPE).stdout.decode('utf-8'))
            subprocess.run(['mv','tmp.gtf', filename.replace(".gz","")])
    else:
        print('Error with ensembl download:', response.status_code)

for index, row in df_sampleSheet.iterrows():
    if pd.isna(row["source"]) == False:
        index_dir = os.path.join(output_index_directory,f"{row['name']}/{row['id']}") # This is to unify the name of folder and all files
        gtf_dir = os.path.join(index_dir,f"annotation/{row['annotation_version']}/gtfs")
        os.makedirs(gtf_dir, exist_ok=True)
        # the first letter must be captialized for the file name of the gtf or fna
        gtf_url = f"https://ftp.ensembl.org/pub/release-{row['annotation_version'].split('_')[-1]}/gtf/{row['name'].lower()}/{row['name'].capitalize()}.{row['id']}.{row['annotation_version'].replace('Ens_','')}.gtf.gz"
        fna_url = f"https://ftp.ensembl.org/pub/release-{row['annotation_version'].split('_')[-1]}/fasta/{row['name'].lower()}/dna/{row['name'].capitalize()}.{row['id']}.dna_sm.toplevel.fa.gz"
        # gtf_filename = os.path.join(gtf_dir,
        #                             f"{row['name'].capitalize()}.{row['id']}.{row['annotation_version'].replace('Ens_','')}.gtf.gz")
        gtf_filename = os.path.join(gtf_dir,
                                    f"{row['name']}.{row['id']}.{row['annotation_version']}.gtf.gz") # This is to unify the name of folder and all files
        download_file(gtf_url, gtf_filename)
        # https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
        # fna_filename = os.path.join(index_dir,
        #                             f"{row['name'].capitalize()}.{row['id']}.dna_sm.toplevel.fa.gz")
        fna_filename = os.path.join(index_dir,
                                    f"{row['name']}.{row['id']}.fa.gz")
        download_file(fna_url, fna_filename)
        df_sampleSheet.at[index, 'gtf_path'] = gtf_filename.replace(".gz","") # only this will replace the values in the dataframe
        df_sampleSheet.at[index, 'fna_path'] = fna_filename.replace(".gz","")
        df_sampleSheet.at[index, 'source'] = ''

df_sampleSheet.to_csv("samplesheet_processed.csv",index=False)