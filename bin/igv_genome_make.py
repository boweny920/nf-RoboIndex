#!/usr/bin/env python

import json
import argparse

parser = argparse.ArgumentParser(description='This script generate igv genome json file')
parser.add_argument('-g', '--gtf_file', help= "Input gtf file associated to genome")
parser.add_argument('-n', '--fna_file' ,help="Output file name")
parser.add_argument('-f', '--fai', help='fai index file of genome')
parser.add_argument('-i', '--id', help='id of genome')
parser.add_argument('-a', '--name', help='name of genome')
parser.add_argument('-v', '--version', help='version of annotation')
args=parser.parse_args()

igenmoe = {
    "id" : args.id,
    "name" : f"{args.name}",
    "fastaURL" : f"https://webfs{args.fna_file}",
    "indexURL" : f"https://webfs{args.fai}",
    "tracks" : [
        {"name" : f"Genes {args.version}",
         "type": "annotation",
         "url" : f"https://webfs{args.gtf_file}",
         "format": "gtf"
         }
    ]
}

outfile_name = f"{args.id}.{args.version}.igv_genome.json"

with open(outfile_name,"w") as newfile:
    json.dump(igenmoe,newfile)

