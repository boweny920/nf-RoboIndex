#!/usr/bin/env python

import  sys
import os
import shutil

infile = sys.argv[1]
destination_file = sys.argv[2]

if os.path.isdir(destination_file):
    destination_dir = destination_file
else:
    destination_dir = os.path.dirname(destination_file)

if os.path.dirname(infile).replace("//","/") == destination_dir.replace("//","/"):
    # if ".gtf" in os.path.basename(infile):
    #     shutil.copy(infile, destination_file)
        # os.remove(infile)
    print(f"{infile} is already in {destination_dir}")

else:
    shutil.copy2(infile, destination_file)
