#!/usr/bin/env python

import  sys
import os
import shutil

infile = sys.argv[1]
destination_file = sys.argv[2]

if infile.replace("//","/") != destination_file.replace("//","/"):
    shutil.move(infile,destination_file)