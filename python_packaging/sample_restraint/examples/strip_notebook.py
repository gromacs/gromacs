#!/usr/bin/env python
# This script takes an ipython notebook as an argument and rewrites the file
# with metadata stripped to make change tracking with git easier.
import sys
import json
import os
import shutil

infile = sys.argv[1]
if not os.path.exists(infile):
    sys.exit("command line argument must be an existing filename")
tempfile = infile + ".tmp"

with open(infile, "r") as fh:
    json_in = json.load(fh)

nb_metadata = json_in["metadata"]


def strip_output_from_cell(cell):
    if "outputs" in cell:
        cell["outputs"] = []
    if "execution_count" in cell:
        cell["execution_count"] = None


for cell in json_in["cells"]:
    strip_output_from_cell(cell)

with open(tempfile, "w") as fh:
    json.dump(json_in, fh, sort_keys=True, indent=1, separators=(",", ": "))
shutil.move(tempfile, infile)
