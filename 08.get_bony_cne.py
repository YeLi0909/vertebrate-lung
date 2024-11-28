#! /usr/bin/env python
import os, re

cne_info = "active_cne.tsv"
bony_cne = set()

with open(cne_info,"r") as f:
    for line in f:
        elements = line.strip().split("\t")
        if elements[1] == "bony":
            bony_cne.add(elements[0])

dis_info = "shared.cnes.2m.plus.distance"
out_file = "bony_cne.distance"

with open(dis_info,"r") as f, open(out_file,"w") as output:
    headline = next(f)
    output.write(headline)
    for line in f:
        elements = line.strip().split("\t")
        cne_id = elements[1]
        total_count = int(elements[2])
        if total_count < 18:
            continue
        paired_count = int(elements[3])
        if paired_count/total_count < 0.8:
            continue
        if cne_id in bony_cne:
            output.write(line)
