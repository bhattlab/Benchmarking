"""
Calculating non-overlapping window coverage and GC content from pile file
Author: Dylan Maghini
Date: Jan 14 2022
"""

import sys

pilefile = sys.argv[1]
windowsize = int(sys.argv[2])

with open(pilefile, "r") as infile:
    # initialize variables
    cur_window = 0
    cur_covg = []
    cur_bases = ""
    cur_contig = ""
    cur_start = ""

    # read in a line
    cur_line = infile.readline().strip()
    while cur_line:
        contig, pos, base, cov, p1, p2 = cur_line.split("\t")
        # initialize if on new contig or previous line hit window size
        if (contig != cur_contig) or (cur_window == windowsize):
            cur_contig = contig
            cur_window = 0
            cur_covg = []
            cur_bases = ""
            cur_start = pos
        # update values with new position
        cur_window += 1
        cur_covg.append(float(cov))
        cur_bases += base
        # if at full window, calculate and output information
        if cur_window == windowsize:
            cur_bases = cur_bases.upper()
            if cur_bases.count("N") == 0:
                avgcovg = sum(cur_covg)/len(cur_covg)
                gccontent = (cur_bases.count("G") + cur_bases.count("C"))/(len(cur_bases)) * 100
                print("\t".join([contig, pos, str(avgcovg), str(gccontent)]))
        cur_line = infile.readline().strip()
