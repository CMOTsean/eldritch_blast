#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:37:03 2020

@author: sean
"""

import sys, getopt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os

pd.options.mode.chained_assignment = None



def repeat_check(plots, infile, outfile, fwd_length, rev_length, n, ODIRA=False):
    df = pd.read_csv(infile, sep="\t", names=["query", "subject", "ID", "length", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bit"])
    df["difference"] = df["send"]-df["sstart"]
    reads_copy_number={}
    no_flank_count=0
    
    ###Loops through each read in blast ouput
    for x in df.subject.unique():
        sub_df = df[df["subject"]==x]
        sub_df.sort_values(by=["send"], inplace=True)
        sub_df = sub_df.reset_index(drop=True)
        
        ###Counts number of reads with no repeat
        if len(sub_df[abs(sub_df["difference"]) > fwd_length+(2*n*0.6)]) > 0:
            reads_copy_number[sub_df.iloc[0,1]] = 1
        
        ###Counts reads which have sequence flanking repeat
        flank_hits_df = sub_df[abs(sub_df["difference"]) > fwd_length+(n/10)]
        if len(flank_hits_df) == 2:
            if ODIRA:
                copyNumb = (min(abs(flank_hits_df.iloc[0,9]-flank_hits_df.iloc[1,8]), abs(flank_hits_df.iloc[1,9]-flank_hits_df.iloc[0,8]))-rev_length)/(fwd_length+rev_length)
            else:
                copyNumb = (min(abs(flank_hits_df.iloc[0,9]-flank_hits_df.iloc[1,8]), abs(flank_hits_df.iloc[1,9]-flank_hits_df.iloc[0,8])))/(fwd_length)
            read_name = flank_hits_df.iloc[0,1]
            reads_copy_number[read_name] = round(copyNumb)+2
            if plots == "flanked":
                repeat_plot(sub_df, x, fwd_length, n, outfile, ODIRA=ODIRA)
        else:
            no_flank_count += 1
        
        ###Plots all reads, regardless of flanking status
        if plots == "all":
            repeat_plot(sub_df, x, fwd_length, n, outfile, ODIRA=ODIRA)
    return reads_copy_number
        
def repeat_plot(sub_df, x, rep_length, n, outdir, flank=True, ODIRA=False):
    plt.figure(x)
    sign_flip=1
    if sum(sub_df["difference"]<(-(rep_length+(n/10))))>0: 
        sign_flip = sign_flip*-1
    labelcount=0
    def len_check(j):
        if j["difference"] < 0:
            return -1
        else:
            return 1
            
    ##Plots the blast hits as horizontal lines, with ABC hits on top, B on bottom
    if ODIRA:
        for i,j  in sub_df.iterrows():
            plt.hlines(len_check(j)*(sign_flip), int(j["sstart"]), int(j["send"]), colors="black")
            if len_check(j)*sign_flip==1:
                labelcount+=1
                plt.text(y=len_check(j)*1.5*sign_flip, x=np.mean([int(j["sstart"]), int(j["send"])]), s=str(labelcount))
    
    ##Plots tandem repeats alternating top and bottom to prevent overlap
    else:
        for i,j  in sub_df.iterrows():
            plt.hlines(sign_flip, int(j["sstart"]), int(j["send"]), colors="black")
            labelcount+=1
            plt.text(y=1.5*sign_flip, x=np.mean([int(j["sstart"]), int(j["send"])]), s=str(labelcount))
            sign_flip = sign_flip * -1
               
    ##Prettify the plot and save it
    plt.ylim(-9, 10)
    plt.tick_params(axis='both',which='both',left=False,labelleft=False)
    plt.ylabel(x)
    plt.savefig(outdir+"/"+x+".png") 
    plt.close()
    



def main(argv):
    
    ###Default options
    ODIRA=False
    n=1000
    ODIRA_fwd=0
    ODIRA_rev=0
    tandem_length=0
    read_names=False
    plots = "None"
    Help_text = """
    This script estimates the copy number of tandem/ODIRA repeats from blast results
    of the repeat + flanking regions against long reads such as minION.
    
    Arguments
    ~~~~~~~~~
    Required:
        -i or --input:          csv containing output of blast search against the reads
                                using -outfmt 6
        -o or --out_dir:        directory destination for plots
        -l or --repeat_length:  the length of repeat (bp) for tandem repeats (not required when using ODIRA option)
    Options:
        -h or --help            displays help message
        -p or --plots           "None", no plotting -Default-
                                "flanked", only plot reads containing the flanking DNA
                                "all", plot all reads
        -n or --flank_length    length of flanking DNA on either side. Default: 1000
        --read_names            prints list of reads with respective copy number to terminal
        --ODIRA                 no argument needed. Specifies ODIRA mode
        --ODIRA_fwd             length of forward repeat in ODIRA mode
        --ODIRA_rev             length of reverse repeat in ODIRA mode
        
    Example cases:
        `python minion_blast.py -i blast_results.csv -o test_dir -l 1500 -p flanked`
        This will give an estimated copy number for a repeat 1500bp long with 1000bp flanking sequence and plot 
        only the reads containing the flanking sequence.
        
        `python minion_blast.py -i blast_results.csv -o test_dir --ODIRA --ODIRA_fwd 800 --ODIRA_rev 400 --read_names`
        This will give an estimated copy number for an ODIRA repeat which has 800bp with 400bp reverse, make no plots,
        and print the reads and their estimated copy number to the terminal.
    """
    ###Get argument options
    try:
        opts, args = getopt.getopt(argv, "hn:l:i:o:p:",["plots=", "help", "input=", "out_dir=", "out_prefix=", "ODIRA_fwd=", "ODIRA_rev=", "ODIRA", "flank_length=", "repeat_length=", "read_names"])
    except getopt.GetoptError:
        print(Help_text)
        sys.exit(2)
        
    ###Assign to variables    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(Help_text)
            sys.exit()
        elif opt in ("-i", "--input"):
            fin = arg
        elif opt in ("-o", "--out_dir"):
            outdir = arg
        elif opt in ("-p", "--plots"):
            if arg == "flanked" or arg == "all":
                plots = arg
            else:
                print(Help_text)
        elif opt in ("-n", "--flank_length"):
            n = int(arg)
        elif opt in ("-l", "--repeat_length"):
            tandem_length = int(arg)
        elif opt in ("--ODIRA"):
            ODIRA = True
        elif opt in ("--ODIRA_fwd"):
            ODIRA_fwd = int(arg)
        elif opt in ("--ODIRA_rev"):
            ODIRA_rev = int(arg)
        elif opt in ("--read_names"):
            read_names = True
    ###Asserts proper repeat lengths are given
    assert fin and outdir, "input file and output directory are required"
    if ODIRA:
        assert ODIRA_fwd > 0, "--ODIRA_fwd required when ODIRA specified"
        assert ODIRA_rev > 0, "--ODIRA_fwd required when ODIRA specified"
        fwd_length = ODIRA_fwd
    else:
        assert tandem_length > 0, "-l or --repeat_length required when ODIRA not specified"
        fwd_length = tandem_length
    ###Makes out directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    ###Call main function
    results = repeat_check(plots=plots, 
                           ODIRA=ODIRA, 
                           infile=fin, 
                           outfile=outdir, 
                           fwd_length=fwd_length, 
                           rev_length=ODIRA_rev, 
                           n=n)
    if read_names == True:
        for i in results:
            print(i, results[i])
    print(Counter(results.values()))

                  
if __name__ == "__main__":
    main(sys.argv[1:])