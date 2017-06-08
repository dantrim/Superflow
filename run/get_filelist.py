#!/usr/bin/env python

import sys

if __name__ == "__main__" :
    group_name = sys.argv[1]
    process_number = sys.argv[2]
    split_dsids = sys.argv[3]
    stored_dir = sys.argv[4]

    samples_to_split = []
    if split_dsids != "X" :
        samples_to_split = split_dsids.split(",")

    injob_filelist = "joblist_%s.txt"%group_name
    samples = []

    lines = open(injob_filelist).readlines()

    split_options = ""

    for line in lines :
        if not line : continue
        line = line.strip()

        split_this_sample = False

        for ds in samples_to_split :
            if ds in line :
                split_this_sample = True

        if split_this_sample :

            split_options = "--sumw --suffix %d"%(int(process_number))

            rname = "./filelists/%s/%s"%(group_name, line) 
            #rname = "./%s/filelists/%s/%s"%(stored_dir, group_name, line) 
            rfiles = open(rname).readlines()
            for rf in rfiles :
                if not rf : continue
                rf = rf.strip()
                samples.append(rf)
        else :
            samples.append(line)

    n_total = len(samples)
    if int(process_number) >= n_total :
        print "get_filelist.py    ERROR input process number %d is out of range, there are only %d samples in process\n"%(int(process_number), n_total)
    else :
        submit_sample = samples[int(process_number)]
        if "susyNt.root" in submit_sample :
            print submit_sample
        else :
            submit_sample = "./%s/filelists/%s/%s"%(stored_dir, group_name, submit_sample)
            print submit_sample

    if split_options != "" :
        print split_options
