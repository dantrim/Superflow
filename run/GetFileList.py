#!/usr/bin/env python
import sys

if __name__=="__main__" :
    group_name = sys.argv[1]
    process_number = sys.argv[2]

    injob_filelist = "joblist_%s.txt"%group_name
    samples = []
    lines = open(injob_filelist).readlines()
    for line in lines :
        if not line : continue
        line = line.strip()
        samples.append(line)
    n_total = len(samples)
    if int(process_number) >= n_total :
        print "GetFileList.py ERROR input process number (%d) is out of range, there are only %d samples in process (max process number allowed: %d)"%(int(process_number), n_total, n_total-1)
        print ""
    else :
        in_job_input = "/filelists/%s/%s"%(group_name, samples[int(process_number)])
        print in_job_input
        #print samples[int(process_number)]
    
