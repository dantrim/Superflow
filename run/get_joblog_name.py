#!/usr/bin/env python

##############################################################################
#
# this script is used to determine the name of the output log for a job
# run on the condor batch system
#
# inputs:
#   1) 'process_group_dir' : name of process being submitted
#   2) 'process_number' : process ID number of condor cluster
#   3) 'split_dsids' : comma separated list (no spaces between) of DSIDS
#                   that should have split submission
#
# daniel.joseph.antrim@cern.ch
# June 2017
#
##############################################################################


import sys
import glob

if __name__ == "__main__" :
    
    process_group_dir = sys.argv[1]
    process_number = sys.argv[2]
    split_dsids = sys.argv[3]

    samples_to_split = []
    if split_dsids != "X" :
        samples_to_split = split_dsids.split(",")

    if not process_group_dir.endswith("/") :
        process_group_dir = process_group_dir + "/"
    dsid_lists = glob.glob(process_group_dir + "*.txt")

    container_names = []
    all_names = []
    containers = {}

    global_process_no = 0
    for ids, ds in enumerate(dsid_lists) :

        containers[ids] = None

        is_split = False
        for sds in samples_to_split :
            if sds in ds :
                is_split = True

        p_list = []
        if is_split :
            lines = open(ds).readlines()
            for line in lines :
                if not line : continue
                line = line.strip()
                all_names.append(line)
                p_list.append(global_process_no)
                global_process_no += 1

        else :
            all_names.append(ds)
            p_list.append(global_process_no)
            global_process_no += 1

        containers[ids] = p_list

    # find which global container this process number belongs to
    global_to_use = -1
    found_it = False
    for global_id, process_id_list in containers.iteritems() :
        if found_it :
            break
        for p_id in process_id_list :
            if int(p_id) == int(process_number) :
                global_to_use = global_id
                found_it = True
                break

    if global_to_use < 0 :
        print "ERROR did not find global idx"
        sys.exit()
    global_name = dsid_lists[global_to_use]
    global_name = global_name.split("/")[-1].replace(".txt","_%d.log"%(int(process_number)))
    print global_name
