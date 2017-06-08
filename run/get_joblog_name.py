#!/usr/bin/env python

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
        #print global_id
        if found_it :
            break
        for p_id in process_id_list :
            #print " > ",p_id
            #print "p_id number : %d %d"%(p_id, process_number)
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
