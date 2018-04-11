#!/usr/bin/env python

import sys

if __name__=="__main__" :
    input_list = sys.argv[1]
    process_no = sys.argv[2]

    log_name = input_list.split("/")[-1]
    log_name = log_name.replace(".txt","_run%s.log"%process_no)
    print log_name
