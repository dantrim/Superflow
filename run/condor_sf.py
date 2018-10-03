#!/bin/env python

import os
import sys
import glob
import argparse
import subprocess
import time

exec_name = "ntupler_nn"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0303/mc/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0303/mc/logs/"

filelist_dir = "/data/uclhc/uci/user/dantrim/n0303val/susynt-read/filelists/"

storage_prefixes = { "mwt2" : "root://fax.mwt2.org:1094/" }

samples = [ "ttbar" ]
dsids_to_split = [ "410472" ]

run_with_systematics = False
tar_name = 'area.tgz'

do_brick = True
do_gp = True
do_uc = True
do_sdsc = False

def check_samples(sample_names, filelist_dir) :

    all_ok = True
    for name in sample_names :
        name_loc = filelist_dir + "/" + name
        if not os.path.isdir(name_loc) :
            print 'ERROR Could not find sample \"%s\" in filelist directory (filelist dir = %s)' % (name, filelist_dir)
            all_ok = False
    if not all_ok :
        print 'Exitting with ERROR'
        sys.exit()

def get_split_samples(txt_files_for_sample) :

    samples_to_split = {}
    samples_to_not_split = []

    for txt in txt_files_for_sample :
        do_split = False
        dsid_for_split = ""
        for dsid in dsids_to_split :
            if dsid in txt :
                do_split = True
                dsid_for_split = dsid
                break
        if do_split :
            samples_to_split[dsid_for_split] = []
            with open(txt, 'r') as txt_file :
                for line in txt_file :
                    line = line.strip()
                    if not line : continue
                    if line.startswith('#') : continue
                    samples_to_split[dsid_for_split].append(line)
        else :
            samples_to_not_split.append(txt)

    return samples_to_split, samples_to_not_split

def bool_string(boolean) :

    return { True : 'true', False : 'false' } [ boolean ]

def make_condor_file(sample, queue_list, condor_filename, exec_name) :

do_brick = True
do_gp = True
do_uc = True
do_sdsc = False

    with open(condor_filename, 'w') as f :

        f.write('universe = vanilla\n')
        f.write('+local=%s\n' % bool_string(do_brick))
        f.write('+site_local=%s\n' % bool_string(do_gp))
        f.write('+uc=%s\n' % bool_string(do_uc))
        f.write('+sdsc=%s\n' % bool_string(do_sdsc))
        f.write('executable = %s\n' % exec_name)
        f.write('should_transfer_files = YES\n')
        f.write('use_x509userproxy = True\n')
        f.write('notification = Never\n')

        for q in queue_list :
            log_base = q.split('/')[-1].replace('.txt','')
            sys_string = '-c'
            if run_with_systematics :
                sys_string = '-a'

    


def submit_sample(sample, verbose = False) :

    txt_files_for_sample = glob.glob(filelist_dir + "/" + sample + "/*.txt")
    n_txt = len(txt_files_for_sample)


    if verbose :
        print 'Found %d text files for sample %s' % (n_txt, sample)

    if n_txt == 0 :
        print 'WARNING Found no text files for sample %s, skipping this sample' % sample
        return

    split_dict, queue_list = get_split_samples(txt_files_for_sample)

    if queue_list :
        condor_filename = 'submit_%s.condor' % sample
        exec_name = 'run_condor_sf_%s.sh' % sample
        condor_file_name, condor_exec_name = make_condor_file(sample, queue_list, condor_filename, exec_name)


    

    

def main() :

    parser = argparse.ArgumentParser( description = 'Launch SuperFlow jobs to the condor batch' )
    parser.add_argument('--sample', default = "",
        help = 'Override samples list by providing a sample name (can be comma-separated-list)')
    parser.add_argument('-v', '--verbose', default = False, action = 'store_true',
        help = 'Turn on verbose mode')
    args = parser.parse_args()

    if args.sample != "" :
        samples = [s.strip() for s in args.sample.split(",")]
    check_samples(samples, filelist_dir)

    n_samples = len(samples)
    for isample, sample in enumerate(samples) :
        print '[%02d/%02d] Submitting %s' % (isample+1, n_samples, sample)
        submit_sample(sample, args.verbose)

    

#____________________________
if __name__ == '__main__' :
    main()
