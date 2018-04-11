#!/bin/env python
import os
import sys
import glob
import subprocess
import time

ana_name = "ntupler_wwbb_novec"
tar_location = "/data/uclhc/uci/user/dantrim/"

#n_split = sys.argv[1]

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug23/mc/ttbar/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug23/mc/ttbar/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar/retry_51/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar/retry_51/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar_pp8/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar_pp8/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/g_oct2_l1topo/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/g_oct2_l1topo/mc/ttbar/logs/"
#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar_amcatnlo/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/ttbar_amcatnlo/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/f_sep28_c1c1ww/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/f_sep28_c1c1ww/mc/ttbar/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/mc/ttbar/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/i_nov22_tight/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/i_nov22_tight/mc/ttbar/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/j_dec7_tight/mc/ttbar_pp8/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/j_dec7_tight/mc/ttbar_pp8/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/k_jan6_bsys/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/k_jan6_bsys/mc/ttbar/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/l_mar21_nobjetcut/mc/ttbar/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/l_mar21_nobjetcut/mc/ttbar/logs/"


filelist_dir = "/data/uclhc/uci/user/dantrim/n0234val/filelists/"
in_job_filelist_dir = "/n0234val/filelists/"
#samples = ["diboson_sherpa_lvlv", "drellyan_sherpa", "wjets_sherpa22", "zjets_sherpa22", "ttV", "singletop"]
#samples = ["ttbar_min"]
#samples = ["Wt"]#, "wwbb_susy2"]
samples = ["ttbar"]
#samples = ["ttbar_pp8"]
#samples = ["retry_mc"]
#samples = ["zjets_and_DY"]
#samples = ["higgs"]
#samples = ["data16_n0234"] #, "data16_n0234"]
#samples = ["data16_n0234"]

samples_to_split = ["410009", "410503", "410505", "410225"]
#samples_to_split = ["341122", "364102", "364114", "364115", "364116", "364118", "364200"]
#splits_to_do = [10, 11, 12, 22, 23]

doBrick = True
doLocal = False  
doSDSC  = False 
doUC    = False 


def get_retrylist() :
    #runlist = ["279867", "280464"]
    runlist = ["364102","364114", "364115", "364116", "364118", "364200"]
    runlist = ["342053"]
    #runlist = ["341122"]
    #retryfile = "/data/uclhc/uci/user/dantrim/n0234val/resub.txt"
    #lines = open(retryfile).readlines()
    #for line in lines :
    #    if not line : continue
    #    line = line.strip()
    #    runlist.append(line)
    return runlist

def main() :
    print "SubmitCondorSF"

    ### retry [begin]
    #retry_lines = open(retry_list).readlines()
    

#    look_for_tarball()
    look_for_condor_script(brick_ = doBrick, local_ = doLocal, sdsc_ = doSDSC, uc_ = doUC)
    look_for_condor_executable()

    for s in samples :
        print "Submtting sample : %s"%s
        suff = ""
        if not s.endswith("/") : suff = "/"
        sample_lists = glob.glob(filelist_dir + s + suff + "*.txt")
        if len(sample_lists) == 0 :
            print "No sample lists in filelist dir!"
            sys.exit()

        #retry_list = get_retrylist()
        #new_lists = []
        #for s_ in sample_lists :
        #    for x in retry_list :
        #        if x in s_ :
        #            new_lists.append(s_)
        #sample_lists = new_lists


        for dataset in sample_lists :

            

            fullname = str(os.path.abspath(dataset))
            dataset_original = dataset
            print "    > %s"%dataset

            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
            print "    >> %s"%dataset

            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
                print "You must call this script from the output directory where the ntuples will be stored!"
                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
                sys.exit()

            do_split_sample = False
            for dsid_split in samples_to_split :
                if dsid_split in dataset :
                    do_split_sample = True

            run_mode = "-c "

            if not do_split_sample :

                # submit the job as usual
                run_cmd = "ARGS="
                run_cmd += '"'
                run_cmd += ' %s '%out_dir
                run_cmd += ' %s '%log_dir
                run_cmd += ' %s '%ana_name
                #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
                run_cmd += ' n0234val '
                run_cmd += ' %s '%dataset
                run_cmd += ' %s '%run_mode # any extra cmd line optino for Superflow executable
                run_cmd += '"'
                run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
                lname = dataset.split("/")[-1].replace(".txt", "")
                run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
                run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
                run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
                run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")

                print run_cmd
                subprocess.call(run_cmd, shell=True)

                time.sleep(0.5)

            elif do_split_sample :
                split_suffix = 0

                split_files = []
                lines = open(dataset_original).readlines()
                for line in lines :
                    if not line : continue
                    line = line.strip()
                    split_files.append(line)

                for split_file in split_files :
                    #if split_suffix not in splits_to_do :
                    #    split_suffix = split_suffix + 1
                    #    continue
                    print "    >>> Sub-file [%d] %s"%(split_suffix, split_file)

                    run_cmd = "ARGS="
                    run_cmd += '"'
                    run_cmd += ' %s '%out_dir
                    run_cmd += ' %s '%log_dir
                    run_cmd += ' %s '%ana_name
                    #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
                    run_cmd += ' n0234val '
                    run_cmd += ' %s '%split_file
                    run_cmd += ' %s --sumw ./n0234val/sumw_file.txt --suffix %d'%(run_mode, split_suffix) # any extra cmd line optino for Superflow executable
                    run_cmd += '"'
                    run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
                    lname = dataset.split("/")[-1].replace(".txt", "")
                    run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
                    run_cmd += ' -append "output = %s%s" '%(log_dir, lname + "_%d.out"%split_suffix)
                    run_cmd += ' -append "log = %s%s" '%(log_dir, lname + "_%d.log"%split_suffix)
                    run_cmd += ' -append "error = %s%s" '%(log_dir, lname + "_%d.err"%split_suffix)

                    print run_cmd
                    subprocess.call(run_cmd, shell=True)

                    time.sleep(0.5)
                    
                    split_suffix = split_suffix + 1
                    

def look_for_tarball() :
    if not os.path.isfile("area.tgz") :
        print "Tarball not found."
        sys.exit()

def look_for_condor_script(brick_ = False, local_ = False, sdsc_ = False, uc_ = False) :

    brick = 'false'
    local = 'false'
    sdsc  = 'false'
    uc    = 'false'
    if brick_ : brick = 'true'
    if local_ : local = 'true'
    if sdsc_  : sdsc = 'true'
    if uc_    : uc = 'true'

    f = open('submitFile_TEMPLATE.condor', 'w')
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick_)
    f.write('+site_local=%s\n'%local_)
    f.write('+sdsc=%s\n'%sdsc_)
    f.write('+uc=%s\n'%uc_)
    #f.write('transfer_input_files = area.tgz.tgz\n')
    f.write('executable = RunCondorSF.sh\n')
    f.write('arguments = $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    #f.write('transfer_output_files = OUTFILE\n')
    #f.write('transfer_output_remaps = OUTFILE_REMAP\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue\n')
    f.close()

def look_for_condor_executable() :
    f = open('RunCondorSF.sh', 'w') 
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " ------- RunCondorSF -------- "\n')
    f.write('output_dir=${1}\n')
    f.write('log_dir=${2}\n')
    f.write('sflow_exec=${3}\n')
    f.write('stored_dir=${4}\n')
    f.write('input=${5}\n')
    f.write('sflow_options=${@:6}\n\n')
    f.write('echo "    output directory   : ${output_dir}"\n')
    f.write('echo "    log directory      : ${log_dir}"\n')
    f.write('echo "    sflow executable   : ${sflow_exec}"\n')
    f.write('echo "    tarred dir         : ${stored_dir}"\n')
    f.write('echo "    sample list        : ${input}"\n')
    f.write('echo "    sflow options      : ${sflow_options}"\n\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xzf area.tgz\n\n')
    f.write('echo "done untarring\n')
    f.write('echo "current directory structure:\n')
    f.write('ls -ltrh\n\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo "Calling : cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n')
    f.write('echo "Directory structure:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('echo "Calling : source RootCore/local_setup.sh"\n')
    f.write('source RootCore/local_setup.sh\n')
    #f.write('echo "Calling : cd SuperRest/"\n')
    #f.write('cd SuperRest/\n')
    #f.write('source setRestFrames.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Calling : ${sflow_exec} -i ${input} ${sflow_options}"\n')
    f.write('${sflow_exec} -i ${input} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    
    


if __name__=="__main__" :
    main()

