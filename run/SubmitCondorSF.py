#!/bin/env python
import os
import sys
import glob
import subprocess

ana_name = "wwlikeAna"
tar_location = "/data/uclhc/uci/user/dantrim/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0216/validation/mc/logs/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0216/validation/mc/Raw/"

filelist_dir = "/data/uclhc/uci/user/dantrim/n0216val/filelists/"
in_job_filelist_dir = "/n0216val/filelists/"
samples = ["sherpaVV"]

#retry_list = "/data/uclhc/uci/user/dantrim/n0216val/Superflow/run/retry.txt"

doBrick = True
doLocal = False
doSDSC  = False
doUC    = False

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
        for dataset in sample_lists :
            fullname = str(os.path.abspath(dataset))
            print "    > %s"%dataset
           # do_this = False
           # for ds in retry_lines :
           #     ds = ds.strip()
           #     ds = ds.replace("./user", "user")
           #     ds = ds.replace(".out" , "")
           #     if ds in fullname :
           #         do_this = True
           # if not do_this :
           #     continue
            

            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
            print "    >> %s"%dataset
            

            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
                print "You must call this script from the output directory where the ntuples will be stored!"
                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
                sys.exit()

            run_cmd = "ARGS="
            run_cmd += '"'
            run_cmd += ' %s '%out_dir
            run_cmd += ' %s '%log_dir
            run_cmd += ' %s '%ana_name
            #run_cmd += ' %s '%(tar_location + "area.tgz")
            run_cmd += ' n0216val '
            run_cmd += ' %s '%dataset
            run_cmd += ' -c '
            run_cmd += '"'
            run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
            lname = dataset.split("/")[-1].replace(".txt", "")
            run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
            run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
            run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
            run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")

            print run_cmd
            subprocess.call(run_cmd, shell=True)

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
    #f.write('transfer_input_files = area.tgz\n')
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
    f.write('echo "Calling : source RestFrames/setup_RestFrames.sh"\n')
    f.write('source RestFrames/setup_RestFrames.sh\n')
    f.write('echo "Calling : source RootCore/local_setup.sh"\n')
    f.write('source RootCore/local_setup.sh\n')
    f.write('echo "Calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('echo "Calling : ${sflow_exec} -i ${input} ${sflow_options}"\n')
    f.write('${sflow_exec} -i ${input} ${sflow_options}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()
    
    


if __name__=="__main__" :
    main()

