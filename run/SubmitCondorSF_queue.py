#!/bin/env python
import os
import sys
import glob
import subprocess
import time

ana_name = "ntupler_wwbb_novec"
tar_location = "/data/uclhc/uci/user/dantrim/"

#n_split = sys.argv[1]

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/b_aug7/data/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/b_aug7/data/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug17/data/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug17/data/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug23/data/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/d_aug23/data/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/mc/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/f_sep28_c1c1ww/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/f_sep28_c1c1ww/mc/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/g_oct2_l1topo/mc/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/g_oct2_l1topo/mc/logs/"

#out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/data/Raw/"
#log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/e_aug31/data/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/mc/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/data/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/h_oct18_sys/data/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/i_nov22_tight/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/i_nov22_tight/mc/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/j_dec7_tight/data/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/j_dec7_tight/data/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/k_jan6_bsys/mc/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/k_jan6_bsys/mc/logs/"

out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/l_mar21_nobjetcut/data/Raw/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0234/l_mar21_nobjetcut/data/logs/"



filelist_dir = "/data/uclhc/uci/user/dantrim/n0234val/filelists/"
in_job_filelist_dir = "/n0234val/filelists/"
samples = ["data_n0234"]
#samples = ["retry_data"]
#samples = ["c1c1_ww", "diboson_sherpa"]
#samples = ["wwbb_susy2"]
#samples = ["singletop", "WtAMC", "drellyan_sherpa", "higgs", "singletop_DS", "ttV", "wjets_sherpa", "zjets_sherpa", "diboson_sherpa"]
#samples = ["singletop", "wwbb_susy2", "drellyan_sherpa", "diboson_sherpa", "higgs", "singletop_DS", "ttV", "wjets_sherpa", "zjets_sherpa"]
#samples = ["wwbb_susy2"]
#samples = ["diboson_sherpa"]
#samples = ["bffN"]
#samples = ["bwn_m1001L20_overlap"]
#samples = ["wwbb_susynt"]
#samples = ["bwn_m1001L20"]
#samples = ["diboson_sherpa"]
#samples = ["zjets_and_DY"]

samples_to_split = ["410009"]

doBrick = True 
doLocal = False 
doSDSC  = False 
doUC    = False 


def get_retrylist() :
    print "GRABBING DSIDS FROM RETRYLIST"
    runlist = ["364110", "364160"]
    #lines = open("data_to_resubmit.txt").readlines()
    #lines = open("dsids_to_test.txt").readlines()
    #for line in lines :
    #    if not line : continue
    #    line = line.strip()
    #    runlist.append(str(line))

    #runlist = ["364102", "364114", "364115", "364116", "364118", "364200"]
    #runlist = []
    #lines = open("retrylist.txt").readlines()
    #for line in lines :
    #    if not line : continue
    #    line = line.strip()
    #    runlist.append(str(line))
    #runlist = ["279867","280464"]
    #retryfile = "/data/uclhc/uci/user/dantrim/n0234val/resub.txt"
    #lines = open(retryfile).readlines()
    #for line in lines :
    #    if not line : continue
    #    line = line.strip()
    #    runlist.append(line)
    return runlist

def main() :
    print "SubmitCondorSF"

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

        number_of_samples = len(sample_lists)
        print "number of samples: ",number_of_samples
        print s
        #print sample_lists
        process_group= s


        script_name = "submitFile_%s.condor"%process_group
        executable_name = "RunCondorSF_%s.sh"%process_group
        build_condor_script(script_name, executable_name, process_group, number_of_samples, doBrick, doLocal, doSDSC, doUC)
        build_job_executable(executable_name, process_group, number_of_samples)


        # Superflow run mode
        run_mode = "-c "

        # build the job command
        run_cmd = "ARGS="
        run_cmd += '"'
        run_cmd += ' %s '%out_dir
        run_cmd += ' %s '%log_dir
        run_cmd += ' %s '%ana_name
        run_cmd += ' n0234val '
        run_cmd += ' %s '%process_group
        run_cmd += ' %s'%run_mode # put here any extra cmd line options for Superflow executable
        run_cmd += '"'
        run_cmd += ' condor_submit %s '%script_name
        run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")

        run_cmd += ' -append "output = %s%s" '%(log_dir, process_group + ".out")
        run_cmd += ' -append "log = %s%s" '%(log_dir, process_group + ".log")
        run_cmd += ' -append "error = %s%s" '%(log_dir, process_group + ".err")

        print run_cmd
        subprocess.call(run_cmd, shell=True)


#        for dataset in sample_lists :
#            print dataset
#            ds_injob = in_job_filelist_dir + dataset.split(in_job_filelist_dir)[1]
#            print ds_injob
#            sys.exit()

            

#            fullname = str(os.path.abspath(dataset))
#            dataset_original = dataset
#            print "    > %s"%dataset
#
#            dataset = "." + dataset[dataset.find(in_job_filelist_dir):]
#            print "    >> %s"%dataset
#
#            if not (str(os.path.abspath(out_dir)) == str(os.environ['PWD'])) :
#                print "You must call this script from the output directory where the ntuples will be stored!"
#                print " >>> Expected submission directory : %s"%os.path.abspath(out_dir)
#                sys.exit()
#
#            do_split_sample = False
#            for dsid_split in samples_to_split :
#                if dsid_split in dataset :
#                    do_split_sample = True
#
#            run_mode = "-c"

#            if not do_split_sample :
#
#                # submit the job as usual
#                run_cmd = "ARGS="
#                run_cmd += '"'
#                run_cmd += ' %s '%out_dir
#                run_cmd += ' %s '%log_dir
#                run_cmd += ' %s '%ana_name
#                #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
#                run_cmd += ' n0234val '
#                run_cmd += ' %s '%dataset
#                run_cmd += ' %s '%run_mode # any extra cmd line optino for Superflow executable
#                run_cmd += '"'
#                run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
#                lname = dataset.split("/")[-1].replace(".txt", "")
#                run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
#                run_cmd += ' -append "output = %s%s" '%(log_dir, lname + ".out")
#                run_cmd += ' -append "log = %s%s" '%(log_dir, lname + ".log")
#                run_cmd += ' -append "error = %s%s" '%(log_dir, lname + ".err")
#
#                print run_cmd
#                subprocess.call(run_cmd, shell=True)
#
#
#            elif do_split_sample :
#                split_suffix = 0
#
#                split_files = []
#                lines = open(dataset_original).readlines()
#                for line in lines :
#                    if not line : continue
#                    line = line.strip()
#                    split_files.append(line)
#
#                for split_file in split_files :
#                    print "    >>> Sub-file [%d] %s"%(split_suffix, split_file)
#
#                    run_cmd = "ARGS="
#                    run_cmd += '"'
#                    run_cmd += ' %s '%out_dir
#                    run_cmd += ' %s '%log_dir
#                    run_cmd += ' %s '%ana_name
#                    #run_cmd += ' %s '%(tar_location + "area.tgz.tgz")
#                    run_cmd += ' n0234val '
#                    run_cmd += ' %s '%split_file
#                    run_cmd += ' %s --sumw --suffix %d'%(run_mode, split_suffix) # any extra cmd line optino for Superflow executable
#                    run_cmd += '"'
#                    run_cmd += ' condor_submit submitFile_TEMPLATE.condor '
#                    lname = dataset.split("/")[-1].replace(".txt", "")
#                    run_cmd += ' -append "transfer_input_files = %s" '%(tar_location + "area.tgz")
#                    run_cmd += ' -append "output = %s%s" '%(log_dir, lname + "_%d.out"%split_suffix)
#                    run_cmd += ' -append "log = %s%s" '%(log_dir, lname + "_%d.log"%split_suffix)
#                    run_cmd += ' -append "error = %s%s" '%(log_dir, lname + "_%d.err"%split_suffix)
#
#                    print run_cmd
#                    subprocess.call(run_cmd, shell=True)
#
#                    split_suffix = split_suffix + 1
                    

def look_for_tarball() :
    if not os.path.isfile("area.tgz") :
        print "Tarball not found."
        sys.exit()

def build_condor_script(script_name, executable_name, proc_group_name, n_samples_in_group, doBrick, doLocal, doSDSC, doUC) :

    brick = 'false'
    local = 'false'
    sdsc  = 'false'
    uc    = 'false'
    if doBrick : brick = 'true'
    if doLocal : local = 'true'
    if doSDSC  : sdsc  = 'true'
    if doUC    : uc    = 'true'

    f = open(script_name, "w")
    f.write('universe = vanilla\n')
    f.write('+local=%s\n'%brick)
    f.write('+site_local=%s\n'%local)
    f.write('+sdsc=%s\n'%sdsc)
    f.write('+uc=%s\n'%uc)
    f.write('executable = %s\n'%executable_name)
    f.write('arguments = $(Process) $ENV(ARGS)\n')
    f.write('should_transfer_files = YES\n')
    f.write('when_to_transfer_output = ON_EXIT\n')
    f.write('use_x509userproxy = True\n')
    f.write('notification = Never\n')
    f.write('queue %d\n'%n_samples_in_group)
    f.close()

def build_job_executable(executable_name, process_group, n_samples_in_group) :
    f = open(executable_name, 'w')
    f.write('#!/bin/bash\n\n\n')
    f.write('echo " -------- RunCondorSF %s -------- "\n'%process_group)
    f.write('process_no=${1}\n')
    f.write('output_dir=${2}\n')
    f.write('log_dir=${3}\n')
    f.write('sflow_exec=${4}\n')
    f.write('stored_dir=${5}\n')
    f.write('group_name=${6}\n')
    f.write('sflow_options=${@:7}\n\n')
    f.write('echo "    process number      : ${process_no}"\n')
    f.write('echo "    output directory    : ${output_dir}"\n')
    f.write('echo "    log directory       : ${log_dir}"\n')
    f.write('echo "    sflow executable    : ${sflow_exec}"\n')
    f.write('echo "    tarred dir          : ${stored_dir}"\n')
    f.write('echo "    process group       : ${group_name}"\n')
    f.write('echo "    sflow options       : ${sflow_options}"\n\n')
    f.write('while (( "$#" )); do\n')
    f.write('    shift\n')
    f.write('done\n\n')
    f.write('work_dir=${PWD}\n')
    f.write('echo "untarring area.tgz"\n')
    f.write('tar -xzf area.tgz\n\n')
    f.write('echo "done untarring"\n')
    f.write('echo " > current directory structure:"\n')
    f.write('ls -ltrh\n\n')
    f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
    f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
    f.write('echo " > calling : cd ${stored_dir}"\n')
    f.write('cd ${stored_dir}\n')
    f.write('echo "Directory structure:"\n')
    f.write('ls -ltrh\n')
    f.write('lsetup fax\n')
    f.write('source susynt-read/bash/setup_root.sh\n')
    f.write('echo " > calling : source RootCore/local_setup.sh"\n')
    f.write('source RootCore/local_setup.sh\n')
    f.write('ls ./filelists/${group_name} > joblist_${group_name}.txt\n')
    f.write('echo "Built in-job filelist for group ${group_name}:"\n')
    f.write('cat joblist_${group_name}.txt\n')
    f.write('echo " > calling : python ./Superflow/run/GetFileList.py ${group_name} ${process_no} > injob_filelist_${group_name}_${process_no}.txt"\n')
    f.write('python ./Superflow/run/GetFileList.py ${group_name} ${process_no} > injob_filelist_${group_name}_${process_no}.txt\n')
    f.write('ls -ltrh\n')
    f.write('input_list_for_process=$(head -1 injob_filelist_${group_name}_${process_no}.txt)\n')
    f.write('echo "input_list_for_process : ${input_list_for_process}"\n')

    f.write('input_list_for_process="./${stored_dir}$input_list_for_process"\n')
    f.write('echo "Found input for group: ${group_name} and process no: ${process_no}"\n')
    f.write('echo "  > ${input_list_for_process}"\n\n')

    f.write('echo " > calling : python ./Superflow/run/GetJobLogName.py ${input_list_for_process} > injob_log_${group_name}_${process_no}.txt"\n')
    f.write('python ./Superflow/run/GetJobLogName.py ${input_list_for_process} ${process_no} > injob_log_${group_name}_${process_no}.txt\n') 
    f.write('ls -ltrh\n')

    f.write('log_for_process=$(head -1 injob_log_${group_name}_${process_no}.txt)\n') 
    f.write('echo "Setting log to: ${log_for_process}"\n')
    #f.write('cd Superflow/\n')
    #f.write('source setRestFrames.sh\n')
    f.write('echo " > calling : cd ${work_dir}"\n')
    f.write('cd ${work_dir}\n')
    f.write('ls -ltrh\n')
    f.write('echo " > calling : ${sflow_exec} -i ${input_list_for_process} ${sflow_options} 2>&1 |tee ${log_for_process}"\n')
    f.write('${sflow_exec} -i ${input_list_for_process} ${sflow_options} 2>&1 |tee ${log_for_process}\n')
    f.write('echo "final directory structure:"\n')
    f.write('ls -ltrh\n')
    f.close()


if __name__=="__main__" :
    main()

