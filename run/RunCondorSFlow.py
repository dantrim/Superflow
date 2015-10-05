#!/bin/env python

import sys
import os
import subprocess
import glob


ana_name = "restframesAna" # name of executable (i.e. your executable running Superflow)
start_dir = "/data/uclhc/uci/user/dantrim/n0216val/"
log_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0216/zee/data/logs/" # output directory for your job and condor logs
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0216/zee/data/raw/" # output directory for your output ntuples

filelist_dir = "/data/uclhc/uci/user/dantrim/n0216val/filelists/data/"
#filelist_dir = "/data/uclhc/uci/user/dantrim/n0216val/filelists/mc/"
samples = ["data"] # name of filelist inside of 'filelist_dir' for samples to run 
#samples = ["zee_sherpa"] # name of filelist inside of 'filelist_dir' for samples to run 

def main() :
    print "RunCondorSFlow"
    print "Looking for sample filelists in main filelist directory : %s"%filelist_dir

    for s in samples :
        print "Submitting sample : %s"%s
        suffix_ = ""
        if not s.endswith("/") : suffix_+= "/"
        sample_lists = glob.glob(filelist_dir + s + suffix_ + "*.txt")
        for dataset in sample_lists :
            print "    > %s"%dataset
            name  = dataset.split("/")
            name  = name[-1]
            name  = name.replace(".txt", "")
            condor_name = name + ".condor"
            condor_sub_file = open(condor_name, "w")
            get_sub_script(ana_name, condor_sub_file, name, dataset, start_dir, out_dir)
            condor_sub_file.close()
            cmd = "condor_submit %s"%condor_name
            subprocess.call(cmd, shell=True)
            cmd = "rm %s"%condor_name
            subprocess.call(cmd, shell=True)

def get_sub_script(ana = "", file_ = None, name = "", dataset = "", startdir = "", outputdir = "") :
    """
    Write a condor submission script with the output logs, etc...
    associated with the currently submitted dataset.

    'ana'       = executable name (your superflow executable)
    'file_'     = the condor submission script that is being written and to be called
    'name'      = dataset descriptive name to be used in the log names 
    'dataset'   = input .txt file to be provided to your executable 'ana' that
                  contains the FAX container names for this sample/dataset
    'startdir'  = main working directory (RootCore dir) for Superflow
    'outputdir' = directory where the output ntuples will end up

    See [1] for documentation on setting up the condor submission scripts.

    [1] http://www.t2.ucsd.edu/twiki2/bin/view/UCLHCWeb/UCIUserDoc#Job_Submission
    """

    if ana == "" :
        print "get_sub_script ERROR   input analysis name (i.e. executable) is empty."
        sys.exit()
    if not file_ :
        print "get_sub_script ERROR    input file is None."
        sys.exit()
    if name == "" :
        print "get_sub_script ERROR    input name is empty."
        sys.exit()
    if dataset == "" :
        print "get_sub_script ERROR    dataset is empty."
        sys.exit()
    if startdir == "" :
        print "get_sub_script ERROR    startdir is empty."
        sys.exit()
    if outputdir == "" :
        print "get_sub_script ERROR    outputdir is empty."
        sys.exit()

    file_.write("universe = vanilla\n")
    file_.write("+local = true\n")      # run on the brick
    file_.write("+site_local = true\n") # run on local batch system of the site
    file_.write("+sdsc = true\n")       # run at SDSC Comet Cluster
    file_.write("+uc = true\n")         # run outside to all UC's
    file_.write("executable = CondorSFlow.sh\n")
    file_.write("arguments = %s %s %s %s\n"%(ana, dataset, startdir, outputdir))
    file_.write("should_transfer_files = YES\n")
    file_.write("when_to_transfer_output = ON_EXIT\n")
    file_.write("log = %s/%s.log\n"%(log_dir, name))
    file_.write("output = %s/%s.out.$(Cluster).$(Process)\n"%(log_dir, name))
    file_.write("error = %s/%s.err.$(Cluster).$(Process)\n"%(log_dir, name))
    file_.write("use_x509userproxy = True\n")
    file_.write("notification = Never\n")
    file_.write("queue\n")

if __name__=="__main__" :
    main()




