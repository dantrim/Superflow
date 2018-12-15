#!/bin/env python

import os
import sys
import glob
import argparse
import subprocess
import time

sf_exec_name = "ntupler_nn"
sf_exec_name = "ntupler_rj_stop2l"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0303/data17/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0303/b_oct7/mc/"
out_dir = "/data/uclhc/uci/user/dantrim/ntuples/n0303/c_oct8/mc16d/"

filelist_dir = "/data/uclhc/uci/user/dantrim/n0303val/susynt-read/filelists/"

#samples = [ 'ttbar_mc16a', 'WtPP8_mc16a', 'WtPP6_mc16a', 'diboson_sherpa_ll_mc16a', 'zjets_sherpa_mc16a', 'drellyan_sherpa_mc16a', 'wjets_sherpa_mc16a', 'rare_top_mc16a', 'single_higgs_mc16a', 'ttX_mc16a' ]
#samples = [ 'ttbar_mc16d', 'WtPP8_mc16d', 'WtPP6_mc16d', 'diboson_sherpa_ll_mc16d', 'zjets_sherpa_mc16d', 'drellyan_sherpa_mc16d', 'wjets_sherpa_mc16d', 'rare_top_mc16d', 'single_higgs_mc16d', 'ttX_mc16d' ]
#samples = [ 'ttbar_mc16a' ]
#samples = [ 'n0303_data1516' ]
samples = [ 'n0303_data17' ]
samples = [ 'ggllvvZZ' ]
samples = [ 'zv_mc16a' ]
samples = [ 'diboson_sherpa_ll_mc16d' ]

dsids_to_split = [ '410472' ] #, '364254' ] #'364254' ]
#dsids_to_split = [ '364254' ] #'364254' ]

run_with_systematics = False
tar_file = '/data/uclhc/uci/user/dantrim/n0303val/area.tgz'
tarred_dir = 'susynt-read'
use_sumw_file = True # True for multi-period running

do_brick = True
do_gp = True
do_uc = True
do_sdsc = False

def special_prefixes() :

    return { #'mwt2' : 'root://fax.mwt2.org:1094/',
            'mwt2' : 'root://fax.mwt2.org:1094//pnfs/uchicago.edu/atlaslocalgroupdisk/rucio/',
              'atlas-xrd' : 'root://atlas-xrd-us.usatlas.org/',
                'griddev' : 'root://griddev03.slac.stanford.edu:2094//xrootd/atlas/atlasscratchdisk/rucio',
                #'griddev' : 'root://griddev03.slac.stanford.edu:2094/',
                'env' : os.environ.get('STORAGEPREFIX','') }

def sflow_exec_arg_string() :

    global use_sumw_file
    global run_with_systematics

    sys_string = '-c'
    if run_with_systematics :
        sys_string = '-a'

    sflow_args = ' %s ' % sys_string
    if use_sumw_file :
        sflow_args += ' --sumw ./susynt-read/data/sumw_file.root '

    return sflow_args

def create_prefixed_lists(base_list_dir, storage_prefix) :

    if storage_prefix in special_prefixes() :
        if special_prefixes()[storage_prefix] :
            storage_prefix = special_prefixes()[storage_prefix]
        else :
            print 'ERROR Requested special prefix (=%s) is invalid' % storage_prefix
            sys.exit()

    os.chdir(base_list_dir)

    filelist_sample_dirs = glob.glob("./filelists_base/*")
    filelist_sample_dirnames = [s.split("/")[-1] for s in filelist_sample_dirs]

    for f in filelist_sample_dirnames :
        cmd = 'mkdir -p ./filelists/%s' % f
        subprocess.call(cmd, shell = True)

    for idir, fdir in enumerate(filelist_sample_dirnames) :
        text_files = glob.glob('./filelists_base/%s/*.txt' % fdir)
        for text_file in text_files :
            new_file = './filelists/%s/%s' % (filelist_sample_dirnames[idir], text_file.split('/')[-1])
            with open(new_file, 'w') as out_file :
                with open(text_file, 'r') as in_file :
                    for line in in_file :
                        line = line.strip()
                        if not line : continue
                        if line.startswith('#') : continue
                        new_did = '%s%s' % (storage_prefix, line)
                        out_file.write(new_did + '\n')

def create_tar(args) :

    print 'Creating tar file'

    global tar_file

    tar_name = args.tar_name
    loc_to_place = '/'.join(tar_file.split("/")[:-1])

    if not os.path.isdir(loc_to_place) :
        print 'ERROR Cannot locate directory to place tar file (attempted dir = %s)' % loc_to_place
        sys.exit()

    full_name = '%s/%s' % (loc_to_place, tar_name)

    if os.path.isfile(full_name) :
        print 'ERROR Tar file already exists, will not continue (tar name = %s)' % full_name
        sys.exit()

    full_dir_to_tar = '%s/%s' % (loc_to_place, tarred_dir)
    if not os.path.isdir(full_dir_to_tar) :
        print 'ERROR Directory to tar (=%s) not found in expected directory (=%s)' % (tarred_dir, loc_to_place)
        sys.exit()



    if args.skip_list_creation :
        if not os.path.isdir('%s/%s/filelists/' % (loc_to_place, tarred_dir) ) :
            print 'ERROR \"filelists\" directory is not found in expected directry (=%s), cannot build tar file' % loc_to_place
            sys.exit()
    elif not os.path.isdir('%s/%s/filelists_base/' % (loc_to_place, tarred_dir)) :
        print 'ERROR \"filelists_base\" directory is not found, cannot create filelists for tar creation'
        sys.exit()


    # move to the dir to dump the tar file
    if args.verbose :
        print 'Moving to dir to create tar: %s' % loc_to_place
    os.chdir(loc_to_place)

    if not args.skip_list_creation :
        base_list_dir = '%s/%s' % (loc_to_place, tarred_dir)
        create_prefixed_lists(base_list_dir, args.prefix)
        os.chdir(loc_to_place)

    things_to_tar = ['build/*'] #, 'filelists/*'] #, 'sumw_file.root']
    things_to_tar.append('data/*')

    global sf_exec_name
    if '_rj_' in sf_exec_name :
        things_to_tar.append('RestFrames/lib/')
        things_to_tar.append('RestFrames/setup_RestFrames_BATCH.sh')
        things_to_tar.append('source/RJTupler/scripts/')

    if not os.path.isdir('%s/%s/filelists/' % (loc_to_place, tarred_dir)) :
        print 'ERROR \"filelists\" directory was not made as expected!'
        sys.exit()

    things_to_tar.append('filelists/*')

    if os.path.isfile('%s/%s/sumw_file.root' % (loc_to_place, tarred_dir)) :
        things_to_tar.append('sumw_file.root')
    things_to_tar = ['./%s/%s' % (tarred_dir, s) for s in things_to_tar]
    tar_cmd = 'tar zcf'
    if args.verbose :
        tar_cmd += 'v'
    tar_cmd += ' %s %s' % (tar_name, ' '.join(things_to_tar))
    subprocess.call(tar_cmd, shell = True)

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

    global tar_file
    global run_with_systematics

    with open(condor_filename, 'w') as f :

        f.write('universe = vanilla\n')
        f.write('+local=%s\n' % bool_string(do_brick))
        f.write('+site_local=%s\n' % bool_string(do_gp))
        f.write('+uc=%s\n' % bool_string(do_uc))
        f.write('+sdsc=%s\n' % bool_string(do_sdsc))
        f.write('executable = %s\n' % exec_name)
        f.write('should_transfer_files = YES\n')
        f.write('transfer_input_files = %s\n' % tar_file)
        f.write('use_x509userproxy = True\n')
        f.write('notification = Never\n')
        if run_with_systematics : 
            f.write('request_memory = 4 GB\n')

        for q in queue_list :

            injob_filelist = './' + filelist_dir[filelist_dir.find(tarred_dir):] + '/'
            injob_filelist += sample
            injob_filelist += '/%s' % q.split('/')[-1]

            log_base = q.split('/')[-1].replace('.txt','')
            sys_string = '-c'
            if run_with_systematics :
                sys_string = '-a'

            arg_string = ' %s %s %s %s ' % (sf_exec_name, tarred_dir, injob_filelist, q)

            sflow_args = sflow_exec_arg_string()

            arg_string += sflow_args

            f.write('\n')
            f.write('arguments = %s\n' % arg_string)
            f.write('output = log_%s.out\n' % log_base)
            f.write('log = log_%s.log\n' % log_base)
            f.write('error = log_%s.err\n' % log_base)
            f.write('queue\n')

def make_executable(exec_name) :

    global sf_exec_name

    with open(exec_name, 'w') as f :
        f.write('#!/bin/bash\n\n\n')
        f.write('echo "--------- %s ----------"\n' % exec_name)
        f.write('hostname\n')
        f.write('echo "start: `date`"\n')
        f.write('echo "input arguments:"\n')
        f.write('sf_exec=${1}\n')
        f.write('tarred_dir=${2}\n')
        f.write('injob_filelist=${3}\n')
        f.write('sample_list=${4}\n')
        f.write('sflow_options=${@:5}\n')
        f.write('echo "   SF executable         : ${sf_exec}"\n')
        f.write('echo "   tarred directory      : ${tarred_dir}"\n')
        f.write('echo "   injob filelist loc    : ${injob_filelist}"\n')
        f.write('echo "   sample list name      : ${sample_list}"\n')
        f.write('echo "   sflow options         : ${sflow_options}"\n\n')
        f.write('while (( "$#" )); do\n')
        f.write('   shift\n')
        f.write('done\n\n')
        f.write('echo "untarring area.tgz"\n')
        f.write('tar -xvf area.tgz\n\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n\n')
        f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
        f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
        f.write('echo "moving"\n') 
        f.write('pushd ${tarred_dir}\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('lsetup fax\n')
        f.write('asetup AnalysisBase,21.2.45\n')
        f.write('source build/x86*/setup.sh\n')
        if '_rj_' in sf_exec_name :
            f.write('source ./source/RJTupler/scripts/setup_restframes.sh --batch\n')
        f.write('echo "moving"\n')
        f.write('popd\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('echo "calling: ${sf_exec} -i ${injob_filelist} ${sflow_options}"\n')
        f.write('${sf_exec} -i ${injob_filelist} ${sflow_options}\n')
        f.write('echo "final directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('echo "finish: `date`"\n')

def make_condor_file_split(dsid, split_files, condor_filename, exec_name) :

    global tar_file
    with open(condor_filename, 'w') as f :
        f.write('universe = vanilla\n')
        f.write('+local=%s\n' % bool_string(do_brick))
        f.write('+site_local=%s\n' % bool_string(do_gp))
        f.write('+uc=%s\n' % bool_string(do_uc))
        f.write('+sdsc=%s\n' % bool_string(do_sdsc))
        f.write('executable = %s\n' % exec_name)
        f.write('should_transfer_files = YES\n')
        f.write('transfer_input_files = %s\n' % tar_file)
        f.write('use_x509userproxy = True\n')
        f.write('notification = Never\n')

        for ifile, split_file in enumerate(split_files) :

            arg_string = ' %d %s %s %s %s ' % (ifile, sf_exec_name, tarred_dir, split_file, dsid)
            sflow_args = sflow_exec_arg_string()
            arg_string += sflow_args

            f.write('\n')
            f.write('arguments = %s\n' % arg_string)
            f.write('output = log_%s_%d.out\n' % (str(dsid), ifile))
            f.write('log = log_%s_%d.log\n' % (str(dsid), ifile))
            f.write('error = log_%s_%d.err\n' % (str(dsid), ifile))
            f.write('queue\n')

def make_executable_split(exec_name) :

    global sf_exec_name

    with open(exec_name, 'w') as f :
        f.write('#!/bin/bash\n\n\n')
        f.write('echo "------------ %s -------------"\n' % exec_name)
        f.write('hostname\n')
        f.write('echo "start: `date`"\n')
        f.write('split_idx=${1}\n')
        f.write('sf_exec=${2}\n')
        f.write('tarred_dir=${3}\n')
        f.write('filename=${4}\n')
        f.write('dsid=${5}\n')
        f.write('sflow_options=${@:6}\n')
        f.write('echo "    split idx        : ${split_idx}"\n')
        f.write('echo "    SF executable    : ${sf_exec}"\n')
        f.write('echo "    tarred directory : ${tarred_dir}"\n')
        f.write('echo "    split file       : ${filename}"\n')
        f.write('echo "    sample dsid      : ${dsid}"\n')
        f.write('echo "    sflow options    : ${sflow_options}"\n')
        f.write('while (( "$#" )); do\n')
        f.write('   shift\n')
        f.write('done\n')
        f.write('echo "untarring area.tgz"\n')
        f.write('tar -xvf area.tgz\n\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n\n')
        f.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
        f.write('source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n')
        f.write('echo "moving"\n')
        f.write('pushd ${tarred_dir}\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('lsetup fax\n')
        f.write('asetup AnalysisBase,21.2.45\n')
        f.write('source build/x86*/setup.sh\n')
        if '_rj_' in sf_exec_name :
            f.write('echo "sourcing RJ:"\n')
            f.write('ls source/RJTupler/scripts/\n')
            f.write('source source/RJTupler/scripts/setup_restframes.sh --batch\n')
        f.write('echo "moving"\n')
        f.write('popd\n')
        f.write('echo "current directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('echo "calling: ${sf_exec} -i ${injob_filelist} ${sflow_options}"\n')
        f.write('${sf_exec} -i ${filename} ${sflow_options} --suffix ${split_idx}\n')
        f.write('echo "final directory structure:"\n')
        f.write('ls -ltrh\n')
        f.write('echo "finish: `date`"\n')

def submit_sample(sample, verbose = False) :

    txt_files_for_sample = glob.glob(filelist_dir + "/" + sample + "/*.txt")
    n_txt = len(txt_files_for_sample)


    if verbose :
        print 'Found %d text files for sample %s' % (n_txt, sample)

    if n_txt == 0 :
        print 'WARNING Found no text files for sample %s, skipping this sample' % sample
        return

    split_dict, queue_list = get_split_samples(txt_files_for_sample)

    start_path = os.getcwd()

    if queue_list :

        if verbose :
            print 'Making \"output\" directory'

        cmd = 'mkdir -p output'
        subprocess.call(cmd, shell = True)
        full_raw = os.path.abspath("./output")
        os.chdir(full_raw)
        if verbose :
            print 'Current dir: %s' % os.getcwd()

        condor_filename = 'submit_%s.condor' % sample
        exec_name = 'run_condor_sf_%s.sh' % sample
        make_condor_file(sample, queue_list, condor_filename, exec_name)
        make_executable(exec_name)

        cmd = 'condor_submit %s' % condor_filename
        subprocess.call(cmd, shell = True)

        os.chdir(start_path)
        if verbose :
            print 'Current dir: %s' % os.getcwd()

    elif len(split_dict.keys()) :

        if verbose :
            print 'Launching split samples'

        os.chdir(start_path)

        for dsid in split_dict :


            split_dir = 'output_%s' % str(dsid)
            cmd = 'mkdir -p %s' % split_dir
            subprocess.call(cmd, shell = True)

            full_split_dir = os.path.abspath(split_dir) 
            os.chdir(full_split_dir)

            if verbose :
                print 'Current dir: %s' % os.getcwd()

            condor_filename = 'submit_split_%s.condor' % str(dsid)
            exec_name = 'run_condor_split_%s.sh' % str(dsid)
            make_condor_file_split(dsid, split_dict[dsid], condor_filename, exec_name)
            make_executable_split(exec_name)

            cmd = 'condor_submit %s' % condor_filename
            subprocess.call(cmd, shell = True)
            os.chdir(start_path)
            if verbose :
                print 'Current dir: %s' % os.getcwd()
            


def get_samples(args) :

    global samples

    if args.sample != "" :
        samples = [s.strip() for s in args.sample.split(",")]
    return samples

def main() :

    parser = argparse.ArgumentParser( description = 'Launch SuperFlow jobs to the condor batch' )
    parser.add_argument('--sample', default = "",
        help = 'Override samples list by providing a sample name (can be comma-separated-list)')
    parser.add_argument('-v', '--verbose', default = False, action = 'store_true',
        help = 'Turn on verbose mode')
    parser.add_argument('-t', '--tar', default = False, action = 'store_true',
        help = 'Create tar file for sending to job location')
    parser.add_argument('--tar-name', default = 'area.tgz',
        help = 'Set name of tar file')
    parser.add_argument('--prefix', default = 'root://fax.mwt2.org:1094/',
        help = 'Provide FAX STORAGEPREFIX to use (used when creating filelists for tar file creation')
    parser.add_argument('--prefixes', action = 'store_true', default = False,
        help = 'Print the list of known/common FAX prefixes')
    parser.add_argument('--skip-list-creation', default = False, action = 'store_true',
        help = 'Do not build filelists during tar file creation, use already existing \"filelist/\" dir')
    args = parser.parse_args()

    if args.prefixes :
        for key, prefix in special_prefixes().items() :
            print ' %s : %s' % (key, prefix)
        return

    pwd = os.getcwd()
    if os.path.abspath(pwd) != os.path.abspath(out_dir) :
        print 'ERROR You must call this script from the designated output directory (out_dir = %s)' % os.path.abspath(out_dir)
        sys.exit()

    if args.tar :
        create_tar(args)

    else :
        samples = get_samples(args)
        check_samples(samples, filelist_dir)

        n_samples = len(samples)
        for isample, sample in enumerate(samples) :
            print '\n[%02d/%02d] Submitting %s' % (isample+1, n_samples, sample)
            submit_sample(sample, args.verbose)

#____________________________
if __name__ == '__main__' :
    main()
