#!/usr/bin/env python
"""
================================================================================
Launch jobs to the condor batch system
Examples
    python submit_to_condor.py path/to/files/*txt -t ./tar_this/ -e superflow_exec -o outputs/

    - If relevant environment variables are set:
    python submit_to_condor.py path/to/files/*txt

Works in general for any executable that takes an input with '-i'. Any
additional executable arguments will require adding 

Author:
    Alex Armstrong <alarmstr@cern.ch>
    with much code borrowed from Daniel Antrim <dantrim@cern.ch>
================================================================================
"""

import sys, os, traceback, argparse
import time
import subprocess
import pdb

################################################################################
# Globals
################################################################################
_include_jigsaw = False
_testing = False
_condor_submit_name = 'submit.condor'
_condor_exec_name = 'run_condor.sh'
_superflow_executables = [
    # Alex executables
    'makeFlatNtuples',
    # Danny executables
    'ntupler_nn',
    'ntupler_rj_stop2l',
]
_grabSumw_executable = 'grabSumw'
# Available sites for condor submissions
_do_brick = True
_do_gp = True
_do_uc = True
_do_sdsc = False # We do not have the necessary permissions, jobs will hang

# User Argument defaults and help information
_help_input_files   = 'Input txt files containing xrootd links'

_df_tar_dir         = os.getenv('TAR_DIR', '$TAR_DIR')
_help_tar_dir       = 'Directory to be tarred and sent with job \
                       [default: %s]' % _df_tar_dir

_df_tar_file        = os.getenv('TAR_FILE', 'area.tgz')
_help_tar_file      = 'Path to tar file to use or create for sending with job \
                       [default: %s]' % _df_tar_file

_df_exec            = os.getenv('CONDOR_EXEC', '$CONDOR_EXEC')
_help_exec          = 'Name of executable to run during job \
                       [default: %s]' % _df_exec

_df_sumw            = os.getenv('SUMW_FILE', '$SUMW_FILE')
_help_sumw          = 'Sum of weights file for multi-period running \
                       [default: %s]' % _df_sumw

_df_output_dir      = os.getenv('BATCH_OUTPUT_DIR', './')
_help_output_dir    = 'Directory for storing job output \
                       [default: %s]' % _df_output_dir

_df_split_dsids     = [] 
_help_split_dsids   = 'DSIDs of input samples that will have one job run per file in the sample'

_help_syst          = 'Run with systematics'

_help_verbose       = 'verbose output'

################################################################################
# MAIN
################################################################################
def main ():
    """
    Main Function

    Assume all unix path input arguments have been expanded into absolute paths
    """

    # Run checks before processing inputs
    global args
    check_environment()
    make_filelists, make_tar = check_inputs(args)

    # Get files to be submitted with jobs
    common_dir = os.path.commonprefix(args.input_files)
    while not os.path.isdir(common_dir):
        common_dir = os.path.dirname(common_dir)
    filelist_dir = common_dir

    # Create tar file if necessary
    if make_tar:
        things_to_tar = get_things_to_tar(args.tar_dir, filelist_dir, args.sumw, _include_jigsaw)
        create_tar(args.tar_dir, args.tar_file, things_to_tar, args.verbose)

    # Submit the jobs

    cwd = os.getcwd()
    os.chdir(args.output_dir)
    submit_jobs(args.input_files,
                args.split_dsids,
                args.executable,
                args.tar_file,
                args.tar_dir,
                args.syst,
                args.verbose)
    os.chdir(cwd)


def create_tar(tar_dir, tar_file, things_to_tar=["./*"], verbose=False) :
    """
    Tar specific files inside a directory
    args:
        tar_dir (str) - directory to be tarred
        things_to_tar (list(str) - things inside tar_dir that will get tarred.
            Default is to tar the whole directory.
        tar_file (str) - output tar file path
        verbose (bool) - run tar command in verbose mode

    """
    tar_dir = os.path.abspath(tar_dir)
    tar_file = os.path.abspath(tar_file)

    # Move to the parent directory of dir to be tarred
    # Save PWD so it can be returned to afterwards
    pwd = os.getcwd()
    par_dir = os.path.dirname(tar_dir)
    tar_dir_name = os.path.basename(tar_dir)
    os.chdir(par_dir)

    # Check directory structure is as expected
    paths_to_tar = [os.path.join(tar_dir_name, x) for x in things_to_tar]
    for x in paths_to_tar:
        x_tmp = x.replace("*","")
        if not os.path.exists(x_tmp):
            print "ERROR :: Requested tar item not found:", x_tmp
        elif os.path.abspath(tar_dir) not in os.path.abspath(x_tmp):
            print "ERROR :: Requested tar item not in tar directory:", x_tmp

    # Check that directory for storing tar file exists
    if not os.path.exists(par_dir):
        print 'ERROR :: Cannot locate directory to place tar file',
        print '(attempted dir = %s)' % par_dir
        sys.exit()

    # Check that directory for storing tar file is not inside the directory to
    # be tarred
    for p in paths_to_tar:
        if os.path.abspath(p) in os.path.abspath(tar_file):
            print "ERROR :: Attempting to tar directory that contains the output",
            print "tar file. Avoid this."
            sys.exit()

    print 'INFO :: Creating tar file'
    print 'INFO :: The following directories from %s will be included:' % tar_dir_name
    for x in things_to_tar:
        print "INFO ::\t\t", x

    tar_cmd = 'tar zcf'
    if verbose :
        tar_cmd += 'v'
    tar_cmd += ' %s %s' % (tar_file, ' '.join(paths_to_tar))
    subprocess.call(tar_cmd, shell = True)
    print "INFO :: Done! Tar file created at %s" % os.path.abspath(tar_file)

    # Return to original directory in case directory was changed
    os.chdir(pwd)

def submit_jobs(input_files, dsids_to_split, exec_name, tar_file, tar_dir, syst, verbose = False) :
    '''

    args:
        tar_dir (list(str)) -
        dsids_to_split (list(str)) -
        exec_name (str) -
        tar_file (str) - path to tarred file for submitting with job
        tar_dir (str) - absolute path to directory that was tarred in tar_file
        syst (bool) - run exectuable with systematics
        verbose (bool) - run with verbose output

    '''
    # Get paths relative to the tar directory for running on job site
    tar_name = os.path.basename(tar_file)
    tar_dir_base = os.path.basename(tar_dir)

    # Build condor file header
    condor_file_str = build_condor_file_header(_condor_exec_name, tar_file, syst)

    # Build condor file queues
    for f in input_files:
        rel_file_path = os.path.relpath(f, os.path.dirname(tar_dir))
        if any(dsid in f for dsid in dsids_to_split):
            condor_file_str += build_condor_file_split_queues(
                f, exec_name, tar_dir_base, syst)
        else:
            condor_file_str += build_condor_file_queues(
                rel_file_path, exec_name, tar_dir_base, syst)

    # Write condor submit file
    with open(_condor_submit_name, 'w') as ofile:
        ofile.write(condor_file_str)

    # Build condor executable
    condor_exec_str = build_condor_executable(exec_name, tar_name, _include_jigsaw)

    # Write condor executable
    with open(_condor_exec_name, 'w') as ofile:
        ofile.write(condor_exec_str)

    # Submit condor jobs
    print "\nSubmitting Condor Jobs!\n"
    cmd = 'condor_submit %s' % _condor_submit_name
    subprocess.call(cmd, shell = True)

def build_condor_file_header(exec_name, tar_file, syst):
    '''
    args:
        exec_name (str) - name of executable
        tar_file (str) - path to tarred file for submitting with job
        syst (bool) - run exectuable with systematics
    '''
    def bool_string(boolean) :
        return { True : 'true', False : 'false' } [ bool(boolean) ]

    header_str = ''
    header_str += 'universe = vanilla\n'
    header_str += '+local=%s\n' % bool_string(_do_brick)
    header_str += '+site_local=%s\n' % bool_string(_do_gp)
    header_str += '+uc=%s\n' % bool_string(_do_uc)
    header_str += '+sdsc=%s\n' % bool_string(_do_sdsc)
    header_str += 'executable = %s\n' % exec_name
    header_str += 'should_transfer_files = YES\n'
    header_str += 'transfer_input_files = %s\n' % tar_file
    header_str += 'use_x509userproxy = True\n'
    header_str += 'notification = Never\n'
    if syst:
        f.write('request_memory = 4 GB\n')

    return header_str

def build_condor_file_split_queues(ifile_path, exec_name, tar_dir, syst):
    '''
    args:
        ifile_path (str) - path to input file
        exec_name (str) - name of executable
        tar_dir (str) - name of tarred directory
        syst (bool) - run exectuable with systematics
    '''
    queue_str = ''

    xrootd_links = []
    with open(ifile_path) as f:
        for line in f:
            line = line.strip()
            if skip_txt_line(line): continue
            xrootd_links.append(line)

    log_base = os.path.basename(ifile_path).replace('.txt','')

    for idx, link in enumerate(xrootd_links):
        exec_args = get_exec_arg_string(exec_name, suffix=idx)
        # Positional arguments for condor executable
        # Order is important. See "build_condor_executable" for expected order
        # They get imported as an environment variable
        arg_string = ' %s %s %s %s' % (exec_name, tar_dir, link, exec_args)

        # Build queue string
        queue_str += '\n'
        queue_str += 'arguments = %s\n' % arg_string
        queue_str += 'output = log_%s_%d.out\n' % (log_base, idx)
        queue_str += 'log = log_%s_%d.log\n' % (log_base, idx)
        queue_str += 'error = log_%s_%d.err\n' % (log_base, idx)
        queue_str += 'queue\n'

    return queue_str

def build_condor_file_queues(ifile_path, exec_name, tar_dir, syst):
    '''
    args:
        ifile_path (str) - path to input file relative to tarred directory
        exec_name (str) - name of executable
        tar_dir (str) - name of tarred directory
        syst (bool) - run exectuable with systematics
    '''
    queue_str = ''

    log_base = os.path.basename(ifile_path).replace('.txt','')
    exec_args = get_exec_arg_string(exec_name)

    # Positional arguments for condor executable
    # Order is important. See "build_condor_executable" for expected order
    # They get imported as an environment variable
    arg_string = ' %s %s %s %s' % (exec_name, tar_dir, ifile_path, exec_args)

    # Build queue string
    queue_str += '\n'
    queue_str += 'arguments = %s\n' % arg_string
    queue_str += 'output = log_%s.out\n' % log_base
    queue_str += 'log = log_%s.log\n' % log_base
    queue_str += 'error = log_%s.err\n' % log_base
    queue_str += 'queue\n'

    return queue_str

def build_condor_executable(exec_name, tar_file, jigsaw=False):
    '''
    '''
    exec_str = ''
    # Environment info
    exec_str += '#!/bin/bash\n\n\n'
    exec_str += 'echo "--------- %s ----------"\n' % exec_name
    exec_str += 'hostname\n'
    exec_str += 'echo "start: `date`"\n'
    exec_str += 'echo "input arguments:"\n'

    # Read in inputs to condor executable from submission script
    # Ordering should agree with what argument order set in submission script
    exec_str += 'executable=${1}\n'
    exec_str += 'tarred_dir=${2}\n'
    exec_str += 'injob_filelist=${3}\n'
    exec_str += 'exec_options=${@:4}\n'
    exec_str += 'echo "   executable            : ${executable}"\n'
    exec_str += 'echo "   tarred directory      : ${tarred_dir}"\n'
    exec_str += 'echo "   injob filelist loc    : ${injob_filelist}"\n'
    exec_str += 'echo "   executable options    : ${exec_options}"\n\n'
    exec_str += 'while (( "$#" )); do\n'
    exec_str += '   shift\n'
    exec_str += 'done\n\n'

    # Untar the tarred file
    exec_str += 'echo "untarring %s"\n' % tar_file
    exec_str += 'tar -xvf %s\n\n' % tar_file

    # Setup environment
    exec_str += 'echo "current directory structure:"\n'
    exec_str += 'ls -ltrh\n\n'
    exec_str += 'export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n'
    exec_str += 'source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n'

    # Setup release inside untarred working directory
    exec_str += 'echo "moving"\n'
    exec_str += 'pushd ${tarred_dir}\n'
    exec_str += 'echo "current directory structure:"\n'
    exec_str += 'ls -ltrh\n'
    exec_str += 'lsetup rucio\n'
    exec_str += 'asetup AnalysisBase,21.2.55\n'
    exec_str += 'source build/x86*/setup.sh\n'
    if jigsaw:
        exec_str += 'echo "sourcing RJ:"\n'
        exec_str += 'ls source/RJTupler/scripts/\n'
        exec_str += 'source ./source/RJTupler/scripts/setup_restframes.sh --batch\n'
    exec_str += 'echo "moving"\n'
    exec_str += 'popd\n'
    exec_str += 'echo "current directory structure:"\n'
    exec_str += 'ls -ltrh\n'

    # Run executable
    exec_str += 'echo "calling: ${executable} -i ${injob_filelist} ${exec_options}"\n'
    exec_str += '${executable} -i ${injob_filelist} ${exec_options}\n'

    # Final check for outputs
    exec_str += 'echo "final directory structure:"\n'
    exec_str += 'ls -ltrh\n'
    exec_str += 'echo "finish: `date`"\n'

    return exec_str

################################################################################
# FORMAT SENSATIVE FUNCTIONS
# functions making non-robust assumptions
################################################################################
def file_name_has_xrootd_prefix(file_name):
    ''' Check that the file name has an acceptable xrootd storage prefix'''
    return file_name.startswith("root://")

def acceptable_input_file_name(file_name):
    return file_name.endswith('txt')

def get_things_to_tar(tar_dir, filelist_dir='', sumw_file='', jigsaw=False):
    '''
    Get path of things to be included in tarring relative to the directory to
    be tarred. Assume all checks have been performed to make sure unix paths
    exist.
    '''
    things_to_tar = []
    things_to_tar.append('build/*')
    things_to_tar.append('data/*')

    if filelist_dir:
        full_path = os.path.abspath(filelist_dir)
        rel_dir = os.path.abspath(tar_dir)
        rel_path = os.path.relpath(full_path, rel_dir)
        things_to_tar.append(rel_path)

    if sumw_file:
        full_path = os.path.abspath(sumw_file)
        rel_dir = os.path.abspath(tar_dir)
        rel_path = os.path.relpath(full_path, rel_dir)
        things_to_tar.append(rel_path)

    if jigsaw:
        things_to_tar.append('source/RestFrames/lib/')
        things_to_tar.append('source/RestFrames/setup_RestFrames_BATCH.sh')
        things_to_tar.append('source/RJTupler/scripts/')

    return things_to_tar

def get_exec_arg_string(exec_name, suffix='', syst='', sumw_file=''):
    '''
    '''
    if exec_name in _superflow_executables:
        return sflow_exec_arg_string(syst, sumw_file, suffix)
    elif exec_name == _grabSumw_executable:
        return ''
    else:
        print "WARNING :: Unknown executable name: %s" % exec_name
        print "INFO :: Running executable with only input file option"
        return ''


def sflow_exec_arg_string(syst=False, sumw_file='', suffix=''):
    '''
    '''
    sflow_args = ''

    if _testing:
        sflow_args += ' -n 1000 '

    # Systematics
    sys_string = '-c'
    if syst :
        sys_string = '-a'
    sflow_args = ' %s ' % sys_string

    # Sum of weights for multi-period processing
    if sumw_file :
        sflow_args += ' --sumw %s ' % sumw_file

    if suffix:
        sflow_args += ' --suffix %s ' % suffix

    return sflow_args

def skip_txt_line(line):
    ''' Check if a line read in from a text file should not be processed'''
    l = line.strip()
    return l.startswith("#") or not l

################################################################################
# SUPPORT FUNCTIONS
################################################################################
def check_environment():
    """ Check if the shell environment is setup as expected """
    assert os.environ['USER'], "USER variable not set"

    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    assert python_ver >= 2.7, ("Running old version of python\n", sys.version)

def check_inputs(args):
    """ Check the input arguments are as expected """

    make_tar_file = False
    make_filelists = False
    ############################################################################
    # Check for incompatable combinations of user options

    ############################################################################
    ## Check that program is being run in correct directory
    #pwd = os.getcwd()
    #if os.path.abspath(pwd) != os.path.abspath(args.output_dir) :
    #    #TODO: Have script move to output directory and then return to cwd
    #    print 'ERROR :: Must call script from the output directory',
    #    print '(= %s)' % os.path.abspath(args.output_dir)
    #    sys.exit()

    # Check that input files were provided
    if not args.input_files:
        print "ERROR :: No input files were provied"

    # Check that either txt files or root files were provided
    if not all(f.endswith(".txt") for f in args.input_files):
        print "ERROR :: Some input files were not an expected format (*.txt)"
        sys.exit()

    # Check that all input files exist if text files are provided
    for f in args.input_files:
        if not os.path.exists(f):
            print "ERROR :: Cannot find input file:", f
            sys.exit()
        if not acceptable_input_file_name(f):
            print "ERROR :: Unpexected input file format: ", f
            print "INFO :: Expecting *.txt or *.root* files"

    print "INFO :: Reading in %d input file(s)" % len(args.input_files)

    # Check a few files to make sure they have prefixes
    import random
    n_files_to_check = min(10, len(args.input_files))
    file_indices = range(len(args.input_files))
    for idx in random.sample(file_indices, n_files_to_check):
        with open(args.input_files[idx], 'r') as f:
            first_file = f.readline().strip()
            if not file_name_has_xrootd_prefix(first_file):
                print "ERROR :: A file was found without proper xrootd",
                print "storage prefixes: %s" % args.input_files[idx] 
                sys.exit()

    # Check that directory for storing tar file exists
    if '/' in args.tar_file: # check only if tar_file is a path
        tar_file_dir = os.path.dirname(args.tar_file)
        if not os.path.exists(tar_file_dir):
            print 'ERROR :: Cannot locate directory to place tar file',
            print '(attempted dir = %s)' % tar_file_dir
            sys.exit()

    # Check if tar file exists.
    # If so, check with user if they want to use it or overwite it.
    # If the user wants to overwrite an old file or no old file exists, then
    # create the new tar file.
    if os.path.exists(args.tar_file):
        usr_msg =  "Tar file already exists: %s\n" % args.tar_file
        usr_msg += "Would you like to [U]se or [O]verwrite it? [U/O] "
        user_op = raw_input(usr_msg)

        # Only accept U or O
        while user_op not in ["U","O"]:
            usr_msg = "Unacceptable answer: %s\n" % user_op
            usr_msg += "Would you like to [U]se or [O]verwrite it? [U/O] "
            user_op = raw_input(usr_msg)

        if user_op == "U":
            make_tar_file = False
        elif user_op == "O":
            make_tar_file = True
    else:
        make_tar_file = True

    # Check that sumw file exists if requested
    if args.sumw and not os.path.exists(args.sumw):
        if args.sumw.startswith("$"):
            args.sumw = ''
        else:
            print "ERROR :: Sum of weights file does not exist:", args.sumw

    return make_filelists, make_tar_file

################################################################################

def get_args():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files',
                        nargs="*",
                        help=_help_input_files)
    parser.add_argument('-t', '--tar-dir',
                        help = _help_tar_dir,
                        default = _df_tar_dir)
    parser.add_argument('--tar-file',
                        help = _help_tar_file,
                        default = _df_tar_file)
    parser.add_argument('-e', '--executable',
                        help = _help_exec,
                        default = _df_exec)
    parser.add_argument('-w', '--sumw',
                        help = _help_sumw,
                        default = _df_sumw)
    parser.add_argument('-o', '--output-dir',
                        help = _help_output_dir,
                        default = _df_output_dir)
    parser.add_argument('--split-dsids',
                        nargs="*",
                        help = _help_split_dsids,
                        default = _df_split_dsids)
    parser.add_argument('--syst',
                        action='store_true',
                        help=_help_syst)
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help=_help_verbose)
    args = parser.parse_args()

    # Change all paths to be absolute
    args.input_files = [os.path.abspath(x) for x in args.input_files]
    args.tar_dir = os.path.abspath(args.tar_dir)
    args.tar_file = os.path.abspath(args.tar_file)
    args.sumw = os.path.abspath(args.sumw) if not args.sumw.startswith("$") else ""
    args.output_dir = os.path.abspath(args.output_dir)

    return args

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        args = get_args()
        if args.verbose:
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()
        main()
        if args.verbose:
            print time.asctime()
            time = (time.time() - start_time)
            print 'TOTAL TIME: %fs'%time,
            print ''
            print '<'*40
    except KeyboardInterrupt, e: # Ctrl-C
        print 'Program ended by keyboard interruption'
        raise e
    except SystemExit, e: # sys.exit()
        print 'Program ended by system exit'
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)


