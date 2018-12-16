#!/usr/bin/env python
"""
================================================================================
Add a prefix to each line of input files
Examples
    python add_prefix_to_filelists.py path/to/files/*txt -p mwt2 -d output/dir/
    python add_prefix_to_filelists.py --show-prefixes

Author:
    Alex Armstrong <alarmstr@cern.ch>
    December 2018
================================================================================
"""

import sys, os, traceback, argparse
import time
import subprocess
import pdb

################################################################################
# GLOBALS
################################################################################
# FAX storage prefixes
PREFIXES = {
    # Note: '/' will be added at the end of the prefix when appending file name
    #'mwt2'         : 'root://fax.mwt2.org:1094/',
    'mwt2_scratch' : 'root://fax.mwt2.org:1094//pnfs/uchicago.edu/atlasscratchdisk/rucio',
    'mwt2_local'   : 'root://fax.mwt2.org:1094//pnfs/uchicago.edu/atlaslocalgroupdisk/rucio',
    #'slac'         : 'root://griddev03.slac.stanford.edu:2094/'
    'slac_scratch' : 'root://griddev03.slac.stanford.edu:2094//xrootd/atlas/atlasscratchdisk/rucio',
    'slac_local'   : 'root://griddev03.slac.stanford.edu:2094//xrootd/atlas/atlaslocalgroupdisk/rucio',
    'atlas-xrd'    : 'root://atlas-xrd-us.usatlas.org/',
    'uk_scratch'   : 'root://xrootd.esc.qmul.ac.uk:1094//atlas/atlasscratchdisk/rucio',
    'env'          : os.environ.get('STORAGEPREFIX','')
}

# User Argument defaults and help information
_help_ifile_names = 'Input file names'

_help_prefix = 'Prefix to prepend to all lines of each input file \
                       [options: %s]' % (', '.join(PREFIXES.keys()))

_df_save_dir = os.getenv("FILELIST_DIR", "$FILELIST_DIR")
_help_save_dir = 'Directory for storing new files with prefixes added \
                  [default: %s]' % _df_save_dir

_help_no_new_dirs = "Store all new files in the save directory without \
                     adding new sub-directories when input files have different \
                     paths. New directories, when added, only include differences \
                     in the absolute paths of input files. This helps prevent files \
                     from being overwritten"
_help_show_prefixes = 'Show prefix options with their storage prefixes'
_help_verbose = 'verbose output'


################################################################################
def main ():
    """ Main Function """

    global args
    check_environment()
    check_inputs(args)

    prefix =  get_storage_prefix(args.prefix)
    print "INFO :: Adding prefix (%s) to %d file(s) and saving at %s" % (
            prefix, len(args.input_files), args.save_dir)

    add_prefix_to_filelists(args.input_files, prefix, args.save_dir)
    
    print "INFO :: Done! New files are saved in", args.save_dir

def add_prefix_to_filelists(input_files, prefix, save_dir, new_directories=True):
    """
    Add a prefix to each line of input files and save in output directory. Some
    lines will get skipped (e.g. blank lines). See skip_txt_line.

    New directories will be created in the output directory if the input files
    do not all share the same parent directory. In this case, the common prefix
    of each path is removed and the remaining paths are used to make the new
    directories. Setting new_directories off will cause all files to be put in
    the save directory which can lead to overwriting if files with different
    paths have the same name

    args:
        input_files (list(str)) - list of files to be processed
        prefix (str) - prefix to be added
        save_dir (str) - directory for storing updated files
        new_directories (bool) - create new directories to prevent overwriting
    
    returns:
        (list(str)) - list of paths to new files 
    """
    # Add needed directory structure if requested
    if new_directories:
        dirs = {os.path.dirname(os.path.abspath(f)) for f in input_files}
        common_dir = os.path.commonprefix(list(dirs))
        while not os.path.exists(common_dir):
            # Common prefix does not always return the common directory
            # as it only returns common characters
            common_dir = os.path.dirname(common_dir)
        dirs_to_make = [os.path.relpath(d, common_dir) for d in dirs]

        for d in dirs_to_make :
            if d == '.': continue # No need to remake current directory
            new_dir = os.path.join(save_dir, d)
            if not os.path.isdir(new_dir):
                cmd = 'mkdir -p %s' % new_dir
                print "INFO :: Making new directory:", new_dir
                subprocess.call(cmd, shell = True)

    # Make new files with prefixes
    overwite_all = False
    first_overwrite_flag = True
    new_files = []
    for f in input_files :
        # Determine new pile path
        if new_directories:
            rel_file_path = os.path.relpath(os.path.abspath(f), common_dir)
        else:
            rel_file_path = f.split('/')[-1]
        new_file = os.path.join(save_dir, rel_file_path)
        new_files.append(new_file)

        # Check if new file already exists
        # If so, check with user if they want to overwrite it or skip
        if os.path.exists(new_file) and not overwite_all:
            if first_overwrite_flag:
                print "INFO :: Attempting to write new files with prefix",
                print "but files were already found to exist:"
                first_overwrite_flag = False
            usr_msg = "File already exists: %s\n" % new_file
            usr_msg += "Would you like to [O]verwrite it, [S]kip and save, "
            usr_msg += "overwrite [A]ll conflicts, or [E]xit program? [O,S,A,E] "
            overwrite_op = raw_input(usr_msg)

            # Only accept O, S, A, or E
            while overwrite_op not in ["O","S","A","E"]:
                usr_msg = "Unacceptable answer: %s\n" % overwrite_op
                usr_msg += "Would you like to [O]verwrite it, [S]kip and save, "
                usr_msg += "overwrite [A]ll conflicts, or [E]xit program? [O,S,A,E] "
                overwrite_op = raw_input(usr_msg)

            if overwrite_op == "O":
                pass
            elif overwrite_op == "S":
                continue
            elif overwrite_op == "A":
                overwite_all = True
            elif overwrite_op == "E":
                sys.exit()

        # Make text for new file
        new_file_str = ''
        with open(f, 'r') as in_file :
            for line in in_file :
                line = line.strip()
                if skip_txt_line(line): continue
                new_line = '%s%s' % (prefix, line)
                new_file_str += new_line + '\n'

        # Write to the new file
        with open(new_file, 'w') as out_file :
            out_file.write(new_file_str)

    # Return paths to new files
    return new_files

################################################################################
# FORMAT SENSATIVE FUNCTIONS
################################################################################
def skip_txt_line(line):
    ''' Check if a line read in from a text file should not be processed'''
    l = line.strip()
    return l.startswith("#") or not l

def acceptable_xrootd_storage_prefix(prefix):
    return prefix.startswith("root://")

def get_storage_prefix(usr_prefix):
    ''' Return full storage prefix '''
    # Check for predefined prefixes or at least an acceptable format of the
    # user provided prefix
    if args.prefix in PREFIXES:
        prefix = PREFIXES[args.prefix]
    elif acceptable_xrootd_storage_prefix(args.prefix):
        prefix = args.prefix
    else:
        print "WARNING :: Unexpected storage prefix: ", args.prefix
        print "INFO :: Either it is not a known prefix shortcut option,",
        print "or the prefix is of an unexpected format."
        print "INFO :: Try running with --show-prefixes"
        prefix = args.prefix

    return prefix

################################################################################
# SUPPORT FUNCTIONS
################################################################################
def check_inputs(args):
    """ Check the input arguments are as expected """
    # Show prefixes, if requested, before doing any other checks
    if args.show_prefixes :
        width = len(max(PREFIXES.keys(), key=len))
        for key, prefix in PREFIXES.items() :
            print '\t%*s : %s' % (width, key, prefix)
        sys.exit()

    # Check that input files were provided
    if not args.input_files:
        print "ERROR :: No input files were provied"

    # Check that all input files exist
    for f in args.input_files:
        if not os.path.exists(f):
            print "ERROR :: Cannot find input file:", f
            sys.exit()
        if not f.endswith('txt'):
            print "ERROR :: Unpexected input file format: ", f
            print "INFO :: Expecting *.txt files"

    # Check the save directory exists
    if not os.path.exists(args.save_dir):
        print "ERROR :: Save directory not found:", args.save_dir

def check_environment():
    """ Check if the shell environment is setup as expected """
    assert os.environ['USER'], "USER variable not set"

    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    assert python_ver >= 2.7, ("Running old version of python\n", sys.version)

################################################################################

def get_args():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files',
                        nargs='*',
                        help=_help_ifile_names)
    parser.add_argument('-p', '--prefix',
                        help = _help_prefix)
    parser.add_argument('-d', '--save-dir',
                        help = _help_save_dir,
                        default = _df_save_dir)
    parser.add_argument('--no-new-dirs',
                        help=_help_no_new_dirs)
    parser.add_argument('--show-prefixes',
                        action='store_true',
                        help = _help_show_prefixes)
    parser.add_argument('-v', '--verbose',
                        help=_help_verbose)

    args = parser.parse_args()
    return args

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        # TODO: Add ability to check standard input so things can be piped
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

