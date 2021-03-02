#!/usr/bin/env python3
from scorpio import __version__

import setuptools
import argparse
import os.path
import snakemake
import sys
import tempfile
import csv
import os
import yaml
from datetime import datetime
from datetime import date
from Bio import SeqIO
import pkg_resources
from . import _program

from reportfunk.funks import report_functions as rfunk
from reportfunk.funks import custom_logger as custom_logger
from reportfunk.funks import log_handler_handle as lh
import scorpiofunks as sfunk

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    usage='''scorpio -i <config.yaml> [options]''')

    io_group = parser.add_argument_group('input output options')
    io_group.add_argument('-i',"--config", action="store",help="Input config file", dest="config")
    io_group.add_argument('-j',"--json", action="store",help="Input json file", dest="json")
    io_group.add_argument('--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('-o','--output-prefix', action="store",help="Output prefix. Default: scorpio",dest="output_prefix")

    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument('--tempdir',action="store",help="Specify where you want the temporary stuff to go Default: $TMPDIR")
    misc_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument("--no-force",action="store_true",help="Dont force run rules",dest="no_force")
    misc_group.add_argument('-t','--threads', action='store',dest="threads",type=int,help="Number of threads")
    misc_group.add_argument("-v","--version", action='version', version=f"scorpio {__version__}")
    
    """
    Exit with help menu if no args supplied
    """
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    """
    Initialising dicts
    """

    # get the default values from scorpiofunks
    config = sfunk.get_defaults()

    """
    Valid inputs are config.yaml/config.yml
    
    """
    # find config file
    if args.config:
        sfunk.add_arg_to_config("config",args.config,config)

    # if a yaml file is detected, add everything in it to the config dict
    if config["config"]:
        sfunk.parse_yaml_file(config["config"], config)
    
    """
    Output directory 
    """
    # default output dir
    sfunk.get_outdir(args.outdir,args.output_prefix,cwd,config)
    figdir  = os.path.join(config["outdir"], "figures")
    if not os.path.exists(figdir):
        os.mkdir(figdir)
    config['figdir'] = figdir

    """
    Get tempdir 
    """
    
    # specifying temp directory, outdir if no_temp (tempdir becomes working dir)
    tempdir = sfunk.get_temp_dir(args.tempdir, args.no_temp,cwd,config)
    config["tempdir"] = tempdir
    """
    Parsing the report_group arguments, 
    config options
    """

    """
    Miscellaneous options parsing

    """
    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["log_string"] = ""
    else:
        quiet_mode = True
        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "

    if args.no_force:
        force_option = False
    else:
        force_option = True
    
    config["forceall"] = force_option

    sfunk.add_arg_to_config("threads",args.threads,config)
    
    try:
        config["threads"]= int(config["threads"])
    except:
        sys.stderr.write(sfunk.cyan('Error: Please specifiy an integer for variable `threads`.\n'))
        sys.exit(-1)
    threads = config["threads"]

    if args.json:
        config["json"] = os.path.join(cwd,args.json)
    else:
        fn = config["filename"]
        fn_unzipped = ".".join(fn.split(".")[:-1])
        if os.path.exists(fn):
            os.system(f"rm {fn}")
        if os.path.exists(fn_unzipped):
            os.system(f"rm {fn_unzipped}")
        os.system(f"wget --password {config['password']} "
                    f"--user {config['username']} {config['url']} "
                    "&& "
                    f"bzip2 -d {fn} ")
        config["json"] = os.path.join(cwd,fn_unzipped)

    config["timestamp"] = sfunk.get_timestamp()
    
    # find the master Snakefile
    snakefile = sfunk.get_snakefile(thisdir)

    if args.verbose:
        
        for k in sorted(config):
            print(sfunk.green(k), config[k])
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=force_option, force_incomplete=True,
                                    workdir=tempdir,config=config, cores=threads,lock=False
                                    )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=force_option,force_incomplete=True,workdir=tempdir,
                                    config=config, cores=threads,lock=False,quiet=True,log_handler=logger.log_handler
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()