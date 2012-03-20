"""
sort file by chromosome, position, with the chromosome
sort order given by the genome module; sorting is done by assuming that column
1 is chromosome and column 2 is a position.  Strand is ignored in sorting.

Output is to stdout
"""

import sys
import logging
import argparse
import shlex
import subprocess
import signal
from gosr.common import arghelpers
from gosr.common import genome
from gosr.common.file import FileOrGzip

def start_unix_sort(mem):
    """start up an external sort process for bed file; returns 2 popen
    objects (the sort and the cut); This is a coroutine! prime and close!"""
    cmdline = shlex.split("sort -S%s -k1,1g -" % mem)
    logging.info("sort call: %s", cmdline)
    sortproc = subprocess.Popen(cmdline, stdin = subprocess.PIPE, 
            stdout = subprocess.PIPE, shell = False)
    cutproc  = subprocess.Popen(shlex.split("cut -f2-"), 
            stdin = sortproc.stdout, shell = False, 
            preexec_fn = sortproc.stdin.close)
    def _exit_nicely(signr, frame):
        logging.warn("Received signal %d; terminating subprocesses and exiting",
                signr)
        sortproc.terminate()
        cutproc.terminate()
        sys.exit(1)
    signal.signal(signal.SIGINT, _exit_nicely)
    try:
        while True:
            data = (yield)
            sortproc.stdin.write(data)
    except GeneratorExit:
        sortproc.stdin.close()
        returncode1 = sortproc.wait()
        returncode2 = cutproc.wait()
        if returncode1 != 0 or returncode2 != 0:
            logging.error("subprocesses exited abnormally (sort: %d, cut: %d)",
                    returncode1, returncode1)
            sys.exit()

def sort_bed(args):
    """sort bed-like file by chrom and start pos"""
    try:
        chroms = getattr(genome, args.genome)
    except AttributeError:
        logging.error("Genome %s not available", args.genome)
        sys.exit(1)
   
    sortproc = start_unix_sort(args.S)
    sortproc.next()

    outlst = []
    n      = 0
    out    = "{0}\t{1}"
    logging.info("Start feeding sort")
    with FileOrGzip(args.infile) as fh:
        for line in fh:
            chrom, pos, _ = line.split("\t", 2)
            outlst.append(out.format(chroms.cpos(chrom, long(pos)), line))
            n += 1
            if n == 100000:
                sortproc.send("".join(outlst))
                outlst = []
                n      = 0
    sortproc.send("".join(outlst))
    sortproc.close()
    logging.info("Done feeding sort; Waiting for sort to finish")

#===============================================================================
# interface
#===============================================================================

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("bed-sort",
            help = """Sort bed file""",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("infile", type = arghelpers.infilename_check,
            help = "Input bed file to be sorted; can be gzip'ed")
    cmdline.add_argument("genome",
            choices = ["mm9", "hg19", "hg18"],
            help = "Input bed file to be sorted")
    cmdline.add_argument("-S", default = "1G",
            help = """memory size; passed on to external unix sort; see 'man
            sort'; [%(default)s]""")
    cmdline.set_defaults(func = sort_bed)

