#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
reads fastq from imput and converts quality scores
from old solexa scores or phred 64 to phred33.
For solexa, this uses the formula

Qp = Qs + 10*log(1+10^(Qs/-10), 10)

For efficiency, this program uses a character table for
string.maketrans that has been previously calculated
to translate old solexa to new solexa character strings

phred64 is a simple left shift

phred33+ makes sure that q values don't exceed 40 (as
it does in recent illumina output; currently works
for input up to Q45).

Output goes to stdout
"""


import sys
import argparse
from string import maketrans

from gosr.common.file import FileOrGzip
from gosr.common import fastq
from gosr.common import arghelpers



solexa   = r""";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
phred33s = r"""""##$$%%&&'()*++,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_"""
solexa_to_phred33 = maketrans(solexa, phred33s)

phred64  = r"""@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
phred33  = r"""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_"""
phred64_to_phred33 = maketrans(phred64, phred33)

phred33plus_to_phred33 = maketrans("JKLMN", "IIIII")

def to_phred33(args):
    if args.score_type == "phred64":
        table = phred64_to_phred33
    elif args.score_type == "solexa":
        table = solexa_to_phred33
    elif args.score_type == "phred33+":
        table = phred33plus_to_phred33
    out = "@{0}\n{1}\n+\n{2}"
    with FileOrGzip(args.fastq) as infile:
        for rid, rseq, rqual in fastq.read(infile):
            print out.format(rid, rseq, rqual.translate(table))

#===============================================================================
# interface
#===============================================================================

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("fastq-to-phred33",
            help = """changes from solexa of phred 64 to phred33""",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("score_type", choices = ["phred64", "solexa", "phred33+"],
            help = "Current score type")
    cmdline.add_argument("fastq", type = arghelpers.infilename_check,
            help = "Fastq file (can be .gz); '-' reads from stdin")
    cmdline.set_defaults(func = to_phred33)
