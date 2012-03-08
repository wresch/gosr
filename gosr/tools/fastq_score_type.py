#! /usr/bin/env python
"""
Determine the quality score type of a fastq file based on the first few
thousand sequences

if fastq file is '-', reads from stdin. Files ending in .gz are 
accepted and uncompressed on the fly.

On stdout returns a single string indicating score type; stderr displays log
"""

import sys
import logging
import collections
import argparse
from gosr.common import fastq
from gosr.common import arghelpers
from gosr.common.file import FileOrGzip

def guess_score_type(letter_freq_dict):
    """returns solexa|phred64|phred33"""
    total = sum(letter_freq_dict.values())
    logging.info("Processed %i quality letters", total)
    letters = [ord(x) for x in letter_freq_dict.keys()]
    maxval = max(letters)
    minval = min(letters)
    logging.info("Score range: %3d - %3d", minval, maxval)
    assert maxval <= 104
    if minval < 59:
        assert maxval <= 73
        logging.info("Score type: phred33")
        return "phred33"
    else:
        if minval < 64:
            logging.info("Score type: solexa")
            return "solexa"
        else:
            logging.info("Score type: phred64")
            return "phred64"

def determine_score_type(args):
    with FileOrGzip(args.fastq) as infile:
        count = 0
        letter_freq = collections.defaultdict(int)
        for rid, seq, qual in fastq.read(infile):
            count += 1
            for l in qual:
                letter_freq[l] += 1
            if count == 5000:
                break
        print guess_score_type(letter_freq)

#===============================================================================
# interface
#===============================================================================

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("fastq-score-type",
            help = """Determine score type for fastq file""",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("fastq", type = arghelpers.infilename_check,
            help = "Fastq file")
    cmdline.set_defaults(func = determine_score_type)

