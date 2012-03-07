#! /usr/bin/env python
# vim: set ft=python :
"""
Calculate density of reads in bins across the genome from bam file.
Output units:  RPKM
Output format: Bedgraph

* Input sort order does not matter
* Output goes to stdout
* Currently ignores chrM and gapped or local alignemts (where the
  aligned length is not the same as the read length).

TODO: add strand-specific binning
TODO: handle gapped/local alignments?
TODO: wiggle format output
TODO: variable size binning
"""

import argparse
import logging
import sys
import numpy
import pysam

from gosr.common import arghelpers


def make_bins(chrominfo, binsize):
    """create a dictionary with one array of bins per chromosome"""
    result = {}
    for name, l in chrominfo.items():
        n_bins = l // binsize  #reads in last bin are discarded
        result[name] = numpy.zeros(n_bins, dtype = numpy.int32)
    return result    

def binbam(bamfile, binsize, fragsize, chrominfo):
    """count aligned reads per bin in bamfile"""
    logging.debug("Making bin data structure")
    bins = make_bins(chrominfo, binsize)
    logging.debug("    ....done")
    n_aln   = 0
    n_unaln = 0
    n_igno  = 0
    shift   = fragsize / 2
    for aln in bamfile:
        if not aln.is_unmapped:
            chrom = bamfile.getrname(aln.tid)
            if chrom == "chrM":
                n_igno += 1
                continue
            if aln.alen != aln.rlen:
                logging.debug("Not a end-to-end alignment, or alignment gapped")
                logging.debug(str(aln))
                n_igno += 1
                continue
            if not aln.is_reverse:
                bin_nr = (aln.pos + shift) // binsize
            else:
                bin_nr = (aln.aend - 1 - shift) // binsize
            n_aln += 1
            try:
                bins[chrom][bin_nr] += 1
            except IndexError:
                logging.debug("BIN OUT OF RANGE: %s: pos[%d] -> bin[%d]", chrom, aln.pos, bin_nr)
        else:
            n_unaln += 1
    # normalization factor
    rpkm_factor = (1e6 / n_aln) * (1000.0 / binsize) 
    logging.info("Aligned reads:   %8d", n_aln)
    logging.info("Unaligned reads: %8d", n_unaln)
    logging.info("Ignored reads:   %8d", n_igno)
    return bins, rpkm_factor

def output_bedgraph(bins, binsize, norm_factor, trackline = ""):
    """write all non-empty bins to bedgraph format strings"""
    bg = "{0}\t{1}\t{2}\t{3:.8f}"
    if trackline != "":
        print "track type=bedGraph alwaysZero=on visibility=full maxHeightPixels=100:80:50 " \
                + trackline
    for chrom in sorted(bins.keys()):
        for i, n in enumerate(bins[chrom]):
            if n > 0:
                print bg.format(chrom, i * binsize, (i + 1) * binsize,
                        norm_factor * n)

################################################################################
# tool interface
################################################################################
def process(args):
    """pipeline driver"""

    bamfile   = pysam.Samfile(args.infile, "rb")
    try:
        chrominfo = dict(zip(bamfile.references, bamfile.lengths))
    except:
        logging.error("Could not extract length information for references from bam file")
        bamfile.close()
        sys.exit(1)
    logging.info("Start binning process")
    try:
        bins, norm_factor = binbam(bamfile, args.binsize, args.frag_size, chrominfo)
        logging.info("DONE")
        output_bedgraph(bins, args.binsize, norm_factor, args.track_line)
    finally:
        bamfile.close()

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("binbam",
            help = "Create a density track of reads in bins from bam file",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("infile", type = arghelpers.infilename_check, 
            help = "Bam file; use '-' for stdin")
    cmdline.add_argument("binsize", type = int,
            help = "Size of bins to use")
    cmdline.add_argument("-s", "--frag-size", type = int, 
            default = 0,
            help = "Size of fragments; reads are shifted by half fragsize [%(default)d]")
    cmdline.add_argument("-t", "--track-line", default = "",
            help = """include trackline 'track type=bedGraph
            alwaysZero=on visibility=full maxHeightPixels=100:80:50 TRACK_LINE' [no trackline]""")
    cmdline.set_defaults(func = process) 
