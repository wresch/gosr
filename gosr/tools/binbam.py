#! /usr/bin/env python
# vim: set ft=python :
"""
Calculate density of reads in bins across the genome from bam file.
Output units:  RPKM
Output format: Bedgraph

* Input sort order does matter
* Output goes to stdout
* Currently ignores chrM and gapped or local alignemts (where the
  aligned length is not the same as the read length).

TODO: handle gapped/local alignments?
TODO: variable size binning
"""

import argparse
import logging
import sys
import itertools
import numpy
import pysam

from gosr.common import arghelpers
from gosr.common import dsp


def make_bins(chrominfo, binsize, by_strand):
    """create a dictionary with one array of bins per chromosome"""
    result = {}
    for name, l in chrominfo.items():
        n_bins = l // binsize  #reads in last bin are discarded
        if not by_strand:
            result[name] = numpy.zeros(n_bins, dtype = numpy.int32)
        else:
            result[name] = [numpy.zeros(n_bins, dtype = numpy.int32),
                            numpy.zeros(n_bins, dtype = numpy.int32)]
    return result

def binbam(bamfile, binsize, fragsize, chrominfo, n_redundancy, by_strand):
    """count aligned reads per bin in bamfile; *bamfile needs to be sorted*"""
    bins = make_bins(chrominfo, binsize, by_strand)
    n_aln   = 0
    n_rmred = 0
    n_igno  = 0
    shift   = fragsize / 2
    for _, alns in itertools.groupby(bamfile, lambda x: (x.tid, x.pos)):
        # split up into plus and minus strand
        all_alns = [x for x in alns if not x.is_unmapped]
        n_aln   += len(all_alns)
        plus     = [x for x in all_alns if not x.is_reverse][0:n_redundancy]
        minus    = [x for x in all_alns if x.is_reverse][0:n_redundancy]
        for aln in itertools.chain(plus, minus):
            chrom = bamfile.getrname(aln.tid)
            if chrom == "chrM":
                n_igno += 1
                continue
            n_rmred += 1
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
                if not by_strand:
                    bins[chrom][bin_nr] += 1
                elif not aln.is_reverse:
                    bins[chrom][0][bin_nr] += 1
                else:
                    bins[chrom][1][bin_nr] += 1
            except IndexError:
                logging.debug("BIN OUT OF RANGE: %s: pos[%d] -> bin[%d]", chrom, aln.pos, bin_nr)
    # normalization factor
    rpkm_factor = (1e6 / n_rmred) * (1000.0 / binsize)
    logging.info("Aligned reads:                %8d", n_aln)
    logging.info(" after removing redundancy:   %8d", n_rmred)
    logging.info(" normalization factor:        %f", rpkm_factor)
    logging.info("Ignored reads:                %8d", n_igno)
    return bins, rpkm_factor

def output_wiggle(bins, binsize, norm_factor, by_strand, name, extra_trackline = ""):
    """write all non-empty bins to bedgraph format strings; always includes
    minimal track line; Output is in 1-based wiggle format."""
    if not by_strand:
        print "track type=wiggle_0 alwaysZero=on visibility=full maxHeightPixels=100:80:50 " \
                + ("name='%s'" % name) + extra_trackline
        for chrom in sorted(bins.keys()):
            print "variableStep chrom=%s span=%d" % (chrom, binsize)
            non_zero_bins = numpy.nonzero(bins[chrom] > 0)
            result = numpy.column_stack((non_zero_bins[0] * binsize + 1,
                bins[chrom][non_zero_bins] * norm_factor))
            numpy.savetxt(sys.stdout, result, "%d\t%.8f")
    else:
        for strand in (0, 1):
            if strand == 0:
                nf = norm_factor
            else:
                nf = -norm_factor
            print "track type=wiggle_0 alwaysZero=on visibility=full maxHeightPixels=100:80:50 " \
                    + ("name='%s[%s]'" % (name, strand and '-' or '+')) + extra_trackline
            for chrom in sorted(bins.keys()):
                print "variableStep chrom=%s span=%d" % (chrom, binsize)
                non_zero_bins = numpy.nonzero(bins[chrom][strand] > 0)
                result = numpy.column_stack((non_zero_bins[0] * binsize + 1,
                    bins[chrom][strand][non_zero_bins] * nf))
                numpy.savetxt(sys.stdout, result, "%d\t%.8f")

def smooth(bins, window_size, by_strand):
    for chrom in bins:
        if not by_strand:
            bins[chrom] = dsp.savitzky_golay_filter(bins[chrom], window_size, 
                    order = 2, deriv = 0)
        else:
            bins[chrom][0] = dsp.savitzky_golay_filter(bins[chrom][0], window_size, 
                    order = 2, deriv = 0)
            bins[chrom][1] = dsp.savitzky_golay_filter(bins[chrom][1], window_size, 
                    order = 2, deriv = 0)

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
    logging.info(" allowing up to %d redundant reads", args.n_redundancy)
    if args.by_strand:
        logging.info("Reporting separate densities for each strand")
    logging.info("Track name: %s", args.name)
    logging.info("Track line extra options: \"%s\"", args.track_line)
    if args.sg > 0:
        logging.info("Smoothing output with savitzky-golay filter, order 2, width %d bins", args.sg)
    try:
        bins, norm_factor = binbam(bamfile, args.binsize, args.frag_size,
                chrominfo, args.n_redundancy, args.by_strand)
    finally:
        bamfile.close()
    
    logging.info("DONE")
    if args.sg > 0:
        smooth(bins, args.sg, args.by_strand)
    output_wiggle(bins, args.binsize, norm_factor, args.by_strand, args.name, args.track_line)

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
    cmdline.add_argument("name",
            help = "Name of track")
    cmdline.add_argument("-f", "--frag-size", type = int,
            default = 0,
            help = "Size of fragments; reads are shifted by half fragsize [%(default)d]")
    cmdline.add_argument("-n", "--n-redundancy", type = int,
            default = 3,
            help = "Number of identical alignments per position allowed [%(default)d]")
    cmdline.add_argument("--sg", type = int,
            default = 0,
            help = """If greater than 0, a Savitzky Golay filter of order 2 with
            the given window size (in terms of bins) is applied before
            outputting [%(default)d]""")
    cmdline.add_argument("-s", "--by-strand", default = False,
            action = "store_true", 
            help = "Output separate densities by strand [%(default)s]")
    cmdline.add_argument("-t", "--track-line", default = "",
            help = """include extra options in track line. 'track type=bedGraph
            alwaysZero=on visibility=full maxHeightPixels=100:80:50' is always
            included""")
    cmdline.set_defaults(func = process)
