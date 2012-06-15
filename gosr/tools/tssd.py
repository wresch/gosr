# vim: set ft=python :
"""
Given reads from bam file, calculate read density around TSSs
listed in GTF file. Bam file does not have to be ordered or indexed.

Note: * For overlapping intervals, one is chosen at random

Also estimates fragment size by calculating the density separately
for plus and minus strand and finding shift that leads to optimal overlap
between the two densities.

Note that this tool calculates read counts, not coverage.
"""

import logging
import argparse
import sys
import numpy
import HTSeq

from gosr.common import arghelpers
from gosr.common import dsp

def overlaps_any(garray, iv):
    steps = list(garray[iv].steps())
    if len(steps) > 1 or steps[0][1] is not None:
        return True
    else:
        return False

def gtf_to_tsspos(gtf, upstream, downstream):
    """extract TSSs from gtf (using 'exon_number' attribute) and
    create a GenomicArray of non-overlapping TSSs plus the upstream
    and downstream region.  Limit is extended by 200 nts on each end
    to allow for shifting to determine fragment size estimate later on"""
    tsspos  = HTSeq.GenomicArray("auto", typecode = 'O', stranded=False)
    n_feat  = 0
    n_tss   = 0
    n_used  = 0
    for feature in gtf:
        n_feat += 1
        if feature.type == "exon" and feature.attr["exon_number"] == "1":
            n_tss += 1
            p      = feature.iv.start_d_as_pos
            if feature.iv.strand == "+":
                window = HTSeq.GenomicInterval(p.chrom, p.pos - upstream,
                        p.pos + downstream + 1, "+")
            elif feature.iv.strand == "-":
                window = HTSeq.GenomicInterval(p.chrom, p.pos - downstream,
                        p.pos + upstream + 1, "-")
            else:
                logging.error("bad strand found in GTF file: %s [line %d]",
                        feature.iv.strand, n_feat)
                sys.exit(1)
            if not overlaps_any(tsspos, window):
                tsspos[window] = window
                n_used += 1
    logging.info("found %d TSSs", n_tss)
    logging.info(" of which %d were used (i.e. non-overlapping)", n_used)
    return tsspos, n_used

def make_density(tsspos, bamfile, up, down):
    """calculate tag density in RPKM around TSSs in tsspos GenomicArray;
    separate densities by whether they represent the left or right side of
    a fragment.  note that this depends on the strand of the feature:
        * for a plus strand feature, a plus strand read is left, a minus
          strand read is right
        * for a minus strand feature, a plus strand read is right, a minus
          strand read is left
    """
    n_reads        = 0
    n_reads_on_tss = 0
    locd           = {'+': {'+': "left", "-": "right"},
                      '-': {"+": "right", "-": "left"}}
    d              = {
            'left'  : numpy.zeros(up + down + 1, dtype = "i"),
            'right' : numpy.zeros(up + down + 1, dtype = "i") }
    for aln in bamfile:
        if aln.aligned:
            n_reads += 1
            alniv = aln.iv
            tss = tsspos[alniv.start_d_as_pos]
            if tss is not None:
                n_reads_on_tss += 1
                pos_in_window = numpy.abs(alniv.start_d - tss.start_d)
                loc = locd[tss.strand][alniv.strand]
                try:
                    d[loc][pos_in_window] += 1
                except IndexError:
                    logging.error("pos_in_window out of bounds: %d",
                            pos_in_window)
                    logging.error("tss: %s", str(tss))
                    logging.error("aln: %s", str(alniv))
                    sys.exit(1)
    logging.info("Reads processed: %9d", n_reads)
    logging.info("Reads on tss:    %9d", n_reads_on_tss)
    return d, n_reads

def dist(x, y):
    assert(len(x) == len(y))
    return numpy.sqrt(numpy.sum(numpy.power(x - y, 2)))

def determine_frag_size(density, extra):
    """take 2 densities (left and right) and shift the right peak to the left
    until the two peaks reach maximal similarity;
    TODO: should data be smoothed before doing this?"""
    x = density["left"]
    y = density["right"]
    peak_dist = []
    shift     = []
    n         = len(x) - 2 * extra - 1

    for i in range(extra + 1):
        b1 = extra - i + 1
        b2 = extra + i + 1
        peak_dist.append(dist(x[b1:(b1 + n)], y[b2:(b2 + n)]))
        shift.append(i)
    optimal_shift = [(s, d) for s, d in zip(shift, peak_dist) if d == min(peak_dist)]
    if len(optimal_shift) > 1:
        logging.warn("more than one possible shift size: %s", optimal_shift)
    logging.info("Inferred fragment size estimate: %d", optimal_shift[0][0] * 2)
    return optimal_shift[0][0] * 2


def output(density, up, down, extra, frag_size, n_tss, n_reads):
    smooth_filter_size = 101
    shift = int(frag_size / 2)
    pos   = range(-up, down + 1)
    n     = len(density["left"]) - 2 * extra - 1
    
    left = (density["left"].astype(float) / n_reads) * 1e9 / n_tss
    left_smooth = dsp.savitzky_golay_filter(left, smooth_filter_size,
            order = 4)
    print "\n".join("{0}|{1}|{2}|left".format(a, b, c) for a, b, c in
            zip(pos[extra:(extra + n)], left[extra:(extra + n)], 
                left_smooth[extra:(extra + n)]))
    
    right = (density["right"].astype(float) / n_reads) * 1e9 / n_tss
    right_smooth = dsp.savitzky_golay_filter(right, smooth_filter_size,
            order = 4)
    print "\n".join("{0}|{1}|{2}|right".format(a, b, c) for a, b, c in
            zip(pos[extra:(extra + n)], right[extra:(extra + n)], 
                right_smooth[extra:(extra + n)]))
    
    combined = numpy.zeros(len(left), dtype = float)
    cs = extra + 1
    left_start = extra - frag_size // 2 + 1
    combined[cs:(cs + n)] += left[left_start:(left_start + n)]
    right_start = extra + frag_size // 2 + 1
    combined[cs:(cs + n)] += right[right_start:(right_start + n)]
    combined_smooth = dsp.savitzky_golay_filter(combined, smooth_filter_size,
            order = 4)
    print "\n".join("{0}|{1}|{2}|combined".format(a, b, c) for a, b, c in
            zip(pos[extra:(extra + n)], combined[extra:(extra + n)], 
                combined_smooth[extra:(extra + n)]))


################################################################################
# tool interface
################################################################################
def process(args):
    """pipeline driver"""
    # additional nts added to allow shifting peaks on top of each other
    # only used for density, not for finding matching reads
    extra      = 200
    up         = args.upstream
    down       = args.downstream
    logging.info("Window: <-- %d --TSS-- %d -->", up, down)

    logging.info("Parsing GTF file [%s]", args.gtffile)
    gtffile = HTSeq.GFF_Reader(args.gtffile)
    tsspos, n_tss_used  = gtf_to_tsspos(gtffile, up + extra, down + extra)

    bamfile = HTSeq.BAM_Reader(args.bamfile)
    density, n_reads = make_density(tsspos, bamfile, up + extra, down + extra)
    if args.frag_size == -1:
        frag_size = determine_frag_size(density, extra)
    else:
        frag_size = args.frag_size
    output(density, up + extra, down + extra, extra,
            frag_size, n_tss_used, n_reads)

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("tssd",
            help = "Calculate read density around TSSs",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("bamfile", type = arghelpers.infilename_check,
            help = "Bam file; use '-' for stdin")
    cmdline.add_argument("gtffile", type = arghelpers.infilename_check,
            help = "GTF annotation file; has to have exon_number attribute.")
    cmdline.add_argument("-u", "--upstream", type = int, default = 2000,
            help = "nts upstream of TSS to include [%(default)s]")
    cmdline.add_argument("-d", "--downstream", type = int, default = 2000,
            help = "nts downstream of TSS to include [%(default)s]")
    cmdline.add_argument("-s", "--frag-size", type = int, default = -1,
            help = "pre determined fragment size; if default [%(default)s] determines size estimate from data")
    cmdline.set_defaults(func = process)
