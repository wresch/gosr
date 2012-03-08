"""
Based on zip file output of fastqc, creates a PDF summary of 
the data file and extracts the data file for later storage in 
the archive

Dependencies:
    * R/ggplot2.
    * Latex
"""

import sys
import os
import re
import string
import tempfile
import logging
import shutil
import zipfile
import subprocess
import argparse
import shlex
from gosr.common import arghelpers

#================================================================================
# DATA
#================================================================================
doc = string.Template(r"""
\documentclass{article}
\usepackage{fullpage}
\usepackage[usenames,dvipsnames]{color}
\usepackage{bookman}
\usepackage{rotating}

\newcommand{\pass}{\large{\textbf{\textcolor{Green}{PASS}}}}
\newcommand{\warn}{\large{\textbf{\textcolor{Orange}{WARN}}}}
\newcommand{\fail}{\large{\textbf{\textcolor{Red}{FAIL}}}}

\title{Run $safe_runid}
\author{gosr fastqc2pdf}

\begin{document}
\renewcommand*\listfigurename{Summary}
\maketitle
\listoffigures
%% ===========================================================================================
<<echo=F>>=
library(ggplot2)
theme.qc                  <- theme_bw(12)
theme.qc$panel.border     <- theme_rect(fill = NA, size = 1, col = "black")
theme.qc$panel.grid.minor <- theme_blank()
theme.qc$panel.grid.major <- theme_blank()

mycol <- list(yellow = rgb(1, 1, 157/255))
@
\section*{Module results}
$figures

\end{document}
""")
#================================================================================
# CODE
#================================================================================

def basic_statistics(status, datastr, tempdir):
    """render the Basic Statistic module for sweave"""
    safe_datastr = datastr.replace("_", "\_")
    lines = safe_datastr.split("\n")
    fields = [c.strip().split("\t") for c in lines \
              if not (c.strip() == "" or c.startswith("#"))]
    table = "\n".join(r"%s & %s \\" % (a, b) for a, b in fields)
    section = string.Template(r"""
\begin{figure}[h!]\centering
\begin{tabular}{rl}
\hline \\
$table
\hline
\end{tabular}
\caption[\$status~~Basic statistics]{\$status~~Basic statistics}
\end{figure}""")
    return section.substitute(locals())

def per_base_sequence_quality(status, datastr, tempdir):
    """render the per base sequence quality module for sweave.
    In some cases, fastqc does not do a boxplot for each base, 
    so only start and end of the read are labeled"""
    datafile = os.path.join(tempdir, "perbase.tsv")
    plotfile = os.path.join(tempdir, "perbase.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "base\tmean\tmedian\tlq\tuq\tp10\tp90"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perBaseQual,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
p<- ggplot(df) +
    geom_boxplot(aes(x = factor(base), lower = lq, upper = uq,
        ymin = p10, ymax = p90, middle = median), 
        stat = "identity", fill = "grey80") +
    geom_line(aes(x = base, y = mean), col = "blue", size = 1) +
    scale_x_discrete("Position", breaks = c(1, max(df$base)),
        labels = c("Start", "End")) +
    ylab("Quality") +
    ylim(0, 40) +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption[\$status~~Per base sequence quality]{\$status~~Per base sequence
quality.  Range of quality values; Box spans 25$^{th}$ to 75$^{th}$ percentiles;
lines span 10$^{th}$ to 90$^{th}$ percentile; blue line is mean}
\end{figure}""")
    return section.safe_substitute(locals())

def per_sequence_quality_scores(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "perseq.tsv")
    plotfile = os.path.join(tempdir, "perseq.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "quality\tn"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perSeqQual,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
p <- ggplot(df) +
    geom_ribbon(aes(quality, ymin = 0, ymax = n), fill = "grey80",
        col = "black") +
    ylab("Frequency") +
    xlab("Quality") +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption{\$status~~Per sequence quality scores.}
\end{figure}""")
    return section.safe_substitute(locals())

def per_base_sequence_content(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "perbasecont.tsv")
    plotfile = os.path.join(tempdir, "perbasecont.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "pos\tG\tA\tT\tC"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perBaseContent,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
df <- cbind(pos = df[, 1], sweep(df[, 2:5], 1, rowSums(df[, 2:5]), "/"))
df <- melt(df, id = "pos")
p <- ggplot(df) +
    geom_line(aes(pos, value, col = variable), size = 1) +
    xlab("Position") +
    scale_y_continuous("Fraction", formatter = "percent",
        limits = c(0, 0.5)) +
    scale_color_manual("", values = c(T = "red", C = "blue", 
        A = "green", G = "black")) +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption{\$status~~Per base sequence content}
\end{figure}""")
    return section.safe_substitute(locals())

def per_base_gc_content(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "perbasegc.tsv")
    plotfile = os.path.join(tempdir, "perbasegc.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "pos\tGC"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perBaseGC,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
p <- ggplot(df) +
    geom_line(aes(pos, GC), size = 1) +
    xlab("Position") +
    ylab("GC content [%]") +
    ylim(0, 100) +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption{\$status~~Per base GC content}
\end{figure}""")
    return section.safe_substitute(locals())

def per_sequence_gc_content(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "perseqgc.tsv")
    plotfile = os.path.join(tempdir, "perseqgc.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "GC\tn"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perSeqGC,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
total <- sum(df$n)
meanGC <- sum(df$n * df$GC) / total
sdGC   <- sqrt(sum((df$GC - meanGC) ^ 2 * df$n / total))
df$theoretical <- dnorm(df$GC, meanGC, sdGC) * total
p <- ggplot(df) +
    geom_bar(aes(GC, n), stat = "identity", fill = "grey80", col = "black") +
    geom_line(aes(GC, theoretical), col = "blue", size = 1) +
    xlab("GC content [%]") +
    ylab("Observations") +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption[\$status~~Per sequence GC content]{\$status~~Per sequence GC content;
Observed: bar graph; theoretical: blue line}
\end{figure}""")
    return section.safe_substitute(locals())

def per_base_n_content(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "perbasen.tsv")
    plotfile = os.path.join(tempdir, "perbasen.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "base\tnn"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<perBaseN,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
p <- ggplot(df) +
    geom_line(aes(base, nn), size = 1) +
    scale_y_continuous("N [%]", limits = c(0, 100)) +
    xlab("Position") +
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption{\$status~~Per base N content}
\end{figure}""")
    return section.safe_substitute(locals())

def sequence_length_distribution(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "seqlen.tsv")
    plotfile = os.path.join(tempdir, "seqlen.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        lines[0] = "len\tn"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<seqlen,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
df <- rbind(c(min(df$len) - 1, 0), df, c(max(df$len) + 1))
p <- ggplot(df) +
    geom_bar(aes(len, n), stat = "identity", fill = "grey80", col = "black") +
    ylab("Observations") +
    xlab("Sequence length") + 
    theme.qc
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption{\$status~~Distribution of sequence lengths in run}
\end{figure}""")
    return section.safe_substitute(locals())

def sequence_duplication_level(status, datastr, tempdir):
    datafile = os.path.join(tempdir, "dup.tsv")
    plotfile = os.path.join(tempdir, "dup.pdf")
    with open(datafile, "w") as df:
        lines = datastr.split("\n")
        total_dup = float(lines.pop(0).split()[3])
        total_dup = "%.2f%%" % total_dup
        lines[0] = "dup\tn"
        for line in lines:
            df.write(line)
            df.write("\n")
    section = string.Template(r"""
\begin{figure}[h!]\centering
<<seqDup,echo=F>>=
df <- read.table("$datafile", sep = "\t", header = T)
df$dupn <- 1:nrow(df)
p <- ggplot(df) +
    geom_line(aes(dupn, n), size = 1) +
    geom_point(aes(dupn, n), shape = 21, fill = "white",
        col = "black", size = 2) +
    ylab("Relative abundance") +
    scale_x_continuous("Duplication level", breaks = df$dupn,
        labels = df$dup) +
    theme.qc +
    opts(title="Sequence duplication level >= $total_dup")
ggsave("$plotfile", p, width = 6, height = 3)
@
\includegraphics[width=6in, height=3in]{$plotfile}
\caption[\$status~~Sequence duplication level]{\$status~~Sequence duplication
level.  This module counts the degree of duplication for every sequence in the
set and creates a plot showing the relative number of sequences with different
degrees of duplication.  Only sequences that occur in the first 200,000
sequences are considered to conserve memory. The last bin contains any sequence
with 10 \emph{or more} duplicates. Reads longer than 75nts are truncated to
50nts.}
\end{figure}""")
    return section.safe_substitute(locals())

def overrepresented_sequences(status, datastr, tempdir):
    safe_datastr = datastr.replace("_", r"\_").replace("%", r"\%")
    lines = safe_datastr.split("\n")
    header = " & ".join(lines.pop(0)[1:].strip().split("\t")) + r"\\ \hline"
    fields = [c.strip().split("\t") for c in lines \
              if not c.strip() == ""]
    tablelist = [header]
    for f in fields:
        f[2] = "%.2f" % float(f[2])
        tablelist.append(" & ".join(f) + r"\\")
    table = "\n".join(tablelist)
    section = string.Template(r"""
\begin{sidewaysfigure}[h!]\centering{\tiny
\begin{tabular}{rrrp{2in}}
\hline \\
$table
\hline
\end{tabular}}
\caption{\$status~~Overrepresented sequences}
\end{sidewaysfigure}""")
    return section.substitute(locals())


def kmer_content(status, datastr, tempdir):
    safe_datastr = datastr.replace("_", r"\_").replace("%", r"\%")
    lines = safe_datastr.split("\n")
    header = " & ".join(lines.pop(0)[1:].strip().split("\t")) + r"\\ \hline"
    fields = [c.strip().split("\t") for c in lines \
              if not c.strip() == ""]
    tablelist = [header]
    for f in fields:
        f[2] = "%.2f" % float(f[2])
        f[3] = "%.2f" % float(f[3])
        tablelist.append(" & ".join(f) + r"\\")
    table = "\n".join(tablelist)
    section = string.Template(r"""
\begin{figure}[h!]\centering{\tiny
\begin{tabular}{rrrrr}
\hline \\
$table
\hline
\end{tabular}}
\caption{\$status~~Kmer content}
\end{figure}""")
    return section.substitute(locals())


modules = {
    "Basic Statistics": basic_statistics,
    "Per base sequence quality": per_base_sequence_quality,
    "Per sequence quality scores": per_sequence_quality_scores,
    "Per base sequence content": per_base_sequence_content,
    "Per base GC content": per_base_gc_content,
    "Per sequence GC content": per_sequence_gc_content,
    "Per base N content": per_base_n_content,
    #"Sequence Length Distribution": sequence_length_distribution,
    "Sequence Duplication Levels": sequence_duplication_level,
    "Overrepresented sequences": overrepresented_sequences,
    "Kmer Content": kmer_content
}

def fastqc2latex(args):
    """parses a fastqc data file into modules and yields
    each module as a (name, pass/fail, data) tuple"""
    # get a file handle to the data file insize the zip archive
    runid      = args.out_prefix
    safe_runid = args.out_prefix.replace("_", r"\_")
    tempdir    = tempfile.mkdtemp()
    sys.exitfunc = lambda: shutil.rmtree(tempdir)
    datafile   = "fastqc_data.txt"
    
    zipf     = zipfile.ZipFile(args.fastqc, "r")
    datafile_zip_path = None
    for zipf_member in zipf.namelist():
        if zipf_member.endswith(datafile):
            datafile_zip_path = zipf_member
            break
    if datafile_zip_path is None:    
        logging.error("zip file did not contain a 'fastqc_data.txt' file")
        zipf.close()
        sys.exit(1)
    with open(os.path.join(tempdir, datafile), "w") as out:
        out.write(zipf.read(datafile_zip_path))
    zipf.close()
    datafile_path = os.path.join(tempdir, datafile)
    data = open(datafile_path, "rU").read()
    
    module = re.compile(r">>(.+?)\s+(pass|fail|warn)\n(.*?)>>END_MODULE", re.DOTALL)
    figures = []
    for module_name, status, data in module.findall(data):
        try:
            figures.append(modules[module_name](status, data, tempdir))
        except KeyError:
            logging.warn("No function for processing module %s registered",
                    module_name)
            continue
    figures       = "\n\n".join(figures)
    rnw_file_path = os.path.join(tempdir, "%s.rnw" % runid)
    with open(rnw_file_path, "w") as out:
        out.write(doc.safe_substitute(locals()))
    
    # create PDF file
    sweave_cmd = "R CMD Sweave --pdf %s.rnw" % runid
    logging.debug(sweave_cmd)
    sweave = subprocess.Popen(shlex.split(sweave_cmd),
            close_fds = True, cwd = tempdir,
            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    sweave_out, sweave_err = sweave.communicate()
    if sweave.returncode != 0:
        logging.error("R CMD Sweave failed with exit code %s", sweave.returncode)
        print >>sys.stderr, sweave_out
        print >>sys.stderr, sweave_err
        sys.exit(1)
    else:
        shutil.move(os.path.join(tempdir, "%s.pdf" % runid), 
            os.path.join(args.dir, "%s.pdf" % runid))
        shutil.move(datafile_path, os.path.join(args.dir, "%s.fastc_data.txt" % runid))



################################################################################
# tool interface
################################################################################

def setup(commands):
    """set up command line parser"""
    cmdline = commands.add_parser("fastqc2pdf",
            help = """Fastqc zip file to PDF converter""",
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description     = __doc__)
    cmdline.add_argument("fastqc", type = arghelpers.infilename_check,
            help = "Fastqc zip file")
    cmdline.add_argument("out_prefix", 
            help = "prefix of output file names")
    cmdline.add_argument("--dir", "-d", type = arghelpers.check_or_make_dir,
            help = "output directory if other than current working directory;\
            will be created if it does not exist [%(default)s]",
            default = "./")
    cmdline.set_defaults(func = fastqc2latex)


