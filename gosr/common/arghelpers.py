import os
import argparse

def infilename_check(s):
    """raises ArgumentTypeError if file s does not exist; handles '-' as special
    case for stdin; returns string unchanged"""
    if s == "-":
        return "-"
    else:
        if not os.path.exists(s):
            msg = "File %s does not exist" % s
            raise argparse.ArgumentTypeError(msg)
        return s
