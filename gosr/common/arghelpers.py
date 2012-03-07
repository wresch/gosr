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

def check_or_make_dir(s):
    """
    check if dir exists; create it if it does not. raise ArgumentTypeError if 
    it can not be created
    """
    if os.path.exists(s):
        if not os.path.isdir(s):
            msg = "%s already exists but is not a directory"
            raise argparse.ArgumentTypeError(msg)
        return s
    else:
        try:
            os.mkdir(s, 0700)
            return s
        except OSError, e:
            raise argparse.ArgumentTypeError("Could not create directory: %s" % e)
