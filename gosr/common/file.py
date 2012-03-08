import subprocess
import logging

class FileOrGzip(object):
    """Wrap regular file, gzipped file, or stdin in a context"""
    def __init__(self, filename):
        self.filetype = None
        if filename.endswith(".gz"):
            self.gz = subprocess.Popen(["zcat", filename], stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE, close_fds = True)
            self.fh = self.gz.stdout
            self.filetype = "gzip"
        elif filename == "-":
            self.fh = sys.stdin
            self.filetype = "stdin"
        else:
            self.fh = open(filename, "r")
            self.filetype = "regular"
    def __enter__(self):
        return self.fh
    def __exit__(self, etype, evalue, traceback):
        if self.filetype in ("gzip", "regular"):
            self.fh.close()
        if self.filetype == "gzip":
            err = self.gz.stderr.read()
            self.fh.close()
            self.gz.wait()
            if self.gz.returncode != 0:
                # ignore complaint about broken pipe from zcat
                if "Broken pipe" not in err:
                    logging.error("gzip exited with return code %d", self.gz.returncode)
                    logging.error(err)
                    sys.exit(1)
        if etype is not None:
            logging.error("An exception occured while in FileOrGzip context:")

