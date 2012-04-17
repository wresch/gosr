"""
basic information about common genomes
"""

import sys
import logging

class Genome(object):
    """Utility class to store and retrieve basic properties of chromosomes of a
    genome"""
    def __init__(self, chromosomes, sizes):
        """takes a list of chromosome names and a list of sizes and creates
        a standard data structure for accessing information about each
        chromosome. chromosomes should be in the order that is the desired sort
        order."""
        sizeL = [long(x) for x in sizes]
        if len(chromosomes) != len(sizes):
            logging.error("list of chromosomes and lengths have to be the same length")
            logging.error(chromosomes)
            logging.error(sizes)
            sys.exit(1)
        self.chromosomes = chromosomes
        self.__size      = dict(zip(chromosomes, sizeL))
        self.__order     = dict(zip(chromosomes, range(len(chromosomes))))
        self.__offsets   = dict(zip(chromosomes, 
            [sum(sizeL[0:i]) for i in range(len(sizeL))]))
        self.__offset_list = sorted(self.__offsets.items(), key = lambda x: x[1])
    
    def size(self, chrom):
        """size of chrom"""
        return self.__size[chrom]
    def order(self, chrom):
        """sort order of chromosome"""
        return self.__order[chrom]
    def offset(self, chrom):
        """0-based offset of start of chrom if chromosomes are concatenated in
        sort order"""
        return self.__offsets[chrom]
    def cpos(self, chrom, pos):
        """Position of base pos (0-based) in chrom if chromosomes were
        concatenated in sort order; result is 0-based"""
        return self.__offsets[chrom] + pos
    def cpos1(self, chrom, pos):
        """Position of base pos (1-based) in chrom if chromosomes were
        concatenated in sort order; result is 0-based"""
        return self.__offsets[chrom] + pos - 1
    def cpos2chrom(self, cpos0):
        """Take a 0-based position in the concatenated genome and
        return a chrom, pos tuple that maps the concatenated genome
        position to a chromosome and a 0-based position"""
        assert cpos0 >= 0
        chrom  = None
        offset = None
        for c, o in self.__offset_list:
            if cpos0 >= o:
                chrom = c
                offset = o
            else:
                break
        return chrom, cpos0 - offset    

#===============================================================================
# Data
#===============================================================================

mm9 = Genome(("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chrX", "chrY", "chrM"),
    (197195432, 181748087, 159599783, 155630120, 152537259, 149517037,
        152524553, 131738871, 124076172, 129993255, 121843856, 121257530,
        120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430,
        166650296, 15902555, 16299))
hg19 = Genome(("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"),
    (249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
        115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983,
        63025520, 48129895, 51304566, 155270560, 59373566, 16571))
hg18 = Genome(("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY",
    "chrM"),
    (247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
        158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
        114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651,
        62435964, 46944323, 49691432, 154913754, 57772954, 16571))
