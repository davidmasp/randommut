"""
This module contains the main class for parsing the genome
This is a class that stores the 4 masks needed for the randomization process

At some point I can include some kind of automatic download from UCSC
"""

import randommut.refseq as rs


def genome_from_path(refseq_path, assembly):
    " File to load a genome directly from files "
    refseq_object = rs.RefSeq(refseq_path, assembly)
    refseq_mask = refseq_object.dress_up_seq()

    chr_list = [] # make this a dictionary?
    for i in refseq_object.get_chr_names():
        chr_tmp = Chromosome(seq_mask=refseq_mask[i], assembly=assembly, chr_id=i)
        chr_list.append(chr_tmp)
    return Genome(chr_list, assembly)

class Genome(object):
    """
    Class to store the 4 masks that will repreesent the available chromosome
    """
    def __init__(self, chr_list, assembly):
        self.chr_list = chr_list
        self.assembly = assembly

    def __len__(self):
        " Get the number of chromosomes stored in the genome "
        return len(self.chr_list)
    def chromosome_list(self):
        " Get the chromosome names "
        return [i for i in self.chr_list]
    def get_assembly(self):
        " Get the inputed assembly "
        return self.assembly
    def chromosome_iterator(self):
        " Yield a chromosome object in the list "
        for i in self.chr_list:
            yield i
    def get_chr_len(self):
        chr_len = {}
        for i in self.chromosome_iterator():
            chr_len[i.chr_id] = len(i)

        return chr_len

class Chromosome(object):
    """
    Class to store the 3 masks that will repreesent the available chromosome
    """
    def __init__(self, seq_mask, assembly, chr_id):
        self.assembly = assembly
        self.seq_mask = seq_mask
        self.chr_id = chr_id
    def n_mask(self):
        " returns the N mask, True is no N, False is N"
        return self.seq_mask[0]
    def purine_mask(self):
        " returns the purine mask, True is A G, False is C T"
        return self.seq_mask[1]
    def strong_mask(self):
        " returns the strong mask, True is C G, False is A T"
        return self.seq_mask[2]
    def __len__(self):
        "returns the length of the chromosome"
        return len(self.seq_mask[2])
