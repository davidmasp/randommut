"""
Functions to store and serialize refseqs
"""

import sys
import re
import numpy as np
from Bio import SeqIO



class RefSeq(object):
    """
    Class to store the original sequence of the Reference Sequence
    """
    def __init__(self, refseq_filepath, assembly):
        self.chr_refseq = SeqIO.to_dict(SeqIO.parse(refseq_filepath, "fasta"))
        self.assembly = assembly
    def __len__(self):
        return len(self.chr_refseq)
    def get_chr_names(self):
        " Return a list with the names of the chromosomes read "
        chrnames = [i for i in self.chr_refseq]
        return chrnames
    def dress_up_seq(self):
        " Mask the genome "
        n_re = re.compile(str("[N]"))
        purine_re = re.compile("[AaGg]")
        strong_re = re.compile("[CcGg]")
        chr_mask = {}

        sys.stderr.write("Starting Serialization\n")

        for chr_id in self.get_chr_names():

            sys.stderr.write("Starting Serialization of {}\n".format(chr_id))

            reference_seq = str(self.chr_refseq[chr_id].seq)

            n_mask = np.empty(len(reference_seq), dtype=bool)
            n_mask.fill(1) # all values to True
            idx = [i.start() for i in re.finditer(n_re, reference_seq)]
            n_mask[idx] = 0 # I put False when the position is masked

            purine_mask = np.empty(len(reference_seq), dtype=bool)
            purine_mask.fill(0) # all values to False
            idx = [i.start() for i in re.finditer(purine_re, reference_seq)]
            purine_mask[idx] = 1 # insert T when regex is matched

            strong_mask = np.empty(len(reference_seq), dtype=bool)
            strong_mask.fill(0) # all values to False
            idx = [i.start() for i in re.finditer(strong_re, reference_seq)]
            strong_mask[idx] = 1

            chr_mask[chr_id] = (n_mask, purine_mask, strong_mask)

        return chr_mask


def maskstoseq(masks):
    " Transform the 3 masks back to a sequence "
    res = np.empty(len(masks[0]), dtype=str)

    npos = np.logical_and(np.logical_not(masks[0]),
                          np.logical_and(np.logical_not(masks[1]),
                                         np.logical_not(masks[2])))

    gpos = np.logical_and(masks[0],
                          np.logical_and(masks[1],
                                         masks[2]))

    cpos = np.logical_and(masks[0],
                          np.logical_and(np.logical_not(masks[1]),
                                         masks[2]))

    apos = np.logical_and(masks[0],
                          np.logical_and(masks[1],
                                         np.logical_not(masks[2])))

    tpos = np.logical_and(masks[0],
                          np.logical_and(np.logical_not(masks[1]),
                                         np.logical_not(masks[2])))

    res[tpos] = "T"
    res[apos] = "A"
    res[cpos] = "C"
    res[gpos] = "G"

    res[npos] = "N"

    return res
