"""
Functions to store and serialize refseqs
"""

import numpy as np
import tqdm
from Bio import SeqIO

def is_set(seq, target):
    """
    Determine if string is in a target set and yield an output corresponding to
    each case
    """
    for i in seq:
        if i in target:
            yield True
        else:
            yield False

def is_set_vector(seq, target):
    """
    Use the numpy library to determine in a large array if values are in
    the set. Very useful to generate the maks.,
    """
    bool_mask = np.array(seq) == target
    return bool_mask


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
        n_set = set(["N"])
        purine_set = set(["A", "G"])
        strong_set = set(["C", "G"])

        chr_mask = {}

        for i in self.get_chr_names():

            #print("Masking Chromosome {}".format(i))

            reference_seq = self.chr_refseq[i]
            # this is the inverse mask, see below
            # maybe a bit dangerous
            n_mask = [i for i in tqdm.tqdm(is_set(reference_seq,
                                                  target=n_set),
                                           total=len(reference_seq))]

            purine_mask = [i for i in tqdm.tqdm(is_set(reference_seq,
                                                       target=purine_set),
                                                total=len(reference_seq))]

            strong_mask = [i for i in tqdm.tqdm(is_set(reference_seq,
                                                       target=strong_set),
                                                total=len(reference_seq))]

            n_mask = np.array(n_mask)
            # here we invert the boolean vector because the positions that are
            # callable are the ones that not have a N and therefore that are
            # not in the set.
            n_mask = np.invert(n_mask)

            purine_mask = np.array(purine_mask)
            strong_mask = np.array(strong_mask)

            chr_mask[i] = (n_mask, purine_mask, strong_mask)

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
