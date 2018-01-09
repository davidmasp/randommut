#!/usr/bin/env python

import os
import pickle
import randommut.genome as gn
import randommut.muts as mt
import randommut.randomize as rnd


def serialize_genome(genome_path, assembly):
    """
    Get a fasta file and transform to a pickle object
    """
    if genome_path.endswith(('.fa', '.fasta')):
        genome = gn.genome_from_path(genome_path, assembly)
        base = os.path.basename(genome_path)

        genome_path_pickle = "{}{}".format(base, ".p")

        pickle.dump(genome, open(genome_path_pickle, "wb"))
    else:
        raise ValueError

def randomize(muts_path, genome_path, assembly, times, winlen):
    """
    perform the randomization
    """

    # genome is a list of chromosomes (iteration through this)
    # mutset is a dictionary of chomosomes

    if genome_path.endswith('.p'):
        genome_path_pickle = genome_path
        genome = pickle.load(open(genome_path_pickle, "rb"))
        if genome.assembly != assembly:
            raise ValueError
    elif genome_path.endswith(('.fa', '.fasta')):
        genome = gn.genome_from_path(genome_path, assembly)

    muts = mt.mutset_from_path(muts_path)

    randomize_output = {}
    for chrom in genome.chr_list:
        chr_id = chrom.chr_id()
        if chr_id in muts:
            mutset = muts[chr_id]
            randomize_output[chr_id] = rnd.rand_single_chr(chrom,
                                                           mutset,
                                                           times,
                                                           winlen)
        else:
            continue

    return randomize_output