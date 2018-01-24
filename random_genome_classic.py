#!/usr/bin/env python

"""
    RANDOM MUT
    ~~~~~~~~~~~~~

    This script uses the package randommut to randomize user set positions.

    :copyright: 2018 by my David Mas @ IRBBarcelona, AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

import os
import sys
import pickle
import pandas as pd
from tqdm import trange
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

def randomize(muts_path, genome_path, assembly, times, winlen, verbose):
    """
    perform  the randomization
    """
    # Optimize this monster please
    # genome is a list of chromosomes (iteration through this)
    # mutset is a dictionary of chomosomes

    if genome_path.endswith('.p'):
        sys.stderr.write("Genome in pickle format detected\n")
        genome_path_pickle = genome_path
        genome = pickle.load(open(genome_path_pickle, "rb"))
        if genome.assembly != assembly:
            raise ValueError
    elif genome_path.endswith(('.fa', '.fasta')):
        genome = gn.genome_from_path(genome_path, assembly)

    sys.stderr.write("Genome read\n")

    muts = mt.mutset_from_path(muts_path)
    sys.stderr.write("Muts read\n")

    randomize_output = {}

    # adding the progress bar
    chromosome_list = genome.chr_list
    progress_chr = trange(len(chromosome_list))
    for chrom_idx in progress_chr:
        chrom = chromosome_list[chrom_idx]
        chr_id = chrom.chr_id
        progress_chr.set_description('Randomizing {}'.format(chr_id))
        if chr_id in muts:
            mutset = muts[chr_id]
            randomize_output[chr_id] = rnd.rand_single_chr(chrom,
                                                           mutset,
                                                           times,
                                                           winlen,
                                                           verbose = verbose)
        else:
            continue
    sys.stderr.write("Rand output generated\n")
    # recover all positions

    full_df = []
    for chrom in genome.chr_list: #opt
        chr_id = chrom.chr_id
        if chr_id in muts:
            mutset = muts[chr_id]
            pos_df = pd.DataFrame(mutset.pos)
            pos_df.columns = ['start', 'end']
            meta_df = pd.DataFrame(mutset.meta)
            meta_df.columns = ['sample', 'ref', 'alt']

            rand_out = randomize_output[chr_id]
            rand_df = pd.DataFrame(rand_out)
            rand_df.columns = ["R{}".format(i+1) for i in range(times)]

            tmp_full = pd.concat([pos_df, meta_df, rand_df], axis=1)
            #import pdb; pdb.set_trace()
            tmp_full["ctx"] = mutset.get_context(chrom) #opt
            tmp_full["chr"] = chr_id
            tmp_full["strand"] = 1
            # this should order the columns
            # chr start end strand ref alt R1 ... Rtimes
            cols = tmp_full.columns.tolist()
            cols = cols[-2:-1] + cols[0:2] + cols[-1:] + cols[2:3] + cols[3:5] + cols[-3:-2] + cols[5:-3]
            tmp_full = tmp_full[cols]
            full_df.append(tmp_full)
        else:
            continue

    #can hapen if only 1 chromosome
    if len(full_df) > 1:
        final_df = full_df[0].append(full_df[1:])
    else:
        final_df = full_df[0]


    return final_df


def remove_clusters():
    """
    Remove the clusters by randomizing one time and getting the randomized
    positions.
    """
    pass


def add_custom_clusters():
    """
    Add clusters based on a defined distributions of distances.
    In particular we need to give mean and standard deviation of the
    normal log10 distribution.
    In order to generate new distances we apply 10**norm(mean,sd)
    """
    pass

def write_randomized_positions(randomize_output, outfilename, compression):
    """
    write the df in a file using pandas
    """
    sys.stderr.write("Writting...\n")
    randomize_output.to_csv(outfilename,
                            sep="\t",
                            header=True,
                            index=False,
                            compression=compression)
    sys.stderr.write("Results file available at {}\n".format(outfilename))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-M", "--mode", type=str,
                        choices=['serialize', 'randomize'],
                        help="The mode to run the software")
    parser.add_argument("-g", "--genome",
                        type=str,
                        help="The path to the reference sequence or serialized genome")
    parser.add_argument("-m", "--muts",
                        type=str,
                        help="The path to the muts file in avinput format",
                        default=None)
    parser.add_argument("-a", "--assembly",
                        type=str,
                        help="Assembly of the genome",
                        default="hg19")
    parser.add_argument("-o", "--outfile", default="rand_positions.tsv")
    parser.add_argument("-t", "--times", type=int, default=10,
                        help="Number of randomized positions")
    parser.add_argument("-w", "--winlen", type=int, default=50000,
                        help="Length of the windows to randomize")
    parser.add_argument("-C", "--compression", type=str,
                        choices=['gzip', 'bz2', None],
                        default=None,
                        help="If the output should be compressed")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="Script returns a more verbose messages.")

    args = parser.parse_args()

    if args.mode == "serialize":
        serialize_genome(genome_path=args.genome, assembly=args.assembly)
    elif args.mode == "randomize":
        if args.muts is None:
            sys.exit('Mut file needed. add with -m tag')
        out_obj = randomize(muts_path=args.muts,
                            genome_path=args.genome,
                            assembly=args.assembly,
                            times=args.times,
                            winlen=args.winlen,
                            verbose = args.verbose)

        write_randomized_positions(randomize_output=out_obj,
                                   outfilename=args.outfile,
                                   compression=args.compression)
    else:
        sys.exit("Wrong mode inputed")
