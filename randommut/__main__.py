"""
    RANDOM MUT
    ~~~~~~~~~~~~~

    This script uses the package randommut to randomize user set positions.

    :copyright: 2018 by my David Mas @ IRBBarcelona, AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

import sys
import argparse
import randommut.randommutworkflow as wf

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
parser.add_argument("-b", "--batch_size", type=int, default=10000,
                    help="Length of the randomization batch. Decrease if memory error.")
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
    wf.serialize_genome(genome_path=args.genome, assembly=args.assembly)
elif args.mode == "randomize":
    if args.muts is None:
        sys.exit('Mut file needed. add with -m tag')
    out_obj = wf.randomize(muts_path=args.muts,
                           genome_path=args.genome,
                           assembly=args.assembly,
                           times=args.times,
                           b_size=args.batch_size,
                           winlen=args.winlen,
                           verbose=args.verbose)

    wf.write_randomized_positions(randomize_output=out_obj,
                                  outfilename=args.outfile,
                                  compression=args.compression)
else:
    sys.exit("Wrong mode inputed")
