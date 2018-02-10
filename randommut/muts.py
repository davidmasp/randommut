"""
This module includes classes and methods to store mutations.
Maybe here the randomiation should occur.
"""

import re
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import randommut.refseq as rs



def mutset_from_path(av_input_path, to0base=True):
    """
    Obtain a dictoniary of mutset onbjects by chromosome
    """

    # colum definition
    # 0. chr
    # 1. start
    # 2. end
    # 3. ref
    # 4. alt
    # 5. strand
    # 6. sample

    table = pd.DataFrame.from_csv(av_input_path,
                                  header=None,
                                  sep="\t",
                                  index_col=False)

    chrom_idx = {}
    chrom_unique = np.unique(table.iloc[:, [0]])

    for chrom in chrom_unique:
        idx = np.array(table.iloc[:, [0]] == chrom)
        chrom_idx[chrom] = [i for i, x in enumerate(idx) if x]

    ms_chr = {}
    for chrom in chrom_idx:
        idx = chrom_idx[chrom]
        #tmp_table = table[idx,]
        ms_chr[chrom] = MutSet(chr_id=chrom,
                               pos_start=table.iloc[idx, 1],
                               pos_end=table.iloc[idx, 2],
                               ref=table.iloc[idx, 3],
                               alt=table.iloc[idx, 4],
                               strand=table.iloc[idx, 5],
                               sample_id=table.iloc[idx, 6],
                               to0base=to0base)

    return ms_chr



class MutSet(object):
    " Class to store a set of mutation "
    def __init__(self,
                 chr_id,
                 pos_start,
                 pos_end,
                 sample_id,
                 ref,
                 alt,
                 strand,
                 to0base=True):

        # each mutset have only one chromosome
        if len(chr_id) != 1:
            if len(np.unique(chr_id)) == 1:
                chr_id = str(np.unique(chr_id))
            else:
                raise ValueError

        self.chr_id = chr_id

        # here I check if all lists have the same length
        common_length = len(pos_start)
        if any(len(lst) != common_length for lst in [pos_end,
                                                     sample_id,
                                                     ref,
                                                     alt,
                                                     strand]):
            raise ValueError

        # if input comes from bed file no need for conversion
        if to0base:
            pos_start = np.array(pos_start,dtype=int) - 1
            pos_end = np.array(pos_end,dtype=int)
            self.pos = np.column_stack((pos_start, pos_end))
        else:
            pos_start = np.array(pos_start)
            pos_end = np.array(pos_end)
            self.pos = np.column_stack((pos_start, pos_end))


        strand_re = re.compile("^[-0]$") # could break not tested
        idx = re.finditer(strand_re, str(strand))

        # this is not super fast but most of the time we will have positive
        for i in idx:
            ref[i] = Seq(ref[i]).reverse_complement()
            alt[i] = Seq(alt[i]).reverse_complement()

        self.meta = np.column_stack((sample_id, ref, alt))
    def __len__(self):
        return len(self.pos)
    def get_chr_id(self):
        "return the chromosome id"
        return self.chr_id
    def get_sample(self):
        " return sample as a list "
        return self.meta[:, 0]
    def get_ref(self):
        " return the references as a list "
        return self.meta[:, 1]
    def get_alt(self):
        " return alternative as a list "
        return self.meta[:, 2]
    def get_tr(self):
        " get the transision in MS6 format"
        def make_tr(arr):
            " small function to generate the tr"
            # this uses the A C format (useful to make generalization here)
            if arr[0] in ["T", "G"]:
                arr[0] = Seq(arr[0]).reverse_complement()
                arr[1] = Seq(arr[1]).reverse_complement()
            return "{}>{}".format(arr[0], arr[1])
        res = np.apply_along_axis(make_tr, axis=1, arr=self.meta[:, [1, 2]])
        return res
    def get_context(self, chromosome_object, k=1):
        """
        Return the context for a position in the mask. Needs the 3 masks
        in the correct order (nmask purine and strong). The k indicates
        how many positions we should move to the upstream and downstream.
        1 per default means the 96 context. It uses by default baseset
        C A, meaning that mutations that have a reference value diferent
        than C or A will be transformed to the other strand by reverese
        complementary
        """
        # Important pos_start != pos_end
        left_end = self.pos[:, 0] - k # this is the start
        rigth_end = self.pos[:, 1] + k  # this is the end

        masks = chromosome_object.seq_mask
        ctx = []

        for i in range(len(left_end)):
            mask_ctx = tuple(j[left_end[i]:rigth_end[i]] for j in masks)
            context = "".join(rs.maskstoseq(mask_ctx))

            ctx.append(context)

        return ctx
    def context_generator(self, chromosome_object, k=1):
        """
        Return the context for a position in the mask. Needs the 3 masks
        in the correct order (nmask purine and strong). The k indicates
        how many positions we should move to the upstream and downstream.
        1 per default means the 96 context. It uses by default baseset
        C A, meaning that mutations that have a reference value diferent
        than C or A will be transformed to the other strand by reverese
        complementary
        """
        # Important pos_start != pos_end
        left_end = self.pos[:, 0] - k # this is the start
        rigth_end = self.pos[:, 1] + k  # this is the end

        masks = chromosome_object.seq_mask

        for i in range(len(left_end)):
            mask_ctx = tuple(j[left_end[i]:rigth_end[i]] for j in masks)
            context = "".join(rs.maskstoseq(mask_ctx))

            yield context

    def get_tr96(self, chromosome_object, k=1):
        " get the mutation subtype in MS96 format"
        ctx = self.get_context(chromosome_object, k=k)
        alt = self.get_alt()
        mut_list = []

        for i, val in enumerate(ctx):
            if val[k] in ["T", "G"]:
                ctx_tr = str(Seq(val).reverse_complement())
                alt_tr = str(Seq(alt[i]).reverse_complement())
            else:
                ctx_tr = val
                alt_tr = alt[i]

            mut = "{}>{}".format(ctx_tr, alt_tr)

            mut_list.append(mut)
        return mut_list
    def divide_batch(self, muts_x_batch):
        " Divide the current mutset object and yield different batches"
        total_muts = len(self)
        iterantions_n = total_muts // muts_x_batch
        iterantions_n = iterantions_n + 1
        
        first_positions = np.multiply(list(range(0, iterantions_n)),
                                      muts_x_batch)
        last_positions = np.multiply(list(range(1, iterantions_n + 1)),
                                     muts_x_batch)
        last_positions[-1] = total_muts
        # I think I do not have to touch this (substracting 1 ) because
        # I use the range functions that already does that.
        print("Number of iterations {}".format(iterantions_n))
        
        for i, s_idx in enumerate(first_positions):
            # s_idx = first_positions[i]
            e_idx = last_positions[i]
            #import pdb; pdb.set_trace()
            length_to_range = e_idx - s_idx
            strand_info = ["+" for i in range(length_to_range)]

            tmp = MutSet(chr_id=self.chr_id,
                         pos_start=self.pos[s_idx:e_idx, 0],
                         pos_end=self.pos[s_idx:e_idx, 1],
                         sample_id=self.meta[s_idx:e_idx, 0],
                         ref=self.meta[s_idx:e_idx, 1],
                         alt=self.meta[s_idx:e_idx, 2],
                         strand=strand_info,
                         to0base=False)

            yield tmp

def num2chr(chr_array):
    " Transform the 37 notation to hg19 "
    res = np.array(["chr{}".format(i) for i in chr_array])
    return res
