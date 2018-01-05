"""
This module includes classes and methods to store mutations.
Maybe here the randomiation should occur.
"""

from tqdm import tqdm
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import randommut.refseq as rs


# this is the key in the iterative vs vectorized war. I could skip the Mut
# object and then get only the mut set !!!

class Mut(object):
    " Class to store a mutation "
    def __init__(self,
                 chr_id,
                 pos_start,
                 pos_end,
                 sample_id,
                 ref,
                 alt,
                 strand,
                 to0base=True): # default because icgc are 1-based
        self.chr_id = chr_id
        if to0base:
            self.pos_start = pos_start - 1 # all should be 0 based
            self.pos_end = pos_end - 1
        else:
            self.pos_start = pos_start  #  coming from a bed file
            self.pos_end = pos_end
        self.sample_id = sample_id
        # now if not added in a explicit way, the software assumes positive
        # strand that should match the reference positions used as a mask
        # by explicit I mean stran being - or 0. Any other thing will be
        # considered positive.
        if strand in ['-', 0]: # should get more consideration
            self.ref = Seq(ref).reverse_complement()
            self.alt = Seq(alt).reverse_complement()
        else:
            self.ref = Seq(ref)
            self.alt = Seq(alt)


    def get_tr(self):
        " get the transision in MS6 format"
        return "{}>{}".format(self.ref, self.alt)
    def get_tr96(self, masks):
        " get the mutation subtype in MS96 format"
        ctx = self.get_context(masks)
        return "{}>{}".format(ctx, self.alt)
    def get_context(self, masks, k=1):
        """
        Return the context for a position in the mask. Needs the 3 masks
        in the correct order (nmask purine and strong). The k indicates
        how many positions we should move to the upstream and downstream.
        1 per default means the 96 context. It uses by default baseset
        C A, meaning that mutations that have a reference value diferent
        than C or A will be transformed to the other strand by reverese
        complementary
        """
        if self.pos_start == self.pos_end:
            pos = self.pos_start
            left_end = pos-k
            rigth_end = pos+k+1
            mask_cont = tuple(i[left_end:rigth_end] for i in masks)
            context = Seq("".join(rs.maskstoseq(mask_cont)))

            if context[k] not in set(["C", "A"]):
                context = context.reverse_complement()

            context = mask_cont
            return context
        else:
            raise ValueError

    def randomize_mutation_same_context(self, masks, wlength, times, k):
        """
        Function to randomize a mutation n times within the available part
        of the window selected by the length in lr. by default 500kb. This
        distance goes in both directions from the center of the original
        mutation therefore windows of 1Mb by default. The available positions
        in this functions are the ones that present the same context

        It needs the corresponding genomic masks. DO CHECK THE CHR
        """
        pos = self.pos_start

        left_end = pos - wlength # this is important to recover coords after
        right_end = pos + wlength + 1 # because it is 0 based

        if left_end < 0:
            left_end = 0

        # end of the common code
        ctx = Seq(self.get_context(masks, k=k))
        ctx_rev = ctx.reverse_complement()

        left_pair = set([ctx[0], ctx_rev[0]])
        center_pair = set([ctx[1], ctx_rev[1]])
        right_pair = set([ctx[2], ctx_rev[2]])

        reg_masks = [i[left_end:right_end] for i in masks] # only 3 times

        left_mask = compute_bimask(reg_masks, biset=left_pair)
        right_mask = compute_bimask(reg_masks, biset=right_pair)
        center_mask = compute_bimask(reg_masks, biset=center_pair)

        left_mask = left_mask[0:-1]
        left_mask.insert(0, 0)

        right_mask = right_mask[1:]
        right_mask.append(0)

        ctx_mask = np.logical_and(left_mask, right_mask)
        ctx_mask = np.logical_and(ctx_mask, center_mask)

        pvector = mask_to_pvector(ctx_mask)
        int_size = len(pvector)

        randomized_positions = np.random.choice(int_size,
                                                p=pvector,
                                                size=times,
                                                replace=True
                                               )
        randomized_positions = randomized_positions + left_end

        return randomized_positions


    def randomize_mutation_same_pair(self, masks, wlength, times):
        """
        Function to randomize a mutation n times within the available part
        of the window selected by the length in lr. by default 500kb. This
        distance goes in both directions from the center of the original
        mutation therefore windows of 1Mb by default. Available positions are
        defined as a set of same pair positions (A:T or C:G)

        It needs the corresponding genomic masks. DO CHECK THE CHR
        """
        pos = self.pos_start

        left_end = pos - wlength # this is important to recover coords after
        right_end = pos + wlength + 1 # because it is 0 based

        if left_end < 0:
            left_end = 0

        # windows goes from pos start to -wlength and +wlength
        pair_mask = masks[2][left_end:right_end] # 3rd is the strong mask
        n_mask = masks[0][left_end:right_end]

        if (self.ref in ["A", "T"]):
            pair_mask = np.logical_not(pair_mask) # transform to weak mask

        master_mask = np.logical_and(n_mask, pair_mask)

        pvector = mask_to_pvector(master_mask)

        #print("length of pvector {}".format(len(pvector)))
        #print("int space size {}".format(int(wlength*2)))
        int_size = len(pvector)

        randomized_positions = np.random.choice(int_size,
                                                p=pvector,
                                                size=times,
                                                replace=True
                                               )
        randomized_positions = randomized_positions + left_end # get coords

        return randomized_positions


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

        common_length = len(pos_start)

        # HOW DO i CHECK THE SAME LENGTH FOR ALL THE ARRAYS?

        self.mutset = {}

        chr_unique = np.unique(chr_id)

        ref_unique = np.unique(ref)

        sample_unique = np.unique(sample_id)

        for i in chr_unique:
            self.mutset[i] = {}
            for j in ref_unique:
                self.mutset[i][j] = {}
                for k in sample_unique:
                    self.mutset[i][j][k] = []

        for i in range(common_length):
            tmp_mut = Mut(chr_id=chr_id[i],
                          pos_start=pos_start[i],
                          pos_end=pos_end[i],
                          sample_id=sample_id[i],
                          ref=ref[i],
                          alt=alt[i],
                          strand=strand[i],
                          to0base=to0base)

            self.mutset[chr_id[i]][ref[i]][sample_id[i]].append(tmp_mut)
    def randomize_mutset_pair(self, times, wlength, genome_object):
        """
        Function to generate random positions conserving the same pair
        """
        rnames = ["r" + str(s) for s in range(1, times+1)]
        rnames.append('ori')
        rnames.append('sample_id')
        rnames.append('ori_context')
        rnames.append('tr')

        rand_df = pd.DataFrame(columns=rnames)

        for i in genome_object.chromosome_iterator():
            if i.chr_id in self.mutset: # into the chr level
                print("Starting chromosome {}".format(i.chr_id))
                for j in ["A", "C", "T", "G"]: #  here I am filtering already rigth?
                    for k in tqdm(self.mutset[i.chr_id][j]): # sample level
                        for m in self.mutset[i.chr_id][j][k]: # mutation level
                            #print(m.get_context(k=1, masks=i.seq_mask))
                            #print(m.get_tr())
                            if i.seq_mask[0][m.pos_start]:
                                # only those that are inside the mask
                                rpos = m.randomize_mutation_same_pair(
                                    masks=i.seq_mask,
                                    times=times,
                                    wlength=wlength).tolist()
                                rpos.append(i.chr_id)
                                rpos.append(m.pos_start)
                                rpos.append(m.sample_id)
                                rpos.append(m.get_context(i.seq_mask))
                                rpos.append(m.get_tr())
                                #zip(rnames,rpos)
                                df2 = pd.DataFrame([rpos], columns=rnames)
                                rand_df = rand_df.append(df2)

        return rand_df

    def randomize_mutset_context(self, times, k, wlength, genome_object):
        """
        to change
        """
        pass





def num2chr(chr_array):
    " Transform the 37 notation to hg19 "
    res = np.array(["chr{}".format(i) for i in chr_array])
    return res


def mask_to_pvector(mask):
    """
    It transforms a mask object (np.boolean array) into a probability vector of
    the same size with an equal value for all the unmasked indexs equivalent
    to the 1 / (number of masked elements). Thus,
    T T T F F F T F F F
    will be
    0.25 0.25 0.25 0 0 0 0.25 0 0 0
    """
    pvector = np.zeros(len(mask)) # the zeros function add 0 as a float

    prob_val = 1 / np.sum(mask)

    pvector[mask] = prob_val

    # necessary? it makes the function crash, choice already handles this
    #if pvector.sum() != 1:
    #    raise ValueError

    return pvector


def compute_bimask(masks, biset):
    """
    Internally transforms inputed masks to desidred bisets of bases. This
    mean that I can transform the 3 masks in a mask that contains all G and Ts.
    """
    # masks reminder
    # index 0 - n mask
    # index 1 - purine (AG)
    # index 2 - strong (CG)

    strong_set = set(["C", "G"])
    weak_set = set(["A", "T"])

    purine_set = set(["A", "G"])
    pyrimidine_set = set(["C", "T"])

    gt_set = set(["G", "T"])
    ca_set = set(["C", "A"])

    # I check which option are we talking about and do the appropiate masks
    if biset == ca_set:
        mask_res = np.logical_xor(masks[1], masks[2])
    elif biset == gt_set:
        mask_res = np.logical_xor(masks[1], masks[2])
        mask_res = np.logical_not(mask_res)
    elif biset == strong_set:
        mask_res = masks[2]
    elif biset == weak_set:
        mask_res = np.logical_not(masks[2])
    elif biset == pyrimidine_set:
        mask_res = np.logical_not(masks[1])
    elif biset == purine_set:
        mask_res = masks[1]

    # important to not forget the n_mask
    mask_res = np.logical_and(mask_res, masks[0])

    return mask_res
