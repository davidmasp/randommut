"""
Instructions to randomize the genome
"""

import sys
import numpy as np
from tqdm import tqdm
from Bio.Seq import Seq
#from memory_profiler import profile

def rand_single_chr_in_batch(chromosome_object,
                             mutset_object_all,
                             times,
                             winlen,
                             verbose,
                             batch_size=1000):
    """
    Calls the rand_single_chr in batch mode
    """
    if len(mutset_object_all) <= batch_size:
        rand_out = rand_single_chr(chromosome_object=chromosome_object,
                                   mutset_object=mutset_object_all,
                                   times=times,
                                   winlen=winlen,
                                   verbose=verbose)
    else:
        out_all = []
        for mutset_object in mutset_object_all.divide_batch(batch_size):
            out_tmp = rand_single_chr(chromosome_object=chromosome_object,
                                      mutset_object=mutset_object,
                                      times=times,
                                      winlen=winlen,
                                      verbose=verbose)
            out_all.append(out_tmp)
        rand_out = np.concatenate(out_all, axis=0)

    return rand_out



def rand_single_chr(chromosome_object, mutset_object, times, winlen, verbose):
    """
    Performs the rendomization process in just one chromosome
    """
    # 1. obtain matrix mask
    # 2. obtain contexts plus ctx_idx = {}
    # 3. determine operation per ctx
    # 4. loop though ctx_idx to obtain final matrix mask
    # 5. randomize
    # 6. get back to mutset idx

    #1
    mask_matrix_raw = generate_mask_matrix(mutset_object=mutset_object,
                                           chromosome_object=chromosome_object,
                                           winlen=winlen)

    #2
    contexts = mutset_object.get_context(chromosome_object)

    ctx_idx = {} # this is important to go back to the mutset ordering

    for idx, ctx in enumerate(contexts):
        ctx_revcomp = Seq(ctx).reverse_complement()
        ctx_revcomp = str(ctx_revcomp)
        if ctx in ctx_idx:
            ctx_idx[ctx].append(idx)
        elif ctx_revcomp in ctx_idx:
            ctx_idx[str(ctx_revcomp)].append(idx)
        else:
            ctx_idx[ctx] = [idx]

    #3 + 4
    ctx_matrix = {}
    for ctx in ctx_idx:
        if verbose:
            tqdm.write("{}\n".format(ctx))
        current_idx = ctx_idx[ctx]
        current_mask_matrix = [i[current_idx,] for i in mask_matrix_raw]

        # all this could go into a function FROM HERE
        # obtain the mask for ctx
        set_left = set(ctx[0])
        set_right = set(ctx[2])
        set_center = set(ctx[1])

        mask_left = compute_bimask(current_mask_matrix, set_left)
        mask_right = compute_bimask(current_mask_matrix, set_right)
        mask_center = compute_bimask(current_mask_matrix, set_center)

        mask_left_shifted = np.apply_along_axis(shift5,
                                                axis=1, # this should be rows
                                                arr=mask_left,
                                                num=1, # this moves to the right
                                                fill_value=False)

        mask_right_shifted = np.apply_along_axis(shift5,
                                                 axis=1, # this should be rows
                                                 arr=mask_right,
                                                 num=-1, # this moves to  left
                                                 fill_value=False)

        mask_ends_context = np.logical_and(mask_left_shifted,
                                           mask_right_shifted)
        mask_context_sense = np.logical_and(mask_center, mask_ends_context)
        # TO HERE

        # obtain the mask for reverse complement
        ctx_revcomp = Seq(ctx).reverse_complement()

        set_left = set(ctx_revcomp[0])
        set_right = set(ctx_revcomp[2])
        set_center = set(ctx_revcomp[1])

        mask_left = compute_bimask(current_mask_matrix, set_left)
        mask_right = compute_bimask(current_mask_matrix, set_right)
        mask_center = compute_bimask(current_mask_matrix, set_center)

        mask_left_shifted = np.apply_along_axis(shift5,
                                                axis=1, # this should be rows
                                                arr=mask_left,
                                                num=1, # this moves to the right
                                                fill_value=False)

        mask_right_shifted = np.apply_along_axis(shift5,
                                                 axis=1, # this should be rows
                                                 arr=mask_right,
                                                 num=-1, # this moves to  left
                                                 fill_value=False)

        mask_ends_context = np.logical_and(mask_left_shifted,
                                           mask_right_shifted)

        mask_context_reverse = np.logical_and(mask_center, mask_ends_context)

        mask_final = np.logical_or(mask_context_reverse, mask_context_sense)

        number_ctx = np.sum(mask_final)
        if verbose:
            tqdm.write("Number of muts in context: {}\n".format(number_ctx))
        ctx_matrix[ctx] = mask_final

    # 5
    rand_res = {}
    for ctx in ctx_matrix:
        rand_result = randomize_mask_matrix(ctx_matrix[ctx], times=times)
        rand_res[ctx] = rand_result

    # 6
    pos = mutset_object.pos
    left_end = pos[:, 0] - winlen # start minus width

    # this is the output, then the cols are max 100 (no memory issue here)
    rand_pos = np.empty((len(pos), times), dtype=int)
    rand_pos.fill(0) # safety first

    for ctx in rand_res:
        idx = ctx_idx[ctx]
        # we go back to oiginal
        original_positions = np.array(left_end[idx]).reshape([len(idx), 1])
        rand_pos[idx, ] = rand_res[ctx] + original_positions
    return rand_pos

# extracted from https://stackoverflow.com/a/42642326/5410410
# preallocate empty array and assign slice by chrisaycock
def shift5(arr, num, fill_value=np.nan):
    """
    function that shifts num elements of a numpy array to the right (positive)
    or left (negative) and incorporates in the missing spots the fill-value
    argument (nan) by default.
    """
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result = arr
    return result


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

    total_pos = np.sum(mask)

    if total_pos != 0:
        prob_val = 1 / total_pos
        pvector[mask] = prob_val
    elif total_pos == 0:
        pvector.fill(-1)
    else:
        raise ValueError

    # necessary? it makes the function crash, choice already handles this
    #if pvector.sum() != 1:
    #    raise ValueError

    return pvector


def randomize_mask_row(mask, times):
    """
    Here we are tranforming a definitive mask into a pvector and then
    chosing randomly the idx that will be selected. These are equivalent
    to the relatives positions in the sequence.
    """
    # why replace = True? The idea is that different randomizations are pure
    # independent experiments, then, this should be rep true.
    pvector = mask_to_pvector(mask)
    int_size = len(pvector) # this should be equibalent to the wl*2

    if pvector[0] != -1:
        rand_idx = np.random.choice(int_size,
                                    p=pvector,
                                    size=times,
                                    replace=True)
    else:
        rand_idx = np.zeros(times)

    return rand_idx

def randomize_mask_matrix(mask_matrix, times):
    """
    It applies the function randomize_mask_row row by row in a
    definitive Mutset Mask
    """
    idx_matrix = np.apply_along_axis(randomize_mask_row,
                                     axis=1, # this is rows
                                     arr=mask_matrix,
                                     times=times) # using the kargs argument

    return idx_matrix

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

    c_set = set(["C"])
    a_set = set(["A"])
    t_set = set(["T"])
    g_set = set(["G"])

    n_set = set(["N"])
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
    elif biset == c_set:
        not_purine = np.logical_not(masks[1])
        mask_res = np.logical_and(not_purine, masks[2])
    elif biset == t_set:
        not_purine = np.logical_not(masks[1])
        not_strong = np.logical_not(masks[2])
        mask_res = np.logical_and(not_purine, not_strong)
    elif biset == g_set:
        mask_res = np.logical_and(masks[2], masks[1])
    elif biset == a_set:
        not_strong = np.logical_not(masks[2])
        mask_res = np.logical_and(masks[1], not_strong)
    elif biset == n_set:
        # see isssue #1 in gitea
        # I am setting the whole mask to F
        # then I add a step in the rand process when if this happens
        # I generate -1 as randomized positions.
        # update: see n. asjasdhajs
        mask_res = np.zeros(np.shape(masks[0]), dtype=bool)
    else:
        sys.stderr.write(biset)
        raise ValueError

    # important to not forget the n_mask
    mask_final = np.logical_and(mask_res, masks[0])

    return mask_final

#@profile
def generate_mask_matrix(mutset_object, chromosome_object, winlen):
    """
    It generates mask matrix
    """
    winlen = winlen + 1 # I think (due to the 0 based stuff)
    total_length = winlen + winlen + 1 # this because of the middle position

    # here is the memory issue
    # NOTE asjasdhajs
    # trying to solve it with this
    # np.ones((2, 2), dtype=bool)
    # from here
    # https://stackoverflow.com/a/21174962/5410410
    # After some test in this approach we reduce the virtual memory generation
    # in the initialitzation phase (reduction of half)
    # Though the RES memory is the same and it behaves similar, with a
    # reduction after a while. I guess this is numpy storing a number vs
    # when it realises in only a boolean array. So it is an improvement I g.
    # Not final solution though.

    n_muts = len(mutset_object.pos)

    mask_matrix = [np.zeros([n_muts, total_length], dtype=bool),
                   np.zeros([n_muts, total_length], dtype=bool),
                   np.zeros([n_muts, total_length], dtype=bool)]

    masks = chromosome_object.seq_mask

    chromosome_length = len(masks[0])

    for j, original_mask in enumerate(masks):
        for i, val in enumerate(mutset_object.pos): #this should enumerate rows
            left_end = val[0] - winlen # this is the start and then i substract
            right_end = val[1] + winlen # this is the end and then i add

            if left_end < 0:
                #tqdm.write("SHORT position found \n")
                false_pos = winlen - val[0]
                # see update in asjasdhajs
                #mask_matrix[j][i, :false_pos] = False
                mask_matrix[j][i, false_pos:] = original_mask[:right_end]
            elif right_end > chromosome_length:
                #tqdm.write("LONG position found \n")
                false_pos = int(right_end) - int(chromosome_length) #this is p
                end_of_world = total_length - false_pos
                mask_matrix[j][i, :end_of_world] = original_mask[left_end:chromosome_length]
                #mask_matrix[j][i, end_of_world:] = False
                # see update in  asjasdhajs
            else:
                mask_matrix[j][i,] = original_mask[left_end:right_end]

    return mask_matrix
