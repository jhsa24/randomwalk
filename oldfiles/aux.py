#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various helper functions for the random walk simulations
Currently a bit empty
"""
import numpy as np

def mean(lst):
    return sum(lst) / len(lst)

"""
Several functions to help mask creation
"""
def combine_masks(mask_list):
    return np.bitwise_or.reduce(mask_list)

def mask_past(index_np, iteration_np, walker_index, walker_iteration):
    return (index_np == walker_index) & (iteration_np >= walker_iteration-1)

def mask_relatives(index_np, iteration_np, walker_index, walker_dict):
    parent = walker_dict[walker_index]["parent"]
    w = walker_dict[walker_index]["walker"]
    if parent == None or w.iteration > 1:
        return np.zeros_like(index_np, dtype=bool)
    pbt = walker_dict[parent]["branch_time"]
    mask1 = (index_np == parent) & (iteration_np >= pbt-1)
    
    sibling = walker_dict[walker_index]["sibling"]
    mask2 = (index_np == sibling) & (iteration_np <= 1)
    
    grandparent = walker_dict[parent]["parent"]
    if (grandparent != None) and (pbt-1 < 0):
        gpbt = walker_dict[grandparent]["branch_time"]
        mask3 = (index_np == grandparent) & (iteration_np >= gpbt-1)
        
        return combine_masks([mask1, mask2, mask3])
    return combine_masks([mask1, mask2])


def mask_parent(index_np, iteration_np, parent_index, walker_iteration, walker_dict):
    if (parent_index == None) or (walker_iteration > 1):
        return np.zeros_like(index_np, dtype=bool)
    grandparent_index = walker_dict[parent_index]["parent"]
    pbt = walker_dict[parent_index]["branch_time"]
    mask1 = (index_np == parent_index) & (iteration_np >= pbt-1)
    if (grandparent_index != None) and (pbt-1 < 0):
        gpbt = walker_dict[grandparent_index]["branch_time"]
        mask2 = (index_np == grandparent_index) & (iteration_np >= gpbt-1)
        return np.bitwise_or(mask1, mask2)
    else:
        return mask1
    
def mask_sibling(index_np, iteration_np, sibling_index, walker_iteration, walker_dict):
    if sibling_index == None or walker_iteration > 1:
        return np.zeros_like(index_np, dtype=bool)
    return (index_np == sibling_index) & (iteration_np <= 1)
    
    