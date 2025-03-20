#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys
import re
#import itertools                                                #for r-length tuples, in sorted order, no repeated elements
#import pandas as pd                                             #pandas tables
#import numpy as np                                              #for calculations
#from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
#from tabulate import tabulate                                   #nice table output
#import matplotlib.pyplot as plt                                 #for molecule display
                     #for molecule display
#from mpl_toolkits.mplot3d import proj3d                         #for fancy arrows in xyz

from lookup_tables import avg_bond_lengths
from helpers import parse_args
from print_tables import PrintTab
from create_plots import CreatePlot

#pd.set_option("display.max_rows", None, "display.max_columns", None)

def classify_bond(pair, bond_length):
    [atom1, atom2] = pair.split("–")

    match = re.search(r"^([a-zA-Z]+)(\d+)", atom1)
    atom1_label = match.group(1)
    match = re.search(r"^([a-zA-Z]+)(\d+)", atom2)
    atom2_label = match.group(1)

    if atom1_label < atom2_label:
        bond_key = atom1_label + "-" + atom2_label
    else:
        bond_key = atom2_label + "-" + atom2_label

    bond_length_dict = avg_bond_lengths[bond_key]
    closest = 'single'
    for (order, length) in bond_length_dict.items():
        length_diff = abs(bond_length - float(length))
        if length_diff < abs(bond_length_dict[closest] - bond_length):
            closest = order

    return [atom1, atom2, bond_length, bond_length_dict[closest], closest]

def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  


    pt = PrintTab(args)
    pt.print_sel_dist_table()
    pt.print_verbose_bond_table()
    pt.print_short_bond_table()
    pt.print_statistics_bond_table()

    bond_table = pt.pr_sel_dist
    classification = [classify_bond(x,y) for (x,y) in zip(bond_table['A-B'], bond_table['distance_calc'])]
    [print(x) for x in classification]

    # cp = CreatePlot(args, pt.xyz_df, pt.sel_dist2)
    

if __name__=="__main__":
    main()