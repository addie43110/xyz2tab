#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys
#import re
#import itertools                                                #for r-length tuples, in sorted order, no repeated elements
#import pandas as pd                                             #pandas tables
#import numpy as np                                              #for calculations
#from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
#from tabulate import tabulate                                   #nice table output
#import matplotlib.pyplot as plt                                 #for molecule display
                     #for molecule display
#from mpl_toolkits.mplot3d import proj3d                         #for fancy arrows in xyz

#from lookup_tables import covalent_radii, atomic_weights, utf_sub_dict
from helpers import parse_args
from print_tables import PrintTab
from create_plots import CreatePlot

#pd.set_option("display.max_rows", None, "display.max_columns", None)


def main():
    args = parse_args()
    
    #for windows console
    sys.stdout.reconfigure(encoding='utf-8')  


    pt = PrintTab(args)

    if args.dihedral:
        print_dihedral_angle(pt.xyz_df, args)

    cp = CreatePlot(args, pt.xyz_df, pt.sel_dist2)
    

if __name__=="__main__":
    main()