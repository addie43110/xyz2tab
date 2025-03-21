import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from tabulate import tabulate
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from helpers import arrow3D

# setattr(Axes3D, 'arrow3D', arrow3D)


class CreatePlot:
    def __init__(self, args, xyz_df, sel_dist2):
        self._args = args
        self._xyz_df = xyz_df
        self._sel_dist2 = sel_dist2

        if args.plane1:
            (atom_names, plane1_df, n1) = self._calc_plane_1()
            self._print_plane_1(atom_names, plane1_df)

        if args.plane2:
            if not args.plane1:
                print('\nWarning! Plane 1 must be defined first. Exit.')
                sys.exit()

            (atom_names, plane_df, n2) = self._calc_plane_2()
            self._print_plane_2(atom_names, plane_df)
            self._print_angle_between_planes(n1, n2)

        if args.show or args.showbl or args.shownl:
            #show the molecule and planes
            self._plot_molecule()


    def _calc_plane_1(self):
        args = self._args
        xyz_df = self._xyz_df
        #Plane No. 1 through selected or all atoms on request
        #get the index of the selected atoms
        a1 = xyz_df.index[xyz_df['atom1_idx'].isin(args.plane1)].tolist()
        #check for range notation symbol ':'
        if ':' in args.plane1:
            #get indices for ':'
            sep_indices = [i for i, x in enumerate(args.plane1) if x == ":"]
            #avoid error message for strange input
            start = xyz_df.index[0]
            #loop over indices for ':'
            for sep_idx in sep_indices: 
                #start is at atom index 0 if 0
                if sep_idx == 0:
                    start = xyz_df.index[0]
                else:
                    try:
                        #start index is the atom index of the atom in input, left from ':'
                        start = xyz_df.index[xyz_df['atom1_idx'] == args.plane1[sep_idx-1]].values.astype(int)[0]
                    except IndexError:
                        #catch malformed inputs
                        print('\nWarning! Malformed input.')
                try: 
                    #stop index is the atom index of the atom in input, right from ':'
                    stop = xyz_df.index[xyz_df['atom1_idx'] == args.plane1[sep_idx+1]].values.astype(int)[0]
                except IndexError:
                    #if there is no atom right from ':' take the last atom
                    stop = xyz_df.index[-1]    
                #atom index range, stop+1 because of numpy 
                atm_range = np.arange(start, stop+1)
                #get unique index, avoid duplicates from range input and singel ato input 
                #e.g. 1 2 4 + range 3 4 5 --> 1 2 3 4 4 5 --> 1 2 3 4 5
                a1 = np.unique(list(a1) + list(atm_range))
                
        #get atom names from the selected atoms
        atom_names1 = xyz_df.iloc[a1,1].tolist()
        #get the x,y,z coordinates from the selected atoms
        xyz_pl1_arr = xyzarr[a1]
            
        #check if atom names from input are in the atom list
        #if not, warning
        atom_names_in_arg = [x for x in args.plane1 if x != ':']
        if not all(elem in atom_names1 for elem in atom_names_in_arg):
            print('\nWarning! One or more atoms could not be found in the input file.')
            
        #no atoms / coordinates for calculation --> exit
        if len(xyz_pl1_arr) == 0:
            #if the list is empty --> exit
            print('\nWarning! No atoms for Plane 1. Exit.')
            sys.exit(1)
            
        #calculate the best fit plane
        c1, n1 = svd_fit(np.asarray(xyz_pl1_arr))
            
        #create data frame for output
        plane1_df = pd.DataFrame()
        plane1_df['Atom'] = atom_names1
        #calculate the distance of each atom to the plane 
        plane1_df['Distance'] = np.dot(np.asarray(xyz_pl1_arr)-c1, n1)
        #sum of squares error
        #sosqf_p1  = plane1_df['Distance'].apply(lambda x: abs(x)**2)

        return (atom_names, plane1_df, n1)

    def _print_plane_1(self, atom_names, plane1_df):
        print('\nBest-fit Plane 1 through', len(atom_names), 'atoms.')
        
        #print some plane related parameters
        #print('')
        #print('Centroid: ', *c1, 'Å')
        #print('Plane normal: ', *n1, 'Å')
        #print('Sum-of-squares error:', f'{sosqf_p1.sum():.4f} Å²')
            
        #print the table with atom names and distances to the plane
        print("\n", tabulate(plane1_df,
            headers=['Atom','Distance to Plane 1 /Å'], 
            tablefmt='github',
            floatfmt=(".4f"),
            showindex=False))
        

    #Plane No. 2 through selected atoms on request
    def _calc_plane_2(self):
        xyz_df = self._xyz_df        
        args = self._args
        
        #get the index of the selected atoms
        a2 = xyz_df.index[xyz_df['atom1_idx'].isin(args.plane2)].tolist()
        
        if ':' in args.plane2:
            
            sep_indices2 = [i for i, x in enumerate(args.plane2) if x == ":"]
            start2 = xyz_df.index[0]
            for sep_idx2 in sep_indices2: 
                
                if sep_idx2 == 0:
                    start2 = xyz_df.index[0]
                else:
                    try:
                        start2 = xyz_df.index[xyz_df['atom1_idx'] == args.plane2[sep_idx2-1]].values.astype(int)[0]
                    except IndexError:
                        print('')
                        print('Warning! Malformed input.')
                try: 
                    stop2 = xyz_df.index[xyz_df['atom1_idx'] == args.plane2[sep_idx2+1]].values.astype(int)[0]
                except IndexError:
                    stop2 = xyz_df.index[-1]    
                    
                atm_range2 = np.arange(start2, stop2+1)
                a2 = np.unique(list(a2) + list(atm_range2))
                
        #get atom names from the selected atoms
        atom_names2 = xyz_df.iloc[a2,1].tolist()
        #get the x,y,z coordinates from the selected atoms
        xyz_pl2_arr = xyzarr[a2]
        
        atom_names_in_arg2 = [x for x in args.plane2 if x != ':']
        if not all(elem in atom_names2 for elem in atom_names_in_arg2):
            print('\nWarning! One or more atoms could not be found in the input file.')

        if len(xyz_pl2_arr) == 0:
            #if the list is empty --> exit
            print('\nWarning! No atoms for Plane 2. Exit.')
            sys.exit(1)
            
        #calculate the best fit plane
        c2, n2 = svd_fit(np.asarray(xyz_pl2_arr))
        
        #create data frame for output
        plane2_df = pd.DataFrame()
        plane2_df['Atom'] = atom_names2
        #calculate the distance of each atom to plane 2 
        plane2_df['DistanceP2'] = np.dot(np.asarray(xyz_pl2_arr)-c2, n2)
        #calculate the distance of each atom to plane 1 
        plane2_df['DistanceP1'] = np.dot(np.asarray(xyz_pl2_arr)-c1, n1)
        #sum of squares error
        #sosqf_p2 = plane2_df['DistanceP2'].apply(lambda x: abs(x)**2)
        return (atom_names2, plane2_df, n2)
        
    def _print_plane_2(self, atom_names, plane_df):
        print('\nBest-fit Plane 2 through', len(atom_names), 'atoms.')
        
        #print some plane related parameters
        #print('')
        #print('Centroid: ', *c2, 'Å')
        #print('Plane normal: ', *n2, 'Å')
        #print('Sum-of-squares error:', f'{sosqf_p2.sum():.4f} Å²')
        
        #print the table with atom names and distances to the plane 2 and plane 1 
        print("\n",tabulate(plane_df,
        headers=['Atom','Distance to Plane 2 /Å','Distance to Plane 1 /Å'], 
        tablefmt='github',
        floatfmt=(".4f"),
        showindex=False))
        
    def _print_angle_between_planes(self, n1, n2):
        #calculate the angle between plane 1 and plane 2
        #no warning if senseless result, e.g. plane 1 = plane 2
        np.seterr(invalid='ignore')
        phi = np.arccos(np.dot(n1,n2))
        print('')
        print('Angle between Plane 1 and Plane 2:', f'{np.degrees(phi):.2f}°')

    def _plot_molecule(self):
        args = self._args
        xyz_df = self._xyz_df
        sel_dist2 = self._sel_dist2
        #show the molecule and planes
        #lists for color asignment
        metals = ['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'Ba', 'Be', \
                    'Bi', 'Ca', 'Cd', 'Ce', 'Cf', 'Cm', 'Co', 'Cr', 'Cs', 'Cu', \
                    'Db', 'Dy', 'Er', 'Es', 'Eu', 'Fe', 'Fm', 'Fr', 'Ga', 'Gd', \
                    'Ge', 'Hf', 'Hg', 'Ho', 'Hs', 'In', 'Ir', 'K', 'La', 'Li', \
                    'Lr', 'Lu', 'Md', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Nd', 'Ni', \
                    'Np', 'Os', 'Pa', 'Pb', 'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Pu', \
                    'Ra', 'Rb', 'Re', 'Rf', 'Rh', 'Rn', 'Ru', 'Sc', 'Sm', 'Sn', \
                    'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'U', \
                    'V', 'W', 'Y', 'Yb', 'Zn', 'Zr']
        green = ['F','Cl']
        brown = ['Br']
        purple = ['P','I']
        orange = ['Si']
        
        #get atom numbers from sel_dist2 data frame, e.g. C11 --> 11
        atom1_num = sel_dist2['atom1_idx'].apply(lambda x: re.sub(r'\D+','', x)).astype(int).tolist()
        atom2_num = sel_dist2['atom2_idx'].apply(lambda x: re.sub(r'\D+','', x)).astype(int).tolist()
        #get atom labels from xyz_df
        atom_label = xyz_df['atom1_idx'].tolist()
        #distance as bond label
        bond_label = sel_dist2['distance_calc'].apply(lambda x: '{:.3f}'.format(x)).tolist()
        
        #if first atom has index 1
        if args.index:
            atom1_num  = [idx - 1 for idx in atom1_num]
            atom2_num  = [idx - 1 for idx in atom2_num]

        #atom1 and atom 2 coordinates
        atom1_coord = xyzarr[atom1_num]
        atom2_coord = xyzarr[atom2_num]
        #put atom1 and atom2 coordinats in an numpy array (for bonds display)
        atom1_2_coord = np.array(list(zip(atom1_coord,atom2_coord)))

        #clumsy but safe
        #for assigning colors to different elemments
        carr = xyz_df.index[xyz_df['element'].isin(['C'])].tolist()
        harr = xyz_df.index[xyz_df['element'].isin(['H'])].tolist()
        narr = xyz_df.index[xyz_df['element'].isin(['N'])].tolist()
        oarr = xyz_df.index[xyz_df['element'].isin(['O'])].tolist()
        sarr = xyz_df.index[xyz_df['element'].isin(['S'])].tolist()
        marr = xyz_df.index[xyz_df['element'].isin(metals)].tolist()
        grarr = xyz_df.index[xyz_df['element'].isin(green)].tolist()
        brarr = xyz_df.index[xyz_df['element'].isin(brown)].tolist()
        parr = xyz_df.index[xyz_df['element'].isin(purple)].tolist()
        orarr = xyz_df.index[xyz_df['element'].isin(orange)].tolist()
        restarr = [item for item in atom1_num if item not in carr]
        restarr = [item for item in restarr if item not in harr]
        restarr = [item for item in restarr if item not in narr]
        restarr = [item for item in restarr if item not in oarr]
        restarr = [item for item in restarr if item not in sarr]
        restarr = [item for item in restarr if item not in marr]
        restarr = [item for item in restarr if item not in grarr]
        restarr = [item for item in restarr if item not in brarr]
        restarr = [item for item in restarr if item not in parr]
        restarr = [item for item in restarr if item not in orarr]
        
        #for atom labeling
        atom_coord_name = zip(xyzarr,atom_label)

        #prepare plot
        fig = plt.figure(figsize=(10,8))
        ax = plt.axes(projection='3d')
        #otherwise molecule looks strange
        ax.set_box_aspect((1, 1, 1))
        
        #clumsy but safe
        #plot atom coordinates with different colors and point sizes
        ax.scatter(*xyzarr[carr].T,s=100/len(xyzarr) * 50,color='black')
        ax.scatter(*xyzarr[harr].T,s=80/len(xyzarr) * 50,color='tan')
        ax.scatter(*xyzarr[narr].T,s=100/len(xyzarr) * 50,color='blue')
        ax.scatter(*xyzarr[oarr].T,s=100/len(xyzarr) * 50,color='green')
        ax.scatter(*xyzarr[sarr].T,s=200/len(xyzarr) * 50,color='yellow')
        ax.scatter(*xyzarr[marr].T,s=300/len(xyzarr) * 50,color='red',alpha=0.85)
        ax.scatter(*xyzarr[grarr].T,s=200/len(xyzarr) * 50,color='green')
        ax.scatter(*xyzarr[brarr].T,s=200/len(xyzarr) * 50,color='brown')
        ax.scatter(*xyzarr[parr].T,s=200/len(xyzarr) * 50,color='purple')
        ax.scatter(*xyzarr[orarr].T,s=200/len(xyzarr) * 50,color='orange')
        ax.scatter(*xyzarr[restarr].T,s=100/len(xyzarr) * 50,color='gray')
        
        #label atoms
        #show no lables if -sn option is activated
        if not args.shownl:
            for coord, label in atom_coord_name:
                ax.text(*(coord+0.12).T, label, fontsize=100/len(xyzarr) + 8 , color='black')
        
        #draw bonds
        for bonds, labels in zip(atom1_2_coord, bond_label):
            ax.plot(*bonds.T, color='gray', linewidth=3.0)
            if args.showbl:
                #show bond labels
                ax.text(*np.average(bonds+0.06,axis=0).T,labels,fontsize=(100/len(xyzarr) + 8)/1.5,color='gray')
        
        ax.scatter(*xyzarr[restarr].T,s=100/len(xyzarr) * 50,color='red')
        #define ranges forplanes
        x_pl=np.sort(xyzarr[:,0])
        y_pl=np.sort(xyzarr[:,1])
        z_pl=np.sort(xyzarr[:,2])

        if args.plane1:
            #plane1 grid
            xx1, yy1 = np.meshgrid((x_pl[0],x_pl[-1]),(y_pl[0],y_pl[-1]))
            if args.plane2:
                #plane2 grid
                xx2, yy2 = np.meshgrid((x_pl[0],x_pl[-1]),(y_pl[0],y_pl[-1]))
            #plane 1 d
            d1 = -c1.dot(n1)
            
            if args.plane2:
                #plane 2 d
                d2 = -c2.dot(n2)
            #plane 1 equation
            z1 = (-n1[0] * xx1 - n1[1] * yy1 - d1) * 1. /n1[2]
            if args.plane2:
                #plane 2 equation
                z2 = (-n2[0] * xx2 - n2[1] * yy2 - d2) * 1. /n2[2]
            #plot plane 1
            surf = ax.plot_surface(xx1, yy1, z1, color='blue', alpha=0.3, label='Plane 1')
            
            if args.plane2:
                #plot plane 2
                surf = ax.plot_surface(xx2, yy2, z2, color='red', alpha=0.3, label='Plane 2')
        
        # show arrows representing the xyz-axes, starting from 0,0,0
        if args.showori:
            arrow_length =sum(abs(i) for i in ax.get_xlim())+sum(abs(i) for i in ax.get_ylim())+sum(abs(i) for i in ax.get_zlim())
            arrow_length = (arrow_length/3)*0.5
            if arrow_length > 3:
                arrow_length = 3
            ax.arrow3D(0,0,0, 0,0,arrow_length,
                        mutation_scale=20,
                        ec ='black',
                        fc='red')
            ax.arrow3D(0,0,0, 0,arrow_length,0,
                        mutation_scale=20,
                        ec ='black',
                        fc='green')
            ax.arrow3D(0,0,0, arrow_length,0,0,
                        mutation_scale=20,
                        ec ='black',
                        fc='blue')
            ax.text(0, 0, arrow_length, 'z',color='red',fontsize=15)
            ax.text(0, arrow_length, 0, 'y',color='green',fontsize=15)
            ax.text(arrow_length, 0, 0, 'x',color='blue',fontsize=15)
            ax.scatter(0,0,0,s=50,color='black',alpha=0.8)

        #no axes
        ax.set_axis_off()
        #tight layout 
        fig.tight_layout()
        #ax.legend(loc='upper left')
        #set z limits for plots, otherwise planes are sometimes very large
        plt.gca().set_zlim(z_pl[0],z_pl[-1])
        #adjust 3d drawing behavior, otherwise molecules are not correctly displayes
        set_axes_equal(ax)
        #show the plot
        plt.show()