import numpy as np
import argparse
import sys
# from matplotlib.patches import FancyArrowPatch
from scipy.spatial.distance import pdist, squareform, cosine

from lookup_tables import covalent_radii, atomic_weights, utf_sub_dict

#numbers to subscript (utf8) numbers
def num_to_subnum(number):
    utf_number=''
    for letter in str(number):
        utf_letter=utf_sub_dict[letter]
        utf_number=utf_number+utf_letter
    return(utf_number)

#removes the upper triangle of the distance matrix and zeros
#e.g., from d(A-B) = 1.234 Å = d(B-A) =1.234 Å, d(B-A) will be removed 
#d(A-B) = 0 Å will be removed as well
def dm_to_series1(df):
    df = df.astype(float) # do not comment this, angle list will be incomplete
    df.values[np.triu_indices_from(df, k=1)] = np.nan
    #replace zeros with nan
    df = df.replace(0, np.nan)
    #return and drop all nan
    return df.unstack().dropna()

#calculate angle from 3 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z) 
#xyzarr is the array of all atomic coordinates
def calc_angle(xyzarr, i, j, k):
    rij = xyzarr[i] - xyzarr[j]
    rkj = xyzarr[k] - xyzarr[j]
    #remove if cosine fails
    #cos_theta = np.dot(rij, rkj)
    #sin_theta = np.linalg.norm(np.cross(rij, rkj))
    #theta = np.arctan2(sin_theta, cos_theta)
    #scipy pdist cosine instead of the 3 lines above
    theta = cosine(rij,rkj)
    theta = np.arccos(1-theta)
    return     np.degrees(theta)

#calculate the dihedral angle from 4 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z); l(x,y,z)
def calc_d_angle(xyzarr, i, j, k, l):
    #no warning if division by zero
    np.seterr(invalid='ignore')
    rji = -1*(xyzarr[j] - xyzarr[i])
    rkj = xyzarr[k] - xyzarr[j]
    rlk = xyzarr[l] - xyzarr[k]
    rkj /= np.linalg.norm(rkj)
    v = rji - np.dot(rji, rkj)*rkj
    w = rlk - np.dot(rlk, rkj)*rkj
    x = np.dot(v, w)
    y = np.dot(np.cross(rkj, v), w)
    return     np.degrees(np.arctan2(y,x))

#calculation of the best-fit plane
#https://gist.github.com/bdrown/a2bc1da0123b142916c2f343a20784b4
def svd_fit(X):
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, V = np.linalg.svd(CX)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = V[-1]
    return C, N

#https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    '''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    # set to 0.38 --> bigger molecules 
    plot_radius = 0.35*max([x_range, y_range, z_range])
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

#draw fancy arrows in x y z
#https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
""" class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)
    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
        xs, ys, zs = proj3d.proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
        
        xs, ys, zs = proj3d.proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)  """
    
""" def arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow) """

def print_dihedral_angle(xyz_df, args):
    #print the dihedral angle on request
    a1=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[0]].tolist()
    a2=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[1]].tolist()
    a3=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[2]].tolist()
    a4=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[3]].tolist()
    try:
        d_angle = calc_d_angle(xyzarr, *a1, *a2, *a3, *a4)
    except TypeError:
        print('\nWarning! Dihedral angle: One or more atoms could not be found in the input file.')
        sys.exit(1)
    print('\nDihedral angle ' + args.dihedral[0] + '-' +args.dihedral[1] + '-' + \
            args.dihedral[2] + '-' + args.dihedral[3] + ':', f'{d_angle:.2f}°')

def parse_args():
    parser = argparse.ArgumentParser( prog='xyz2tab', 
        description = "Print bond, lengths angles and more from xyz files.")
    #parser.add_argument("filename", 
    #    help = "filename, xyz; e.g. mymolecule.xyz")
    parser.add_argument("allfrags_dir")
    parser.add_argument('-n', '--name', type=str, help="name of the parent molecule (the one being fragmented)", default="unnamed")
    #exclude atoms
    parser.add_argument('-ea','--excludeAt', nargs="+", type=str,
        help='exclude bonds and angles to specified atoms; e.g. -ea N1 or -ea N1 N2')
    #exclude elements
    parser.add_argument('-ee','--excludeEl', nargs="+", type=str,
        help='exclude bonds and angles to specified elements; e.g. -ee C or -ee C N')
    #sort by value
    parser.add_argument('-sa','--sortasc', default=0, action='store_true',
        help='sort values for bond lengths and angles ascending')
    #sort by value
    parser.add_argument('-sd','--sortdes', default=0, action='store_true',
        help='sort values for bond lengths and angles descending')
    #sort by name
    parser.add_argument('-sae','--sortascEl', default=0, action='store_true',
        help='ascending alphabetical sort of elements')
    #sort by name
    parser.add_argument('-sde','--sortdesEl', default=0, action='store_true',
        help='descending alphabetical sort of elements')
    #include contacts
    parser.add_argument('-ic','--includeCon', nargs="+", type=str,
        help='include contacts, e.g. -ic C0 C11 or -ic C0 C11 N14')
    #calculate dihedral angle of selected atoms
    parser.add_argument('-d','--dihedral', nargs=4, type=str,
        help='calculate the dihedral angle of 4 atoms, e.g. -d C0 C11 C12 N13')  
    #calculate the best plane 1
    parser.add_argument('-p1','--plane1', nargs='+', type=str,
        help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
            and show the distances to plane 1, \
            e.g. -p1 C0 C11 C12 N13 or -p1 C0 : N13 or -p1 :')
    #calculate the best plane 2
    parser.add_argument('-p2','--plane2', nargs='+', type=str,
        help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
            and show the distances to plane 1 and 2 and the angle to plane 1, \
            e.g. -p2 C21 C22 C23 N23 or -p2 C21 : N23 or -p2 :')
    #add +x% to radius
    parser.add_argument('-r','--radius', default=8, type=float,
        help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 8 %%')
    #verbose
    parser.add_argument('-v','--verbose', default=0, action='store_true',
        help='verbose print, includes 2 more tables')
    #index
    parser.add_argument('-i','--index', default=0, action='store_true',
        help='the index for the first atom is 1 instead of 0')
    #plot
    parser.add_argument('-s','--show', default=0, action='store_true',
        help='plot xyz coordinates, bonds and planes')
    #plot with bond lengths
    parser.add_argument('-sb','--showbl', default=0, action='store_true',
        help='same as -s with bond lengths')
    #plot with no labels
    parser.add_argument('-sn','--shownl', default=0, action='store_true',
        help='same as -s with no labels')
    #show orientation
    parser.add_argument('-so','--showori', default=0, action='store_true',
        help='plot three arrows along the xyz-axes at 0,0,0 to show the orientation of the molecule')

    return parser.parse_args(sys.argv[1:])