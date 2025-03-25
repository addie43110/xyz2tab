#+x % to covalence radius
#radii_ext = 8 / 100 

#covalent radii from Alvarez (2008)
#DOI: 10.1039/b801115j
covalent_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28,
    'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
    'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
    'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
    'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
    'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
    'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
    'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
    'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
    'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
    'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
    'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
    'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
    'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
    'Am': 1.80, 'Cm': 1.69
}

#atomic weights
atomic_weights = {
    'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
    'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
    'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
    'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
    'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
    'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
    'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
    'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
    'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
    'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
    'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
    'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
    'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
    'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
    'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
    'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
    'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
    'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
    'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
    'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
    'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
    'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
    'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
    'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
}

#dict for numbers to subscript numbers
utf_sub_dict = {
    "0" : "₀",
    "1" : "₁",
    "2" : "₂",
    "3" : "₃",
    "4" : "₄",
    "5" : "₅",
    "6" : "₆",
    "7" : "₇",
    "8" : "₈",
    "9" : "₉",
}

avg_bond_lengths = {
    'C-C' : { 'single' : 1.535, 'double' : 1.339, 'triple' : 1.203},
    'C-H' : { 'single' : 1.090},
    'C-O' : { 'single' : 1.430, 'double' : 1.210},
    'C-N' : { 'single' : 1.430, 'double' : 1.380, 'triple' : 1.160},
    'H-O' : { 'single' : 0.960},
    'N-O' : { 'single' : 1.360},
    'H-N' : { 'single' : 1.020}
}