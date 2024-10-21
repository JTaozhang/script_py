import os
import sys

import numpy as np
from numpy import ndarray
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter


structure = Structure.from_file('GeSe6x6_90.0(20).cif')

coords_z = structure.cart_coords[:,2]
print(coords_z[0:3])
bottom_coords_z = [coord_z for coord_z in coords_z if 12.0 < coord_z < 13.0]
top_coords_z = [coord_z for coord_z in coords_z if 15.0 < coord_z <16.0]
aver_bottom_z = np.mean(bottom_coords_z)
aver_top_z =np.mean(top_coords_z)

print("average bottom coords:",aver_bottom_z,"the number of coords:",len(bottom_coords_z))
print("average top coords:",aver_top_z,"the number of coords:",len(top_coords_z))
print("average interlayer distance:",aver_top_z-aver_bottom_z)

###############################translate the top layer###########################
add_dist=3.005750000000001-(aver_top_z-aver_bottom_z)
coords= structure.cart_coords

indices=[i for i,z in enumerate(coords_z) if z >15]
for indice in indices:
    coords[indice,2] += add_dist
"""
for i, coord in enumerate(coords):
    if coord[2] >3.5:
        coords[i,2] += add_dist

frac_coords=structure.lattice.get_fractional_coords(coords)
test_array: ndarray=np.concatenate((frac_coords[:,0].reshape(-1,1),structure.frac_coords[:,0].reshape(-1,1)),axis=1)
#print(test_array)
assert (frac_coords[:,0]-structure.frac_coords[:,0] < 1e-10).all()
"""
frac_coords=structure.lattice.get_fractional_coords(coords)
test_array: ndarray=np.concatenate((frac_coords[:,0].reshape(-1,1),structure.frac_coords[:,0].reshape(-1,1)),axis=1)
assert (frac_coords[:,0]-structure.frac_coords[:,0] < 1e-10).all()
new_struct=Structure(structure.lattice,structure.species,coords,coords_are_cartesian=True)
print(add_dist)
##################################translate ending####################################
coords_z = new_struct.cart_coords[:,2]
bottom_coords_z = [coord_z for coord_z in coords_z if 12.0 < coord_z < 13.0]
top_coords_z = [coord_z for coord_z in coords_z if 15.0 < coord_z <16.0]
aver_bottom_z = np.mean(bottom_coords_z)
aver_top_z =np.mean(top_coords_z)

print("average bottom coords:",aver_bottom_z,"the number of coords:",len(bottom_coords_z))
print("average top coords:",aver_top_z,"the number of coords:",len(top_coords_z))
print("average interlayer distance:",aver_top_z-aver_bottom_z)

cif_writer = CifWriter(new_struct)
cif_writer.write_file('90.0_6x6.cif')