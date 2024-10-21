import os
import sys
from tqdm import tqdm
import numpy as np
import time
from pymatgen.core.structure import Structure

start_time = time.time()
assert len(sys.argv) == 1
openmx_overlap_path = "/share/home/zhangtao/software/openmx-overlap/openmx3.9/source/openmx"                          #sys.argv[2]
pot_path = "/share/home/zhangtao/software/openmx-pure/openmx3.9/DFT_DATA19"                                           #sys.argv[3]


struct_nums=[90.1]
def cif_to_openmx(struct,struc_index):
    global openmx_overlap_path, pot_path
    stru_unit=struct
    stru = Structure(
        stru_unit.lattice,
        stru_unit.species,
        stru_unit.frac_coords,
        coords_are_cartesian=False,
        to_unit_cell=True)
    species_list = [str(species) for species in stru_unit.species]
    species_Se_list = [species for species in species_list if species == 'Se']
    num_Se = int(len(species_Se_list))
    num_atom=int(len(stru))
    frac_coords = stru.frac_coords
    frac_coords_str = ''
    for i in range(len(stru)):
        if i < num_Se:
            frac_coords_str += f' {i + 1} Ge {frac_coords[i, 0]} {frac_coords[i, 1]} {frac_coords[i, 2]} 2.0  2.0  0.0  0.0  0.0  0.0  0  '  # need set for openmx format
        else:
            frac_coords_str += f' {i + 1} Se {frac_coords[i, 0]} {frac_coords[i, 1]} {frac_coords[i, 2]} 3.0  3.0  0.0  0.0  0.0  0.0  0  '  # need set for openmx format
        if i != len(stru) - 1:
            frac_coords_str += '\n'

    save_str = f"""System.Name                       openmx
DATA.PATH                         {pot_path}
HS.fileout                        on

Species.Number                     2
<Definition.of.Atomic.Species
  Ge    Ge7.0-s3p2d2      Ge_PBE19
  Se    Se7.0-s3p2d2      Se_PBE19
Definition.of.Atomic.Species>
Atoms.Number                      {num_atom}
Atoms.SpeciesAndCoordinates.Unit  FRAC
<Atoms.SpeciesAndCoordinates
{frac_coords_str}
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit            Ang
<Atoms.UnitVectors
{stru_unit.lattice}
Atoms.UnitVectors>

scf.XcType                    GGA-PBE
scf.SpinPolarization          off
scf.ElectronicTemperature     300.0
scf.energycutoff              300.0
scf.maxIter                   500
scf.EigenvalueSolver          band
scf.Kgrid                     1  1  1
scf.Mixing.Type               rmm-diisk
scf.Init.Mixing.Weight        0.05
scf.Min.Mixing.Weight         0.01
scf.Max.Mixing.Weight         0.30
scf.Mixing.History            25
scf.Mixing.StartPulay         15
scf.criterion                 1.0e-7

################vdW correction#################
scf.dftD                     on            # on|off, default=off
version.dftD                  3            # 2|3, default=2
DFTD3.damp                   bj            # zero|bj, default=bj
DFTD.Unit                    AU            # Ang|AU
DFTD.rcut_dftD            100.0            # default=100 (DFTD.Unit)
DFTD.cncut_dftD              40            # default=40 (DFTD.Unit)
DFTD.IntDirection         0 0 0            # default=1 1 1 (1:on 0:off),direction x,y,z;DFTD.IntDirection control the vdW interaction between the central cell and the periodic images.


MD.Type                       nomd
MD.maxIter                    1
MD.TimeStep                   1.0
MD.Opt.criterion              0.0003

Dos.fileout                   off
Dos.Erange                    -25.0  20.0
Dos.Kgrid                     2  2  2

Band.dispersion               off
Band.Nkpath                   1
<Band.kpath
   1  0.000000 0.000000 0.000000   0.000000 0.000000 0.000000   G G
Band.kpath>
### END ###
"""  # openmx, openmx_in.dat
    os.makedirs(f'./{struc_index}', exist_ok=True)
    with open(f'{struc_index}/openmx_in.dat',
              'w') as save_f:
        save_f.write(save_str)
    stru.to(fmt='poscar',
            filename=f'/share/home/zhangtao/work/GeSe/deeph/GeSe/work_dir/overlap/New/{struc_index}/POSCAR')
###main execution##########################################################
for struct_num in tqdm(struct_nums,desc="Processing structures"):
    stru_unit = Structure.from_file(f"./GeSe.cif")
    cif_to_openmx(stru_unit,struct_num)



end_time = time.time()
elapsed_time = end_time - start_time
print('spent time:', elapsed_time)
