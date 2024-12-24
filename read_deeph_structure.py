import numpy as np
import os,sys
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import argparse, warnings
from tqdm import tqdm

periodic_table = {'Ac': 89, 'Ag': 47, 'Al': 13, 'Am': 95, 'Ar': 18, 'As': 33, 'At': 85, 'Au': 79, 'B': 5, 'Ba': 56,
                  'Be': 4, 'Bi': 83, 'Bk': 97, 'Br': 35, 'C': 6, 'Ca': 20, 'Cd': 48, 'Ce': 58, 'Cf': 98, 'Cl': 17,
                  'Cm': 96, 'Co': 27, 'Cr': 24, 'Cs': 55, 'Cu': 29, 'Dy': 66, 'Er': 68, 'Es': 99, 'Eu': 63, 'F': 9,
                  'Fe': 26, 'Fm': 100, 'Fr': 87, 'Ga': 31, 'Gd': 64, 'Ge': 32, 'H': 1, 'He': 2, 'Hf': 72, 'Hg': 80,
                  'Ho': 67, 'I': 53, 'In': 49, 'Ir': 77, 'K': 19, 'Kr': 36, 'La': 57, 'Li': 3, 'Lr': 103, 'Lu': 71,
                  'Md': 101, 'Mg': 12, 'Mn': 25, 'Mo': 42, 'N': 7, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10, 'Ni': 28,
                  'No': 102, 'Np': 93, 'O': 8, 'Os': 76, 'P': 15, 'Pa': 91, 'Pb': 82, 'Pd': 46, 'Pm': 61, 'Po': 84,
                  'Pr': 59, 'Pt': 78, 'Pu': 94, 'Ra': 88, 'Rb': 37, 'Re': 75, 'Rh': 45, 'Rn': 86, 'Ru': 44, 'S': 16,
                  'Sb': 51, 'Sc': 21, 'Se': 34, 'Si': 14, 'Sm': 62, 'Sn': 50, 'Sr': 38, 'Ta': 73, 'Tb': 65, 'Tc': 43,
                  'Te': 52, 'Th': 90, 'Ti': 22, 'Tl': 81, 'Tm': 69, 'U': 92, 'V': 23, 'W': 74, 'Xe': 54, 'Y': 39,
                  'Yb': 70, 'Zn': 30, 'Zr': 40, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
                  'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
reversed_periodic_table={v:k for k,v in periodic_table.items()}

def read_deeph_data(dirpath:str,savepath:str):
    element_path=os.path.join(dirpath,"element.dat")
    position_path=os.path.join(dirpath,"site_positions.dat")
    lat_path=os.path.join(dirpath,"lat.dat")
    assert os.path.isfile(element_path), f"File {element_path} does not exist."
    assert os.path.isfile(position_path), f"File {position_path} does not exist."
    assert os.path.isfile(lat_path), f"File {lat_path} does not exist."
    elements=np.loadtxt(element_path)
    positions=np.loadtxt(position_path)
    positions=np.transpose(positions)
    latt=np.loadtxt(lat_path)
    elements=elements.astype(int)
    # print(elements)

    # print("element:",elements[0:2],"number of atom",len(elements),"position:",np.shape(positions),"lattice:",np.shape(latt))
    # element_labels=[]
    # for element in elements:
    #     label=reversed_periodic_table(element)
    #     element_labels.append(label)
    # # element_type,element_counts=np.unique(element_labels,return_counts=True)
    # element_counts={}
    # for element in element_labels:
    #     if element in element_counts:
    #         element_counts[type]+=1
    #     else:
    #         element_type[type]=1
    stru = Structure(
        latt,
        elements,
        positions,
        coords_are_cartesian=True,
        to_unit_cell=True)
    stru.to(fmt="poscar",filename=f"{savepath}/POSCAR")

def main():
    workdir = parsed_args.input_dir
    savedir = parsed_args.output_dir

    dir_name=10
    for i in tqdm(range(576)):
        workpath=os.path.join(workdir,str(i))
        savepath=os.path.join(savedir,str(i))
        os.makedirs(savepath, exist_ok=True)
        os.makedirs(workpath, exist_ok=True)
        read_deeph_data(workpath,savepath)

    print("finish structure transformation!!!")

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir", "-o", type=str, default="./",
        help="dir of output "
    )
    parser.add_argument(
        "--input_dir", "-i", type=str, default="./",
        help="dir of output"
    )
    return parser.parse_args()
if __name__ == '__main__':
    parsed_args = parse_commandline()
    main()