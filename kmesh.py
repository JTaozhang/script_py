from Brillouin_zone import get_brillouin_zone_2d,get_reciprocal_lattice,plot2D_BZ
import numpy as np
from scipy.spatial import Voronoi, ConvexHull
from matplotlib.path import Path
from typing import List, Dict, Tuple, Optional, Union, Any
import matplotlib.pyplot as plt
import sys
from ase.io import read

def points_in_voronoi_region(points_to_check:np.ndarray, vor:Voronoi, region_index:int) -> np.ndarray:
    """
    判断多个点是否在指定 Voronoi 区域内（包括边界）

    参数：
        points_to_check: (N, D) ndarray，待判断的点（二维或三维）
        vor: scipy.spatial.Voronoi 对象
        region_index: int，目标 Voronoi 区域编号（如原点对应区域）

    返回：
        mask: (N,) bool 数组，True 表示该点在区域内
    """
    region = vor.regions[region_index]
    if not region or -1 in region:
        raise ValueError("Voronoi 区域无效或为开放区域，无法判断包含关系")

    vertices = vor.vertices[region]
    dim = vertices.shape[1]
    points_to_check = np.atleast_2d(points_to_check)
    # print("vertices",vertices)

    if dim == 2:
        hull = ConvexHull(vertices)
        vertices=vertices[hull.vertices]
        path = Path(vertices)
        return np.array([path.contains_point(p, radius=1e-6) for p in points_to_check])
    
    elif dim == 3:
        try:
            hull = ConvexHull(vertices)
            A, b = hull.equations[:, :-1], hull.equations[:, -1]
            return np.all(A @ points_to_check.T + b[:, None] <= 1e-6, axis=0)
        except:
            raise ValueError("该 Voronoi 区域无法构成有效凸包")
    else:
        raise ValueError("仅支持二维或三维 Voronoi 区域")


def write_kpoint(kpoints:List)->str:
    kpnum = len(kpoints)
    print("number of kpoint:",kpnum)
    file = open('KPOINTS','w+')
    file.write(f'k-points in BZ \n {kpnum}\n Reciprocal\n')
    for k in kpoints:
        file.write('{}   {}   {}   1.0\n'.format(k[0],k[1],k[2]))
    file.close()


def produce_kp(kxend:int,kyend:int,kzend:int,kmeshx:int,kmeshy:int,kmeshz:int,kxstart:int=0,kystart:int=0,kzstart:int=0,dim=2)->np.ndarray:

    if dim ==2:
        knum = kmeshx * kmeshy
        kxpath = np.linspace(kxstart,kxend,kmeshx)
        kypath = np.linspace(kystart,kyend,kmeshy)
        kzpath = np.zeros((knum))
        kx_grid, ky_grid = np.meshgrid(kxpath, kypath)
        extra_points = np.stack([kx_grid.ravel(), ky_grid.ravel(),kzpath], axis=-1)
    if dim ==3:
        knum = kmeshx * kmeshy* kmeshz
        kxpath = np.linspace(kxstart,kxend,kmeshx)
        kypath = np.linspace(kystart,kyend,kmeshy)
        kzpath = np.linspace(kzstart,kzend,kmeshz)
        extra_points = np.stack([kx_grid.ravel(), ky_grid.ravel(),kzpath.ravel()], axis=-1)
    return extra_points

def main(poscar_path,kmesh):

    kmeshx,kmeshy,kmeshz=kmesh,kmesh,kmesh
    dim=2   ###DIMENSION OF SYSTEM
    kxstart = -1.0
    kxend   = 1.0
    kystart = -1.0
    kyend   = 1.0
    kzstart = 0.0
    kzend = 1.0
    extra_points=produce_kp(kxend,kyend,kzend,kmeshx,kmeshy,kmeshz,kxstart,kystart,kzstart)

    ###retrive the reciprocal lattice matrix
    atoms=read(f"{poscar_path}/POSCAR")
    cell=atoms.cell.array
    rec_cell=get_reciprocal_lattice(cell)
    print("reciprocal lattice matrix:",rec_cell)
    kpoints_carts = np.einsum('ik,kl->il', extra_points, rec_cell)


    if dim==2:
        kpoints_carts=kpoints_carts[:,:2]
        assert np.shape(kpoints_carts)[-1]==2
    elif dim==3:
        kpoints_carts=kpoints_carts
    else:
        raise ValueError("仅支持二维或三维 kpoints")

    region_index,vor,edge=get_brillouin_zone_2d(rec_cell)
    filetered_points=points_in_voronoi_region(kpoints_carts,vor,region_index)
    # print(extra_points[filetered_points])
    write_kpoint(extra_points[filetered_points])

    ######################################plottting BZ############################################
    fig=plt.figure(figsize=(6,6))
    ax=fig.add_subplot(111)
    plot2D_BZ(edge,rec_cell,ax)
    ax.scatter(kpoints_carts[filetered_points][:,0],kpoints_carts[filetered_points][:,1],color='b',s=5)
    plt.quiver([0,0], [0,0], rec_cell[:2,0], rec_cell[:2,1], angles='xy', scale_units='xy', scale=1.0, color='k')
    plt.savefig("kpoints@2D.png",dpi=450)
    plt.show()
    plt.close()

if __name__=="__main__":
    if len(sys.argv) != 3:
        print("用法: python kmesh.py  poscar_path kmesh")
        sys.exit(1)
    poscar_path=sys.argv[1]
    kmesh = int(sys.argv[2])
    print("process on!!!")
    main(poscar_path,kmesh)
    print("process finishing!!!")