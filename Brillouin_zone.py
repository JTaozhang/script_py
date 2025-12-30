import numpy as np
from scipy.spatial import Voronoi
from typing import List, Dict, Tuple, Optional, Union, Any
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_reciprocal_lattice(cell:np.ndarray):
    return np.linalg.inv(cell).T*2*np.pi 


def get_brillouin_zone_3d(cell:np.ndarray):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """

    rec_cell = np.asarray(cell, dtype=float)
    assert rec_cell.shape == (3, 3)

    px, py, pz = np.tensordot(rec_cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]


    vor = Voronoi(points)
    # 找到最接近原点的点的索引
    origin_index = np.argmin(np.linalg.norm(points, axis=1))
    # 获取该点对应的 Voronoi 区域编号
    region_index = vor.point_region[origin_index]
    # 获取该区域的顶点索引
    region = vor.regions[region_index]
    if -1 in region or len(region) == 0:
        raise ValueError("Voronoi region around origin is unbounded or invalid.")
    bz_facets = []
    bz_ridges = []
    bz_vertices = region

    # 3D: collect ridges and facets connected to origin
    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        if origin_index in pid and -1 not in rid:
            ridge = vor.vertices[np.r_[rid, [rid[0]]]]
            bz_ridges.append(ridge)
            bz_facets.append(vor.vertices[rid])

    bz_vertices = list(set(bz_vertices))

    return region_index, vor, bz_ridges
def get_brillouin_zone_2d(rec_cell:np.ndarray):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """
    rec_cell2d = np.asarray(rec_cell, dtype=float)[:2,:2]
    assert rec_cell2d.shape == (2, 2)
    px, py = np.tensordot(rec_cell2d, np.mgrid[-1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel()]
    vor = Voronoi(points)
    # 找到最接近原点的点
    origin_index = np.argmin(np.linalg.norm(points, axis=1))
    region_index = vor.point_region[origin_index]
    region = vor.regions[region_index]

    if -1 in region or len(region) == 0:
        raise ValueError("Voronoi region around origin is unbounded or invalid.")
    bz_facets = []
    bz_ridges = []
    bz_vertices = region
    # 2D: ridges = edges of polygon
    polygon = vor.vertices[region]
    for i in range(len(polygon)):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % len(polygon)]
        bz_ridges.append(np.array([p1, p2]))
        bz_facets.append(np.array([p1, p2]))

    return region_index, vor, bz_ridges

def plot3D_BZ(data,rec_cell,ax):


    for xx in data:
        ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], color='k', lw=1.0)
    b1,b2,b3=np.linalg.norm(rec_cell, axis=1)
    ax.set_xlim(-b1, b1)
    ax.set_ylim(-b2, b2)
    ax.set_zlim(-b3, b3)


    # plt.savefig("Brillouin3D.png",dpi=450)
    # plt.show()
    # plt.close()
    return ax

def plot2D_BZ(data,rec_cell,ax):
    rec_cell2d=rec_cell[:2,:2]

    for xx in data:
        ax.plot(xx[:, 0], xx[:, 1], color='k', lw=1.0)
    print("shape of rec_cell2d",np.shape(rec_cell2d))
    b1,b2=np.linalg.norm(rec_cell2d, axis=1)
    ax.set_xlim(-b1, b1)
    ax.set_ylim(-b2, b2)
    return ax


def main():
    cell= np.array([[32.2932854426, 0.0000000000,        0.0000000000],
        [-16.1466427212,       27.9668055649,        0.0000000000],
            [0.0000000000,        0.0000000000,       35.0000000000]])
    print(get_reciprocal_lattice(cell))

    rec_cell=get_reciprocal_lattice(cell)
    region_index,vor,e=get_brillouin_zone_2d(rec_cell)
    fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_subplot(111, projection='3d')
    ax =fig.add_subplot(1,1,1)
    print("region_index",region_index,type(region_index))
    plot2D_BZ(e,rec_cell,ax)
    plt.savefig("Brillouin2D.png",dpi=450)
    plt.show()
    plt.close()

if __name__=="__main__":
    main()
