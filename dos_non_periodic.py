import numpy as np
import os,sys
import argparse
import matplotlib.pyplot as plt

def dos_plot(dat,fermi):
    fig, ax = plt.subplots(1, 1)
    ax.plot(dat[:, 0] - fermi, dat[:, 1])

    ax.set_ylabel('density of states')
    ax.set_xlabel('Energy (eV)')
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-0.5, 80.0)
    # 在指定的 x 坐标上绘制一条虚线
    x_coordinate = 0  # 你想要绘制虚线的 x 坐标
    ax.axvline(x=x_coordinate, color='r', linestyle='--', linewidth=1)
    # 保存图像
    plt.savefig('dos_plot.png')  # 你可以更改文件名和格式，如 'dos_plot.pdf'

def func_gaussian(E, En, sigma):
    occu_probability = np.exp(-((E - En) / sigma) ** 2)
    return occu_probability


def func_dos(E: np.array, En: np.array, sigma):
    m, n = len(E), len(En)

    dos = []
    for i in range(m):
        states = 0
        for j in range(n):
            state = func_gaussian(E[i], En[j], sigma)
            states = states + state
            # print(E[i],En[j],state,states)

        dos_line = [E[i], states]
        dos.append(dos_line)
    return dos



def read_data(file: str):
    assert os.path.isfile(file), f"File does not exist {file}"
    print(f"File exists:{file}")

    with open(file, "r") as f:
        lines = f.readlines()
        count = 0
        eigen_values = []
        print(len(lines))
        for line in lines:
            if not line.isspace():

                parts = line.strip().split()
                try:
                    eigen_value = float(parts[1])
                    eigen_values.append(eigen_value)
                    count += 1
                except ValueError:
                    pass
            else:
                pass
        print("the infomation of eigenvalues:", np.shape(eigen_values), type(eigen_values), "number of eigenvalues:",
              count)
    return eigen_values

def main(fermi:float,sigma:float,point_density:float,input:str,output:str):

    energy_window=[-5.0,-3.0]
    sampling_point=int((energy_window[1]-energy_window[0])/point_density)
    #eigenvalues=np.array(read_data('./GeSe_90.0-kpoint-G.BANDDAT1'))+fermi
    #eigenvalues=np.array(read_data('./GeSe_90.2.BANDDAT1'))
    eigenvalues=np.array(read_data(input))+fermi
    eigen_min,eigen_max=np.min(eigenvalues),np.max(eigenvalues)
    #energy=np.linspace(eigen_min,eigen_max,sampling_point)
    energy=np.linspace(energy_window[0],energy_window[1],sampling_point)
    print(f'energy range from {eigen_min} to {eigen_max}')
    dat=np.array(func_dos(energy,eigenvalues,sigma))
    dos_plot(dat,fermi)
    #print(energy[0],eigenvalues[0])
    print(np.shape(dat))
    np.savetxt(output,dat)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate the DOS for non-periodic system!!")
    parser.add_argument(
        "--fermi",
        "-f",
        type=float,
        required=True,
        default= -4.048,
        help="the fermi energy of the system, default -4.048 eV",
    )
    parser.add_argument(
        "--sigma",
        "-s",
        type=float,
        default=0.01,
        required=True,
        help="the smearing width of Gaussian method default 0.01 eV",
    )
    parser.add_argument(
        "--point_density",
        "-pd",
        type=float,
        default=0.001,
        help="the sampling density of energy window (-5.0,-3.0)!!",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default="openmx.BANDDAT1",
        help="the eigenvalue file of system!!",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="dos.dat",
        help="the output data file name!!",
    )
    args = parser.parse_args()
    main(**vars(args))