# input working path and number of frames
import os
import sys
from ase.io import read,write

def main():
    if len(sys.argv) < 3:
        print("Usage: python movie.py base_directory frames_numbers")
        sys.exit()

    base_dir = sys.argv[1]                  # working path
    frames_number=int(sys.argv[2])               #number of frames
    movie_file = open('movie.xyz', 'w')
    #print(base_dir)
    trajs=[]

    for i in range(frames_number):
        dir_name = os.path.join(base_dir, str(i))
        #print('this is the work path', dir_name)
        if not os.path.isdir(dir_name):
            print(f"this is not a dir {os.path.isdir(dir_name)}")
            continue
        os.chdir(dir_name)
        poscars = [f for f in os.listdir() if f.endswith('POSCAR')]
        for f in os.listdir():
            if f.endswith('POSCAR'):
                trajs1 = read(f, ':', 'vasp')
                #print(trajs1)
                trajs=trajs+trajs1
        os.chdir(base_dir)
    os.system(f" cd .")
    write(movie_file, trajs, 'extxyz')


if __name__ == '__main__':
    main()




