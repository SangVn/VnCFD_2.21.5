# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Hàm xuất lưới
def save_mesh(file_name, nodes_list):
    print('Writing mesh: ', file_name)
    f = open(file_name, 'w')
    f.write('TITLE = "vncfd python"\n')
    f.write('VARIABLES = "X", "Y"\n')

    for n, nodes in enumerate(nodes_list):
        Ni, Nj = nodes.shape[0], nodes.shape[1]
        f.write('ZONE T="%d", I= %d, J= %d\n' % (n+1, Ni, Nj))
        for j in range(Nj):
            for i in range(Ni):
                f.write('%f %f\n' % (nodes[i, j, 0], nodes[i, j, 1]))
    f.close()
    print("Done!")
    
def generate_mesh_mach2(Ni=101, Nj=41):
    # Kích thước vùng tính toán
    ly, lx = 4.0, 10.0
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới 
    nodes = np.zeros((Ni, Nj, 2))
    # tọa độ x tại các điểm lưới
    dx = lx / Ni
    x = np.linspace(0, lx, Ni)
    # tọa độ y của biên dưới
    y0 = np.zeros(Ni)
    # index i tương ứng vị trí x = 2, 4 trên biên dưới
    i2 = int(2./dx)
    i4 = int(4./dx)

    y0[i2:i4] = (x[i2:i4]-2.)*np.tan(np.pi/12)
    y0[i4:] = 2.0*np.tan(np.pi/12)

    # khoảng cách dy giữa hai điểm cùng cột
    dy = np.array([(ly-y)/(Nj-1) for y in y0])
    
    # xác định tọa độ (x, y) của từng điểm 
    for i in range(Ni):
        for j in range(Nj):
            nodes[i, j, 0] = x[i]
            nodes[i, j, 1] = y0[i]+j*dy[i]
            
    np.savez("../mach2/mesh.npz", nodes)
    save_mesh("../mach2/mesh.dat", [nodes])
    return



# Tạo lưới bài toán dòng chảy quanh hình trụ
# x^n - 1 = (x-1)*(1+x+x^2+x^3+...+x^(n-1))
# Rmax-Rmin =  dr*(1+s+s^2+s^3+...+s^(N-1)) = dr*(s^N - 1)/(s-1)
# (Rmax-Rmin)- dr*(S^N-1)/(s-1) = 0
# Task 1: Re=40, Ni=Nj=101, dr0=1.25e-3, ratio=1.1
# Task 2: Re=150, Ni=Nj=201, dr0=1e-3, ratio=fsolve(f(s)) # Rmax=200*(2*Rmin) or #ratio = 1.035

# def generate_mesh_cylinder(Ni=201, Nj=201, dr0=1.0e-3, ratio=None):
def generate_mesh_cylinder(Ni=101, Nj=101, dr0=1.25e-3, ratio=1.1):
    # chia góc 2Pi ra thành Ni điểm lưới
    alpha = np.linspace(0.0, -2 * np.pi, Ni)
    nodes = np.zeros((Ni, Nj, 2))
    Rmin = 0.05
    Rmax = 200*(2*Rmin)
    if ratio is None:
        f = lambda s: (Rmax-Rmin) - dr0*(s**(Nj-1)-1)/(s-1)
        ratio = fsolve(f,2)
        print("ratio: ", ratio, f(ratio))

    for i in range(Ni):
        r = Rmin
        dr = dr0
        for j in range(Nj):
            nodes[i, j, 0] = r * np.cos(alpha[i])
            nodes[i, j, 1] = r * np.sin(alpha[i])
            r += dr
            dr *= ratio
    print("Rmax: ", r)

    np.savez("../cylinder/mesh.npz", nodes)
    save_mesh("../cylinder/mesh.dat", [nodes])
    return

def generate_mesh_cavity(Ni=65, Nj=65):
    # Kích thước vùng tính toán
    ly, lx = 0.0254, 0.0254
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới
    nodes = np.zeros((Ni, Nj, 2))
    # tọa độ x tại các điểm lưới
    x = np.linspace(0, lx, Ni)
    y = np.linspace(0, ly, Nj)
    # xác định tọa độ (x, y) của từng điểm
    for i in range(Ni):
        for j in range(Nj):
            nodes[i, j, 0] = x[i]
            nodes[i, j, 1] = y[j]
    np.savez("../lid_driven_cavity/mesh.npz", nodes)
    save_mesh("../lid_driven_cavity/mesh.dat", [nodes])
    return

def generate_mesh_backstep(Ni=211, Nj=91):
    x = np.linspace(-11,10,Ni)
    y = np.linspace(0,9,Nj)
    nodes = np.zeros((Ni,Nj,2))
    for i in range(Ni):
      for j in range(Nj):
        nodes[i,j,0] = x[i]
        nodes[i,j,1] = y[j]
        
    nodes1 = nodes[0:111,10:]
    nodes2 = nodes[110:]
    np.savez("../backstep/mesh.npz", nodes1, nodes2)
    save_mesh("../backstep/mesh.dat", [nodes1, nodes2])
    return

# example: mesh_file = "../backstep/mesh.npz"
def split_mesh():
    # source_file = "../backstep/mesh.npz"
    # dist_file = "../backstep/mesh_4blocks.npz"
    source_file = "../cylinder/mesh.npz"
    dist_file = "../cylinder/mesh_4blocks.npz"

    n = 4
    mesh = np.load(source_file)
    nb = len(mesh.files)
    print("Block number = ", nb, "become -> ", nb*n)
    blocks = []
    for file in mesh.files:
        nodes = mesh[file]
        size_i = nodes.shape[0]
        node_id = int(size_i/n)
        for i in range(n):
            nodes_i = nodes[i*node_id : (i+1)*node_id+1]
            blocks.append(nodes_i)

    np.savez(dist_file, blocks[0], blocks[1], blocks[2], blocks[3])

    #test
    mesh = np.load(dist_file)
    for file in mesh.files:
        nodes = mesh[file]
        plt.plot(nodes[:,:,0], nodes[:,:,1], '+')
    plt.show()

if __name__ == '__main__':
    generate_mesh_cavity()
    # generate_mesh_mach2()
    # generate_mesh_cylinder()
    # generate_mesh_backstep()
    # split_mesh()

