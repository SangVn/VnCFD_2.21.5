# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

import numpy as np
from mesh_generator import save_mesh
import matplotlib.pyplot as plt

def convert_nasa_mesh(infile, outfile1, outfile2, opt=0):
	f = open(infile, 'r')
	size = f.readline()
	data = f.readline()
	f.close()
	size = size.split()
	data = data.split()
	print(size)

	size = np.array([int(s) for s in size])
	data = np.array([float(d) for d in data])

	N = size[0]*size[1]
	x = data[0:N].reshape(size[1], size[0])
	y = data[N:2*N].reshape(size[1], size[0])
	nodes = np.zeros((size[1], size[0], 2))

	r = 1.0
	if opt==1: r = 0.0254 #inch -> m
	if opt==2: r = 0.3048 #ft -> m

	nodes[:, :, 0] = x*r
	nodes[:, :, 1] = y*r

	save_mesh(outfile1, [nodes])
	
	#chuyển ma trận về dạng chuẩn: nodes(nx, ny)
	nodes = nodes.transpose((1,0,2))
	print("Save mesh ", outfile2)
	np.savez(outfile2, nodes)
	return

# VnCFD_2.21.5 uses only structured mesh. If you have a CGNS mesh,
#     * open it with paraview, press "Apply",
#     * click "Information" to know mesh size (dimension)
#     * (if it's 3D mesh, make a "Slice" to get 2D mesh),
#     * select "File" -> "Save Data" -> type "File name" (e.g. "mesh.txt") -> Enter
#     * to convert it to VnCFD_2.21.5 mesh format, modify and run "convert_cgns_mesh_()" function

# for one block mesh
def convert_cgns_mesh_1():
	path = "naca0012/mesh0.txt" # mesh path
	size = (47, 103) # mesh size
	z0 = np.loadtxt(path, skiprows=1, delimiter=",", usecols=(0, 1)) # (0, 1) ~ (x, y)
	z0 = z0.reshape(size[0],size[1],2).transpose(1,0,2)
	# np.savez("../naca0012/mesh.npz", z0)
	# save_mesh("../naca0012/mesh.dat", [z0])
	# test
	plt.plot(z0[:, :, 0], z0[:, :, 1], "g+")
	plt.show()


# https://turbmodels.larc.nasa.gov/backstep_grids.html
def convert_cgns_mesh_2():
	z0 = np.loadtxt("nasa_backstep/zone0.txt", skiprows=1, delimiter=",", usecols=(0, 2)) # (0, 2) ~ (x, y)
	z1 = np.loadtxt("nasa_backstep/zone1.txt", skiprows=1, delimiter=",", usecols=(0, 2))
	z2 = np.loadtxt("nasa_backstep/zone2.txt", skiprows=1, delimiter=",", usecols=(0, 2))
	z3 = np.loadtxt("nasa_backstep/zone3.txt", skiprows=1, delimiter=",", usecols=(0, 2))

	z0 = z0.reshape(65,65,2).transpose(1,0,2)
	z1 = z1.reshape(65,25,2).transpose(1,0,2)
	z2 = z2.reshape(113,97,2).transpose(1,0,2)
	z3 = z3.reshape(113,33,2).transpose(1,0,2)

	# joint 2 zones in 1, except the first column of the second zone
	z01 = np.concatenate((z0, z1[1:]))
	z23 = np.concatenate((z2, z3[1:]))

	np.savez("../backstep/mesh.npz", z01, z23)
	save_mesh("../backstep/mesh.dat", [z01, z23])

	# test
	plt.plot(z01[:, :, 0], z01[:, :, 1], "g+")
	plt.plot(z23[:, :, 0], z23[:, :, 1], "b+")
	plt.show()

if __name__ == "__main__":
	convert_cgns_mesh_1()
	# convert_cgns_mesh_2()
	# path = "../nozzle/"
	# convert_nasa_mesh("nasa_cdnozzle/cdnoz.txt", path+"mesh.dat", path+"mesh.npz", opt=1)
	# path = "../flat_plate/"
	# convert_nasa_mesh("nasa_plate/fplam.txt", path+"mesh.dat", path+"mesh.npz", opt=2)

