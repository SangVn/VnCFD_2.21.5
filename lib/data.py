# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

import matplotlib.pyplot as plt
from numpy import zeros, zeros_like, save, load, log10, concatenate
from .fortranlib import cell_topo, face_topo, p2u, u2p, misclosure_rho, reconstruction_vector, reconstruction_component
from .fortranlib import eu_time_step, bound_flux_roe, bound_flux_p2f, inner_flux_reconstr_godunov, inner_flux_reconstr_tvd
from .fortranlib import cell_jacobian, face_jacobian, gradient_cell, gradient_face, calc_face_p, flux_diff
from .fortranlib import ns_time_step, gradient_inner_face
from .service import Temperature, Mach
from setting import P_init_field, boco_list, joint_list, CFL, path_dir, reconstr, joint, boco_auto_split, mpi_send_recv
import setting as sett
from mpi4py import MPI
from sys import exit, stdout


class Bound:
    def __init__(self, face_normal, face_area, icell_P, icell_G, icell_half_size, gcell_P, gcell_half_size, side_P,
                 side_G, icell_dPn, icell_res, sign_ic, nd):
        self.size = len(face_area)
        self.face_normal = face_normal
        self.face_area = face_area
        self.side_P = side_P
        self.side_G = side_G

        self.icell_P = icell_P
        self.icell_G = icell_G
        self.icell_dPn = icell_dPn
        self.icell_res = icell_res
        self.icell_half_size = icell_half_size

        self.gcell_half_size = gcell_half_size
        self.jcell_P = None # joint cell
        self.jcell_G = zeros_like(side_G)

        self.sign_ic = sign_ic  # if inner_cell is left_cell: sign_ice = -1; esle: sign_ice = +1
        self.nd = nd # iface -> nd=1, jface -> nd=2

        if sign_ic == 1: 
            self.gcell_P = gcell_P
            self.side_Pl = self.gcell_P # left cell is ghost cell
            self.side_Pr = icell_P      # right cell is inner cell
        else:
            self.gcell_P = gcell_P
            self.side_Pl = icell_P      # left cell is inner cell
            self.side_Pr = self.gcell_P # right cell is ghost cell

    def split(self, s, e):
        new_bound = Bound(self.face_normal[s:e], self.face_area[s:e], self.icell_P[s:e], self.icell_G[:, s:e],
                          self.icell_half_size[s:e], self.gcell_P[s:e], self.gcell_half_size[s:e], self.side_P[s:e],
                          self.side_G[:, s:e], self.icell_dPn[s:e], self.icell_res[s:e], self.sign_ic, self.nd)
        return  new_bound

    def reconstr(self):
        if self.sign_ic==1:
            self.side_Pr = self.icell_P - self.icell_dPn
        else:
            self.side_Pl = self.icell_P + self.icell_dPn

    def flux_roe(self):
        self.icell_res += self.sign_ic * bound_flux_roe(self.face_normal, self.face_area, self.side_Pl, self.side_Pr)

    def flux_p2f(self):
        self.icell_res += self.sign_ic * bound_flux_p2f(self.face_normal, self.face_area, self.side_P)

    def joint(self, other_icell_P, other_jcell_G, other_icell_half_size):
        self.jcell_P = other_icell_P
        self.jcell_G = other_jcell_G
        self.gcell_half_size[:] = other_icell_half_size[:]
        if self.sign_ic == 1:
            self.side_Pl = self.jcell_P
        else:
            self.side_Pr = self.jcell_P

    def joint_mpi(self):
        sendbuf = self.icell_half_size.copy()
        recvbuf = mpi_send_recv(sendbuf, self.sign_ic)
        self.gcell_half_size[:] = recvbuf[:]

    def joint_gradient_face(self):
        if MPI.COMM_WORLD.Get_size()>1:
            for i in range(2):
                sendbuf = self.icell_G[i].copy()
                recvbuf = mpi_send_recv(sendbuf, self.sign_ic)
                self.jcell_G[i,:] = recvbuf[:]

        if self.sign_ic == 1: #left cell is ghost cell (joint cell)
            self.side_G[:] = gradient_inner_face(self.gcell_P, self.gcell_half_size, self.jcell_G, self.icell_P,
                                     self.icell_half_size, self.icell_G, self.nd)[:]
        else: # left cell is inner cell
            self.side_G[:] = gradient_inner_face(self.icell_P, self.icell_half_size, self.icell_G, self.gcell_P,
                                     self.gcell_half_size, self.jcell_G, self.nd)[:]



class euBlock:
    def __init__(self, id, node_X, boco):
        self.id = id
        self.file_field = path_dir+"field"+str(id)+".npy"
        self.size = (node_X.shape[0]-1, node_X.shape[1]-1) # size of inner cells

        self.cell_number = self.size[0]*self.size[1]
        self.cell_X, self.all_cell_half_size, self.cell_volume = cell_topo(node_X)
        self.cell_half_size = self.all_cell_half_size[1:-1, 1:-1]
        self.cell_J, self.cell_J_inv = cell_jacobian(node_X)
        # all_cell = inner_cell + ghost_cell
        self.all_cell_P = zeros((self.size[0]+2, self.size[1]+2, 4), order='F')
        self.cell_P = self.all_cell_P[1:-1, 1:-1]
        self.cell_U = zeros((self.size[0], self.size[1], 4), order='F')
        self.cell_U_prev = zeros((self.size[0], self.size[1], 4), order='F')
        self.cell_res = zeros((self.size[0], self.size[1], 4), order='F')
        self.cell_dP = zeros((2, self.size[0], self.size[1], 4), order='F')
        self.cell_G = zeros((2, self.size[0], self.size[1], 4), order='F')

        self.iface_X, self.iface_normal, self.iface_area, self.jface_X, self.jface_normal, self.jface_area = face_topo(node_X)
        self.iface_P = zeros((self.size[0] + 1, self.size[1], 4), order='F')
        self.iface_G = zeros((2, self.size[0] + 1, self.size[1], 4), order='F')
        self.jface_P = zeros((self.size[0], self.size[1] + 1, 4), order='F')
        self.jface_G = zeros((2, self.size[0], self.size[1] + 1, 4), order='F')

        bound0 = Bound(self.iface_normal[0], self.iface_area[0], self.cell_P[0], self.cell_G[:, 0],
                       self.all_cell_half_size[1, 1:-1], self.all_cell_P[0, 1:-1], self.all_cell_half_size[0, 1:-1],
                       self.iface_P[0], self.iface_G[:, 0], self.cell_dP[0, 0], self.cell_res[0], 1, 1)
        bound1 = Bound(self.iface_normal[-1], self.iface_area[-1], self.cell_P[-1], self.cell_G[:, -1],
                       self.all_cell_half_size[-2, 1:-1], self.all_cell_P[-1, 1:-1], self.all_cell_half_size[-1, 1:-1],
                       self.iface_P[-1], self.iface_G[:, -1], self.cell_dP[0, -1], self.cell_res[-1], -1, 1)
        bound2 = Bound(self.jface_normal[:, 0], self.jface_area[:, 0], self.cell_P[:, 0], self.cell_G[:, :, 0],
                       self.all_cell_half_size[1:-1, 1], self.all_cell_P[1:-1, 0], self.all_cell_half_size[1:-1, 0],
                       self.jface_P[:, 0], self.jface_G[:, :, 0], self.cell_dP[1, :, 0], self.cell_res[:, 0], 1, 2)
        bound3 = Bound(self.jface_normal[:, -1], self.jface_area[:, -1], self.cell_P[:, -1], self.cell_G[:, :, -1],
                       self.all_cell_half_size[1:-1, -2], self.all_cell_P[1:-1, -1], self.all_cell_half_size[1:-1, -1],
                       self.jface_P[:, -1], self.jface_G[:, :, -1], self.cell_dP[1, :, -1], self.cell_res[:, -1], -1, 2)
        bounds = [bound0, bound1, bound2, bound3]


        # boco_list = (boco_0, boco_1, boco_2, boco_3)
        # example: boco_i = [(no_slip, None, 10), (joint, 10, None)]
        self.boco_list = boco
        self.bounds = []
        for i,boco_i in enumerate(self.boco_list):
            bound_list = []
            if len(boco_i)==1:
                bound_list.append(bounds[i])
            else: # >1
                for bc in boco_i:
                    new_bound = bounds[i].split(bc[1], bc[2])
                    bound_list.append(new_bound)
            self.bounds.append(bound_list)

        comm = MPI.COMM_WORLD
        if comm.Get_size()>1:
            rank = comm.Get_rank()
            if rank==0 and self.id==1:
                self.boco_list[0], self.boco_list[1] = self.boco_list[1], self.boco_list[0]
                self.bounds[0], self.bounds[1] = self.bounds[1], self.bounds[0]

    def set_joint_mpi(self):
        for i, bound in enumerate(self.bounds):
            for j, boco in enumerate(self.boco_list[i]):
                if boco[0] is joint: bound[j].joint_mpi()

    def write_field(self):
        save(self.file_field, self.cell_P)
        
    def read_field(self, nf=1):
        cell_P = load(self.file_field)
        for i in range(1,nf):
            if (self.cell_P.shape == cell_P.shape): break
            next_file_field = path_dir+"field"+str(self.id+i)+".npy"
            next_cell_P = load(next_file_field)
            cell_P = concatenate((cell_P, next_cell_P))  
              
        if (self.cell_P.shape != cell_P.shape):
            print("Error: block%d.cell_P.shape != field%d.shape" % (self.id, self.id))
            print("Did you want to export or plot data? Please check 'export_plot' option!")
            exit()
            
        self.cell_P[:] = cell_P[:]
        p2u(self.cell_P, self.cell_U)

    def init_field(self, P0):
        self.cell_P[:, :] = P0
        p2u(self.cell_P, self.cell_U)
        self.write_field()

    def set_boco(self):
        for i, bound in enumerate(self.bounds):
            for j, boco in enumerate(self.boco_list[i]): boco[0](bound[j], 1)

    def reconstruction(self):
        if reconstr=="godunov": return
        elif reconstr=="tvd_vector":
            reconstruction_vector(self.all_cell_P, self.all_cell_half_size, self.cell_dP)
        else: #tvd_component
            reconstruction_component(self.all_cell_P, self.all_cell_half_size, self.cell_J, self.cell_J_inv, self.cell_dP)

    def flux_convection(self):
        if reconstr=="godunov":
            inner_flux_reconstr_godunov(self.iface_normal, self.iface_area, self.jface_normal, self.jface_area,
                                        self.cell_P, self.cell_res)
        else: # reconstr=="tvd"
            inner_flux_reconstr_tvd(self.iface_normal, self.iface_area, self.jface_normal, self.jface_area, self.cell_P,
                                    self.cell_dP, self.cell_res)
        # boundary_flux
        for i, bound in enumerate(self.bounds):
            for j, boco in enumerate(self.boco_list[i]): boco[0](bound[j], 2)

    def new_U(self, dt):
        self.cell_U += self.cell_res*dt/self.cell_volume[:, :, None]
        self.cell_res[:] = 0.0

    def new_P(self):
        u2p(self.cell_U, self.all_cell_P)

    def calculate_misclosure(self):
        return misclosure_rho(self.cell_U_prev[:,:,0], self.cell_P[:,:,0])

    def time_step_global(self):
        dt = eu_time_step(self.cell_P, self.cell_half_size)
        return dt

    def iteration(self, dt, rkn=0):
        if rkn==0: self.cell_U_prev[:] = self.cell_U[:]
        self.set_boco()
        self.reconstruction()
        self.flux_convection()
        self.new_U(dt)
        if rkn==1:
            self.cell_U = 0.5*(self.cell_U + self.cell_U_prev)
        self.new_P()



class nsBlock(euBlock):
    def __init__(self, id, node_X, boco):
        super().__init__(id, node_X, boco)
        self.iface_J, self.jface_J = face_jacobian(node_X, self.cell_X)

    def gradient_cell(self):
        gradient_cell(self.all_cell_P, self.all_cell_half_size, self.cell_J, self.cell_J_inv, self.cell_G)

    def gradient_face(self):
        gradient_face(self.all_cell_P, self.all_cell_half_size, self.cell_G, self.iface_G, self.jface_G)
        # boundary_gradient
        for i, bound in enumerate(self.bounds):
            for j, boco in enumerate(self.boco_list[i]): boco[0](bound[j], 3)

    def calc_face_P(self):
        calc_face_p(self.cell_P, self.cell_half_size, self.iface_P, self.jface_P)

    def time_step_global(self):
        dt = ns_time_step(self.cell_P, self.cell_half_size)
        return dt

    def flux_diffusion(self):
        flux_diff(self.iface_normal, self.iface_area, self.iface_J, self.iface_G, self.iface_P,\
                self.jface_normal, self.jface_area, self.jface_J, self.jface_G, self.jface_P, self.cell_res)

    def iteration(self, dt, rkn=0):
        if rkn==0: self.cell_U_prev[:] = self.cell_U[:]
        self.set_boco()
        self.reconstruction()
        self.flux_convection()
        self.gradient_cell()
        self.gradient_face()
        self.calc_face_P()
        self.flux_diffusion()
        self.new_U(dt)
        if rkn==1:
            # new_u_rkn1(self.cell_U, self.cell_U_prev)
            self.cell_U = 0.5*(self.cell_U + self.cell_U_prev)
        self.new_P()



class Blocks():
    def __init__(self, btype="eu", mesh_file=None, boco=boco_list, joint=joint_list, export_plot=False):
        self.state_file = path_dir+'solver.state'
        if mesh_file is None: self.mesh_file = path_dir+"mesh.npz"
        else: self.mesh_file = mesh_file
        self.time_step_global = 0.0
        self.boco_list = boco
        self.joint_list = joint
        self.runge_kutta = 1
        if reconstr[:3]=="tvd": self.runge_kutta = 2
        self.node_X = []
        self.blocks = []
        self.cell_number = 0.0
        self.num_files_per_block = 1

        comm = MPI.COMM_WORLD
        mpi_size = comm.Get_size()
        rank = comm.Get_rank()
        if mpi_size==1:
            mesh = load(self.mesh_file)
            self.len = len(mesh.files)
            sn = 1
            if export_plot==True:
                with open(self.state_file, 'r') as fs:
                    first_line = fs.readline().split()
                    mpirun_n = int(first_line[1])
                    sn = int(mpirun_n/self.len)
                    self.num_files_per_block = sn+(int(first_line[1])%2)

            for i in range(self.len):
                node_X = mesh[mesh.files[i]]
                self.node_X.append(node_X)
                if btype=="eu":
                    block = euBlock(i*sn+1, node_X, boco[i])
                else: #btype=="ns"
                    block = nsBlock(i*sn+1, node_X, boco[i])
                self.blocks.append(block)
                self.cell_number += block.cell_number
            self.set_joint(joint)
        else: # size > 1, one process - one block
            node_X = None
            blocks_size = None
            new_boco = boco
            if rank==0:
                mesh = load(self.mesh_file)
                node_X = [mesh[name] for name in mesh.files]
                blocks_size = len(mesh.files)
                if (mpi_size > blocks_size):
                    node_X, new_boco = self.mpi_auto_split(node_X, boco, mpi_size, blocks_size)

            node_X = comm.scatter(node_X, root=0)
            blocks_size = comm.bcast(blocks_size, root=0)
            if (mpi_size > blocks_size):
                new_boco = comm.bcast(new_boco, root=0)

            self.len = 1
            self.node_X.append(node_X)
            if btype=="eu":
                block = euBlock(rank+1, node_X, new_boco[rank])
            else:# btype=="ns"
                block = nsBlock(rank+1, node_X, new_boco[rank])
            self.blocks.append(block)
            self.cell_number += block.cell_number
            self.set_joint(joint)


    def __getitem__(self, item):
        '''Lấy phần tử của dãy Blocks.'''
        return self.blocks[item]

    def mpi_auto_split(self, node_X_list, boco_list, mpi_size, blocks_size):
        if mpi_size<=blocks_size: return
        print("Warning: auto_split was called, split_size = %d, blocks_size = %d\n" %(mpi_size, blocks_size))
        N = int(mpi_size / blocks_size)
        new_node_X_list = []
        new_boco_list = []
        for bi,node_X in enumerate(node_X_list):
            if(bi==blocks_size-1): N = mpi_size-len(new_node_X_list)
            if(N==1):
                new_node_X_list.append(node_X)
                new_boco_list.append(boco_list[bi])
            else:
                new_size = int(node_X.shape[0]/N)
                for i in range(N):
                    id_start = i * new_size
                    id_end = (i + 1) * new_size
                    if i!=N-1: new_node_X = node_X[id_start: id_end+1]
                    else: new_node_X = node_X[id_start:]
                    new_node_X_list.append(new_node_X)

                    new_boco_i = []
                    if(i==0): new_bc_0 = boco_list[bi][0]
                    else: new_bc_0 = [(joint, None, None)]
                    new_boco_i.append(new_bc_0)
                    if(i==N-1): new_bc_1 = boco_list[bi][1]
                    else: new_bc_1 = [(joint, None, None)]
                    new_boco_i.append(new_bc_1)

                    for j in range(2,4):
                        bc_j = boco_list[bi][j]
                        if len(bc_j)==1: #[(boco, None, None)]
                            new_boco_i.append(bc_j)
                        else:
                            temp = boco_auto_split(bc_j, id_start, id_end)
                            new_boco_i.append(temp)
                    new_boco_list.append(new_boco_i)
        return new_node_X_list, new_boco_list

    def set_joint(self, joint_list):
        if MPI.COMM_WORLD.Get_size() > 1:
            for block in self.blocks:
                block.set_joint_mpi()
            return
        if joint_list is None: return
        for joint in joint_list:
            bound1 = self.blocks[joint[0]].bounds[joint[1]][joint[2]]
            bound2 = self.blocks[joint[3]].bounds[joint[4]][joint[5]]
            bound1.joint(bound2.icell_P, bound2.icell_G, bound2.icell_half_size)
            bound2.joint(bound1.icell_P, bound1.icell_G, bound1.icell_half_size)

    def read_field(self):
        '''Đọc trường khí động.'''
        for block in self.blocks: block.read_field(self.num_files_per_block)

    def write_field(self):
        '''Ghi trường khí động.'''
        for block in self.blocks: block.write_field()

    def init_field(self, P0=P_init_field):
        '''Thiết lập điều kiện ban đầu, ghi trường khí động.'''
        for block in self.blocks:
            block.init_field(P0)
            block.write_field()
        # self.blocks[0].init_field(P0)
        # self.blocks[1].init_field([P0[0], 0.0, 0.0, P0[3]])
        if MPI.COMM_WORLD.Get_rank()==0:
            with open(self.state_file, 'w') as f:
                mpi_size = MPI.COMM_WORLD.Get_size()
                f.write("mpirun: %d processors\n" % mpi_size)
                f.write('iter  time  misclosure_rho:\n')
                f.write("{0:6d} {1:6f} {2:4.2f}\n".format(0, 0.0, 0.0))

    def set_time_step(self, cfl=CFL):
        '''Xác định bước thời gian trong toàn bộ vùng tính.'''
        self.time_step_global = 1.0e6
        for block in self.blocks:
            dt = block.time_step_global()
            self.time_step_global = min(self.time_step_global, dt)
        self.time_step_global *= cfl
        return self.time_step_global

    def iteration(self,rkn=0):
        '''Thực hiện bước lặp thời gian.'''
        for block in self.blocks: block.iteration(self.time_step_global,rkn)

    def calculate_misclosure(self):
        sum_dr = 0.0
        for block in self.blocks:
            sum_dr += block.calculate_misclosure()
        return sum_dr

    def run(self):
        comm = MPI.COMM_WORLD
        mpi_size = comm.Get_size()
        rank = comm.Get_rank()

        # thiết lập các thông số:
        time_target, iter_target, write_field_frequency_time, write_field_frequency_iter, print_frequency_iter = \
        sett.time_target, sett.iter_target, sett.write_field_frequency_time, sett.write_field_frequency_iter, sett.print_frequency_iter
        export_field_frequency_iter, misclosure_target = sett.export_field_frequency_iter, sett.misclosure_target

        if time_target is None: time_target = 1.0e10
        if iter_target is None: iter_target = 1e10
        if write_field_frequency_time is None: write_field_frequency_time = 1.0e10
        if write_field_frequency_iter is None: write_field_frequency_iter = 1e10
        if print_frequency_iter is None: print_frequency_iter = 10
        if misclosure_target is None: misclosure_target = -16.0
        # đọc iter và time từ file state

        with open(self.state_file, 'r') as f:
            for line in f: pass
            line = line.split()
            loop_iter = int(line[0])
            loop_time = float(line[1])
            mscl_rho = float(line[2])

        stop = (loop_time >= time_target or loop_iter >= iter_target or mscl_rho<=misclosure_target)
        if (mscl_rho <= misclosure_target): stop = True

        if rank==0:
            print("iteration, dt, time, misclosure_rho")
            print("{0:6d}, {1:8.3e}, {2:6f}, {3:4.2f}".format(loop_iter, 0.0, loop_time, mscl_rho))
            stdout.flush()

        sum_cell_number = self.cell_number

        if mpi_size>1: sum_cell_number = comm.reduce(sum_cell_number, op=MPI.SUM, root=0)

        self.read_field()
        while (stop is False):
            # runge kutta loop
            for rkn in range(self.runge_kutta):
                dt = self.set_time_step()
                if mpi_size > 1:
                    dt = comm.allreduce(dt, MPI.MIN)
                    self.time_step_global = dt
                self.iteration(rkn)

            loop_iter += 1
            loop_time += dt

            if (loop_time >= time_target or loop_iter >= iter_target): stop = True

            # hiển thị các thông số cơ bản lên terminal
            if (not(loop_iter%print_frequency_iter) or (stop == True)):
                mscl_rho = self.calculate_misclosure()
                if mpi_size > 1:
                    mscl_rho = comm.reduce(mscl_rho, op=MPI.SUM, root=0)
                if rank == 0:
                    mscl_rho = log10(mscl_rho / sum_cell_number)
                    if (mscl_rho <= misclosure_target): stop = True
                    print("{0:6d}, {1:8.3e}, {2:6f}, {3:4.2f}".format(loop_iter, dt, loop_time, mscl_rho))
                    stdout.flush()
                    with open(self.state_file, 'a') as f:
                        f.write("{0:6d} {1:6f} {2:4.2f}\n".format(loop_iter, loop_time, mscl_rho))
                stop = comm.bcast(stop, root=0)

            # ghi lại kết quả giữa chừng
            period_i = int(loop_iter / write_field_frequency_iter)
            period_t = int(loop_time / write_field_frequency_time)
            if (stop == True) or (loop_iter == period_i * write_field_frequency_iter)\
                or (loop_time >= period_t * write_field_frequency_time > loop_time - dt):
                self.write_field()
                if rank == 0:
                    print('\nwrite_field at iteration: %d, time: %f\n' % (loop_iter, loop_time))
                    stdout.flush()

            # xuất dữ liệu giữa chừng ở định dạng tecplot data
            if (mpi_size == 1) and (export_field_frequency_iter != None) and not(loop_iter%export_field_frequency_iter):
                self.export_block_data(str(loop_iter)+"field.dat")

    def export_block_data(self, filename=None):
        if MPI.COMM_WORLD.Get_size()>1:
            if MPI.COMM_WORLD.Get_rank() == 0: print("export_block_data: do it with only one process!")
            exit()

        self.read_field()

        '''Xuất trường khí động vào file định dạng Tecplot block data.'''
        if filename is None: filename ="field.dat"
        filename = path_dir + filename
        print('Write block data to: %s' % filename)
        f = open(filename, 'w')
        f.write('TITLE = "vncfd field"\n')
        f.write('VARIABLES = "X", "Y", "rho", "u", "v", "p", "Mach", "T"\n')

        n = 0
        for block in self.blocks:
            nodes = self.node_X[n]
            n += 1
            f.write('ZONE T="%d", I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)\n' % (
            n, nodes.shape[0], nodes.shape[1]))

            X_p, Y_p = nodes[:, :, 0].T.ravel(), nodes[:, :, 1].T.ravel()
            for x in X_p: f.write('%.8f ' % x)
            f.write('\n')
            for y in Y_p: f.write('%.8f ' % y)
            f.write('\n')

            cells_P = block.cell_P.transpose(1,0,2).reshape(-1, 4)
            for i in range(4):
                for cell in cells_P: f.write('%.8f ' % cell[i])
                f.write('\n')
            for cell in cells_P:
                M = Mach(cell)
                f.write('%.8f ' % M)
            for cell in cells_P:
                T = Temperature(cell)
                f.write('%.8f ' % T)
            f.write('\n')
        f.close()
        
    def get_scale(self, axis_limit):
        x_min, x_max = 1e9, -1e9
        y_min, y_max = 1e9, -1e9
        for node in self.node_X:
            x_min = min(x_min, node[:, :, 0].min())
            x_max = max(x_max, node[:, :, 0].max())
            y_min = min(y_min, node[:, :, 1].min())
            y_max = max(y_max, node[:, :, 1].max())
        r = [x_min, x_max, y_min, y_max]
        for i in range(len(axis_limit)):
            if axis_limit[i]!=None: r[i] = axis_limit[i]
        s = (r[1]-r[0])/(r[3]-r[2])
        return s

    # cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'spring', 'summer', 'autumn', 'winter', 'cool',...]
    def plot_field(self, field='rho', opt='pcolor', axis_limit=[None, None, None, None],
                   cmap='viridis',clbar=False, savefig=False):
        self.read_field()
        fid = {'rho': 0, 'u': 1, 'v': 2, 'p': 3}
        fname = {'rho': 'Density', 'u': 'Velocity_X', 'v':'Velocity_Y', 'p':'Pressure'}
        s = self.get_scale(axis_limit)
        fig, ax = plt.subplots(figsize=(6*s,6))
        for n,block in enumerate(self.blocks):
            if opt == 'pcolor':
                cs = ax.pcolor(self.node_X[n][:,:,0], self.node_X[n][:,:,1], block.cell_P[:,:,fid[field]], cmap=cmap)
            else:
                cs = plt.contourf(block.cell_X[:,:,0], block.cell_X[:,:,1], block.cell_P[:,:,fid[field]], cmap=cmap)
        if clbar is True: fig.colorbar(cs)
        plt.title("VnCFD: "+fname[field])

        if savefig is True:
            iname = path_dir + field + ".png"
            if MPI.COMM_WORLD.Get_size() > 1:
                iname = path_dir + "r" + str(MPI.COMM_WORLD.Get_rank()) +"_"+ field+".png"
            fig.savefig(iname, bbox_inches="tight")
        plt.show()

    def plot_mesh(self, opt="bound", axis_limit=[None, None, None, None], savefig=False):
        s = self.get_scale(axis_limit)
        fig = plt.figure(figsize=(6*s,6))
        for n,block in enumerate(self.blocks):
            if opt=="node": plt.plot(self.node_X[n][:, :, 0], self.node_X[n][:, :, 1], "k+")
            if opt=="cell": plt.plot(block.cell_X[:, :, 0], block.cell_X[:, :, 1], "k+")
            plt.plot(block.iface_X[0, :, 0], block.iface_X[0, :, 1], "r-")
            plt.plot(block.iface_X[-1, :, 0], block.iface_X[-1, :, 1], "g-")
            plt.plot(block.jface_X[:, 0, 0], block.jface_X[:, 0, 1], "m-")
            plt.plot(block.jface_X[:, -1, 0], block.jface_X[:, -1, 1], "b-")
        plt.title("VnCFD: bound(0, 1, 2, 3) ~ (red, green, magenta, blue)")
        if savefig is True:
            iname = path_dir + "mesh.png"
            if MPI.COMM_WORLD.Get_size()>1:
                iname = path_dir + "r"+ str(MPI.COMM_WORLD.Get_rank()) + "_mesh.png"
            plt.savefig(iname, bbox_inches="tight")
        plt.show()
