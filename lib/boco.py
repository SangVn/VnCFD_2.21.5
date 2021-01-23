# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from .fortranlib import p2f, bc_symmetry, bc_outflow, bc_inflow
from numpy import empty_like
from mpi4py import MPI

P_freestream = None
U_moving_wall = None
pressure_exit = None
outflow_opt = None

def set_boco_const(P_fstream, U_wall=None, outlet_pressure=None, outlet_opt="fixed"): #or "mean"
    global P_freestream, U_moving_wall, pressure_exit, outflow_opt
    P_freestream = P_fstream
    U_moving_wall = U_wall
    if outlet_pressure is None: pressure_exit = P_freestream[3]
    else: pressure_exit = outlet_pressure
    if outlet_opt=="fixed": outflow_opt = 1
    else: outflow_opt = 2 # "mean"

# boco_func(bound, opt):
# opt==1: set side_P, gcell_P
# opt==2: reconstruction, flux calc.
# opt==3: bound gradient calc.
#

def null(bound,opt):
    pass

# ghost_cell_P = P_freestream
def farfield(bound, opt):
    if opt==1:
        bound.side_P[:] = 0.5*(bound.icell_P[:]+P_freestream)
        bound.gcell_P[:] = P_freestream
    elif opt==2:
        bound.reconstr()
        bound.flux_roe()
    else:
        pass

# ghost_cell_P = inner_cell_P
def neumann(bound, opt):
    if opt==1:
        bound.side_P[:] = bound.icell_P[:]
        bound.gcell_P[:] = bound.side_P[:]
    elif opt==2:
        bound.flux_p2f()
    else:
        pass

def symmetry(bound, opt):
    if opt==1:
        bound.side_P[:] = bc_symmetry(bound.face_normal, bound.icell_P)[:]
        bound.gcell_P[:] = bound.side_P[:]
    elif opt==2:
        bound.flux_p2f()
    else:
        pass

def no_slip(bound, opt):
    if opt==1:
        bound.side_P[:, (0,3)] = bound.icell_P[:, (0,3)]
        bound.side_P[:, (1,2)] = 0.0
        bound.gcell_P[:] = bound.side_P[:]
    elif opt==2:
        bound.flux_p2f()
    else:
        bound.side_G[bound.nd%2,:, (1,2)] = 0.0

def moving_wall(bound, opt):
    if opt==1:
        bound.side_P[:, (0,3)] = bound.icell_P[:, (0,3)]
        bound.side_P[:, (1,2)] = U_moving_wall
        bound.gcell_P[:] = bound.side_P[:]
    elif opt==2:
        bound.flux_p2f()
    else:
        bound.side_G[bound.nd%2,:, (1,2)] = 0.0

def inflow(bound, opt):
    if opt==1:
        bound.side_P[:] = bc_inflow(bound.face_normal, bound.icell_P, P_freestream)[:]
        # bound.gcell_P[:] = bound.side_P[:]
        bound.gcell_P[:] = P_freestream
    elif opt==2:
        bound.flux_p2f()
        # bound.flux_roe()
    else:
        pass

def outflow(bound, opt):
    if opt==1:
        bound.side_P[:] = bc_outflow(bound.face_normal, bound.icell_P, pressure_exit, outflow_opt)
        bound.gcell_P[:] = bound.side_P[:]
    elif opt==2:
        bound.flux_p2f()
    else:
        pass

# for mpich, Microsoft MPI
# try it with openmpi latest version
def mpi_send_recv(sendbuf, bound_sign_ic):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    if size > 1:
        rank = comm.Get_rank()
        other_rank = rank - bound_sign_ic
        if other_rank == -1: other_rank = size - 1
        elif other_rank == size: other_rank = 0
        recvbuf = empty_like(sendbuf)

        comm.Send(sendbuf, dest=other_rank, tag=77)
        comm.Recv(recvbuf, source=other_rank, tag=77)

        return recvbuf

# If you use openmpi, maybe you have a deadlock problem with big data,
# this function splits it to sub-data and send-recv them.
# To use it, change it's name to mpi_send_recv() and
# change above function's name to other
def mpi_send_recv_sub(sendbuf, bound_sign_ic):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    if size > 1:
        rank = comm.Get_rank()
        other_rank = rank - bound_sign_ic
        if other_rank == -1: other_rank = size - 1
        elif other_rank == size: other_rank = 0
        recvbuf = empty_like(sendbuf)

        lenbuf = sendbuf.shape[0]
        if lenbuf>125:
            s, e = 0, 125
            while True:
                comm.Send(sendbuf[s:e], dest=other_rank, tag=77)
                comm.Recv(recvbuf[s:e], source=other_rank, tag=77)
                if e is None: break
                s = e
                e = e+125
                if e>=lenbuf: e = None
        else:
            comm.Send(sendbuf, dest=other_rank, tag=77)
            comm.Recv(recvbuf, source=other_rank, tag=77)

        return recvbuf


def joint(bound, opt):
    if opt==1:
        if MPI.COMM_WORLD.Get_size()==1:
            bound.gcell_P[:] = bound.jcell_P[:]
        else:
            sendbuf = bound.icell_P.copy()
            recvbuf = mpi_send_recv(sendbuf, bound.sign_ic)
            bound.gcell_P[:] = recvbuf[:]
        bound.side_P[:] = 0.5*(bound.icell_P[:]+bound.gcell_P[:])
        # dg = bound.gcell_half_size[:, bound.nd-1]
        # di = bound.icell_half_size[:, bound.nd-1]
        # ds = dg+di
        # bound.side_P[:] = ((bound.icell_P*dg[:,None]+bound.gcell_P*di[:,None])/ds[:,None])[:]
    elif opt==2:
        bound.reconstr()
        if MPI.COMM_WORLD.Get_size()>1:
            if bound.sign_ic==1: sendbuf = bound.side_Pr.copy()
            else: sendbuf = bound.side_Pl.copy()
            recvbuf = mpi_send_recv(sendbuf, bound.sign_ic)
            if bound.sign_ic==1: bound.side_Pl[:] = recvbuf[:]
            else: bound.side_Pr[:] = recvbuf[:]
        bound.flux_roe()
    else:
        pass
        # bound.joint_gradient_face()


def boco_auto_split(bc_j, id_start, id_end):
    new_boco_i = []
    for bc_k in bc_j:
        if bc_k[1] == None:
            if bc_k[2] == None or bc_k[2] > id_end:
                new_boco_i.append((bc_k[0], None, None))
            elif bc_k[2] < id_start:
                pass
            elif bc_k[2] < id_end:
                new_boco_i.append((bc_k[0], None, bc_k[2] - id_start))
            else:
                print("case ???")
        else:
            if bc_k[1] > id_end:
                pass
            elif bc_k[1] < id_start:
                if bc_k[2] == None or bc_k[2] > id_end:
                    new_boco_i.append((bc_k[0], None, None))
                else:
                    new_boco_i.append((bc_k[0], None, bc_k[2] - id_start))
                pass
            else:
                if bc_k[2] == None or bc_k[2] > id_end:
                    new_boco_i.append((bc_k[0], bc_k[1] - id_start, None))
                else:
                    new_boco_i.append((bc_k[0], bc_k[1] - id_start, bc_k[2] - id_start))

    return new_boco_i
