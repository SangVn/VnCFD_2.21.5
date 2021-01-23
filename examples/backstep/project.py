# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from lib.service import Pm2P, PRe_US2SI, reynolds_number


# P_freestream = Pm2P(M=0.25, T=293.15, p=101325.0, alf=0.0)
#nasa case: https://turbmodels.larc.nasa.gov/backstep_val.html
P_freestream =  PRe_US2SI(M=0.128, T_rankine=537.0, Re=36000.0, L=1.0)
# print("P_freestream: ", P_freestream)
# print("Reynolds: ", reynolds_number(P_freestream, L=1.0))
P_init_field = P_freestream
p_exit = P_freestream[3]

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
#   _________________
#   |___1___|   2    |
#           |________|
def set_boco():
    set_boco_const(P_freestream, outlet_pressure=p_exit, outlet_opt="fixed") #or "mean"
    bc_0 = [(inflow, None, None)]
    bc_1 = [(joint, None, None)]
    bc_2 = [(symmetry, None, 6), (no_slip, 6, None)]
    bc_3 = [(symmetry, None, 6), (no_slip, 6, None)]
    blk1_bc_list  = [bc_0, bc_1, bc_2, bc_3]

    # bc_0 = [(symmetry, None, 10), (joint, 10, None)] # for simple mesh
    bc_0 = [(no_slip, None, 48), (joint, 48, None)] # for nasa mesh

    bc_1 = [(outflow, None, None)]
    bc_2 = [(no_slip, None, None)]
    bc_3 = [(no_slip, None, None)]
    blk2_bc_list  = [bc_0, bc_1, bc_2, bc_3]

    boco_list = [blk1_bc_list, blk2_bc_list]

    #[(left block, left bound, edge_id
    # right block, right bound, edge_id), (...)]
    joint_list = [(0,1,0, 1,0,1)]

    return boco_list, joint_list

#--------------------------------------------------------------------------#
# case: 4 blocks mesh
#   _________________________
#   |__1__|__2__|  3  |  4  |
#               |_____|_____|
# joint_list = [(0, 1, 0, 1, 0, 0), (1, 1, 0, 2, 0, 1), (2, 1, 0, 3, 0, 0)]
#--------------------------------------------------------------------------#

boco_list, joint_list = set_boco()


# các thông số khác
# solver = "eu"
solver = "ns"
# reconstr = "godunov"
#reconstr = "tvd_vector"
reconstr = "tvd_component"

# CFL = 0.45
CFL = 0.85
time_target = None
iter_target = 10000
misclosure_target = None

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = 2000

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# export field ở dạng tecplot data (only for serial run!)
export_field_frequency_iter = None
