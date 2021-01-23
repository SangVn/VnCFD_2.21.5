# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from lib.service import  P_US2SI, reynolds_number


# Bài toán lid-driven Cavity
# https://www.grc.nasa.gov/WWW/wind/valid/cavity/cavity.html

P_freestream = P_US2SI(M=0.05, T_rankine=460.0, p_psi=0.0425, alf=0.0)
P_init_field = [P_freestream[0], 0.0, 0.0, P_freestream[3]]
U_wall = [P_freestream[1], P_freestream[2]]

# print("P_freestream: ", P_freestream)
# print("Reynolds: ", reynolds_number(P_freestream, 0.0254))

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def set_boco():
    set_boco_const(P_freestream, U_wall=U_wall)
    bc_0 = [(no_slip, None, None)]
    bc_1 = [(no_slip, None, None)]
    bc_2 = [(no_slip, None, None)]
    bc_3 = [(moving_wall, None, None)]
    blk1_bc_list = [bc_0, bc_1, bc_2, bc_3]

    boco_list = [blk1_bc_list]
    joint_list = None
    return boco_list, joint_list

boco_list, joint_list = set_boco()

# các thông số khác
# solver = "eu"
solver = "ns"
# reconstr = "godunov"
#reconstr = "tvd_vector"
reconstr = "tvd_component"

CFL = 0.5
time_target = None
iter_target = 25000
misclosure_target = None

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = 10000

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# export field ở dạng tecplot data (only for serial run!)
export_field_frequency_iter = None
