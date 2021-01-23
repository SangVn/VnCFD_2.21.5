# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from lib.service import Pm2P

P_freestream = Pm2P(M=2.0, T=293.15, p=101325.0)
# P_freestream = Pm2P(M=1.8, T=293.15, p=101325.0)
P_init_field = P_freestream

# print("P_freestream: ", P_freestream)

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def set_boco():
    set_boco_const(P_freestream)
    bc_0 = [(farfield, None, None)]
    bc_1 = [(neumann, None, None)]
    bc_2 = [(no_slip, None, None)]
    bc_3 = [(symmetry, None, None)]
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

CFL = 0.45
time_target = None
iter_target = 20000
misclosure_target = -12.0

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = None

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# export field ở dạng tecplot data (only for serial run!)
export_field_frequency_iter = None
