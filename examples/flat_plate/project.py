# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from lib.service import P_US2SI, reynolds_number


# https://www.grc.nasa.gov/WWW/wind/valid/fplam/fplam01/fplam01.html

# Case A
P_freestream = P_US2SI(M=0.1, T_rankine=700.0, p_psi=6.0, alf=0.0)
P_init_field = P_freestream
p_exit = P_freestream[3]

# print("P_freestream: ", P_freestream)
# print("Reynolds:", reynolds_number(P_freestream, L=0.3048))

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def set_boco():
    set_boco_const(P_freestream, outlet_pressure=p_exit, outlet_opt="fixed") #or mean
    # (..., None, None) nghĩa là bắt đầu từ side đầu tiên tới side cuối cùng trên biên
    bc_0 = [(inflow, None, None)]
    bc_1 = [(outflow, None, None)]
    bc_2 = [(symmetry, None, 10), (no_slip, 10, None)]
    bc_3 = [(farfield, None, None)]
    blk1_bc_list = [bc_0, bc_1, bc_2, bc_3]

    boco_list = [blk1_bc_list]
    joint_list = None
    return boco_list, joint_list


boco_list, joint_list = set_boco()

# các thông số khác
# solver = "eu"
solver = "ns"
reconstr = "godunov"
#reconstr = "tvd_vector"
# reconstr = "tvd_component"

CFL = 0.85
time_target = None
iter_target = 2000
misclosure_target = None

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = 5000

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# export field ở dạng tecplot data (only for serial run!)
export_field_frequency_iter = None
