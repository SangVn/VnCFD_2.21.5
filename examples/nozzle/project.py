# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from lib.service import Pm2P, P2Pt


P_freestream = Pm2P(M=0.2, T=293.15, p=101325.0, alf=0.0)
P_init_field = P_freestream

# áp suất toàn phần: pt
# xét 3 trường hợp p_exit/pt = 0.89, 0.75, 0.4
parm_total = P2Pt(P_freestream)
p_exit = parm_total[1]*0.75

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def set_boco():
    set_boco_const(P_freestream, outlet_pressure=p_exit, outlet_opt="fixed")
    bc_0 = [(inflow, None, None)]
    bc_1 = [(outflow, None, None)]
    bc_2 = [(symmetry, None, None)]
    bc_3 = [(symmetry, None, None)]
    blk1_bc_list = [bc_0, bc_1, bc_2, bc_3]

    boco_list = [blk1_bc_list]
    joint_list = None
    return boco_list, joint_list

boco_list, joint_list = set_boco()

# các thông số khác
solver = "eu"
# solver = "ns"
reconstr = "godunov"
#reconstr = "tvd_vector"
# reconstr = "tvd_component"

CFL = 0.85
time_target = None
iter_target = 10000
misclosure_target = None

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = None

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# export field ở dạng tecplot data (only for serial run!)
export_field_frequency_iter = None
