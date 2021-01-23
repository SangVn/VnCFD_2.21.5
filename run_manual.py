# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

# !***************************************************************!
# !                 https://www.facebook.com/VnCFD                !
# !               https://vncfdgroup.wordpress.com                !
# !***************************************************************!

# Run Script in Terminal or Command-line
# python3 run_manual.py
# mpirun -n N python3 run_manual.py (N=1,2,3..)

# import sys
import time
from mpi4py import  MPI
from lib.data import Blocks
from lib.service import  plot_state
from setting import  path_dir, solver

def run():
    blocks = Blocks(btype=solver, export_plot=False)
    blocks.init_field()
    blocks.run()
    # blocks.export_block_data("eu_godunov.dat")
    blocks.plot_field(field="u", opt="pcolor", savefig=False, clbar=False)
    # blocks.plot_mesh(opt="bound", savefig=False)
    # plot_state([path_dir + "solver.state"], legend=["ns_tvd_component"], save_path=path_dir+"stat.png")

if __name__ == '__main__':
    t0 = time.time()
    run()
    if MPI.COMM_WORLD.Get_rank() == 0: print("run-time: {0:6f} s".format(time.time() - t0))