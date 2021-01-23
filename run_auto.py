# coding: utf-8
# Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

# !***************************************************************!
# !                 https://www.facebook.com/VnCFD                !
# !               https://vncfdgroup.wordpress.com                !
# !***************************************************************!

# Run Script in Terminal or Command-line
# python3 run_auto.py
# mpirun -n N python3 run_auto.py (N=1,2,3..)

import sys
import time
from mpi4py import  MPI
from lib.data import Blocks
from lib.service import  plot_state
from setting import  path_dir, solver

def run(opt):
    if opt[0] == "init":
        blocks = Blocks(btype=solver)
        blocks.init_field()
    elif opt[0] == "run":
        blocks = Blocks(btype=solver)
        blocks.run()
    elif opt[0] == "export":
        if MPI.COMM_WORLD.Get_size()>1:
            if MPI.COMM_WORLD.Get_rank() == 0: print("Export data: do it with only one process!")
            exit()
        filename ="field.dat"
        if len(opt) > 1: filename = opt[1]
        blocks = Blocks(btype=solver, export_plot=True)
        blocks.export_block_data(filename)
    elif opt[0] == "plot_field":
        if MPI.COMM_WORLD.Get_size()>1:
            if MPI.COMM_WORLD.Get_rank() == 0: print("plot_field: do it with only one process!")
            exit()
        blocks = Blocks(btype=solver, export_plot=True)
        plt_opt = ["p", "pcolor"]
        if len(opt) > 1: plt_opt[0] = opt[1]
        if len(opt) > 2: plt_opt[1] = opt[2]
        blocks.plot_field(field=plt_opt[0], opt=plt_opt[1], savefig=False, clbar=False)
    elif opt[0] == "plot_mesh":
        blocks = Blocks(btype=solver)
        plt_opt = "cell"
        if len(opt) > 1: plt_opt = opt[1]
        blocks.plot_mesh(opt=plt_opt, savefig=True)
    elif opt[0] == "plot_state":
        if MPI.COMM_WORLD.Get_rank() == 0:
            plot_state([path_dir+"solver.state"], legend=None, save_path=None)
    elif opt[0] == "help":
        if MPI.COMM_WORLD.Get_rank() == 0:
            print("Options: ")
            print("         - help          display this information")
            print("         - init          init field")
            print("         - run           run calculation")
            print("         - export        export field")
            print("                         - filename    (default: field.dat)")
            print("         - plot_field    post-processing")
            print("                         - field_name  (default: rho | u, v, p)")
            print("                         - plot_option (default: pcolor | contourf)")
            print("         - plot_mesh     check mesh")
            print("                         - type        (default: bound | node, cell)")
            print("         - plot_state    plot state file")
    else:
        if MPI.COMM_WORLD.Get_rank() == 0:
            print("Invalid option! Run again and try option: help")
        sys.exit()

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    if comm.Get_rank() == 0:
        if len(sys.argv) < 2:
            opt = input("Option --> ")
            opt = opt.split()
        else: opt = sys.argv[1:]
    else:
        opt = None
    if comm.Get_size()>1:
        opt = comm.bcast(opt, root=0)

    t0 = time.time()
    run(opt)
    if MPI.COMM_WORLD.Get_rank() == 0:
        print("run-time: {0:6f} s".format(time.time() - t0))
        sys.stdout.flush()