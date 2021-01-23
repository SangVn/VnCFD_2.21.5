Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

### About VnCFD_2.21.5
   VnCFD_2.21.5 is a free, open source computational fluid dynamics (CFD) software 
   released for Educational Purposes. It's written in a combination of Python, Fortran and MPI.
   The core of calculation, written in Fortran, and the MPI help to speed up the simulations, 
   while Python is comfortable to control Classes and Data. The VnCFD_2.21.5 Code, which is short, clear and 
   easy to read, helps you save a lot of time spent on learning the basic of CFD coding.

### Installation Requirements
   - Python3, numpy, matplotlib, mpi4py, scipy
   - gcc, g++, gfortran (Linux) or MinGW-w64 (Windows)
   - mpich (Linux) or Microsoft MPI (Windows)
   - paraview

**NOTE: READ The Installation Guide to Install All The Requirements**

### Compiling
1. In the "lib" directory, open Terminal (Command-line, PowerShell),
   compile the "fortranlib" by one of these commands:
   - ./compile_fortranlib
   - python -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib
   - (use "python3" instead of "python" in Linux)
   
### Running Examples

**Note: Read more information in Tutorial Guide**

1. Open file "setting.py"
   - Choose any example and uncomment 2 lines *from ..., path_dir ...*, then save the file

2. Open file "project.py" in the chosen example directory (e.g. *examples/name/project.py*)
   - Check or set new parameters, save the file

3. To start computation, run script "run_auto.py" by one of these commands:
   - python run_auto.py option
   - mpiexec -n N python run_auto.py option
   - mpirun -n N python run_auto.py option  (N=1,2,3..)
   
4. Options
   
     Task                | Option
   ----------------------|------------
   - to view help        | help
   - to init field       | init
   - to run calculation  | run
   - to export data      | export
   - to plot field       | plot_field
   - to plot mesh        | plot_mesh
   - to plot convergence | plot_state

*Note: Press "ctrl+C" to interrupt the calculation*

5. Post-processing
   - Use option **plot_...**
   - Use **run_manual.py** to have more flexibility
   - Recommend: use **paraview** to view **field.dat** in the chosen example directory
   
<center>VnCFD_2.21.5 = Python + MPI + Fortran</center>
<img src="doc/VnCFD_2.21.5.png">
