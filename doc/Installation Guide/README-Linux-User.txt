NOTE. I'm so glad that you are linux-user like me!

1. Installation Requirements
	- Check version and install gcc, g++, gfortran if it's necessary
		+ sudo apt-get install gfortran
		+ sudo apt-get install gcc
		+ sudo apt-get install g++

    - Install mpich
		+ sudo apt-get install mpich
        + (if you use "openmpi", maybe there is deadlock problem with big data, then
          read the solution - "mpi_send_recv_sub() function" in "boco.py" module)

    - Open Terminal, check python and pip version:
		+ python -V
		+ pip -V

	- Install Python 3 if Python version is 2
	    + sudo apt-get install python3

	- Install pip if it does not exist
	    + sudo apt-get install python3-pip

	- Install site-packages
	    + pip3 install matplotlib
		+ pip3 install numpy
		+ pip3 install sicpy
		+ pip3 install mpi4py (read more information on https://mpi4py.readthedocs.io/en/stable/install.html)

    - (To reinstall site-packages, for example mpi4py
        + pip3 uninstall mpi4py
        + pip3 install mpi4py --no-cache-dir)

2. Test F2PY
	- To wrap the Fortran subroutine, open Terminal, execute the command:
		+ python3 -m numpy.f2py -c fortranlib.f95 -m fortranlib
	- Run test_f2py.py module:
		+ python3 test_f2py.py
		
3. Test MPI4PY
	- Run test_mpi4py.py module
		+ Execute the command: 
		  - python3 test_mpi4py.py
		-> Output: size: 1, rank: 0
		+ Execute one of these commands:
		  - mpiexec -n 2 python3 test_mpi4py.py
		  - mpirun -n 2 python3 test_mpi4py.py
		-> Output: size: 2, rank: 0 | size: 2, rank: 1		
	
4. Install Paraview for CFD Post-Processing
	- https://www.paraview.org/

5. To edit Python code, download and install Pycharm Community
    - https://www.jetbrains.com/pycharm/

6. Compiling VnCFD_2.21.5
    In the folder "lib", open Terminal, compile the "fortranlib" by one of these commands:
    - ./compile_fortranlib
    - python3 -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib
