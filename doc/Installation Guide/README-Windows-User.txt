NOTE. So many troubles in one installation guide. So I'd recommend you to use Linux distribution(e.g. Ubuntu).

1. Install Mingw-w64
	- Download MinGW-W64 "GCC-8.1.0 x86_64-posix-sjlj" from sourceforge:
	    (Note: Make sure it is "x86_64-posix-sjlj", there are some other packages)
		https://sourceforge.net/projects/mingw-w64/files/mingw-w64/
	- Unzip and copy the "mingw64" folder into C: drive (C:\mingw64)
	- You may need to install "7zip" to unzip file: https://www.7-zip.org/
	- Add "C:\mingw64\bin" to user Path:
		+ "Right-click" on "This Computer", click "Properties"
		+ Click "Advanced system settings" on the left
		+ Select "Environment Variables..." at the bottom
		+ Select "Path" in the "System variables" section
		+ Click "Edit", then click "New" and add path "C:\mingw64\bin", finally, click "Ok"
	- Check gcc, g++, gfortran version:
		+ Open "Command-line"
		+ Execute command: "gcc -v"
		+ Execute command: "g++ -v"
		+ Execute command: "gfortran -v"

2. Install Microsoft MPI
	- Download Microsoft MPI installer "msmpisetup.ext" and "msmpisdk.msi" from
	  https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi#ms-mpi-downloads
	- Run both files
	- Add paths, where you installed the MS-MPI (in my case: C:\Program Files (x86)\Microsoft SDKs\MPI\ and
	  C:\Program Files\Microsoft MPI\Bin), to user Path (as add "C:\mingw64\bin" to user Path).

3. Install Latest Python 3 Release
	- Download and run installer from https://www.python.org/downloads/windows/
	- In Setup window, check "Add Python to PATH" 
	- Check Python version in "Command-line": 
		+ Press "Windows+R" to open "Run" box, then type "cmd" and press "Enter"
		+ Execute command: type "python -- version" and press "Enter"

4. Install site-packages:
	- Open "Command-line", check the version of "pip":
		+ Execute command: "pip -V"
	- To install matplotlib, numpy, scipy, mpi4py, run:
		+ pip install matplotlib
		+ pip install numpy
		+ pip install scipy
		+ pip install mpi4py

5. Test F2PY | Install Visual Studio Community
	- In this folder, press "Shift+Right-click", select open "PowerShell"
	- Or press "ctrl+L" then "ctrl+C" to copy folder path, open "Command-line", 
	  type "cd" then press "ctrl+V" to paste the path and press "Enter"
	- Wrap the Fortran subroutine to Python by one of 2 commands
		+ f2py -c fortranlib.f95 -m fortranlib  
		+ python -m numpy.f2py -c fortranlib.f95 -m fortranlib
	- If "vscvarsall.bat" not found, download and run Visual Studio Community installer from
		+ https://visualstudio.microsoft.com/
		+ select "Python development"
		+ on the right in "Summary" deselect all, then select "Python native development tools"
		+ (may be you need "Desktop development with C++" too)
	- (Read more about this error on stackoverflow: 
		https://stackoverflow.com/questions/48826283/compile-fortran-module-with-f2py-and-python-3-6-on-windows-10)
	- Run f2py again
	- If error again, find where "vscvarsall.bat" is located 
	  (In my case: "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build") and add this path
    to user Path (as add "C:\mingw64\bin" path)
	- Run test_f2py.py module:
		+ python test_f2py.py
	- If "ImportError: DLL load failed..." appear, copy and paste the "libfortran...dll" (generated in "fortranlib\.libs") 
	  into this directory (where there's "fortranlib...amd64.pyd")
	- Run test_f2py again

6. Test MPI4PY
	- Run test_mpi4py.py module
		+ Execute command: python test_mpi4py.py
		-> Output: size: 1, rank: 0
		+ Execute command:
		    mpiexec -n 2 python test_mpi4py.py
			or
			mpirun -n 2 python test_mpi4py.py
		-> Output: size: 2, rank: 0 | size: 2, rank: 1

	- "ImportError: DLL load failed..." maybe appear if you used "Anaconda environment" to  install "mpi4py",
	  then you need to remove and reinstall it
		+ conda remove mpi4py (or pip uninstall mpi4py)
		+ pip install mpi4py  (or pip install mpi4py --no-cache-dir)
		+ run test_mpi4py.py again

7. Install Paraview for CFD Post-Processing
	- Download and install Paraview from: https://www.paraview.org/
	
8. If you don't have any Code Editors, download and install Pycharm Comunity and Notepad++:
	- https://www.jetbrains.com/pycharm/
	- https://notepad-plus-plus.org/

9. Compiling VnCFD_2.21.5
    - In folder "lib", open the Command-line (PowerShell), compile the "fortranlib" by one of these commands:
        + ./compile_fortranlib
        + python -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib
    - copy and paste the "libfortran...dll" (generated in "fortranlib\.libs") into this directory (if it's necessary, as in section 5.)
