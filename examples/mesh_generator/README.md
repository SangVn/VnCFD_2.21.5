Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

### Run py-module in terminal with the command
  - python module.py

### Run fortran module
  - gfortran name.f 
  - ./a.out

### CDV Nozzle
  [Source: NASA](https://www.grc.nasa.gov/WWW/wind/valid/cdv/cdv01/cdv01.html)

### Laminar Flat Plate
  [Source: NASA](https://www.grc.nasa.gov/WWW/wind/valid/fplam/fplam01/fplam01.html)
  
### naca0012
  [Source: NASA](https://turbmodels.larc.nasa.gov/naca0012_grids.html)  

### Note
  VnCFD_2.21.5 uses only structured mesh. If you have a CGNS mesh,
     * open it with paraview, press "Apply",
     * click "Information" to know mesh size (dimension)
     * (if it's 3D mesh, make a "Slice" to get 2D mesh),
     * select "File" -> "Save Data" -> type "File name" (e.g. "mesh.txt") -> Enter
     * to convert it to *VnCFD_2.21.5 mesh format*, modify and run *convert_cgns_mesh_1()* function in module *convert_mesh.py*

