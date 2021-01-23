@echo off
Rem Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

Rem Terminal: ./compile_fortranlib

python -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib
