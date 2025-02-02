#!/bin/bash
ulimit -s unlimited
ifort -c flash_drought_def.f90
ifort -c slow_drought_def.f90
ifort -g -traceback -fpe0 -CB -mcmodel=large FD_def-8days_identification.f90 -L/data/software/netcdf-4.4.1/lib -lnetcdff -lnetcdf -I/data/software/netcdf-4.4.1/include -o FD_def-8days_identification.exe flash_drought_def.o slow_drought_def.o
./FD_def-8days_identification.exe
rm FD_def-8days_identification.exe
rm slow_drought_def.o
rm flash_drought_def.o
