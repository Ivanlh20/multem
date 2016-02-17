Introduction
============

MULTEM is a collection of routines written in C++ with CUDA to simulate TEM experiments.
Currently, there are three supported ways to use MULTEM:
 - Matlab: using the mex interface
 - C++: using the library itself
 - GUI: using the user graphical interface 
 
The library is under heavy development and subject to change.
The Matlab interface is the recommended way for researchers.

Using MULTEM with Matlab on Linux
=================================
The following steps work using Matlab R2015a and CUDA 7.5. It assume that Visual studio 2013, g++4.7 or Clang compiler is installed in your operating system:

-First of all, you have to add to the matlab path the following folders: Crystalline_Materials, Matlab_functions and mex_executables.

-Then you need to modified the `mex_CUDA_xxxx` file locate in the 'mex_files' folder, which corresponds to your operating system:
for Windows `mex_CUDA_win64`
for Linux `mex_CUDA_glnxa64`
for Mac `mex_CUDA_maci64`

-Go the the line which contains the following definition`NVCCFLAGS="$xxxxx"` and replace 'xxxxx' by 'SINGLE_CARD' or 'MULTI_CARD'. 'SINGLE_CARD'/'MULTI_CARD' defines the compute capability of your Nvidia graphic card.

- run the `compile_multem.m` script. This will create the required executables files, which allow to run the examples.