#MULTEM

##Introduction
============

**MULTEM** is a collection of routines written in C++ with CUDA to perform accurate and fast multislice simulations for different TEM experiments as: CTEM, STEM, ISTEM, ED, PED, CBED, ADF-TEM, ABF-HC, EFTEM and EELS.

Currently, there are three supported ways to use MULTEM:
- C++: using the library itself
- Matlab: using the mex interface
- GUI: using the user graphical interface 

The library is under heavy development and subject to change.
The Matlab interface is the recommended way for researchers.

Using MULTEM with Matlab
=================================
The following steps work using Matlab R2015a and CUDA 7.5. It assumes that Visual studio 2013, g++4.7 or Clang compiler is installed in your operating system:

- First of all, you have to add to the Matlab path the following folders: Crystalline_Materials, Matlab_functions and mex_executables.

- Then you need to modify the `mex_CUDA_xxxx` file located in the 'mex_files' folder, which corresponds to your operating system:
  1. for Windows `mex_CUDA_win64`
  2. for Linux `mex_CUDA_glnxa64`
  3. for Mac `mex_CUDA_maci64`

- Go to the line which contains the following definition`NVCCFLAGS="$xxxxx"` and replace 'xxxxx' by 'SINGLE_CARD' or 'MULTI_CARD'. 'SINGLE_CARD'/'MULTI_CARD' defines the compute capability of your Nvidia graphic card.

- Run the `compile_multem.m` script. This will create the required executable files to run the examples.