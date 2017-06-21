# MULTEM

Introduction
============

**MULTEM** is a collection of routines written in C++ with CUDA to perform accurate and fast multislice simulations for different TEM experiments as: HRTEM, STEM, ISTEM, ED, PED, CBED, ADF-TEM, ABF-HC, EFTEM and EELS.

Currently, there are three supported ways to use MULTEM:
- C++: using the library itself
- Matlab: using the mex interface
- GUI: using the user graphical interface 

The library is under heavy development and subject to change.
The Matlab interface is the recommended way for researchers.

Remarks
=================================
In order to use the GPU capability of MULTEM, you need a Nvidia Graphic card with **compute capability greater than 2.0** and **CUDA 8.0** installed in your operating system. You can check the compute capability of your graphic card using the following [nvidia website](https://developer.nvidia.com/cuda-gpus)

Using GUI interface
=================================
The precompile GUI interface is only available for Windows operating system.

- Go to [https://github.com/Ivanlh20/MULTEM/releases](https://github.com/Ivanlh20/MULTEM/releases) and download `MULTEM_binary.7z`.
- Execute `vc_redist.x64.exe` located in `gui_bin` folder.
- Execute `multem.exe`.

Using precompile mexfiles for Matlab
=================================
The precompile mexfiles are only available for Windows operating system.

- Go to [https://github.com/Ivanlh20/MULTEM/releases](https://github.com/Ivanlh20/MULTEM/releases) and download `MULTEM_binary.7z`.
- Execute `vc_redist.x64.exe` located in `mex_bin` folder.
- Add to the Matlab path the following folders: crystalline_materials, matlab_functions and mex_bin.
- Run the examples located in 'mex_examples_multem'.

Compiling MULTEM for Matlab
=================================
The following steps work using Matlab R2017a and CUDA 8.0. It assumes that Visual studio 2015, g++4.9 or Clang compiler is installed in your operating system:

- First of all, you have to add to the Matlab path the following folders: crystalline_materials, matlab_functions and mex_bin.

- Then you need to modify the `mex_CUDA_xxxx` file located in the `mex_files_multem` folder, which corresponds to your operating system:
  1. for Windows `mex_CUDA_win64`
  2. for Linux `mex_CUDA_glnxa64`
  3. for Mac `mex_CUDA_maci64`

- Go to the line which contains the following definition`NVCCFLAGS="$xxxxx ..."` and replace 'xxxxx' by 'SINGLE_CARD' or 'MULTI_CARD'. 'SINGLE_CARD'/'MULTI_CARD' defines the compute capability of your graphic card.

- Run the `compile_mex_multem.m` script. This will create the required executable files to run the examples.
- Run the examples located in 'mex_examples_multem'.
