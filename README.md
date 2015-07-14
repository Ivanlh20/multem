Introduction
============

MULTEM is a collection of routines written in C++ with CUDA (and soon also
OpenCL) to simulate (S)TEM experiments.
Currently, there are two supported ways to use MULTEM:
 - Matlab: using the mex interface
 - C++: using the library itself
The library is under heavy development and subject to change.
The Matlab interface is the recommended way

Using MULTEM with Matlab on Linux
=================================
The following works using Matlab R2015a, GCC 4.9, and CUDA 7.0, Linux 64-bit, assuming nvcc/gcc is in PATH:
 - modify the file `mex_functions/mex_CUDA_glnxa64.xml` so that the `LINKLIBS` line points to you CUDA SDK libraries instead of MATLAB's own versions:
   As an example:

        LINKLIBS="-Wl,-rpath-link,$MATLABROOT/bin/$ARCH -L&quot;$MATLABROOT/bin/$ARCH&quot; -lmx -lmex -lmat -lm -lmwgpu -L/opt/cuda/lib64 -lcufft -lcudart"
   
 - run the `compile_mex_files.m` script from within MATLAB. This will create all the MULTEM mex files in the `Mex_Executables` directory.
 
Using MULTEM with Matlab on Windows 64-bit
==========================================

The following works using Matlab R2015a, Visual Studio 2013 (Community), CUDA 7.0, and Windows 64-bit.
 - run the `compile_mex_files.m` script from within MATLAB. 
