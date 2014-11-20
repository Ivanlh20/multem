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
The following steps work using Matlab R2014a, GCC 4.9, and CUDA 6.5, Linux 64-bit, assuming nvcc/gcc is in PATH:
 - copy the `mex_CUDA_glnxa64.xml` file from `/usr/local/MATLAB/R2014a/toolbox/distcomp/gpu/extern/src/mex/glnxa64/` to `~/.matlab/R2014a`
 - Edit the file with a text editor so that the LINKLIBS line points to you CUDA SDK libraries instead of MATLAB's own versions:
   As an example:

        LINKLIBS="-Wl,-rpath-link,$MATLABROOT/bin/$ARCH -L&quot;$MATLABROOT/bin/$ARCH&quot; -lmx -lmex -lmat -lm -lmwgpu -L/opt/cuda/lib64 -lcufft -lcudart"
   
 - run the `compileMULTEM.m` script. This will create `MULTEMMat.mexa64` and `CreateCrystalbyLayers.mexa64`, which need to be copied alongside the files using the functions.
 - To run the examples, copy `MULTEMMat.mexa64` to the `Multislice examples` directory and `CreateCrystalbyLayers.mexa64` to the `Build Crystal` directory so the scripts in those directories can call the respective functions.
