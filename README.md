# MULTEM

Introduction
============

**MULTEM** is a collection of routines written in C++ with CUDA to perform accurate and fast multislice simulations for different TEM experiments as: HRTEM, STEM, ISTEM, ED, PED, CBED, ADF-TEM, ABF-HC, EFTEM and EELS. It is developed by Ivan Lobato (Ivanlh20@gmail.com).

Currently, there are three supported ways to use MULTEM:
- C++: using the library itself
- Matlab: using the mex interface
- GUI: using the user graphical interface 

The library is under heavy development and subject to change.
The Matlab interface is the recommended way for researchers.

Remarks
=================================
In order to use the GPU capability of MULTEM, you need a Nvidia Graphic card with **compute capability greater than 2.0** and **CUDA 8.0** installed in your operating system. You can check the compute capability of your graphic card using the following nvidia website: https://developer.nvidia.com/cuda-gpus.

Using precompile GUI interface
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

Building MULTEM for Matlab
=================================
The following steps work using Matlab R2017a and CUDA 8.0. It assumes that Visual studio 2015 professional, g++4.9 or Clang(Xcode 8.x) compiler is installed in your operating system. Additionally, Multem also requires fftw3, blas and lapack libraries.

- First of all, you have to set a C++ compiler to Matlab by executing the following comand: `mex -setup cpp`. Be aware that Matlab 2017a only support the above compilers.
- Then add to the Matlab path the following folders: crystalline_materials, matlab_functions and mex_bin.
- Run the `compile_mex_multem.m` script. This will create the required executable files to run the examples.
- Run the examples located in `mex_examples_multem`.

Troubleshooting
=================================
- If MULTEM do not compile with the above procedures, one of the following procedures might fix it

  **for Windows:**
  
  	- Verify the installation of Visual studio 2015 professional.
  	- Verify the installation of Cuda 8.0 (https://developer.nvidia.com/cuda-downloads).
  	
  **for Linux:**
  
  	- Verify that gcc-4.9 and g++4.9 are the default compilers installed in your operating system. In Ubuntu, it can be installed by executing the following commands:
  	
  		* `sudo add-apt-repository ppa:ubuntu-toolchain-r/test`
		* `sudo apt-get update`
		* `sudo apt-get install gcc-4.9 g++-4.9`
		* `sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9`

  	- Verify the correct installation of Cuda 8.0 (https://developer.nvidia.com/cuda-downloads).
  	
    - Verify the installation of fftw3 libraries. In Ubuntu, it can be installed by executing the following command: 
    	* `sudo apt-get install libfftw3-dev libfftw3-doc`
    
    - Verify the installation of blas and lapack libraries. In Ubuntu, it can be installed by executing the following command: 
    	* `sudo apt-get install libblas-dev liblapack-dev`
    
- Verify the intallation path of cuda 8.0, fftw3, blas and lapack. Their installation path should be specified in the `MEX.m` file located at `matlab_functions`.

**Please cite MULTEM in your publications if it helps your research:**

    @article{LVAV16_1,
      Author = {I.Lobato and S.Van Aert and J.Verbeeck},
      Journal = {Ultramicroscopy},
      Title = {Progress and new advances in simulating electron microscopy datasets using MULTEM},
      Year = {2016},
  	  volume  = {168},
      pages   = {17-27}      
    }
    
     @article{LD15_2,
      Author = {I. Lobato and D. Van Dyck},
      Journal = {Ultramicroscopy},
      Title = {MULTEM: A new multislice program to perform accurate and fast electron diffraction and imaging simulations using Graphics Processing Units with CUDA},
      Year = {2015},
  	  volume  = {156},
      pages   = {9-17}      
    } 
    
**if you use our parameterization of the electronscattering factors, please cite the following article:** 

	@Article{LD14_1,
  	Title = {{An accurate parameterization for the scattering factors, electron densities and electrostatic potentials for neutral atoms that obey all physical constraints}},
  	Author = {I. Lobato and D. Van Dyck},
  	Journal = {Acta Crystallographica Section A},
  	Year = {2014},
  	Pages = {636-649},
  	Volume = {70}
	}