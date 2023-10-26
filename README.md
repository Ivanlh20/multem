# MULTEM

## Introduction
**MULTEM** is a powerful and advanced software package designed to provide researchers with a versatile tool for simulating a wide range of electron microscopy experiments. Developed by Ivan Lobato (Ivanlh20@gmail.com), MULTEM is built on a collection of C++ routines with CUDA support, enabling it to perform efficient and accurate multislice simulations for various TEM experiments.

MULTEM uses the widely adopted multislice method to simulate electron scattering and wave propagation in a crystal. This method involves dividing the crystal into thin slices and calculating the electron scattering and wave propagation in each slice. This allows for accurate and efficient simulations of various electron microscopy experiments such as high-resolution TEM (HRTEM), scanning TEM (STEM), imaging STEM (ISTEM), electron diffraction (ED), precession electron diffraction (PED), convergent beam electron diffraction (CBED), annular dark field-TEM (ADF-TEM), annular bright field Hollow Cone (ABF-HC), energy filtered TEM (EFTEM), and electron energy loss spectroscopy (EELS).

MULTEM's implementation is further enhanced by its support for CUDA, a parallel computing platform developed by NVIDIA. This feature enables MULTEM to use graphics processing units (GPUs) for simulations, greatly reducing computation time and increasing simulation speed.

Currently, there are three ways to use MULTEM::
- C++: directly using the library
- Matlab: using the provided mex interface
- GUI: using the user-friendly graphical interface 

Please note that the library is under active development and subject to change. The Matlab interface is the recommended way for researchers to use MULTEM.

## Remarks

In order to use the GPU capability of MULTEM, you need a Nvidia Graphic card with **compute capability greater than 3.5** and **CUDA 11.8** installed in your operating system. You can check the compute capability of your graphic card using the following nvidia website: https://developer.nvidia.com/cuda-gpus.

### Using precompiled GUI interface

The precompiled GUI interface is only available for Windows operating system.

- Go to [https://github.com/Ivanlh20/MULTEM/releases](https://github.com/Ivanlh20/MULTEM/releases) and download `MULTEM_binary.7z`.
- Execute `vc_redist.x64.exe` located in `gui_bin` folder.
- Execute `multem.exe`.

### Using precompiled mexfiles for Matlab

The precompiled mexfiles are only available for Windows operating system and Ubuntu 18.04-based Linux distributions.

- Go to [https://github.com/Ivanlh20/MULTEM/releases](https://github.com/Ivanlh20/MULTEM/releases) and download `MULTEM.zip`.
- Execute `vc_redist.x64.exe` located in `mex_bin` folder. (Windows only)
- Add the following folders to the Matlab path: crystalline_materials, matlab_functions and mex_bin.
- Run the examples located in 'mex_examples_multem'.

### Building MULTEM for Matlab

The following steps have been tested and found to work with Matlab 2022b and CUDA 11.8. It is assumed that a C++ compiler such as Visual Studio 2019 Community, g++11.3 or Clang (Xcode 10.x) is installed on your operating system. Additionally, MULTEM also requires the fftw3, BLAS, and LAPACK libraries to be installed.The following steps have been tested and found to work with Matlab 2022b and CUDA 11.8. It is assumed that a C++ compiler such as Visual Studio 2019 Community, g++11.3 or Clang (Xcode 10.x) is installed on your operating system. Additionally, MULTEM also requires the fftw3, BLAS, and LAPACK libraries to be installed.

- Firstly, a C++ compiler must be set for Matlab by executing the following command: `mex -setup cpp`. It is important to note that Matlab 2022b only supports the compilers listed above.
- Next, add the following folders to the Matlab path: crystalline_materials, matlab_functions, and mex_bin.
- Run the script `compile_mex_multem.m`. This will create the necessary executable files to run the examples.
- Finally, run the examples located in the `mex_examples_multem folder`.

### Troubleshooting

- If MULTEM does not compile with the above procedures, one of the following procedures might fix it

  **for Windows:**
  
  	- Verify the installation of Visual studio 2019 community.
  	- Verify the installation of Cuda 11.8 (https://developer.nvidia.com/cuda-downloads).
  	
  **for Linux:**
  
  	- Verify that gcc-11.3 and g++11.3 are the default compilers installed in your operating system. In Ubuntu, it can be installed by executing the following commands:
  	  ```bash
      sudo apt-get update
      sudo apt-get install gcc-11.3 g++-11.3
      ```

  	- Verify the correct installation of Cuda 11.8 (https://developer.nvidia.com/cuda-downloads).
  	
    - Verify the installation of fftw3 libraries. In Ubuntu, it can be installed by executing the following command: 
      ```bash
      sudo apt-get install libfftw3-dev libfftw3-doc
      ```
    
    - Verify the installation of blas and lapack libraries. In Ubuntu, it can be installed by executing the following command: 
      ```bash
      sudo apt-get install libblas-dev liblapack-dev
      ```

- Verify the installation path of cuda 11.8, fftw3, blas and lapack. Their installation paths should be specified in the [ilm_mex.m](./matlab_functions/ilm_mex.m).

**Please cite MULTEM in your publications if it helps your research:**
```bibtex
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
 ```
**if you use our parameterization of the electronscattering factors, please cite the following article:** 
```bibtex
	@Article{LD14_1,
  	Title = {{An accurate parameterization for the scattering factors, electron densities and electrostatic potentials for neutral atoms that obey all physical constraints}},
  	Author = {I. Lobato and D. Van Dyck},
  	Journal = {Acta Crystallographica Section A},
  	Year = {2014},
  	Pages = {636-649},
  	Volume = {70}
  }
```