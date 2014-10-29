@echo off
rem
rem MEX OPTS for use with nvmex (nvcc CUDA compiler) on 64-bit windows.   
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%

rem *****These two lines are all that should need to be tweaked (and target GPU below):
rem *****(Make sure you have sintalled the MS SDK and VS express 10)

set VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio 10.0
set SDKDIR=C:\Program Files (x86)\Microsoft SDKs\Windows\v7.1A

set VCINSTALLDIR=%VSINSTALLDIR%\VC
set PATH=%VCINSTALLDIR%\BIN\;%VSINSTALLDIR%\Common7\IDE;%SDKDIR%\bin\;%VSINSTALLDIR%\Common7\Tools;%VCINSTALLDIR%\VCPackages;%MATLAB_BIN%;%PATH%

set INCLUDE=%VCINSTALLDIR%\INCLUDE;%SDKDIR%\include;%INCLUDE%
set LIB=%SDKDIR%\Lib\x64;%MATLAB%\extern\lib\win64;%VCINSTALLDIR%\LIB\amd64;%LIB%

set MW_TARGET_ARCH=win64

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=nvcc
set COMPFLAGS=-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_35,code=compute_35 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_30,code=compute_30 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20 -c -Xcompiler "/c /Zp8 /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD"
set OPTIMFLAGS=-Xcompiler "/O2 /Oy- /DNDEBUG"
set DEBUGFLAGS=-Xcompiler "/Zi /Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb""
set NAME_OBJECT= 

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win64\microsoft
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /MAP /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x /MACHINE:X64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib

set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/DEBUG /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
set POSTLINK_CMDS1=del %LIB_NAME%.x
set POSTLINK_CMDS2=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%";2 -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
rem *****set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest" 
