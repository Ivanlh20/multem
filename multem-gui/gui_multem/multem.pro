CONFIG += qt thread
QT       += core

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

DESTDIR = "executable"
TEMPLATE = app
TARGET = multem

multem_path = "$$dirname(PWD)"
multem_src_path = "$${multem_path}/src"
multem_lib_path = "$${multem_path}/mex_bin"
multem_gui_path = "$${multem_path}/gui_multem"

HEADERS  = \
    mainwindow.h \
    q_types.h \
    q_colormap.h \
    q_widget.h \
    q_load_specimen.h \
    q_image.h \
    q_data_viewer.h \
    q_colormap.h

SOURCES = \
    main.cpp \
    mainwindow.cpp \
    q_widget.cpp

CUDA_SOURCES = multem.cu

OTHER_FILES  += $${CUDA_SOURCES}

LIBS_EXT = -L'"$${multem_src_path}"' -llibfftw3f-3 -llibfftw3-3 -llibblas -lliblapack

# Include source code path
INCLUDEPATH += $${multem_src_path}
DEPENDPATH += $${multem_src_path}
win32::LIBS += $${LIBS_EXT}

#################################################################
#-------------------------Cuda setup-----------------------------
#################################################################

# System machine
SYSTEM_TYPE = 64   # '32' or '64', depending on OS bitness

# Compute architecture= 20 30 32 35 37 50 60
CUDA_COMPUTE_ARCH = 60
for(_a, CUDA_COMPUTE_ARCH):{
    formatted_arch = $$join(_a,' ','-gencode=arch=compute_$${_a},code=sm_', '')
    CUDA_ARCH += $$formatted_arch
}

# Path to cuda toolkit install
macx:CUDA_DIR = "/Developer/NVIDIA/CUDA-8.0"
linux:CUDA_DIR = "/usr/local/cuda-8.0"
win32:CUDA_DIR = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0"
CUDA_DIR = $$clean_path($${CUDA_DIR})

# Path to cuda SDK install
macx:CUDA_SDK = "/Developer/NVIDIA/CUDA-8.0/samples"
linux:CUDA_SDK = "/usr/local/cuda-8.0/samples"
win32:CUDA_SDK = "C:\ProgramData\NVIDIA Corporation\CUDA Samples\v8.0"
CUDA_SDK = $$clean_path($${CUDA_SDK})

# Cuda libraries
macx:CUDA_LIBS += "$${CUDA_DIR}/lib"
linux:CUDA_LIBS += "$${CUDA_DIR}/lib$${SYSTEM_TYPE}"
win32:CUDA_LIBS += "$${CUDA_DIR}/lib/x$${SYSTEM_TYPE}"
CUDA_LIBS = $$clean_path($${CUDA_LIBS})

QMAKE_LIBDIR += $${CUDA_LIBS}

# Includes cuda
INCLUDEPATH += $${CUDA_DIR}/include
DEPENDPATH += $${CUDA_DIR}/include
LIBS += -L$${CUDA_LIBS} -lcuda -lcudart -lcufft
win32:INCLUDEPATH += $${CUDA_SDK}\common\inc

# Join the includes and libraries in a line
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')
CUDA_LIBS = -L'"$${CUDA_LIBS}"' -lcuda -lcudart -lcufft

# nvcc flags (ptxas option verbose is always useful)
CONFIG(debug, debug|release){
    win32:MSVCRT_LINK_FLAG_DEBUG = /MDd
    win32:NVCC_FLAGS = "-g --compiler-options=/FS,/bigobj,/GR,/W3,/wd4819,/EHs,/nologo,/Od,/Zi,/RTC1,/MDd"
    #win32:NVCCFLAGS += -D_DEBUG -Xcompiler $${MSVCRT_LINK_FLAG_DEBUG}
    OBJECTS_DIR = debug
    CUDA_OBJ_DIR = debug
}
else{
    win32:MSVCRT_LINK_FLAG_RELEASE = /MD
    win32:NVCC_FLAGS = "--compiler-options=/FS,/bigobj,/GR,/W3,/wd4819,/EHs,/nologo,/O2,/MD,/Oy-,/DNDEBUG"
    OBJECTS_DIR = release
    CUDA_OBJ_DIR = release
}

CUDA_DEFINES += "--compiler-options=/D_CRT_SECURE_NO_DEPRECATE,/D_SCL_SECURE_NO_DEPRECATE,/D_SECURE_SCL=0"

NVCC_FLAGS += $${CUDA_ARCH} $${CUDA_DEFINES}

QMAKE_CFLAGS_DEBUG += /Zi /MDd /MP2
QMAKE_CXXFLAGS_DEBUG += /MDd
QMAKE_CFLAGS_RELEASE += /O2 /MD /MP2
QMAKE_CXXFLAGS_RELEASE += /MD

macx:CUDA_OBJ_EXT = .o
linux:CUDA_OBJ_EXT = .o
win32:CUDA_OBJ_EXT = .obj

cuda.input = OTHER_FILES
cuda.output = $${CUDA_OBJ_DIR}/${QMAKE_FILE_BASE}_cuda$${CUDA_OBJ_EXT}
cuda.commands = '"$${CUDA_DIR}/bin/nvcc.exe"' -c --machine $${SYSTEM_TYPE} $${NVCC_FLAGS} $${CUDA_INC} $${CUDA_LIBS} $${LIBS_EXT} -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
cuda.dependency_type = TYPE_C
QMAKE_EXTRA_COMPILERS += cuda

message("$${cuda.commands}")

RESOURCES   = multem.qrc
