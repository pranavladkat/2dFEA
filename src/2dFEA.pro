TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

APP_DIR = /home/pranavpr/Google Drive/Classes/wind-turbine/panel_method/code/3dPanel/3dPanel
PETSC_DIR = /home/pranavpr/petsc
PETSC_ARCH = arch-linux2-c-debug
PETSC_BINS = $$PETSC_DIR/$$PETSC_ARCH/bin
MPI_BIN = /usr/local/openmpi/bin

INCLUDEPATH += -I /home/pranavpr/petsc/include -I /home/pranavpr/petsc/arch-linux2-c-debug/include -I /usr/local/openmpi/include \
							 /usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/include \

LIBS += -Wl,-rpath,/home/pranavpr/petsc/arch-linux2-c-debug/lib -Wl,-rpath,/home/pranavpr/petsc/arch-linux2-c-debug/lib -L/home/pranavpr/petsc/arch-linux2-c-debug/lib -lpetsc -llapack -lblas -lX11 -lssl -lcrypto -lpthread -lm -Wl,-rpath,/usr/local/openmpi/lib -L/usr/local/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lmpi_cxx -lstdc++ -ldl -lmpi -lgcc_s -lpthread -ldl \
				-L/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64  -leng -lmx

unix:QMAKE_RPATHDIR += /usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXX = mpicxx

HEADERS += \
		mesh.hpp \
		material.hpp \
		preprocessor.hpp \
		quadrature.hpp \
		element.hpp \
		stiffelement.hpp \
		boundarycondition.hpp \
    functions.h

OTHER_FILES += \
		lapac_example.txt



