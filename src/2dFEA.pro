TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -llapack -lblas

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXX = mpicxx

HEADERS += \
		mesh.hpp \
		material.hpp

OTHER_FILES += \
		lapac_example.txt



