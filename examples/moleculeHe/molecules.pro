 TEMPLATE = app
CONFIG += console
CONFIG -= qt


LIBS += -llapack -lblas -larmadillo -L/usr/lib64/mpich2/lib -L/usr/lib/ -lmpi #-lmpich
INCLUDEPATH += /usr/include/mpich2-x86_64/ /usr/include/mpi/



SOURCES += main.cpp \
    vmcsolver.cpp \
    WaveFunction.cpp \
    hamiltonian.cpp \
    lib.cpp

HEADERS += \
    vmcsolver.h \
    WaveFunction.h \
    hamiltonian.h \
    lib.h

#Maximize (-O3)
release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

#C++ 11
COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

#MPI Parallelization
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
