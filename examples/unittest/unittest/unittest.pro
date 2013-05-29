#-------------------------------------------------
#
# Project created by QtCreator 2013-04-25T10:38:13
#
#-------------------------------------------------

CONFIG   += console
CONFIG   -= qt
TEMPLATE = app

LIBS += -lUnitTest++ -llapack -lblas -larmadillo -L/usr/lib64/mpich2/lib -L/usr/lib/ -lmpi #-lmpich
INCLUDEPATH += /usr/include/unittest++ /usr/include/mpich2-x86_64/ /usr/include/mpi/


SOURCES += \
    main.cpp \

include(../../defaults.pri)

SOURCES += $$SRC_DIR/WaveFunction.cpp \
            $$SRC_DIR/hamiltonian.cpp \
            $$SRC_DIR/slaterDeterminant.cpp \
            $$SRC_DIR/correlation.cpp \
            #$$SRC_DIR/main.cpp \
            $$SRC_DIR/vmcsolver.cpp \
            $$SRC_DIR/lib.cpp

HEADERS += \
    $$SRC_DIR/vmcsolver.h \
    $$SRC_DIR/lib.h \
    $$SRC_DIR/WaveFunction.h \
    $$SRC_DIR/hamiltonian.h \
    $$SRC_DIR/slaterDeterminant.h \
    $$SRC_DIR/correlation.h \
    $$SRC_DIR/lib.h

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

