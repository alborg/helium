TEMPLATE = app
CONFIG += console
CONFIG -= qt

win32 {

LIBS += -LC:\Libs\armadillo\include\ -larmadillo \
    -LC:\Libs\ -llapack \

INCLUDEPATH = C:\Libs\armadillo\include "C:/Program Files (x86)/MPICH2/include"
}

unix {

LIBS += -llapack -lblas -larmadillo "C:/Program Files (x86)/MPICH2/lib/mpi.lib"

}

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    ../WaveFunction.cpp \
    WaveFunction.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    ../WaveFunction.h \
    WaveFunction.h

