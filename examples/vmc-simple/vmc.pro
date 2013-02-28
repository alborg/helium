TEMPLATE = app
CONFIG += console
CONFIG -= qt

win32 {

LIBS += -LC:\Libs\armadillo\include\ -larmadillo -LC:\Libs\ -llapack -LC:\Libs\ -lblas

    #"C:/Program Files (x86)/MPICH2/lib/mpi.lib"

INCLUDEPATH += C:\Libs C:\Libs\armadillo\include #"C:/Program Files (x86)\MPICH2\include"

DEPENDPATH += C:\Libs

}

unix {

LIBS += -llapack -lblas -larmadillo -L/usr/lib64/mpich2/lib -lmpich
INCLUDEPATH += /usr/include/mpich2-x86_64/

}

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    WaveFunction.cpp \
    hamiltonian.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    WaveFunction.h \
    hamiltonian.h

