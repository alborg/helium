TEMPLATE = app
CONFIG += console
CONFIG -= qt

win32 {

LIBS += -LC:\Libs\armadillo\include\ -larmadillo \
    -LC:\Libs\ -llapack \

INCLUDEPATH = C:\Libs\armadillo\include \
}

unix {

LIBS += -llapack -lblas -larmadillo

}

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp

HEADERS += \
    vmcsolver.h \
    lib.h

