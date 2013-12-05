TEMPLATE = app
CONFIG += console debug
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    colorimage.cpp \
    utils.cpp

QMAKE_CXXFLAGS += -std=c++0x

## MPI Settings
#QMAKE_CXX = mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = mpicc
#
#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

## C++11 macros
#Q_COMPILER_RVALUE_REFS
#Q_COMPILER_DECLTYPE
#Q_COMPILER_VARIADIC_TEMPLATES
#Q_COMPILER_AUTO_TYPE
#Q_COMPILER_EXTERN_TEMPLATES
#Q_COMPILER_DEFAULT_DELETE_MEMBERS
#Q_COMPILER_CLASS_ENUM
#Q_COMPILER_INITIALIZER_LISTS
#Q_COMPILER_LAMBDA
#Q_COMPILER_UNICODE_STRINGS

HEADERS += \
    settings.h \
    colorimage.h \
    utils.h \
    ControlDict.h
