####################################################################
## Makefile
####################################################################
## 1. Need to install gsl (http://www.gnu.org/software/gsl/)
## 	  This is necessary

## 2. You will need to set the CXX environment variable to c++ compiler
# CXX=g++
CXX=mpic++
CCOMPILER=$(CXX)
# mypath=$HOME/workroom/0prog/Stackel
mypath=../../../dependencies/
#gjy note: The main makefile is the one in TACT/aa/, so the relative path is about TACT/aa/ (i.e. a more ../ for TACT/)

## 3. You will need to set the path to gsl

GSLPATH=${mypath}

####################################################################
## Tests:
## To run the tests you need google test and to set the path:
GTESTPATH=/usr/local/
# GTESTPATH=$mypath/googletest-master
# <!-- GTESTPATH=/Users/jls/work/code/other/gtest-1.7.0/ -->



####################################################################
## Optional
## 1. Use Torus (https://github.com/PaulMcMillan-Astro/Torus)
### 1.1. Install from https://github.com/PaulMcMillan-Astro/Torus along
###		with ebf library (http://sourceforge.net/projects/ebfformat/files/libebf/c-cpp/)
### 1.2. Set Torus and ebf paths

TORUSPATH=../../Torus-master/
# TORUSPATH=/home/darkgaia/workroom/0prog/Stackel/Torus-master

EBFPATH=${mypath}
# EBFPATH=/home/darkgaia/workroom/0prog/Stackel/libebf_c_cpp-0.0.3

### 1.3. Compile with "make TORUS=1"

## 2. Use LAPACK (http://www.netlib.org/lapack/)
### 2.1. Install LAPACK from http://www.netlib.org/lapack/
### 2.2. Set LAPACK path

LAPACKPATH=/usr/lib/x86_64-linux-gnu/lapack/

### 2.3. Compile with "make LAPACK=1"



####################################################################
## Python module
####################################################################
## These paths are only important if you want to make python modules
## Relies on boost
## compile the shared library as
## cd aa/; make lib/aa_py.so; cd -

BOOSTINCPATH =/usr/local/include/boost
BOOSTLIBPATH =-L/usr/local/lib
PYTHONINCPATH=/usr/include/python2.7/
PYTHONLIBPATH=-L/usr/lib/python2.7/



####################################################################
## Compiler flags
DEBUG = -ggdb3
# -ffast-math -Ofast
# OPTS = $(CXXFLAGS) -std=c++0x -Wall -O2 -fopenmp -fPIC
# OPTS = $(CXXFLAGS) -std=c++14 -Wall -O2 -fopenmp -fPIC
OPTS = $(CXXFLAGS) -std=c++17 -Wall -O2 -fopenmp -fPIC -ggdb3



####################################################################
##
ifeq ($(TORUS),1)
	TORUSFLAG=-DTORUS
else
	TORUSFLAG=
endif

ifeq ($(LAPACK),1)
	LAPACKFLAG=-DLAPACK
else
	LAPACKFLAG=
endif

export TORUS
export TORUSFLAG
export LAPACK
export LAPACKFLAG

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBNAME =soname
    LIBFLAG=-shared
endif
ifeq ($(UNAME_S),Darwin)
    LIBNAME=install_name
    LIBFLAG=-dynamiclib -undefined suppress -flat_namespace
endif

export LIBNAME
export LIBFLAG
export UNAME_S
