
default: clean scfcoef so exe

runprog: run

some: clean so

clean:
	rm -rf *.o *.mod *.mod0 *.so *.exe
	rm -rf a2b2c2.txt delta_i.txt record.txt some_lmn_foci_Pot.txt Integrate.status
	# rm -rf orbit/*.dat orbit/*.png



##-----------------------------------------------
## compiler and options
CCOMPILER0 = gcc
CCOMPILER = g++ #mpi++
# CFLAG = -Wall -O2 -fopenmp -fPIC
CFLAG = -std=c++17 -Wall -O2 -fopenmp -fPIC
#[learn code] use extern "C"{include"your_C_lib_code.c"}, otherwise "undefined reference to ..." because Fortran~C instead of C++
LIBFLAG = -shared #enable this options if make shared objects; do not enable this if make executablefile
LIBNAME = # -soname
STATIC = -static-libstdc++

FORTRANCOMPILER = gfortran #mpi++
FORTRANFLAG = -Wall -O3 -fopenmp -fPIC #-ffixed-line-length-0 #-ffree-line-length-0

INC = -I./
LIBSO = -L./ -lORB
ADDLIBPATH = -Wl,--rpath=./
LIBOTHER = -lgsl -lgslcblas -lm -lgfortran -lpthread



##----------------------------------------------
SRC_ALL:=$(wildcard *.c *.cpp *.f *.f90)
#[learn code] use <:=> and <+=> instead of = because <:=> is to change value in order
#\ <=> or < = > is same
FILTEROUT=CoefSCF.f out.cpp
SRC_ALL:=$(filter-out $(FILTEROUT), $(SRC_ALL))

OBJ_ALL:=$(SRC_ALL)
OBJ_ALL:=$(OBJ_ALL:%.c=%.o)
OBJ_ALL:=$(OBJ_ALL:%.cpp=%.o)
OBJ_ALL:=$(OBJ_ALL:%.f=%.o)
OBJ_ALL:=$(OBJ_ALL:%.f90=%.o)



##-----------------------------------------------
%.o: %.f
	$(FORTRANCOMPILER) $(FORTRANFLAG) -c $< -o $@ $(INC)
%.o: %.f90
	$(FORTRANCOMPILER) $(FORTRANFLAG) -c $< -o $@ $(INC)
%.o: %.c
	$(CCOMPILER0) $(CFLAG) -c $< -o $@ $(INC)
%.o: %.cpp
	$(CCOMPILER) $(CFLAG) -c $< -o $@ $(INC)

so: $(OBJ_ALL)
	# echo "$(SRC_ALL)"
	# echo "$(OBJ_ALL)"
	$(CCOMPILER) $(CFLAG) $(LIBFLAG) $(SRC_ALL) -o libORB.so $(INC) $(ADDLIBPATH) $(LIBOTHER)

exe:
	rm -rf out.o
	$(CCOMPILER) $(CFLAG) -c out.cpp -o out.o $(INC)
	$(CCOMPILER) $(CFLAG) out.o -o out.exe $(ADDLIBPATH) $(LIBSO) $(LIBOTHER) $(INC)



##-----------------------------------------------
#to prepare
scfcoef:
	# cp ../your_path_to_SCF_setting/scfmod ./
	# cp ../your_path_to_SCF_setting/scfpar ./
	# cp ../your_path_to_SCF_setting/scficoef ./
	# cp ../your_path_to_SCF_setting/scfocoef ./
	# gfortran -O3 CoefSCF.f -o CoefSCF.exe #-ffree-form
	cp ../../SCF_coeff_pot/scfmod ./
	cp ../../SCF_coeff_pot/scfpar ./
	cp ../../SCF_coeff_pot/scficoef ./
	cp ../../SCF_coeff_pot/scfocoef ./
	gfortran -O3 CoefSCF.f -o CoefSCF.exe #-ffree-form

#to run
run:
	# ./CoefSCF.exe
	# cp scfocoef scficoef
	# rm -rf Orbit.dat NaN-record
	./out.exe
