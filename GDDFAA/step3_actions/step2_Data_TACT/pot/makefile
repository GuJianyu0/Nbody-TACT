#### 只有源文件.cpp出目标文件.o, 头文件.h没有; 因为data_interface.h文件里都有了, 没用到其.cpp, 所以之前编译通过
include ../Makefile.inc

default: clean all

all:  libpot_js.so
forceall:	clean all

#------------------------------------------------------
SRC = src/
POT = pot/
OBJ = obj/

OBJECTS = $(OBJ)potential.o $(OBJ)Multipole.o $(OBJ)orbit.o
#gjy add data_interface.o

HEADERS = $(wildcard inc/*.h)

#------------------------------------------------------

GENERALPATH = ../general/
CTPATH = ../general/coordtransforms/
GNUPLOTLIBPATH = $(GENERALPATH)gnuplot

ifeq ($(TORUS), 1)
	TORUSINCPATH = -I$(TORUSPATH)/src/ -I$(TORUSPATH)/src/utils/ -I$(TORUSPATH)/src/pot
	TORUSLIBPATH = -L$(TORUSPATH)/obj/ -L$(TORUSPATH)/WDlib/obj/
else
	TORUSINCPATH =
	TORUSLIBPATH =
endif

COORDINCPATH = $(CTPATH)inc/
COORDLIBPATH = $(CTPATH)

INC = -I../../ -Iinc/ $(TORUSINCPATH) -I$(GENERALPATH) -I$(GSLPATH)include/ -I$(COORDINCPATH)

LIBWLRPATH = -Wl,-rpath,$(TORUSPATH)/obj/ -Wl,-rpath,$(TORUSPATH)/WDlib/obj/  -Wl,-rpath,$(EBFPATH)/lib/

LIB = $(LIBWLRPATH) -L$(GENERALPATH) -L$(GNUPLOTLIBPATH) $(TORUSLIBPATH) -L$(GSLPATH)lib/ \
	  -L$(COORDLIBPATH) -lcoords $(TORUSLIBS) -lplot -lgsl -lgslcblas \
	  -lm $(EXTRACLIB)

BOOSTLIB = -lboost_python -lpython2.7
LIBPY = #$(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)

#gjy add: DataInterface
LIBDATASO = -L../DataInterface/ -lDataInterface -lgfortran
#------------------------------------------------------

$(OBJ)%.o : $(SRC)%.cpp ${HEADERS}
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) $(TORUSFLAG)

multipole_test.exe: $(OBJECTS)
	$(CCOMPILER) $(OPTS) -o multipole_test.exe $(OBJECTS) $(LIB)

libpot_js.so: $(OBJECTS)
	$(CCOMPILER) $(LIBFLAG) $(OPTS) -o libpot_js.so -Wl,-$(LIBNAME),$(CURDIR)/libpot_js.so $(OPTS) $(OBJECTS) $(INC) $(LIB) $(LIBDATASO) $(TORUSFLAG)
#gjy add: $(LIBDATASO)

libpot_js_nt.so: $(OBJECTS)
	$(CCOMPILER) $(LIBFLAG) $(OPTS) -o libpot_js_nt.so -Wl,-$(LIBNAME),$(CURDIR)/libpot_js_nt.so $(OPTS) $(OBJECTS) $(INC) $(LIB)

clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe
