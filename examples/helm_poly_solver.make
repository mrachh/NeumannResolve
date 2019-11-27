
EXEC = int3

FC = gfortran
FFLAGS = -O2 -c -w 
FLINK = gfortran -w -o $(EXEC) 
FEND = -L/usr/local/opt/openblas/lib -lopenblas 
#FEND = -lopenblas

ifeq ($(OMP),ON)
FFLAGS = -O2 -c -w --openmp
FLINK = gfortran -w -o $(EXEC) --openmp

endif

SRC = ../src

EFOL = ../data


.PHONY: all clean list

SOURCES =  helm_poly_solver.f \
  $(SRC)/prini.f \
  $(SRC)/legeexps.f \
  $(SRC)/quaplot.f \
  $(SRC)/second_f90.f \
  $(SRC)/pplot2.f \
  $(SRC)/hkrand.f \
  $(SRC)/hank103.f \
  $(SRC)/dlaran.f \
  $(SRC)/adapgaus_new.f \
  $(SRC)/cadapgaum.f \
  $(SRC)/cadapgau.f \
  $(SRC)/orthom.f \
  $(SRC)/pnpoly.f \




OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	mkdir -p $(EFOL)
	$(FLINK) $(OBJECTS) $(FEND)
	mv $(EXEC) ./$(EFOL)/
	cd $(EFOL) && ./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



