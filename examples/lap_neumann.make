
EXEC = int2

FC = gfortran
FFLAGS = -O2 -c -w 
FLINK = gfortran -w -o $(EXEC)
FEND = -L/usr/local/opt/openblas/lib -lopenblas

SRC = ../src

EFOL = ../helm_neuamnn_Data


.PHONY: all clean list

SOURCES =  lap_neumann.f \
  $(SRC)/prini.f \
  $(SRC)/legeexps.f \
  $(SRC)/quaplot.f \
  $(SRC)/second_f90.f \
  $(SRC)/pplot2.f \
  $(SRC)/hkrand.f \
  $(SRC)/dlaran.f \
  $(SRC)/adapgaus_new.f \
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



