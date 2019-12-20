
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

SOURCES =  lap_poly_solver.f \
  $(SRC)/prini.f \
  $(SRC)/legeexps.f \
  $(SRC)/quaplot.f \
  $(SRC)/second_f90.f \
  $(SRC)/pplot2.f \
  $(SRC)/hkrand.f \
  $(SRC)/hank103cc.f \
  $(SRC)/dlaran.f \
  $(SRC)/adapgaus_new.f \
  $(SRC)/dgmres_blas.f \
  $(SRC)/dgmres_fmm.f \
  $(SRC)/cadapgaum.f \
  $(SRC)/cadapgau.f \
  $(SRC)/orthom.f \
  $(SRC)/pnpoly.f \
  $(SRC)/levrtree2d.f \
  $(SRC)/cfmm2dpart.f \
  $(SRC)/lfmm2drouts.f \
  $(SRC)/d2mtreeplot.f \
  $(SRC)/d2tstrcr_omp.f \
  $(SRC)/laprouts2d.f \
  $(SRC)/l2dterms.f \
  $(SRC)/vol_plot_routs.f \




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



