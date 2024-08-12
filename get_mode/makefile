HOME = .
BIN = $(HOME)/bin
LIB = ./lib
#LIB = /opt/Modes/QmXD/lib

SRC = .
OBJ = $(HOME)/obj
LAPACK = $(HOME)/lib/LAPACK/SRC

#F90 = ifort
#FFLAGS =  -extend_source #-B108 -f -W # -extend_source 
F90 = gfortran 
FFLAGS =  -funroll-all-loops -fno-f2c -O3 #-llapack
#FFLAGS = -ffixed-line-length-132 -g
LFLAGS = -L$(LIB)


forward_new: coupling.h \
	$(OBJ)/decipher_cluster.o \
	$(OBJ)/futil.o \
	$(OBJ)/get_events.o \
	$(OBJ)/tdiffer.o \
	$(OBJ)/forward_new.o \
	$(OBJ)/mantle.o \
	$(OBJ)/mantle_model.o \
	$(OBJ)/model_vector.o \
	$(OBJ)/multiplet2.o \
	$(OBJ)/receiver.o \
	$(OBJ)/rot_ell.o \
	$(OBJ)/singlet.o \
	$(OBJ)/source.o
	$(F90) $(FFLAGS) $(LFLAGS) -o $(BIN)/forward_new \
	$(OBJ)/decipher_cluster.o \
	$(OBJ)/futil.o \
	$(OBJ)/get_events.o \
	$(OBJ)/tdiffer.o \
	$(OBJ)/forward_new.o \
	$(OBJ)/mantle.o \
	$(OBJ)/mantle_model.o \
	$(OBJ)/model_vector.o \
	$(OBJ)/multiplet2.o \
	$(OBJ)/receiver.o \
	$(OBJ)/rot_ell.o \
	$(OBJ)/singlet.o \
	$(OBJ)/source.o \
	-lpack -lblas -ltmg

make_setup:  coupling.h \
       $(OBJ)/make_setup.o
	$(F90) $(FFLAGS) $(LFLAGS) -o $(BIN)/make_setup \
       $(OBJ)/make_setup.o

$(OBJ)/convert.o: coupling.h $(SRC)/convert.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/convert.o $(SRC)/convert.f

$(OBJ)/get_Q.o: coupling.h $(SRC)/get_Q.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/get_Q.o $(SRC)/get_Q.f

$(OBJ)/decipher_cluster.o: coupling.h $(SRC)/decipher_cluster.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/decipher_cluster.o $(SRC)/decipher_cluster.f

$(OBJ)/dls.o: coupling.h $(SRC)/dls.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/dls.o $(SRC)/dls.f

$(OBJ)/futil.o: coupling.h $(SRC)/futil.f
	$(F90) $(FFLAGS) -c -o   $(OBJ)/futil.o $(SRC)/futil.f

$(OBJ)/get_cmt.o: coupling.h $(SRC)/get_cmt.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/get_cmt.o $(SRC)/get_cmt.f

$(OBJ)/get_events.o: coupling.h $(SRC)/get_events.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/get_events.o $(SRC)/get_events.f

$(OBJ)/get_spectrum.o: coupling.h $(SRC)/get_spectrum.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/get_spectrum.o $(SRC)/get_spectrum.f

$(OBJ)/tdiffer.o: coupling.h $(SRC)/tdiffer.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/tdiffer.o $(SRC)/tdiffer.f

$(OBJ)/deglitch.o: coupling.h $(SRC)/deglitch.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/deglitch.o $(SRC)/deglitch.f

$(OBJ)/invert.o: coupling.h $(SRC)/invert.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/invert.o $(SRC)/invert.f

$(OBJ)/inv_par.o: coupling.h $(SRC)/inv_par.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/inv_par.o $(SRC)/inv_par.f

$(OBJ)/inv_par_fd.o: coupling.h $(SRC)/inv_par_fd.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/inv_par_fd.o $(SRC)/inv_par_fd.f

$(OBJ)/make_setup.o: coupling.h $(SRC)/make_setup.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/make_setup.o $(SRC)/make_setup.f

$(OBJ)/mantle.o: coupling.h $(SRC)/mantle.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/mantle.o $(SRC)/mantle.f

$(OBJ)/mantle_model.o: coupling.h $(SRC)/mantle_model.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/mantle_model.o $(SRC)/mantle_model.f

$(OBJ)/model_vector.o: coupling.h $(SRC)/model_vector.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/model_vector.o $(SRC)/model_vector.f

$(OBJ)/model_vector_fd.o: coupling.h $(SRC)/model_vector_fd.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/model_vector_fd.o $(SRC)/model_vector_fd.f

$(OBJ)/multiplet.o: coupling.h $(SRC)/multiplet.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/multiplet.o $(SRC)/multiplet.f

$(OBJ)/multiplet2.o: coupling.h $(SRC)/multiplet2.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/multiplet2.o $(SRC)/multiplet2.f

$(OBJ)/QL6.o: coupling.h $(SRC)/QL6.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/QL6.o $(SRC)/QL6.f

$(OBJ)/receiver.o: coupling.h $(SRC)/receiver.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/receiver.o $(SRC)/receiver.f

$(OBJ)/reinvert.o: coupling.h $(SRC)/reinvert.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/reinvert.o $(SRC)/reinvert.f

$(OBJ)/rot_ell.o: coupling.h $(SRC)/rot_ell.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/rot_ell.o $(SRC)/rot_ell.f

$(OBJ)/forward.o: coupling.h $(SRC)/forward.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/forward.o $(SRC)/forward.f

$(OBJ)/forward_new.o: coupling.h $(SRC)/forward_new.f90
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/forward_new.o $(SRC)/forward_new.f90

$(OBJ)/comp_gmatrix.o: coupling.h $(SRC)/comp_gmatrix.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/comp_gmatrix.o $(SRC)/comp_gmatrix.f

$(OBJ)/singlet.o: coupling.h $(SRC)/singlet.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/singlet.o $(SRC)/singlet.f

$(OBJ)/source.o: coupling.h $(SRC)/source.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/source.o $(SRC)/source.f

$(OBJ)/test.o: coupling.h $(SRC)/test.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/test.o $(SRC)/test.f

$(OBJ)/variance.o: coupling.h $(SRC)/variance.f
	$(F90) $(FFLAGS)  -c -o   $(OBJ)/variance.o $(SRC)/variance.f


clean:
	rm $(OBJ)/*.o
