###########################################################
# Note: 
# 1) The compliation is seperated to three parts
# 2) Firt part for NMA
#    Second part for PWave
#    Third part for WPpt
###########################################################

###########################################################
# some flags
###########################################################
#for cpu
CC1 = mpiifort
CFLAGS1 = 
INC1 = -I/path-to-fftw/include/
LDFLAGS1 = 
LIBS1 = /path-to-fftw/lib/libfftw3.a

#for gpu
CC2 = pgfortran
CFLAGS2 = -O3 -cuda
INC2 =
LDFLAGS2 =
LIBS2 =

###########################################################
# source files
###########################################################
SOURCES_PWAVE =                     \
	$(wildcard PWave/cor_1D.f90)      \
	$(wildcard PWave/wl_temp.f90)     \
	$(wildcard PWave/temporal_coh.f90)  \
	$(wildcard PWave/wl_spat.f90)  \
	$(wildcard PWave/spatial_coh.f90) \
	$(wildcard PWave/pwave.f90)    
SOURCES_NMA =                   \
	$(wildcard NMA/*.f90)

###########################################################
# object files
###########################################################
OBJ_PWAVE = $(SOURCES_PWAVE:.f90=.o)
OBJ_NMA = $(SOURCES_NMA:.f90=.o)

###########################################################
# executables
###########################################################
all: NMa PWave_GPU run
NMa: $(OBJ_NMA)
	$(CC1) $(LDFLAGS1) $^ -o $@ $(LIBS1)
PWave_GPU: $(OBJ_PWAVE)
	$(CC2) $(CFLAGS2) $(LDFLAGS2) $^ -o $@ $(LIBS2)
run:
	@mkdir ../bin; cp NMa PWave_GPU ../bin; pyinstaller --onefile --distpath=../bin WPpt/WPpt.py
	
###########################################################
# rules for building object files
###########################################################

NMA/%.o: NMA/%.f90
	$(CC1) $(CFLAGS1) $(INC1) -c $< -o $@
	
PWave/%.o: PWave/%.f90
	$(CC2) $(CFLAGS2) $(INC2) -c $< -o $@
	
###########################################################
#clean up
###########################################################
clean:
	rm -fr */*.o *.mod PWave_GPU NMa *.spec ../bin

###########################################################
