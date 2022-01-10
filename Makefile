include arch.make
ODIR = ../bin
dummy_build_folder := $(shell mkdir -p $(ODIR))

objects1:=NMA.o
all: ../bin/NMA
../bin/NMA: $(objects1)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(INFLAGS) $(LDFLAGS) $(LIBS)

objects2:=PWave.o
all: ../bin/PWave
../bin/PWave: $(objects2)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(INFLAGS) $(LDFLAGS) $(LIBS)
	
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $< $(INFLAGS)
	touch $*.mod
clean:
	rm -f ../NMA ../PWave *.mod *.o
