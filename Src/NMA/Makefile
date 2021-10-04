include arch.make

objects:=NMA.o
all: ../NMA
../NMA: $(objects)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(INFLAGS) $(LDFLAGS) $(LIBS)

%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../NMA *.mod *.o
