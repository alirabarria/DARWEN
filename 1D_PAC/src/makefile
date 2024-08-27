# Makefile for PAC program

FC      = gfortran
FFLAGS  = -O3 -march=native -fno-second-underscore -Wall -fbackslash

BINDIR  = .
EXEC    = pac

OBJECTS = pac_data.o \
           pac.o \
           read_input.o \
           solve_continuity_equation.o \
           compute_newtonraphson.o \
           compute_dlsodes.o \
           compute_cranknicolson.o \
           get_coefficients.o \
           evaluate_derjac.o \
           write_results.o \
           solid_body_rotation.o \
           dlsodes.o

COMMON  = pac.common \
           dimension.common

# Add -fallow-argument-mismatch only for dlsodes.o
dlsodes.o: FFLAGS += -fallow-argument-mismatch

all: $(EXEC)
	echo "             Done!"

$(EXEC) : $(OBJECTS) $(COMMON)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(BINDIR)/$(EXEC)

clean:
	rm *.o
	rm *.mod
	rm $(EXEC)
