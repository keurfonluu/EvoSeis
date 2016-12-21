#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!  Copyright (c) 2015 by
#!  Magnitude, France  and  MINES ParisTech, France
#!  All rights reserved.
#!
#!  This software is furnished under a license and may be used and copied
#!  only in  accordance with  the  terms  of  such  license  and with the
#!  inclusion of the above copyright notice. This software or  any  other
#!  copies thereof may not be provided or otherwise made available to any
#!  other person.  No title to and ownership of  the  software is  hereby
#!  transferred.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#         Makefile for Fortran programs
#
# ======================================================================
# Declarations for compiler (comment or decomment as necessary)
# ======================================================================
#
# ---- Gnu complier (gcc)
FC      :=  gfortran
#
# ---- Gnu MPI compiler
# FC      :=  mpif90
#
# flags for debugging or for maximum performance or MPI, comment as necessary
#
# ---- Option for gfortran compiler
#FFLAGS  :=  -Ofast
#FFLAGS  :=  -O3 -ffree-line-length-none
#FFLAGS  :=  -O3 -ffree-line-length-none -Wall -Wextra -fbounds-check
#FFLAGS  :=  -O3 -ffree-line-length-none -cpp -Ddo_mpi
#FFLAGS  :=  -O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -flto
FFLAGS  :=  -O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -flto -fbackslash -fopenmp
#
# ======================================================================
# Declarations of executables to be compiled and various dependances
# ======================================================================
# Name of executable
TARGET1 :=  evoloc.exe
TARGET2	:=	evowave.exe

# Directories
SRCDIR  :=  src
OBJDIR  :=  obj
BINDIR	:=	bin
MAIN    :=  $(SRCDIR)
EIKDIR  :=  $(SRCDIR)/fteik3d
FFTDIR	:=	$(SRCDIR)/fftpack

# Link objects to create executable (tab required on second line)
OBJS1   :=  $(OBJDIR)/FTeik3d_2.0.o \
						$(OBJDIR)/forlab.o \
						$(OBJDIR)/fftpack5.o \
            $(OBJDIR)/fftpack.o \
						$(OBJDIR)/optimizers.o \
            $(OBJDIR)/evoseis.o \
            $(OBJDIR)/evoloc.o \

OBJS2   :=  $(OBJDIR)/FTeik3d_2.0.o \
						$(OBJDIR)/forlab.o \
						$(OBJDIR)/fftpack5.o \
						$(OBJDIR)/fftpack.o \
						$(OBJDIR)/optimizers.o \
						$(OBJDIR)/evoseis.o \
						$(OBJDIR)/evowave.o \

# These routines depend on include file - recompile if include modified
ALL		  :=	$(TARGET1) $(TARGET2)
all: $(ALL)
evoloc: $(TARGET1)
evowave: $(TARGET2)

# ======================================================================
# General rules, these should not require modification
# General rules for building ".o" objects from fortran programs or subroutines
# ======================================================================

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)
	mkdir -p $(BINDIR)/gui
	cp $(SRCDIR)/*.py $(BINDIR)/gui/

$(OBJDIR)/%.o: $(EIKDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $^ -o $@ -J$(OBJDIR)

$(OBJDIR)/%.o: $(FFTDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -w -c $^ -o $@ -J$(OBJDIR)

$(OBJDIR)/%.o: $(MAIN)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $^ -o $@ -J$(OBJDIR)

$(TARGET1): $(OBJS1) | $(BINDIR)
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJS1)

$(TARGET2): $(OBJS2) | $(BINDIR)
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJS2)

# Utilities

.PHONY: all evoloc evowave clean veryclean

clean:
	rm -rf $(ALL) $(OBJDIR) $(BINDIR)

veryclean:
	rm -rf $(ALL) $(OBJDIR) $(BINDIR) output
# ======================================================================
# That's all
# ======================================================================
