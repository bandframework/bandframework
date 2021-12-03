#***********************************************************************
#
#    Copyright 2018, I.J. Thompson
#
#    This file is part of FRESCOX.
#
#    FRESCOX is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    FRESCOX is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FRESCOX. If not, see <http://www.gnu.org/licenses/>.
#
#    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
#    LICENSE
#
#    The precise terms and conditions for copying, distribution and
#    modification are contained in the file COPYING.
#
#***********************************************************************
#	frescox: for Fortran90   Version frxy 7a
#
#	Ian J. Thompson,     I-Thompson@llnl.gov
#	September 2019
#
#
#MACH = $(MACHTYPE)
#MACH = `arch`
####------------ Choose appropriate architecture
####------------ by removing '# ' comment prefix on one line,
####------------ or else 'setenv MACH ...' in your shell startup
####------------ (if uncertain: choose std)
# MACH=alpha		# std alpha OS
# MACH=cray		# Cray with scilib
# MACH=ev5		# alpha ev5 with dxml
# MACH=ev6		# alpha ev6 with dxml
# MACH=hp		# hp workstation
# MACH=i386		# intel mac g95
# MACH=i686		# linux pgf90 with lapack and blas
# MACH=i686-g		# linux pgf90 with lapack and blas and -g
# MACH=i86pc		# Solaris pc on intel (if exists?)
# MACH=ipsc		# Intel iPSC/860 hypercube
# MACH=intel		# linux intel ifc
# MACH=linux		# linux generic f90
# MACH=origin		# SGI origin
# MACH=rs6k		# IBM rs6000
# MACH=salf		# Windows Salford f90
# MACH=std		# standard unix
# MACH=sun4-g		# Solaris sparc with -g
# MACH=sun4-spl		# Solaris sparc with sun performance lib
# MACH=sun4		# Solaris sparc
# MACH=ultra2		# Solaris ultra2 sparc
# MACH=vax		# Vax vms
# MACH=xlf		# IBM xlf90, on mac

# MACH=x86_64           # AMD 64 linux pgf90 without lapack and blas
# MACH=g95              # Gnu fortran g95 without lapack and blas
# MACH=altix            # SGI Altix 3700 IA64 arch with ifort v8.1

####------------ Choose precompiled blas or blas source
#                with LBLAS and BLAS in your architecture's *.def file

FC= gfortran -Wall
F90C = $(FC)
FORM =  -w  -cpp
PARALLEL = frxx15
BINSI = frxx4.$(O) frxx7a.$(O)
BINS = frxx0.$(O) $(BINSI)


#### Without nag library:
#
NAG = nagstub.$(O)

#### Without blas library:
BLAS = lapack.$(O)
#LBLAS = -llapack -lblas -lg2c -lm

FFLAGS = -O2
FFLAGS1 = -O1
# FFLAGS = -g -O2
# FFLAGS1 = -g -O1
LFLAGS  =

O = o
E =
FGO = $(F90C) $(FORM) $(FFLAGS) -o fr2nl
LINK = $(F90C) $(FFLAGS) -o frescox $(OBJ) $(LIBS) $(LFLAGS) $(LBLAS)
LINKS = $(F90C) $(FFLAGS) -o sfrescox $(OBJSRCH) $(LIBS) $(LFLAGS) $(LBLAS)
# INST = install sfrescox frescox $(FRESCOXLIB)
INST = cp  sfrescox frescox fr2nl /content/Frescoxex ; chmod a+rx /content/Frescoxex  /content/Frescoxex 
RM = rm -f

# (or edit in the required file if your 'make' does not do 'include's)

LOCAL = f14arm64
SEARCH = minuit-cern.$(O)

OBJS= globx7.$(O) version_number.$(O) \
	frxx1.$(O) frxx2.$(O) frxx3.$(O) frxx3b.$(O) $(BINS) frxx5.$(O) frxx6.$(O)  \
	erwinrc4.$(O) frxx7b.$(O) frxx8.$(O) frxx9.$(O) frxx10.$(O) frxx11.$(O) bincc.$(O)\
	frxx13.$(O) frxx16.$(O) frxx17.$(O) frxx21.$(O) lagmesh.$(O) \
	frxx27.$(O) frxx22.$(O) \
	erwinmpi6.$(O) cdc.$(O) $(SEARCH) \
	$(PARALLEL).$(O) $(NAG) $(BLAS) $(LOCAL).$(O)

OBJ= $(OBJS) fresco.$(O)
OBJSRCH= $(OBJS) sfresco.$(O)

target: $(HELPERS) frescox$(E) sfrescox$(E)  fr2nl$(E)

fr2nl$(E): fr2nl.f
	$(F90C) $(FORM) -o fr2nl -O1 fr2nl.f

frescox$(E): $(OBJ)
	$(LINK)

sfrescox$(E): $(OBJSRCH)
	$(LINKS)

install:  target
	$(INST)

objs: $(OBJ)

# f14$(MACH).$(O): globx7.$(O) f14$(MACH).f
version_number.$(O): version_number.f version.f
fresco.$(O): globx7.$(O) fresco.f
sfresco.$(O): globx7.$(O) sfresco.f
frxx0.$(O): globx7.$(O) frxx0.f $(LOCAL).$(O)
frxx1.$(O): globx7.$(O) frxx1.f coulomb.f veffpot.f makeset.f bpmfus.f solveset.f usescatwf.f
frxx2.$(O): globx7.$(O) frxx2.f
frxx3.$(O): globx7.$(O) frxx3.f
frxx3b.$(O): globx7.$(O) frxx3b.f
erwinrc4.$(O): globx7.$(O) erwinrc4.f
erwinmpi6.$(O): globx7.$(O) erwinmpi6.f pzladistf.f
frxx4.$(O): globx7.$(O) frxx4.f
frxx5.$(O): globx7.$(O) frxx5.f
bincc.$(O): globx7.$(O) bincc.f
frxx6.$(O): globx7.$(O) frxx6.f
frxx7a.$(O): globx7.$(O) frxx7a.f
frxx7b.$(O): globx7.$(O) frxx7b.f
frxx8.$(O): globx7.$(O) frxx8.f
frxx9.$(O): globx7.$(O) frxx9.f
frxx10.$(O): globx7.$(O) frxx10.f
frxx11.$(O): globx7.$(O) frxx11.f
$(LOCAL).$(O): globx7.$(O) $(LOCAL).f
frxx16.$(O): globx7.$(O) frxx16.f
frxx17.$(O): globx7.$(O) frxx17.f boundstate.f usescatwf.f
frxx21.$(O): globx7.$(O) frxx21.f boundstate.f usescatwf.f
frxx22.$(O): globx7.$(O) frxx22.f
frxx27.$(O): globx7.$(O) frxx27.f  usescatr.f
cdc.$(O): globx7.$(O) cdc.f
lagmesh.$(O): lagmesh.f frxx21.$(O)
minuit-cern.$(O): minuit-cern.f # d506cm.inc

F90FLAGS = $(FFLAGS)

.f.obj:
	echo Compile $<
	$(F90C) $(FORM) $(FFLAGS)  $<
.f.o:
	$(F90C) $(FORM) -c $(FFLAGS)  $<
#  #inserted the next two lines for mac source files.
#  .F.o:
#  	$(F90C) $(FORM) -c $(FFLAGS)  $<

frxx1.o: frxx1.f
	$(F90C) $(FORM) $(FFLAGS1) -c frxx1.f
minuit-cern.o: minuit-cern.f
	$(F90C) $(FORM) $(FFLAGS1) -c minuit-cern.f


#########
clean :
	$(RM) *.$(O)
	$(RM) core
	$(RM) frescox$(E)
	$(RM) sfrescox$(E)
	$(RM) fr2nl$(E)
	$(RM) *.M
	$(RM) *.mod
	$(RM) *.MIX
