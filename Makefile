# CODE VERSION ######################################################

VERSION = 10.0

# executables to make

nsga2Rule=n2
amosaRule=mo
nsga2Folder=nsga2
amosaFolder=amosa

ifeq ($(METHOD),-DAMOSA)
	execFile = fly_amosa
	FLYEXECS = unfold printscore execFile scramble 
	methodRule=$(amosaRule)
	methodFolder=$(amosaFolder)
	methodINCLUDES=-I../amosa
else ifeq ($(METHOD), -DNSGA2)
	execFile = fly_nsga2
	FLYEXECS = unfold printscore execFile scramble
	methodRule=$(nsga2Rule)
	methodFolder=$(nsga2Folder)
	methodINCLUDES=-I../nsga2
endif



# FLAGS FOR -v FOR ALL EXECUTABLES ##################################
# this passes user and host name, compiler and version to the com-
# pile so it can be printed in the -v message of each executable

USRFLAG  = -DUSR=\"$(USER)\"
HOSTFLAG = -DMACHINE=\"$(HOST)\"
COMPFLAG = -DCOMPILER=\"$(CC)\"
FLAGFLAG = -DFLAGS=\"optimized\"
VERSFLAG = -DVERS=\"$(VERSION)\"

VFLAGS = $(USRFLAG) $(HOSTFLAG) $(COMPFLAG) $(FLAGFLAG) $(VERSFLAG) 

# find out about which architecture we're on and set compiler 
# accordingly

AMOSAFLAGS = -DAMOSA
NSGA2FLAGS = -DNSGA2

OSTYPE=linux-gnu
#ifeq ($(OSTYPE),linux-gnu)

# This should recognize most linux machines by looking for the string "linux" in $OSTYPE 
# The gcc version specified manually due to the miss configuration of my Mac
ifneq (,$(findstring linux,$(OSTYPE)))
	MPICC = /usr/local/bin/mpicc
	CC = gcc
	MPIFLAGS = $(CCFLAGS) -DOMPI_CC
	DEBUGFLAGS = $(DEBUGFLAGS) -DOMPI_CC
	PROFILEFLAGS = $(PROFILEFLAGS) -DOMPI_CC
	FLYEXECS = unfold printscore $(execFile) scramble ##fly_sa.mpi
	SUNDIALS = /usr/local
endif

ifeq ($(OSTYPE),osf1)
	CC = cc
endif

#echo $(CC)

# find the compiler and set corresponding flags

ifeq ($(CC),cc)
	CCFLAGS = -std1 -fast -DALPHA_DU -DNDEBUG 
	DEBUGFLAGS = -std1 -O0 -g
	PROFILEFLAGS = -O2 -g1 
	LIBS = -lm -ldxml -lgsl -lgslcblas
	FLIBS = $(LIBS)
	KCC = /bin/kcc
	KFLAGS = -ckapargs=' -scalaropt=3 -optimize=5 -roundoff=3 -arl=4 '
# uncomment 2 lines below if you don't want kcc
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
endif

ifeq ($(CC),icc)
# lucas flags for profiling and debugging 2-6-03
# 	CCFLAGS = -O3 -DNDEBUG
        CCFLAGS = -O3 -xW -tpp7 -ipo 
        PRECFLAGS = -mp -prec_div -pc80 
#       DEBUGFLAGS = -g  -inline_debug_info -O1  -xW -tpp7
        DEBUGFLAGS = -g  -inline_debug_info -O0
#	DEBUGFLAGS = -g
#       PROFILEFLAGS = -prof_dir profile_data -prof_gen -O2  -xW -tpp7
	PROFILEFLAGS = -p -qp -O2 -xW -tpp7
#       USEPROFFLAGS = -prof_use  -prof_dir profile_data  -O3 -xW -tpp7 -ipo -opt_report
	LIBS = -limf -lgsl -lgslcblas
r/local	LIBS = -lm
	FLIBS = -limf -lgsl -lgslcblas -static
	KCC = $(CC)
	KFLAGS = $(CCFLAGS)
        export ICC = "yes"
endif

ifeq ($(CC),gcc)
  	CCFLAGS = -Wall -m64 -O2 -std=gnu99 -DHAVE_SSE2 $(METHOD) 
   	PROFILEFLAGS = -g -pg -O2 -DHAVE_SSE2
	LIBS = -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L$(SUNDIALS)/lib
	FLIBS = -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L$(SUNDIALS)/lib
	KCC = $(CC)
	KFLAGS = $(CCFLAGS)
endif

# debugging?

ifdef DEBUG
	CCFLAGS = $(DEBUGFLAGS)
	FLAGFLAG = -DFLAGS=\"debugging\"
else
	DEBUG = "Off"
endif

ifdef PROFILE
	CCFLAGS = $(PROFILEFLAGS)
	FLAGFLAG = -DFLAGS=\"profiling\"
	KCC = $(CC)
	KFLAGS =
else
	PROFILE = "Off"
endif

ifdef USEPROFILE
        CCFLAGS = $(USEPROFFLAGS)
endif

# export all variables that Makefiles in subdirs need
# 2012 july 25, I needed to add the location of my sundials libraries - A. Crombach

export INCLUDES = -I. -I../lam -I/usr/local/include -I$(SUNDIALS)/include $(methodINCLUDES) -I/usr/include/malloc
export CFLAGS = -std=gnu99 $(CCFLAGS) $(INCLUDES)
export VFLAGS
export CC
export KCC
export MPICC
export KFLAGS
export LIBS
export FLIBS
export MPIFLAGS
export FLYEXECS

export NSGA2FLAGS
export AMOSAFLAGS
export METHOD
export execFile

#define targets

fly: lsa 
	cd fly && $(MAKE)

deps: 
	cd fly && $(MAKE) -f basic.mk Makefile && chmod +w Makefile

lsa: $(methodRule)
	cd lam && make

$(methodRule):
	cd $(methodFolder) && make

clean:
	rm -f core* *.o *.il
	rm -f */core* */*.o */*.il
	rm -f fly/unfold fly/printscore fly/scramble
	rm -f fly/fly_amosa fly/fly_nsga2 #fly/fly_sa.mpi
	#rm -f amosa/*.o

veryclean:
	rm -f core* *.o *.il
	rm -f */core* */*.o */*.il */*.slog */*.pout */*.uout
	rm -f fly/unfold fly/printscore fly/scramble
	rm -f fly/fly_amosa fly/fly_nsga2 #fly/fly_sa.mpi
	rm -f lam/gen_deviates
	rm -f fly/Makefile
	rm -f fly/zygotic.cmp.c

help:
	@echo "make: this is the Makefile for fly code"
	@echo "      always 'make deps' first after a 'make veryclean'"
	@echo ""
	@echo "      the following targets are available:"
	@echo "      lsa:       make object files in the lam directory only"
	@echo "      fly:       compile the fly code (which is in 'fly')"
	@echo "      clean:     gets rid of cores and object files"
	@echo "      veryclean: gets rid of executables and dependencies too"
	@echo ""
	@echo "      your current settings are:"   
	@echo "      compiler:  $(CC)"
	@echo "      flags:     $(CFLAGS)"
	@echo "      debugging: $(DEBUG)"
	@echo "      profiling: $(PROFILE)"
	@echo "      os type:   $(OSTYPE)"
	@echo ""


