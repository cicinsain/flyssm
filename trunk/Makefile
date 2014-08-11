# CODE VERSION ######################################################

VERSION = 10.0

# executables to make

FLYEXECS = unfold printscore fly_sa scramble 

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

#OSTYPE=linux-gnu
#ifeq ($(OSTYPE),linux-gnu)

# To make 32 or 64 bits version, (I think SSm can be run only with 32 bit version of Matlab) you can set it here, or just override it by executing 
# make BITS=32


#BITS = 32
BITS = 64

# This should recognize most linux machines by looking for the string "linux" in $OSTYPE 
ifneq (,$(findstring linux,$(OSTYPE)))
	MPICC = /usr/lib64/openmpi/bin/mpicc
	CC = gcc
	MPIFLAGS = $(CCFLAGS) -DMPI
	DEBUGFLAGS = $(DEBUGFLAGS) -DMPI
	PROFILEFLAGS = $(PROFILEFLAGS) -DMPI
	FLYEXECS = unfold printscore fly_sa scramble

	ifeq ($(BITS),32)
		#32 bit
		SUNDIALS = $(HOME)/My/lib32/sundials-2.5.0
		GSL = $(HOME)/My/lib32/gsl-1.13
		LGSL = -L$(GSL)/lib
		MATLAB_HOME = $(HOME)/soft/MATLAB/R2012a
		M = m32
		ARCH = glnx86
	endif

	ifeq ($(BITS),64)
		#64 bit 
		SUNDIALS = $(HOME)/My/lib64/sundials-2.5.0
		MATLAB_HOME = $(HOME)/soft/MATLAB/R2013a
		M = m64
		ARCH = glnxa64
		LGSL = 
	endif

	LSUNDIALS = -L$(SUNDIALS)/lib 
endif

#ifeq ($(OSTYPE),osf1)
#	CC = cc
#endif

# find the compiler and set corresponding flags

#ifeq ($(CC),cc)
#	CCFLAGS = -std1 -fast -DALPHA_DU -DNDEBUG
#	DEBUGFLAGS = -std1 -O0 -g
#	PROFILEFLAGS = -O2 -g1 
#	LIBS = -lm -ldxml -lgsl -lgslcblas
#	FLIBS = $(LIBS)
#	KCC = /bin/kcc
#	KFLAGS = -ckapargs=' -scalaropt=3 -optimize=5 -roundoff=3 -arl=4 '
# uncomment 2 lines below if you don't want kcc
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
#endif

#ifeq ($(CC),icc)
# lucas flags for profiling and debugging 2-6-03
# 	CCFLAGS = -O3 -DNDEBUG
#        CCFLAGS = -O3 -xW -tpp7 -ipo 
#        PRECFLAGS = -mp -prec_div -pc80 
#       DEBUGFLAGS = -g  -inline_debug_info -O1  -xW -tpp7
#        DEBUGFLAGS = -g  -inline_debug_info -O0
#	DEBUGFLAGS = -g
#       PROFILEFLAGS = -prof_dir profile_data -prof_gen -O2  -xW -tpp7
#	PROFILEFLAGS = -p -qp -O2 -xW -tpp7
#       USEPROFFLAGS = -prof_use  -prof_dir profile_data  -O3 -xW -tpp7 -ipo -opt_report
#	LIBS = -limf -lgsl -lgslcblas
#	LIBS = -lm
#	FLIBS = -limf -lgsl -lgslcblas -static
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
 #       export ICC = "yes"
#endif

ifeq ($(CC),gcc)
#  	CCFLAGS = -Wall -m64 -O2 -std=gnu99 -DHAVE_SSE2 
  	CCFLAGS = -g3 -O0 -std=gnu99 -DHAVE_SSE2 -$(M) -fPIC -DPIC 
#  	CCFLAGS = -g3 -O0 -std=gnu99 -DHAVE_SSE2 -fPIC -DPIC #64 bit
   	PROFILEFLAGS = -g -pg -O2 -DHAVE_SSE2
	LIBS = -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial $(LSUNDIALS) $(LGSL)
	FLIBS = $(LIBS)
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

export INCLUDES = -I. -I/users/jjaeger/dcicin/NetBeansProjects/SSm/trunk/util -I/users/jjaeger/dcicin/NetBeansProjects/SSm/trunk/fly -I/usr/local/include -I$(MATLAB_HOME)/extern/include -I$(SUNDIALS)/include
#-I$(HOME)/local/include
export CFLAGS = $(CCFLAGS) $(INCLUDES) $(CLFLAGS)
#export CLFLAGS = 
#-L$(HOME)/local/lib
export VFLAGS
export CC
export KCC
export MPICC
export KFLAGS
export LIBS
export FLIBS
export MPIFLAGS
export FLYEXECS

#define targets
.PHONY: fly util

ssm:	objects
	$(MATLAB_HOME)/bin/mex -$(ARCH) -g $(INCLUDES) -output ggn $(LIBS) RootOfAllEvol.c ./selected/*.o 	

objects: util
	 cd fly && $(MAKE) objects

fly:	util
	cd fly && $(MAKE)

util:
	cd util && $(MAKE)

deps: 
	cd fly && $(MAKE) -f basic.mk Makefile && chmod +w Makefile

clean:
	rm -f core* *.o *.il
	rm -f */core* util/*.o fly/*.o */*.il
	rm -f fly/unfold fly/printscore fly/scramble util/gen_deviates
	rm -f fly/fly_sa fly/fly_sa.mpi

veryclean:	clean
	rm -f */*.slog */*.pout */*.uout
	rm -f fly/Makefile
	rm -f fly/zygotic.cmp.c

help:
	@echo "make: this is the Makefile for fly code"
	@echo "      always 'make deps' first after a 'make veryclean'"
	@echo ""
	@echo "      the following targets are available:"
	@echo "      util:      make object files in the util directory only"
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

