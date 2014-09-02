# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 


#fly_sa objects
FOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
         ../lam/error.o ../lam/distributions.o ../lam/random.o ../lam/ioTools.o solvers.o score.o ../lam/dSFMT.o ../lam/dSFMT_str_state.o
# serial code
FSOBJ = moves.o ../lam/lsa.o savestate.o flyMOP.o

NSGA2OBJ = ../nsga2/allocate.o ../nsga2/decode.o ../nsga2/fillnds.o ../nsga2/mutation.o ../nsga2/rank.o ../nsga2/auxiliary.o ../nsga2/display.o ../nsga2/initialize.o ../nsga2/nsga2.o ../nsga2/report.o ../nsga2/crossover.o ../nsga2/dominance.o ../nsga2/list.o ../nsga2/problemdef.o ../nsga2/sort.o ../nsga2/crowddist.o ../nsga2/eval.o ../nsga2/merge.o ../nsga2/rand.o ../nsga2/tourselect.o
AMOSAOBJ = ../amosa/amosa_real.o
SSOBJ = ../ss/allocate.o ../ss/evaluate.o ../ss/init.o ../ss/input.o ../ss/linkedlist.o ../ss/local_search.o ../ss/problemdef.o ../ss/recombine.o ../ss/refine.o ../ss/report.o ../ss/sort.o ../ss/ss.o ../ss/ssTools.o ../ss/stats.o ../ss/update.o 

ESSOBJ = ../ess/ess.o ../ess/essAllocate.o ../ess/essEvaluate.o ../ess/essGoBeyond.o ../ess/essIO.o ../ess/essInit.o ../ess/essLocalSearch.o ../ess/essProblem.o ../ess/essRand.o ../ess/essRecombine.o ../ess/essSort.o ../ess/essStats.o ../ess/essTools.o

ifeq ($(METHOD), -DAMOSA)
	METHODOBJ = $(AMOSAOBJ)
else ifeq ($(METHOD), -DNSGA2)
	METHODOBJ = $(NSGA2OBJ)
else ifeq ($(METHOD), -DSS)
	METHODOBJ = $(SSOBJ)
else ifeq ($(METHOD), -DESS)
	METHODOBJ = $(ESSOBJ)
endif

# parallel code
#FPOBJ = moves-mpi.o fly_sa-mpi.o ../lam/lsa-mpi.o savestate-mpi.o


#printscore objects
POBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
       ../lam/error.o ../lam/distributions.o ../lam/random.o ../lam/ioTools.o solvers.o score.o printscore.o ../lam/dSFMT.o ../lam/dSFMT_str_state.o

#unfold objects
UOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
	 ../lam/error.o ../lam/distributions.o ../lam/random.o ../lam/ioTools.o solvers.o score.o unfold.o ../lam/dSFMT.o ../lam/dSFMT_str_state.o

#scramble objects
SOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
	 ../lam/error.o ../lam/distributions.o ../lam/random.o ../lam/ioTools.o solvers.o score.o scramble.o ../lam/dSFMT.o ../lam/dSFMT_str_state.o

SOURCES = `ls *.c`

#Below here are the rules for building things

all: $(FLYEXECS)

# special cases: dependencies and flags for individual .c files

flyMOP.o: flyMOP.c
	$(CC) -c $(CFLAGS) $(VFLAGS) flyMOP.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS) printscore.c

scramble.o: scramble.c
	$(CC) -c $(CFLAGS) $(VFLAGS) scramble.c

unfold.o: unfold.c
	$(CC) -c $(CFLAGS) $(VFLAGS) unfold.c

zygotic.o: zygotic.c
	$(CC) -c $(CFLAGS) $(KFLAGS) zygotic.c


# parallel stuff

#fly_sa-mpi.o: fly_sa.c
#	$(MPICC) -c -o fly_sa-mpi.o $(MPIFLAGS) $(CFLAGS) $(VFLAGS) fly_sa.c

# lsa-mpi.o: ../lam/lsa.c
# 	$(MPICC) -c -o ../lam/lsa-mpi.o $(MPIFLAGS) $(CFLAGS) ../lam/lsa.c

# moves-mpi.o: moves.c 
# 	$(MPICC) -c -o moves-mpi.o $(MPIFLAGS) $(CFLAGS) moves.c

# savestate-mpi.o: savestate.c
# 	$(MPICC) -c -o savestate-mpi.o $(MPIFLAGS) $(CFLAGS) savestate.c

# executable targets: serial ...

$(execFile): $(FOBJ) $(FSOBJ)
	$(CC) -o $(execFile) $(CFLAGS) $(LDFLAGS) $(FOBJ) $(FSOBJ) $(METHODOBJ) $(FLIBS) 

printscore: $(POBJ)
	$(CC) -o printscore $(CFLAGS) $(LDFLAGS) $(POBJ) $(LIBS) 

unfold: $(UOBJ)
	$(CC) -o unfold $(CFLAGS) $(LDFLAGS) $(UOBJ) $(LIBS) 

scramble: $(SOBJ)
	$(CC) -o scramble $(CFLAGS) $(LDFLAGS) $(SOBJ) $(LIBS) 

# ... and parallel

#fly_sa.mpi: $(FOBJ) $(FPOBJ)
#	$(MPICC) -o fly_sa.mpi $(CFLAGS) $(LDFLAGS) $(FOBJ) $(FPOBJ) $(FLIBS)

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(INCLUDES) -M ${SOURCES} >> $@
	chmod -w $@

#Automatically generated dependencies list#
fly_io.o: fly_io.c /usr/include/string.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/xlocale.h /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 ../lam/mathLib.h fly_io.h /usr/include/stdio.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/ctype.h /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h ../lam/ioTools.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/sys/time.h ../lam/global.h \
 ../lam/error.h ../lam/distributions.h moves.h ../lam/sa.h maternal.h \
 score.h maternal.h
flyMOP.o: flyMOP.c /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/limits.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/syslimits.h \
 /usr/include/limits.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/posix1_lim.h \
 /usr/include/bits/local_lim.h /usr/include/linux/limits.h \
 /usr/include/bits/posix2_lim.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h ../lam/distributions.h integrate.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h maternal.h ../lam/random.h ../lam/sa.h \
 moves.h maternal.h score.h fly_io.h /usr/include/ctype.h \
 ../lam/ioTools.h /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/sys/time.h ../lam/global.h ../lam/error.h solvers.h \
 ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
integrate.o: integrate.c /usr/include/stdio.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/sys/times.h \
 ../lam/error.h ../lam/global.h maternal.h ../lam/distributions.h \
 integrate.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h \
 ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h score.h \
 fly_io.h /usr/include/ctype.h ../lam/ioTools.h /usr/include/sys/time.h \
 ../lam/global.h ../lam/error.h moves.h ../lam/sa.h maternal.h
maternal.o: maternal.c \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/ctype.h /usr/include/xlocale.h /usr/include/string.h \
 ../lam/error.h ../lam/global.h maternal.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 ../lam/distributions.h integrate.h /usr/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_spline.h \
 /usr/include/gsl/gsl_interp.h /usr/include/gsl/gsl_inline.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h fly_io.h \
 ../lam/ioTools.h /usr/include/sys/time.h ../lam/global.h ../lam/error.h \
 moves.h ../lam/sa.h maternal.h score.h
moves.o: moves.c /usr/include/math.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 moves.h ../lam/sa.h ../lam/distributions.h maternal.h ../lam/random.h \
 fly_io.h /usr/include/ctype.h ../lam/ioTools.h /usr/include/sys/time.h \
 ../lam/global.h ../lam/error.h moves.h score.h maternal.h fly_io.h
printscore.o: printscore.c /usr/include/math.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 integrate.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h maternal.h ../lam/distributions.h \
 maternal.h score.h fly_io.h /usr/include/ctype.h ../lam/ioTools.h \
 /usr/include/sys/time.h ../lam/global.h ../lam/error.h moves.h \
 ../lam/sa.h score.h solvers.h ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
savestate.o: savestate.c /usr/include/math.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/stdlib.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 moves.h ../lam/sa.h ../lam/distributions.h maternal.h ../lam/random.h
score.o: score.c /usr/include/math.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/ctype.h /usr/include/xlocale.h /usr/include/string.h \
 score.h maternal.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 ../lam/distributions.h fly_io.h ../lam/ioTools.h /usr/include/sys/time.h \
 ../lam/global.h ../lam/error.h moves.h ../lam/sa.h maternal.h \
 integrate.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h \
 ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h
scramble.o: scramble.c /usr/include/stdlib.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h moves.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 ../lam/sa.h ../lam/distributions.h maternal.h ../lam/error.h score.h \
 maternal.h fly_io.h /usr/include/ctype.h ../lam/ioTools.h \
 /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/sys/time.h ../lam/global.h ../lam/error.h moves.h score.h \
 zygotic.h /usr/include/sys/resource.h /usr/include/bits/resource.h \
 ../lam/random.h scramble.h
solvers.o: solvers.c \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 /usr/include/math.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/limits.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/syslimits.h \
 /usr/include/limits.h /usr/include/bits/posix1_lim.h \
 /usr/include/bits/local_lim.h /usr/include/linux/limits.h \
 /usr/include/bits/posix2_lim.h ../lam/error.h ../lam/global.h solvers.h \
 ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h maternal.h \
 ../lam/distributions.h maternal.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h integrate.h /usr/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_spline.h \
 /usr/include/gsl/gsl_interp.h /usr/include/gsl/gsl_inline.h
translate.o: translate.c /usr/include/ctype.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/xlocale.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/string.h ../lam/error.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h maternal.h \
 /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h \
 /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 ../lam/distributions.h score.h fly_io.h ../lam/ioTools.h \
 /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/sys/time.h ../lam/global.h ../lam/error.h moves.h \
 ../lam/sa.h maternal.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h
unfold.o: unfold.c /usr/include/stdlib.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/sys/types.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h \
 integrate.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h maternal.h ../lam/distributions.h \
 maternal.h solvers.h ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h zygotic.h score.h fly_io.h \
 /usr/include/ctype.h ../lam/ioTools.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/sys/time.h ../lam/global.h \
 ../lam/error.h moves.h ../lam/sa.h score.h fly_io.h moves.h
zygotic.o: zygotic.c /usr/include/ctype.h /usr/include/features.h \
 /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h \
 /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/xlocale.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stddef.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/string.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h ../lam/error.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/float.h maternal.h \
 ../lam/distributions.h solvers.h ../sundials/include/cvode/cvode.h \
 ../sundials/include/sundials/sundials_nvector.h \
 ../sundials/include/sundials/sundials_types.h \
 ../sundials/include/sundials/sundials_config.h \
 ../sundials/include/cvode/cvode_band.h \
 ../sundials/include/cvode/cvode_direct.h \
 ../sundials/include/sundials/sundials_direct.h \
 ../sundials/include/sundials/sundials_band.h \
 ../sundials/include/cvode/cvode_spgmr.h \
 ../sundials/include/cvode/cvode_spils.h \
 ../sundials/include/sundials/sundials_iterative.h \
 ../sundials/include/sundials/sundials_spgmr.h \
 ../sundials/include/nvector/nvector_serial.h \
 ../sundials/include/sundials/sundials_dense.h \
 ../sundials/include/sundials/sundials_math.h integrate.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h fly_io.h ../lam/ioTools.h \
 /usr/include/sys/time.h ../lam/global.h ../lam/error.h moves.h \
 ../lam/sa.h maternal.h score.h
