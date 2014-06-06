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
flyMOP.o: flyMOP.c \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/limits.h \
  /usr/include/limits.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/machine/limits.h \
  /usr/include/i386/limits.h /usr/include/i386/_limits.h \
  /usr/include/sys/syslimits.h /usr/include/stdlib.h \
  /usr/include/Availability.h /usr/include/AvailabilityInternal.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/unistd.h /usr/include/sys/unistd.h \
  /usr/include/sys/_types/_posix_vdisable.h \
  /usr/include/sys/_types/_seek_set.h /usr/include/sys/_types/_gid_t.h \
  /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h \
  /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h \
  ../lam/error.h ../lam/global.h ../lam/distributions.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h maternal.h ../lam/random.h \
  ../lam/sa.h moves.h score.h fly_io.h /usr/include/ctype.h \
  /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h \
  ../lam/ioTools.h /usr/include/math.h /usr/include/sys/time.h \
  /usr/include/time.h /usr/include/_structs.h \
  /usr/include/sys/_structs.h /usr/include/sys/_types/_clock_t.h \
  solvers.h /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h zygotic.h
fly_io.o: fly_io.c /usr/include/string.h /usr/include/_types.h \
  /usr/include/sys/_types.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h \
  /usr/include/i386/_types.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_null.h /usr/include/sys/_types/_rsize_t.h \
  /usr/include/sys/_types/_errno_t.h /usr/include/sys/_types/_ssize_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/secure/_common.h /usr/include/math.h ../lam/mathLib.h \
  fly_io.h /usr/include/stdio.h /usr/include/sys/_types/_va_list.h \
  /usr/include/sys/_types/_off_t.h /usr/include/secure/_stdio.h \
  /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_wint_t.h \
  ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  ../lam/ioTools.h /usr/include/stdlib.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/sys/time.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_clock_t.h /usr/include/sys/_select.h \
  ../lam/error.h ../lam/distributions.h moves.h ../lam/sa.h maternal.h \
  score.h
integrate.o: integrate.c /usr/include/stdio.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/_types.h \
  /usr/include/sys/_types.h /usr/include/machine/_types.h \
  /usr/include/i386/_types.h /usr/include/sys/_types/_va_list.h \
  /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_null.h \
  /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h \
  /usr/include/secure/_stdio.h /usr/include/secure/_common.h \
  /usr/include/stdlib.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h \
  /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h \
  /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h \
  /usr/include/math.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  /usr/include/string.h /usr/include/sys/_types/_rsize_t.h \
  /usr/include/sys/_types/_errno_t.h /usr/include/strings.h \
  /usr/include/secure/_string.h /usr/include/sys/times.h \
  /usr/include/sys/_types/_clock_t.h ../lam/error.h ../lam/global.h \
  maternal.h ../lam/distributions.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h solvers.h \
  /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h zygotic.h \
  /usr/include/time.h /usr/include/_structs.h \
  /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h score.h fly_io.h \
  /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h moves.h \
  ../lam/sa.h
maternal.o: maternal.c \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  /usr/include/stdlib.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/_types.h \
  /usr/include/sys/_types.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h \
  /usr/include/i386/_types.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/math.h \
  /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/strings.h \
  /usr/include/secure/_string.h /usr/include/secure/_common.h \
  ../lam/error.h ../lam/global.h maternal.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/secure/_stdio.h ../lam/distributions.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h zygotic.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_types/_clock_t.h \
  /usr/include/sys/_types/_time_t.h fly_io.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h moves.h \
  ../lam/sa.h score.h
moves.o: moves.c /usr/include/math.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/stdlib.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h ../lam/error.h \
  ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  moves.h ../lam/sa.h ../lam/distributions.h maternal.h ../lam/random.h \
  fly_io.h /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_clock_t.h /usr/include/sys/_select.h score.h
printscore.o: printscore.c /usr/include/math.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/stdio.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_null.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/stdlib.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h \
  /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h \
  /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h \
  /usr/include/string.h /usr/include/sys/_types/_rsize_t.h \
  /usr/include/sys/_types/_errno_t.h /usr/include/strings.h \
  /usr/include/secure/_string.h /usr/include/unistd.h \
  /usr/include/sys/unistd.h /usr/include/sys/_types/_posix_vdisable.h \
  /usr/include/sys/_types/_seek_set.h /usr/include/sys/_types/_gid_t.h \
  /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h \
  /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h \
  ../lam/error.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  integrate.h /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h maternal.h ../lam/distributions.h \
  score.h fly_io.h /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/time.h /usr/include/_structs.h \
  /usr/include/sys/_structs.h /usr/include/sys/_types/_clock_t.h moves.h \
  ../lam/sa.h solvers.h /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h zygotic.h
savestate.o: savestate.c /usr/include/math.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/stdio.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_null.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/stdlib.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/machine/types.h \
  /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h \
  /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h \
  /usr/include/unistd.h /usr/include/sys/unistd.h \
  /usr/include/sys/_types/_posix_vdisable.h \
  /usr/include/sys/_types/_seek_set.h /usr/include/sys/_types/_gid_t.h \
  /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h \
  /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h \
  ../lam/error.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  moves.h ../lam/sa.h ../lam/distributions.h maternal.h ../lam/random.h
score.o: score.c /usr/include/math.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/stdlib.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/ctype.h \
  /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h \
  /usr/include/string.h /usr/include/sys/_types/_rsize_t.h \
  /usr/include/sys/_types/_errno_t.h /usr/include/strings.h \
  /usr/include/secure/_string.h score.h maternal.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  ../lam/distributions.h fly_io.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_clock_t.h /usr/include/sys/_select.h \
  ../lam/error.h moves.h ../lam/sa.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h solvers.h \
  /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h
scramble.o: scramble.c /usr/include/stdlib.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/_types.h \
  /usr/include/sys/_types.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h \
  /usr/include/i386/_types.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/sys/types.h /usr/include/sys/_types/_blkcnt_t.h \
  /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_gid_t.h \
  /usr/include/sys/_types/_in_addr_t.h \
  /usr/include/sys/_types/_in_port_t.h /usr/include/sys/_types/_ino_t.h \
  /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_key_t.h \
  /usr/include/sys/_types/_nlink_t.h /usr/include/sys/_types/_clock_t.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_useconds_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_isset.h \
  /usr/include/sys/_types/_fd_copy.h \
  /usr/include/sys/_types/_pthread_cond_t.h \
  /usr/include/sys/_types/_pthread_condattr_t.h \
  /usr/include/sys/_types/_pthread_mutex_t.h \
  /usr/include/sys/_types/_pthread_mutexattr_t.h \
  /usr/include/sys/_types/_pthread_once_t.h \
  /usr/include/sys/_types/_pthread_rwlock_t.h \
  /usr/include/sys/_types/_pthread_rwlockattr_t.h \
  /usr/include/sys/_types/_pthread_t.h \
  /usr/include/sys/_types/_pthread_key_t.h \
  /usr/include/sys/_types/_fsblkcnt_t.h \
  /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/unistd.h \
  /usr/include/sys/unistd.h /usr/include/sys/_types/_posix_vdisable.h \
  /usr/include/sys/_types/_seek_set.h /usr/include/sys/select.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_select.h \
  /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h moves.h \
  ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  ../lam/sa.h ../lam/distributions.h maternal.h ../lam/error.h score.h \
  fly_io.h /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h ../lam/ioTools.h /usr/include/math.h \
  /usr/include/sys/time.h /usr/include/time.h /usr/include/_structs.h \
  /usr/include/sys/_structs.h zygotic.h ../lam/random.h scramble.h
solvers.o: solvers.c \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  /usr/include/math.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/stdlib.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/limits.h \
  /usr/include/limits.h /usr/include/machine/limits.h \
  /usr/include/i386/limits.h /usr/include/i386/_limits.h \
  /usr/include/sys/syslimits.h ../lam/error.h ../lam/global.h solvers.h \
  /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h maternal.h \
  ../lam/distributions.h zygotic.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_types/_clock_t.h \
  /usr/include/sys/_types/_time_t.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h
translate.o: translate.c /usr/include/ctype.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/runetype.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_ct_rune_t.h \
  /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h \
  /usr/include/sys/_types/_wint_t.h /usr/include/stdlib.h \
  /usr/include/Availability.h /usr/include/AvailabilityInternal.h \
  /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h \
  /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h \
  /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h \
  /usr/include/i386/signal.h /usr/include/machine/_mcontext.h \
  /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_null.h /usr/include/machine/types.h \
  /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h \
  /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h \
  /usr/include/string.h /usr/include/sys/_types/_rsize_t.h \
  /usr/include/sys/_types/_errno_t.h /usr/include/sys/_types/_ssize_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/secure/_common.h ../lam/error.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  maternal.h /usr/include/stdio.h /usr/include/sys/_types/_va_list.h \
  /usr/include/sys/_types/_off_t.h /usr/include/secure/_stdio.h \
  ../lam/distributions.h score.h fly_io.h ../lam/ioTools.h \
  /usr/include/math.h /usr/include/sys/time.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_clock_t.h /usr/include/sys/_select.h moves.h \
  ../lam/sa.h zygotic.h
unfold.o: unfold.c /usr/include/stdlib.h /usr/include/Availability.h \
  /usr/include/AvailabilityInternal.h /usr/include/_types.h \
  /usr/include/sys/_types.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h \
  /usr/include/i386/_types.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h \
  /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/sys/_types/_ct_rune_t.h /usr/include/sys/_types/_rune_t.h \
  /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_off_t.h \
  /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h \
  /usr/include/secure/_common.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/unistd.h /usr/include/sys/unistd.h \
  /usr/include/sys/_types/_posix_vdisable.h \
  /usr/include/sys/_types/_seek_set.h /usr/include/sys/_types/_gid_t.h \
  /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h \
  /usr/include/sys/_types/_fd_def.h /usr/include/sys/_types/_timespec.h \
  /usr/include/sys/_types/_time_t.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h \
  /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h \
  /usr/include/time.h /usr/include/_structs.h \
  /usr/include/sys/_structs.h /usr/include/sys/_types/_clock_t.h \
  ../lam/error.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  integrate.h /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h maternal.h ../lam/distributions.h \
  solvers.h /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h zygotic.h score.h fly_io.h \
  /usr/include/ctype.h /usr/include/runetype.h \
  /usr/include/sys/_types/_wint_t.h ../lam/ioTools.h /usr/include/math.h \
  /usr/include/sys/time.h moves.h ../lam/sa.h
zygotic.o: zygotic.c /usr/include/ctype.h /usr/include/sys/cdefs.h \
  /usr/include/sys/_symbol_aliasing.h \
  /usr/include/sys/_posix_availability.h /usr/include/runetype.h \
  /usr/include/_types.h /usr/include/sys/_types.h \
  /usr/include/machine/_types.h /usr/include/i386/_types.h \
  /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_ct_rune_t.h \
  /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h \
  /usr/include/sys/_types/_wint_t.h /usr/include/stdio.h \
  /usr/include/Availability.h /usr/include/AvailabilityInternal.h \
  /usr/include/sys/_types/_va_list.h /usr/include/sys/_types/_null.h \
  /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h \
  /usr/include/secure/_stdio.h /usr/include/secure/_common.h \
  /usr/include/stdlib.h /usr/include/sys/wait.h \
  /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h \
  /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h \
  /usr/include/machine/signal.h /usr/include/i386/signal.h \
  /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h \
  /usr/include/mach/i386/_structs.h \
  /usr/include/sys/_types/_sigaltstack.h \
  /usr/include/sys/_types/_ucontext.h \
  /usr/include/sys/_types/_pthread_attr_t.h \
  /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_uid_t.h \
  /usr/include/sys/resource.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/stdint.h \
  /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h \
  /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h \
  /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h \
  /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h \
  /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h \
  /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h \
  /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h \
  /usr/include/machine/endian.h /usr/include/i386/endian.h \
  /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h \
  /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h \
  /usr/include/machine/types.h /usr/include/i386/types.h \
  /usr/include/sys/_types/___offsetof.h /usr/include/sys/_types/_dev_t.h \
  /usr/include/sys/_types/_mode_t.h /usr/include/string.h \
  /usr/include/sys/_types/_rsize_t.h /usr/include/sys/_types/_errno_t.h \
  /usr/include/strings.h /usr/include/secure/_string.h \
  /usr/include/math.h ../lam/error.h ../lam/global.h \
  /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/5.1/include/float.h \
  maternal.h ../lam/distributions.h solvers.h \
  /usr/local/include/cvode/cvode.h \
  /usr/local/include/sundials/sundials_nvector.h \
  /usr/local/include/sundials/sundials_types.h \
  /usr/local/include/sundials/sundials_config.h \
  /usr/local/include/cvode/cvode_band.h \
  /usr/local/include/cvode/cvode_direct.h \
  /usr/local/include/sundials/sundials_direct.h \
  /usr/local/include/sundials/sundials_band.h \
  /usr/local/include/cvode/cvode_spgmr.h \
  /usr/local/include/cvode/cvode_spils.h \
  /usr/local/include/sundials/sundials_iterative.h \
  /usr/local/include/sundials/sundials_spgmr.h \
  /usr/local/include/nvector/nvector_serial.h \
  /usr/local/include/sundials/sundials_dense.h \
  /usr/local/include/sundials/sundials_math.h integrate.h \
  /usr/local/include/gsl/gsl_errno.h /usr/include/errno.h \
  /usr/include/sys/errno.h /usr/local/include/gsl/gsl_types.h \
  /usr/local/include/gsl/gsl_spline.h \
  /usr/local/include/gsl/gsl_interp.h \
  /usr/local/include/gsl/gsl_inline.h zygotic.h /usr/include/time.h \
  /usr/include/_structs.h /usr/include/sys/_structs.h \
  /usr/include/sys/_types/_timespec.h /usr/include/sys/_types/_clock_t.h \
  /usr/include/sys/_types/_time_t.h fly_io.h ../lam/ioTools.h \
  /usr/include/sys/time.h /usr/include/sys/_types/_fd_def.h \
  /usr/include/sys/_types/_suseconds_t.h \
  /usr/include/sys/_types/_fd_setsize.h \
  /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h \
  /usr/include/sys/_types/_fd_isset.h /usr/include/sys/_types/_fd_zero.h \
  /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h moves.h \
  ../lam/sa.h score.h
