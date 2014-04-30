
RM=rm -f
# CC=gcc-4.8
# CC=clang
# CFLAGS=-Wall -g -std=c99 -O2
# -fsanitize=address -fno-omit-frame-pointer
# -O1 -fsanitize=address -fno-omit-frame-pointer
# -O2 -fsanitize=address -fno-omit-frame-pointer
OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
SS_INCLUDES=-I../fly -I../lam -I../ss
MAIN=_ss.o
all:$(MAIN)
$(MAIN):$(OBJS)
	# $(CC) $(CFLAGS) $(OBJS) -o $(MAIN) -lm
	# $(CC) $(CCFLAGS) $(SS_INCLUDES)  -c $<
%.o: %.c ss.h rand.h
	$(CC) $(CFLAGS) $(SSFLAGS) $(SS_INCLUDES) -c $<
clean:
	$(RM) $(OBJS)
