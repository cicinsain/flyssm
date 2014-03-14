# Makefile for compiling NSGA-II source code
#CC=/usr/local/bin/gcc-4.8
#LD=/usr/local/bin/gcc-4.8
RM=rm -f
# Amir:
# CFLAGS=-Wall -g -O2 -std=gnu99
# CFLAGS=-Wall -ansi -pedantic -g
NSGA_INCLUDES=-I../fly -I../lam -I../nsga2
OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
MAIN=nsga2_.o
all:$(MAIN)
$(MAIN):$(OBJS)
#	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm
%.o: %.c nsga2.h rand.h
	$(CC) $(CCFLAGS) $(NSGA_INCLUDES)  -c $<
clean:
	$(RM) $(OBJS)
