# Makefile for compiling NSGA-II source code

# CFLAGS=-I/usr/include/malloc
# Amir:
# CFLAGS=-Wall -g -O2 -std=gnu99
# CFLAGS=-Wall -ansi -pedantic -g
AMOSA_INCLUDES=-I../fly -I../lam -I../amosa -I/usr/include/malloc
OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
AMOSA_FLAGS=-Wno-unused-variable
# MAIN=amosa_rea.o
RM=rm -f

all:$(OBJS)
# $(MAIN):$(OBJS)

%.o: %.c amosa.h
	$(CC) $(CFLAGS) $(AMOSA_FLAGS) $(AMOSA_INCLUDES) -c $<
clean:
	$(RM) $(OBJS)


# all:
# 	gcc $(CFLAGS) -O2 -lm amosa_real.c -o amosa

# clean:
# 	rm amosa