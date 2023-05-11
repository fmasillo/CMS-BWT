VERSION = 1.0
SHELL = /bin/sh

CC = g++
CFLAGS = -std=c++20 -O3 -funroll-loops -march=native -DNDEBUG #-pg
#CFLAGS = -std=c++11 -g -Wall -Wextra -pedantic -funroll-loops -DNDEBUG
OBJS = libsais/src/libsais.o CMS-BWT-functions.o main.o
#LIBS = -ldl

EXEC = clean cms-bwt

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: $(EXEC)

cms-bwt: $(OBJS)
	$(CC) $(CFLAGS) -o cms_bwt $(OBJS)

clean:
	/bin/rm -f *.o libsais/src/*.o
	
