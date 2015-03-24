ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)

all: p_resolution
# testsample

p_resolution:p_resolution.C PipPimNpi.o
	$(CC) -o $@ $^


%.o : %.C
	$(CC) -c $^ -o $@


.PHONY:clean
clean:
	rm -f analysis *.o src/*.o

cleanall:
	rm -f analysis *.o src/*.o *.eps *.pdf *.ps
