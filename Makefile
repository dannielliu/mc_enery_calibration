ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)

all: CreateSample CreateSample2 CreateSample_Ks
# testsample



CreateSample:CreateSample.C
	$(CC) $^ -o $@

CreateSample2:CreateSample2.C
	$(CC) $^ -o $@

testsample:testsample.C src/function.o
	$(CC) -lRooFitCore -lRooFit $^ -o $@

% : %.C
	$(CC) $^ -o $@


.PHONY:clean
clean:
	rm -f analysis *.o src/*.o

cleanall:
	rm -f analysis *.o src/*.o *.eps *.pdf *.ps
