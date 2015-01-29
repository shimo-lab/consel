CC = gcc
CFLAGS = -g -Wall
#CFLAGS = -O3 -Wall
LDLIBS = -lm


targets= makermt catpv catci catrep treeass consel randrep catmt seqmt makerep catass
objects= misc.o rand.o tree.o freadmat.o opt.o mt19937.o
bindir=../bin

all: $(targets)

makermt: misc.o rand.o freadmat.o mt19937.o
catpv: misc.o
treeass: tree.o misc.o
consel: misc.o rand.o opt.o mt19937.o
catmt: misc.o rand.o mt19937.o
seqmt: freadmat.o misc.o

makerep: misc.o rand.o freadmat.o  mt19937.o
catci: misc.o
catrep: misc.o
catass: misc.o
randrep: misc.o rand.o  mt19937.o

test: rand.o misc.o mt19937.o
tabpv1: misc.o
tabpv2: misc.o
tabpv4: misc.o

clean:
	rm -f $(targets)
	rm -f $(objects)

install:
	mkdir -p $(bindir)
	cp $(targets) $(bindir)
