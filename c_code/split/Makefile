CC = gcc
CFLAGS=-O3
#ZPATH=/nfs/1000g-work/G1K/work
ZPATH=/nfs/production/reseq-info/work
LDFLAGS=-L$(ZPATH)/lib -Wl,-rpath,$(ZPATH)/lib

all: split

split: split.o
	$(CC) $(CFLAGS) $(LDFLAGS) split.o -lz -o split

split.o: split.c
	$(CC) $(CFLAGS) -I$(ZPATH)/include -c split.c

clean:
	\rm -f *.o split
