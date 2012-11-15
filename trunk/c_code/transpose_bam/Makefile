CC = gcc
CFLAGS=-O3 -D_GNU_SOURCE
SAMPATH=/nfs/1000g-work/G1K/work/bin/samtools_latest/
LDFLAGS=-L$(SAMPATH) -Wl,-rpath,$(SAMPATH)

all: transpose_bam

transpose_bam: transpose_bam.o
	$(CC) $(CFLAGS) $(LDFLAGS) transpose_bam.o -lbam -lm -lz -o transpose_bam

transpose_bam.o: transpose_bam.c
	$(CC) $(CFLAGS) -I$(SAMPATH) -c transpose_bam.c

clean:
	\rm -f *.o transpose_bam
