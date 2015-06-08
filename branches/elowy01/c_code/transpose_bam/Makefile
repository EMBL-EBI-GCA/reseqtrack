CC = gcc
CFLAGS=-O3 -D_GNU_SOURCE
HTSPATH=/nfs/1000g-work/G1K/work/bin/htslib/
SAMPATH=/nfs/1000g-work/G1K/work/bin/samtools/
LDFLAGS=-L$(HTSPATH) -L$(SAMPATH) -Wl,-rpath,$(HTSPATH),-rpath,$(SAMPATH)

all: transpose_bam

transpose_bam: transpose_bam.o
	$(CC) $(CFLAGS) $(LDFLAGS) transpose_bam.o -lhts -lbam -lm -lz -lpthread -o transpose_bam

transpose_bam.o: transpose_bam.c
	$(CC) $(CFLAGS) -I$(SAMPATH) -c transpose_bam.c

clean:
	\rm -f *.o transpose_bam
