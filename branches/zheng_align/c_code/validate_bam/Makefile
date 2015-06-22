CC = gcc
CFLAGS=-O3 -D_GNU_SOURCE
#SAMPATH=/nfs/1000g-work/G1K/work/bin/samtools_dev/samtools
SAMPATH=/nfs/production/reseq-info/work/bin/samtools
LDFLAGS=-L$(SAMPATH) -Wl,-rpath,$(SAMPATH)

all: validate_bam

validate_bam: validate_bam.o
	$(CC) $(CFLAGS) $(LDFLAGS) validate_bam.o -lssl -lbam -lm -lz -o validate_bam

validate_bam.o: validate_bam.c
	$(CC) $(CFLAGS) -I$(SAMPATH) -c validate_bam.c

clean:
	\rm -f *.o validate_bam
