CC = gcc
CFLAGS=-O3 -Wall
SAMPATH=/nfs/1000g-work/G1K/work/bin/git_clone_samtools
LDFLAGS=-L$(SAMPATH) -Wl,-rpath,$(SAMPATH)

all: validate_bam

validate_bam: validate_bam.o
	$(CC) $(CFLAGS) $(LDFLAGS) validate_bam.o -lssl -lbam -lm -lz -o validate_bam

validate_bam.o: validate_bam.c
	$(CC) $(CFLAGS) -I$(SAMPATH) -c validate_bam.c

clean:
	\rm -f *.o validate_bam
