CC = gcc
CFLAGS=-O3 -D_GNU_SOURCE
HTSPATH=/nfs/production/reseq-info/work/bin/htslib-1.3/
SAMPATH=/nfs/production/reseq-info/work/bin/samtools-1.3/
LDFLAGS=-L$(HTSPATH) -L$(SAMPATH) -Wl,-rpath,$(HTSPATH),-rpath,$(SAMPATH)

all: validate_bam

validate_bam: validate_bam.o
	$(CC) $(CFLAGS) $(LDFLAGS) validate_bam.o -lssl -lcrypto -lbam -lhts -lm -lz -o validate_bam

validate_bam.o: validate_bam.c
	$(CC) $(CFLAGS) -I$(SAMPATH) -I$(HTSPATH) -c validate_bam.c

clean:
	\rm -f *.o validate_bam
