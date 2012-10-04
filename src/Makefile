all: pairend

pairend: pairend.c fastq.c
	gcc -Wall -fopenmp -fomit-frame-pointer -O2 $+ -o $@

clean:
	-rm pairend

.PHONY: all clean
