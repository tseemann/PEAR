all: pairend

pairend: pairend.o fastq.o
	gcc -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

pairend.o: pairend.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

fastq.o: fastq.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

clean:
	-rm pairend pairend.o fastq.o

.PHONY: all clean
