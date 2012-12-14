all: pairend

pairend: pairend.o fastq.o bintest.o args.o emp.c statistics.o
	gcc -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

pairend.o: pairend.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

fastq.o: fastq.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

bintest.o: bintest.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

args.o: args.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

emp.o: emp.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

statistics.o: statistics.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

clean:
	-rm pairend pairend.o fastq.o args.o bintest.o statistics.o

.PHONY: all clean
