all: pear

pear: pear.o fastq.o args.o emp.c statistics.o
	gcc -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

pear.o: pear.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

fastq.o: fastq.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

args.o: args.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

emp.o: emp.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

statistics.o: statistics.c
	gcc -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

clean:
	-rm pear pear.o fastq.o args.o statistics.o

.PHONY: all clean
