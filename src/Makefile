all: pear

pear: pear.o args.o emp.c statistics.o reader.o  # fastq.o
	gcc -g -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

pear.o: pear.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -O2 $+ -lm -o $@

fastq.o: fastq.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

args.o: args.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

emp.o: emp.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

statistics.o: statistics.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

reader.o: reader.c
	gcc -g -c -Wall -fopenmp -fomit-frame-pointer -lm -O2 $+ -o $@

clean:
	-rm pear pear.o fastq.o args.o statistics.o

.PHONY: all clean
