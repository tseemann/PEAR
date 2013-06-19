CC     = gcc
CFLAGS = -Wall -fomit-frame-pointer -O3 -lm -march=native
INSTALL_PATH ?= /usr/local/bin
all: pear

pear: pear-pt.o args.o statistics.o reader.o  # fastq.o
	$(CC) $(CFLAGS) $+ -lm -lpthread -o $@

pear-pt.o: pear-pt.c
	$(CC) -c $(CFLAGS) $+ -o $@

args.o: args.c
	$(CC) -c $(CFLAGS) $+ -o $@

statistics.o: statistics.c
	$(CC) -c $(CFLAGS) $+ -o $@

reader.o: reader.c
	$(CC) -c $(CFLAGS) $+ -o $@

clean:
	-rm pear pear-pt.o args.o statistics.o reader.o

install:
	-cp pear $(INSTALL_PATH)

.PHONY: all clean
