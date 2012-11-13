#ifndef ARGS_H
#define ARGS_H
#include <stdio.h>
#include <getopt.h>

struct user_args 
 {
   char       * fastq_left;
   char       * fastq_right;
   int          min_asm_len;
   int          max_asm_len;
   int          qual_thres;
   int          score_method;
   int          min_overlap;
 };

void usage (void);
int decode_switches (int argc, char * argv[], struct user_args * sw);
#endif
