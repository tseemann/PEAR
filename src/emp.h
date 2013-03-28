#ifndef EMP_H
#define EMP_H

#include "fastq.h"

struct emp_freq
 {
   int freqa;
   int freqc;
   int freqg;
   int freqt;
   int freqn;

   double total;
   
   double pa;
   double pc;
   double pg;
   double pt;

   double q;
 };

struct emp_freq * get_emp_freq (int nReads, int Len, struct reads_info ** ri_left, struct reads_info ** ri_right);

#endif
