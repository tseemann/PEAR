#include <stdio.h>
#include <stdlib.h>
#include "emp.h"
#include "fastq.h"

struct emp_freq *
get_emp_freq (int nReads, int Len, struct reads_info ** ri_left, struct reads_info ** ri_right)
{
  struct emp_freq * ef;
  int i, j, k;
  char c;

  ef = (struct emp_freq *) malloc (sizeof (struct emp_freq));

  ef->freqa = ef->freqc = ef->freqg = ef->freqt = 0;

  for (i = 0; i < nReads; ++ i)
   {
     for (j = 0; j < Len; ++ j)
      {
        for (k = 1, c = ri_left[i]->data[j];  k >= 0; --k)
         {
           switch (c)
            {
              case 'A':
              case 'a':
                ++ ef->freqa;
                break;
              
              case 'C':
              case 'c':
                ++ ef->freqc;
                break;

              case 'G':
              case 'g':
                ++ ef->freqg;
                break;

              case 'T':
              case 't':
                ++ ef->freqt;
                break;
            }
           c = ri_right[i]->data[j];
         }
      }
   }
  printf ("A: %d C: %d G: %d T: %d\n", ef->freqa, ef->freqc, ef->freqg, ef->freqt );
  
  ef->total = ef->freqa + ef->freqc + ef->freqg + ef->freqt;
  ef->pa = ef->freqa / (double)ef->total;
  ef->pc = ef->freqc / (double)ef->total;
  ef->pg = ef->freqg / (double)ef->total;
  ef->pt = ef->freqt / (double)ef->total;
  
  ef->q = ef->pa * ef->pa + ef->pc * ef->pc + ef->pg * ef->pg + ef->pt * ef->pt;

  printf ("Total: %f\n", ef->total);
  printf ("pa: %f\n", ef->pa);
  printf ("pc: %f\n", ef->pc);
  printf ("pg: %f\n", ef->pg);
  printf ("pt: %f\n", ef->pt);

  printf ("q: %f\n", ef->q);

  return ef;

}
