#include <stdio.h>
/** @file statistics.h
  */

#include "pvalue.h"
extern int match_score;
extern int mismatch_score;



/** @brief  test
  *
  *  Statistical test
  *
  * @param pval
  *   The given maximum p-value
  *
  * @param f
  *  The expected number of matches for the current assembled sequences
  *
  * @param min_overlap
  *  The minimum overlap as input by the user in the command-line arguments
  *
  * @param q
  *  The probability of random base match
  */
int
stat_test (double pval, double f, int min_overlap, double q)
{
  double * table_ptr;
  double cutoff;


  if (pval == 0.01)
   {
     table_ptr = precomp_01[min_overlap];
     //cutoff = precomp_01[min_overlap - 1][(int)(q * 100 - 1)];
   }
  else if (pval == 0.05)
   {
     table_ptr = precomp_05[min_overlap];
   }
  else if (pval == 0.001)
   {
     table_ptr = precomp_001[min_overlap];
   }
  else
   {
     table_ptr = precomp_0001[min_overlap];
   }

  cutoff = table_ptr[(int)(q * 100)];

//  printf ("%f %f %d\n", cutoff, f, min_overlap );

  if (f >= cutoff) return (1);

  return (0);
}
