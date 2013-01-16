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
  if (min_overlap > 99){min_overlap = 99;}	
	
  if (pval == 1.0)
   {
	   return (1);
   }   
  else if (pval == 0.01)
   {
     table_ptr = precomp_01[min_overlap];
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

  if (f >= cutoff) return (1);

  return (0);
}

/*
int stat_test2 (double pval, double oes, int min_overlap, double q)
{

    if (pval == 1.0)
   {
	   return (1);
   } 

  if (oes > 33) return (1);

  return (0);
}
*/


int stat_test2 (double pval, double oes, int min_overlap, double q)
{
  double * table_ptr;
  double cutoff;
  if (min_overlap > 99){min_overlap = 99;}
  if (q > 0.49){q = 0.49;}		
	
  if (pval == 1.0)
   {
	   return (1);
   }   
  else if (pval == 0.01)
   {
     table_ptr = precomp2_01[min_overlap];
   }
  else if (pval == 0.05)
   {
     table_ptr = precomp2_05[min_overlap];
   }
  else if (pval == 0.001)
   {
     table_ptr = precomp2_001[min_overlap];
   }
  else
   {
     table_ptr = precomp2_0001[min_overlap];
   }

  cutoff = table_ptr[(int)(q * 100)];

  if (oes > cutoff) return (1);

  return (0);
}

