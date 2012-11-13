#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI                      3.14159265358979323846L
#define SMOOTHING               1000


long double * precomp;


int
init_precomp (int read_len, double freqa, double freqc, double freqg, double freqt, int min_asm_len)
{
  long double   c, r, prod;
  double        pa, pc, pg, pt, freqall;
  double        q, p;
  int           k, i, t;
  int           dbg = 0;

  freqall = freqa + freqc + freqg + freqt;
  pa = freqa / freqall;
  pc = freqc / freqall;
  pg = freqg / freqall;
  pt = freqt / freqall;

  printf ("=> Precomp START\n");

  printf ("%lf %lf %lf %lf\n", pa, pc, pg, pt);
  /* probability of a match */
  q = pa * pa + pc * pc + pg * pg + pt * pt;

  /* allocate mem */
  precomp = (long double *) calloc ((SMOOTHING + 1), sizeof (long double));
  if (!precomp)
   {
     return (0);
   }

  /* compute the terms of the binomial test */
  for (t = 0; t <= SMOOTHING; ++t)
   {
     prod = 1;
     for (i = min_asm_len; i <= read_len; ++i)
      {
        r = 0;
        for (k = 0; k <= (int)(((double)t / SMOOTHING) * i); ++ k)
         {
           /* compute stirling's approximation */
           p   = (double)k / (double)i;
           c   = (
                  (1 / sqrtl (2 * PI)) * 
                  (sqrtl (i) / (sqrtl ( k * (i - k)))) * 
                  powl(p, -k) * 
                  powl (1 - p, - (i - k))
                 );
           if (k == 0 || k == i ) c = 1;
           if (t == SMOOTHING) printf ("i: %d  k: %d (i,k): %Lf\n", i, k, c);
           r   += c * powl (q,k) * powl (1 - q, i - k);
           if (r > 1) r = 1;
           if (t == SMOOTHING) {printf ("%Lf\n", r); ++dbg; if (dbg == 100) exit(1);}
         }
        prod *= r;
      }
     prod *= prod;
     precomp[t] = prod;
   }
     exit(1);

  /* store the cummulative probabilities */
/*  for ( i = 1; i <= read_len; ++ i)
   {
     to   = ((i + 1) * i) / 2 - 1;
     from = to + i - 1;
     for (k = from; k >= to; --k)
      {
        precomp[k] += precomp[k + 1];
      }
   }*/

  printf ("long double precomp[i] = \n {");
  for (i = 0; i <= SMOOTHING; ++ i)
   {
     if (i % 20 == 0) printf ("\n  "); 
     printf (" %.10LG", precomp[i]);
   }
  printf ("\n };\n");

  printf ("=> Precomp END\n");
  printf ("t=500: %LG\n", precomp[500]);
  printf ("t=700: %LG\n", precomp[700]);
  exit(1);
  return (1);
}

