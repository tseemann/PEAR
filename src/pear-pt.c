#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "args.h"
#include "emp.h"
#include "reader.h"
#include "async.h"

/** @file pear-pt.c
    @brief Main file containing scoring and assembly related functions (pthreads version)
*/
#define         PHRED_INIT                       33
#define         THREAD_MIN_PACKET_SIZE           500
#define         THREAD_PACKET_SIZE_DELTA         20

static pthread_mutex_t cs_mutex_wnd  = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t cs_mutex_io   = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t cs_mutex_out  = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  cs_mutex_cond = PTHREAD_COND_INITIALIZER;

struct thread_global_t thr_global;


//int stat_test (double, double, int, double);
int stat_test2 (double, double, int, double);

double
assemble_overlap (struct read_t * left, struct read_t * right, int base_left, int base_right, int ol_size, struct read_t * ai);

/*
double
assemble_overlap_ef (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai, struct emp_freq  * ef);
*/

/* TODO: Dynamically allocate them */
double      sc_eq[256][256];
double     sc_neq[256][256];
double     sc_eqA[256][256];
double     sc_eqC[256][256];
double     sc_eqG[256][256];
double     sc_eqT[256][256];
double   sc_neqAC[256][256];
double   sc_neqAG[256][256];
double   sc_neqAT[256][256];
double   sc_neqCA[256][256];
double   sc_neqCG[256][256];
double   sc_neqCT[256][256];
double   sc_neqGA[256][256];
double   sc_neqGC[256][256];
double   sc_neqGT[256][256];
double   sc_neqTA[256][256];
double   sc_neqTC[256][256];
double   sc_neqTG[256][256];
double     qs_mul[256][256];

int match_score    = 1;
int mismatch_score = 1;

/** @brief Trimming of forward part of unassembled sequence
  *
  * Finds two adjacent quality scores that are smaller than the minimum quality score threshold \a min_quality.
  * It then \e trims the part starting from the rightmost quality score position.
  *
  * @param read
  *   Forward sequence
  *
  * @param min_quality
  *   Minimum quality score threshold
  *
  * @param uncalled
  *   Variable to store the ratio of uncalled bases
  *
  * @return
  *   Returns the length of the trimmed sequence
  */
int
trim (struct read_t * read, int min_quality, double * uncalled)
{
  int                   i;
  char                * qscore;
  char                * data;

  qscore = read->qscore + 1;
  data   = read->data   + 1;
  *uncalled = 0;

  i = 1;
  if (!*data)
   {
     if (*(data - 1) == 'N' || *(data - 1) == 'n') ++ *uncalled;
     return (i);
   }

  while (*data)
   {
     if (*(data - 1) == 'N' || *(data - 1) == 'n') ++ (*uncalled);
     if (*(data - 1) - PHRED_INIT < min_quality && *data - PHRED_INIT < min_quality)
      {
        *data   = 0;
        *qscore = 0;
        *uncalled = (*uncalled) / (i);
        return (i);
      }
     ++ i;
     ++ data;
   }
  if (*(data - 1) == 'N' || *(data - 1) == 'n') ++ (*uncalled);
  *uncalled = (*uncalled) / i;
  return (i);
}

/** @brief Trimming of reverse part (reversed + complemented) of unassembled sequence
  *
  * Finds two adjacent quality scores that are smaller than the minimum quality score threshold \a min_quality.
  * It then \e trims the part starting from the rightmost quality score position.
  *
  * @param read
  *   Forward sequence
  *
  * @param min_quality
  *   Minimum quality score threshold
  *
  * @param uncalled
  *   Variable to store the ratio of uncalled bases
  *
  * @return
  *   Returns the length of the trimmed sequence
  */
int
trim_cpl (struct read_t * read, int min_quality, double * uncalled)
{
  int                   i;
  char                * qscore;
  char                * data;
  int                   len;

  qscore = read->qscore;
  data   = read->data;
  *uncalled = 0;
  len = 0;

  while (*data)
   {
     ++len;
     ++data;
   }

  data = read->data;
  for (i = len - 1; i > 0; -- i)
   {
     if (data[i] == 'N' || data[i] == 'n') ++ *uncalled;
     if (qscore[i] - PHRED_INIT < min_quality && qscore[i - 1] - PHRED_INIT < min_quality)
      {
        qscore[i - 1] = 0;
        data[i - 1]   = 0;
        memmove (data,   data + i,   len - i + 1);
        memmove (qscore, qscore + i, len - i + 1);
        *uncalled = (*uncalled) / (len - i);
        return (len - i);
      }
   }
  if (*data == 'N' || *data == 'n') ++ *uncalled;
  *uncalled = (*uncalled) / len;
  return (len);
}

/** @brief Initialize table of precomputed scores
    
    This function computes a table of scores for every possible combination of
    basepair and quality scores. The element \a sc_eqX[i][j] is the score
    for a basepair \a X with quality scores \a i and \j in the two sequences.
    The element \a sc_neqXY[i][j] is the score of basepairs \a X and \a Y with
    quality scores \a i and \j, respectively.

    @param match
      The match weight to be used in scoring

    @param mismatch
      The mismatch weight to be used in scoring

    @param ef
      Structure containing empirical frequencies
*/
void init_scores (int match, int mismatch, struct emp_freq * ef)
{
  int           i, j;
  double        ex, ey;
  double        pa2, pc2, pg2, pt2, pagt2, pcgt2, pact2, pacg2, pacg, pact, pagt, pcgt;
  double        pac2, pag2, pat2, pcg2, pct2, pgt2;
  double        p2acg, p2act, p2agt, p2cgt;

  pa2 = ef->pa * ef->pa;
  pc2 = ef->pc * ef->pc;
  pg2 = ef->pg * ef->pg;
  pt2 = ef->pt * ef->pt;

  pacg2 = pa2 + pc2 + pg2;  pacg = ef->pa + ef->pc + ef->pg; p2acg = pacg * pacg;
  pact2 = pa2 + pc2 + pt2;  pact = ef->pa + ef->pc + ef->pt; p2act = pact * pact;
  pagt2 = pa2 + pg2 + pt2;  pagt = ef->pa + ef->pg + ef->pt; p2agt = pagt * pagt;
  pcgt2 = pc2 + pg2 + pt2;  pcgt = ef->pc + ef->pg + ef->pt; p2cgt = pcgt * pcgt;

  pac2 = (ef->pa + ef->pc) * (ef->pa + ef->pc);
  pag2 = (ef->pa + ef->pg) * (ef->pa + ef->pg);
  pat2 = (ef->pa + ef->pt) * (ef->pa + ef->pt);
  pcg2 = (ef->pc + ef->pg) * (ef->pc + ef->pg);
  pct2 = (ef->pc + ef->pt) * (ef->pc + ef->pt);
  pgt2 = (ef->pg + ef->pt) * (ef->pg + ef->pt);


  for (i = 0; i < 256; ++ i)
   {
     for (j = 0; j < 256; ++ j) 
      {
        ex = pow (10.0, - (i - PHRED_INIT) / 10.0);
        ey = pow (10.0, - (j - PHRED_INIT) / 10.0);

        sc_eq[i][j]  =    match * ((1 - ex) * (1 - ey) + (ex * ey) / 3.0);
        sc_neq[i][j] = mismatch * (1 - (1.0 / 3.0) * (1 - ex) * ey - (1.0 / 3.0) * (1 - ey) * ex - (2.0 / 9.0) * ex * ey);
        qs_mul[i][j] = ex * ey;

        sc_eqA[i][j] = match * (1 - ex) * (1 - ey) + (ex * ey) * pcgt2 / p2cgt;
        sc_eqC[i][j] = match * (1 - ex) * (1 - ey) + (ex * ey) * pagt2 / p2agt;
        sc_eqG[i][j] = match * (1 - ex) * (1 - ey) + (ex * ey) * pact2 / p2act;
        sc_eqT[i][j] = match * (1 - ex) * (1 - ey) + (ex * ey) * pacg2 / p2acg;

        sc_neqAC[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pc / pcgt) - (1 - ex) * ey * (ef->pa / pagt) - ex * ey * (pg2 + pt2) / pgt2);
        sc_neqCA[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pa / pagt) - (1 - ex) * ey * (ef->pc / pcgt) - ex * ey * (pg2 + pt2) / pgt2);

        sc_neqAG[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pg / pcgt) - (1 - ex) * ey * (ef->pa / pact) - ex * ey * (pc2 + pt2) / pct2);
        sc_neqGA[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pa / pact) - (1 - ex) * ey * (ef->pg / pcgt) - ex * ey * (pc2 + pt2) / pct2);

        sc_neqAT[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pt / pcgt) - (1 - ex) * ey * (ef->pa / pacg) - ex * ey * (pc2 + pg2) / pcg2);
        sc_neqTA[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pa / pacg) - (1 - ex) * ey * (ef->pt / pcgt) - ex * ey * (pc2 + pg2) / pcg2);

        sc_neqCG[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pg / pagt) - (1 - ex) * ey * (ef->pc / pact) - ex * ey * (pa2 + pt2) / pat2);
        sc_neqGC[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pc / pact) - (1 - ex) * ey * (ef->pg / pagt) - ex * ey * (pa2 + pt2) / pat2);

        sc_neqCT[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pt / pagt) - (1 - ex) * ey * (ef->pc / pacg) - ex * ey * (pa2 + pg2) / pag2);
        sc_neqTC[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pc / pacg) - (1 - ex) * ey * (ef->pt / pagt) - ex * ey * (pa2 + pg2) / pag2);

        sc_neqGT[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pt / pact) - (1 - ex) * ey * (ef->pg / pacg) - ex * ey * (pa2 + pc2) / pac2);
        sc_neqTG[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pg / pacg) - (1 - ex) * ey * (ef->pt / pact) - ex * ey * (pa2 + pc2) / pac2);

     }
   }
}

inline void
scoring_ef (char dleft, char dright, char qleft, char qright, int score_method, double * score, double * oes, int match, int mismatch, struct emp_freq * ef)
{
  double tmp;

  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += (ef->q * match - (1 - ef->q) * mismatch);
          *oes    = *score;
          break;
        case 2:
          tmp     = (1 - ef->q) * mismatch;
          *oes   += (ef->q * match - tmp);
          *score -= tmp; 
          break;
        case 3:
          *score -= mismatch;
          *oes   += (ef->q * match - (1 - ef->q) * mismatch);
          break;
      }
   }
  else if (dleft == dright)     /* equal */
   {
     switch (score_method)
      {
        case 1:
          switch (dleft)
           {
             case 'A':
               *score += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'C':
               *score += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'G':
               *score += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'T':
               *score += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft] / match) * mismatch);
               break;
           }
          *oes    = *score;
          break;
        case 2:
          switch (dleft)
           {
             case 'A':
               *score += sc_eqA[(int)qright][(int)qleft];
               *oes += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'C':
               *score += sc_eqC[(int)qright][(int)qleft];
               *oes   += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'G':
               *score += sc_eqG[(int)qright][(int)qleft];
               *oes   += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'T':
               *score += sc_eqT[(int)qright][(int)qleft];
               *oes   += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft] / match) * mismatch);
               break;
           }
          break;
        case 3:
          switch (dleft)
           {
             case 'A':
               *oes += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'C':
               *oes += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'G':
               *oes += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft] / match) * mismatch);
               break;
             case 'T':
               *oes += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft] / match) * mismatch);
               break;
           }
          *score += match;
          break;
      }
   }
  else          /* not equal */
   {
     switch (score_method)
      {
        case 1:
           switch  (dleft)
           {
             case 'A':
               switch (dright)
                {
                  case 'C':
                    *score = *score - (sc_neqAC[(int)qleft][(int)qright] - (1 - sc_neqAC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score = *score - (sc_neqAG[(int)qleft][(int)qright] - (1 - sc_neqAG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score = *score - (sc_neqAT[(int)qleft][(int)qright] - (1 - sc_neqAT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    *score = *score - (sc_neqCA[(int)qleft][(int)qright] - (1 - sc_neqCA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score = *score - (sc_neqCG[(int)qleft][(int)qright] - (1 - sc_neqCG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score = *score - (sc_neqCT[(int)qleft][(int)qright] - (1 - sc_neqCT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    *score = *score - (sc_neqGA[(int)qleft][(int)qright] - (1 - sc_neqGA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *score = *score - (sc_neqGC[(int)qleft][(int)qright] - (1 - sc_neqGC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score = *score - (sc_neqGT[(int)qleft][(int)qright] - (1 - sc_neqGT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    *score = *score - (sc_neqTA[(int)qleft][(int)qright] - (1 - sc_neqTA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *score = *score - (sc_neqTC[(int)qleft][(int)qright] - (1 - sc_neqTC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score = *score - (sc_neqTG[(int)qleft][(int)qright] - (1 - sc_neqTG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
           }
          *oes = *score;
          break;
        case 2:
          switch  (dleft)
           {
             case 'A':
               switch (dright)
                {
                  case 'C':
                    *score -= sc_neqAC[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqAC[(int)qleft][(int)qright] - (1 - sc_neqAC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score -= sc_neqAG[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqAG[(int)qleft][(int)qright] - (1 - sc_neqAG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score -= sc_neqAT[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqAT[(int)qleft][(int)qright] - (1 - sc_neqAT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqCA[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqCA[(int)qleft][(int)qright] - (1 - sc_neqCA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score -= sc_neqCG[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqCG[(int)qleft][(int)qright] - (1 - sc_neqCG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score -= sc_neqCT[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqCT[(int)qleft][(int)qright] - (1 - sc_neqCT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqGA[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqGA[(int)qleft][(int)qright] - (1 - sc_neqGA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *score -= sc_neqGC[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqGC[(int)qleft][(int)qright] - (1 - sc_neqGC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *score -= sc_neqGT[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqGT[(int)qleft][(int)qright] - (1 - sc_neqGT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqTA[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqTA[(int)qleft][(int)qright] - (1 - sc_neqTA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *score -= sc_neqTC[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqTC[(int)qleft][(int)qright] - (1 - sc_neqTC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *score -= sc_neqTG[(int)qleft][(int)qright];
                    *oes    = *oes - (sc_neqTG[(int)qleft][(int)qright] - (1 - sc_neqTG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
           }
          break;
        case 3:
          *score -= mismatch;
           switch  (dleft)
           {
             case 'A':
               switch (dright)
                {
                  case 'C':
                    *oes = *oes - (sc_neqAC[(int)qleft][(int)qright] - (1 - sc_neqAC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *oes = *oes - (sc_neqAG[(int)qleft][(int)qright] - (1 - sc_neqAG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqAT[(int)qleft][(int)qright] - (1 - sc_neqAT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    *oes = *oes - (sc_neqCA[(int)qleft][(int)qright] - (1 - sc_neqCA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *oes = *oes - (sc_neqCG[(int)qleft][(int)qright] - (1 - sc_neqCG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqCT[(int)qleft][(int)qright] - (1 - sc_neqCT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    *oes = *oes - (sc_neqGA[(int)qleft][(int)qright] - (1 - sc_neqGA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *oes = *oes - (sc_neqGC[(int)qleft][(int)qright] - (1 - sc_neqGC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqGT[(int)qleft][(int)qright] - (1 - sc_neqGT[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    *oes = *oes - (sc_neqTA[(int)qleft][(int)qright] - (1 - sc_neqTA[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'C':
                    *oes = *oes - (sc_neqTC[(int)qleft][(int)qright] - (1 - sc_neqTC[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                  case 'G':
                    *oes = *oes - (sc_neqTG[(int)qleft][(int)qright] - (1 - sc_neqTG[(int)qleft][(int)qright] / mismatch) * match);
                    break;
                }
               break;
           }
          break;
      }
   }
}

/* TODO: Remember to speed up this function by doing something with the multiplication and division of match/mismatch */
inline void 
scoring (char dleft, char dright, char qleft, char qright, int score_method, double * score, double * oes, int match, int mismatch)
{
  double tmp;

  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += (0.25 * match - (1 - 0.25) * mismatch);
          *oes    = *score;
          break;
        case 2:
          tmp     = (1 - 0.25) * mismatch;
          *oes   += (0.25 * match - tmp);
          *score -= tmp; 
          break;
        case 3:
          *score -= mismatch;
          *oes += (0.25 * match - (1 - 0.25) * mismatch);
          break;
      }
   }
  else if (dleft == dright)     /* equal */
   {
     switch (score_method)
      {
        case 1:
          *score += (sc_eq[(int)qright][(int)qleft] - (1 - sc_eq[(int)qright][(int)qleft] / match) * mismatch);
          *oes    = *score;
          break;
        case 2:
          tmp     = sc_eq[(int)qright][(int)qleft];
          *oes   += (tmp - (1 - sc_eq[(int)qright][(int)qleft] / match) * mismatch);
          *score += tmp;
          break;
        case 3:
          *score += match;
          *oes   += (sc_eq[(int)qright][(int)qleft] - (1 - sc_eq[(int)qright][(int)qleft] / match) * mismatch);
          break;
      }
   }
  else          /* not equal */
   {
     switch (score_method)
      {
        case 1:
          *score = *score - (sc_neq[(int)qright][(int)qleft] - (1 - sc_neq[(int)qright][(int)qleft] / mismatch) * match);
          *oes    = *score;
          break;
        case 2:
          tmp     = sc_neq[(int)qright][(int)qleft];
          *oes    = *oes - (tmp - (1 - sc_neq[(int)qright][(int)qleft] / mismatch) * match);
          *score -= tmp;
          break;
        case 3:
          *score -= mismatch;
          *oes    = *score - (sc_neq[(int)qright][(int)qleft] - (1 - sc_neq[(int)qright][(int)qleft] / mismatch) * match);
          break;
      }
   }
}

inline int
assembly_ef (struct read_t * left, struct read_t * right, int match_score, int mismatch_score, struct emp_freq * ef, struct user_args  * sw)
{
  int                   i,j;
  int                   n;
  double                score;
  double                oes;
  double                best_score = 0;
  double                best_oes = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  int                   nMatch;
  int                   asm_len = 0;
  int                   st_pass;
  double                uncalled = 0;
  
  n = strlen (left->data);
     
  /* compute score for every overlap */
  score = 0;
  oes   = 0;
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {     
     nMatch = 0;
     score = 0;
     oes   = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef (left->data[n - i + j], right->data[j], left->qscore[n - i + j], right->qscore[j], sw->score_method, &score, &oes, match_score, mismatch_score, ef);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }
     if (score > best_score)
      {
        best_overlap = i;
        best_score   = score;
        best_oes     = oes;
      }
   }

  /* compute for score for runthrough case */
  for (i = n - 1; i > 0; --i)
   {
     score  = 0;
     oes    = 0;
     nMatch = 0;
     for (j = 0; j < i; ++j)
      {
        scoring_ef (left->data[j], right->data[n - i + j], left->qscore[j], right->qscore[n - i + j], sw->score_method, &score, &oes, match_score, mismatch_score, ef);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        run_through  = 1;
        best_overlap = i;
        best_score   = score;
        best_oes     = oes;
      }
   }


  if (sw->test == 1)
   {
     //st_pass = stat_test2 (sw->p_value, best_score, sw->min_overlap, ef->q);
     st_pass = stat_test2 (sw->p_value, best_oes, sw->min_overlap, ef->q);
   }
  else
   {
     //st_pass = stat_test2 (sw->p_value, best_score, best_overlap, ef->q);
     st_pass = stat_test2 (sw->p_value, best_oes, best_overlap, ef->q);
   }

  if (!st_pass) return (0);


  /* do the assembly!!!! */
  if (!run_through)
   {
     if (best_overlap == 0)
      {
        asm_len = 2 * n;

        for (j = 0; j < n; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = 0; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 1;
         }
        else
         {
           return (0);
         }
      }
     else if (best_overlap == n )
      {
        asm_len         = n;

        for (j = 0; j < asm_len; ++ j)
          if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[j] == 'N' || right->data[j] == 'n')) ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 0;
           assemble_overlap (left, right, 0, 0, n, left);
           left->data[n]   = 0;
           left->qscore[n] = 0;
         }
        else
         {
           return (0);
         }
      }
     else
      {
        asm_len = 2 * n - best_overlap;
        for (j = 0; j < n - best_overlap; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = n - best_overlap; j < n; ++ j)
          if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[j - n + best_overlap] == 'N' || right->data[j - n + best_overlap] == 'n'))  ++uncalled;
        for (j = best_overlap; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 1;
           
           assemble_overlap (left, right, n - best_overlap, 0, best_overlap, left);
           memmove (right->data,   right->data   + best_overlap,  n - best_overlap);
           memmove (right->qscore, right->qscore + best_overlap,  n - best_overlap);
           /* THIS IS WRONG */
           //memcpy (right->data,   right->data   + best_overlap,  n - best_overlap);
           //memcpy (right->qscore, right->qscore + best_overlap,  n - best_overlap);

           right->data[n   - best_overlap] = 0;
           right->qscore[n - best_overlap] = 0;
         }
        else
         {
           return (0);
         }
      }
   }     /* run-through case */
  else
   {
     asm_len = best_overlap;
     
     /* compute uncalled */
     for (j = 0; j < asm_len; ++ j)
       if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[n - best_overlap + j] == 'N' || right->data[n - best_overlap + j] == 'n')) ++uncalled;
     uncalled /= asm_len;

     if (asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
      {
        assemble_overlap (left, right, 0, n - best_overlap, best_overlap, left);
        
        left->data[best_overlap]   = 0;
        left->qscore[best_overlap] = 0;
        *(left->data - 1) = 0;   /* flag that it's one piece */
      }
     else
      {
        return (0);
      }
   }

  /* TODO: Remember to reset *(let->data - 1) */

  /* TODO: Optimize this in assemble_overlap? */

  return (1);
}

inline int
assembly (struct read_t * left, struct read_t * right, int match_score, int mismatch_score, struct user_args  * sw)
{
  int                   i,j;
  int                   n;
  double                score;
  double                oes;
  double                best_score = 0;
  double                best_oes   = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  int                   nMatch;
  int                   asm_len = 0;
  int                   st_pass;
  double                uncalled = 0;
  
  n = strlen (left->data);
     
  /* compute score for every overlap */
  score = 0;
  oes   = 0;
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {     
     nMatch = 0;
     score = 0;
     oes   = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring (left->data[n - i + j], right->data[j], left->qscore[n - i + j], right->qscore[j], sw->score_method, &score, &oes, match_score, mismatch_score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }
     if (score > best_score)
      {
        best_overlap = i;
        best_score   = score;
        best_oes     = oes;
      }
   }

  /* compute for score for runthrough case */
  for (i = n - 1; i > 0; --i)
   {
     score  = 0;
     oes    = 0;
     nMatch = 0;
     for (j = 0; j < i; ++j)
      {
        scoring (left->data[j], right->data[n - i + j], left->qscore[j], right->qscore[n - i + j], sw->score_method, &score, &oes, match_score, mismatch_score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        run_through  = 1;
        best_overlap = i;
        best_score   = score;
        best_oes     = oes;
      }
   }


  if (sw->test == 1)
   {
     //st_pass = stat_test2 (sw->p_value, best_score, sw->min_overlap, 0.25);
     st_pass = stat_test2 (sw->p_value, best_oes, sw->min_overlap, 0.25);
   }
  else
   {
     //st_pass = stat_test2 (sw->p_value, best_score, best_overlap, 0.25);
     st_pass = stat_test2 (sw->p_value, best_oes, best_overlap, 0.25);
   }

  if (!st_pass) return (0);


  /* do the assembly!!!! */
  if (!run_through)
   {
     if (best_overlap == 0)
      {
        asm_len = 2 * n;

        for (j = 0; j < n; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = 0; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 1;
         }
        else
         {
           return (0);
         }
      }
     else if (best_overlap == n )
      {
        asm_len         = n;

        for (j = 0; j < asm_len; ++ j)
          if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[j] == 'N' || right->data[j] == 'n')) ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 0;
           assemble_overlap (left, right, 0, 0, n, left);
           left->data[n]   = 0;
           left->qscore[n] = 0;
         }
        else
         {
           return (0);
         }
      }
     else
      {
        asm_len = 2 * n - best_overlap;
        for (j = 0; j < n - best_overlap; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = n - best_overlap; j < n; ++ j)
          if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[j - n + best_overlap] == 'N' || right->data[j - n + best_overlap] == 'n'))  ++uncalled;
        for (j = best_overlap; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           *(left->data - 1) = 1;
           
           assemble_overlap (left, right, n - best_overlap, 0, best_overlap, left);
           memmove (right->data,   right->data   + best_overlap,  n - best_overlap);
           memmove (right->qscore, right->qscore + best_overlap,  n - best_overlap);
           /* THIS IS WRONG */
           //memcpy (right->data,   right->data   + best_overlap,  n - best_overlap);
           //memcpy (right->qscore, right->qscore + best_overlap,  n - best_overlap);

           right->data[n   - best_overlap] = 0;
           right->qscore[n - best_overlap] = 0;
         }
        else
         {
           return (0);
         }
      }
   }     /* run-through case */
  else
   {
     asm_len = best_overlap;
     
     /* compute uncalled */
     for (j = 0; j < asm_len; ++ j)
       if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[n - best_overlap + j] == 'N' || right->data[n - best_overlap + j] == 'n')) ++uncalled;
     uncalled /= asm_len;

     if (asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
      {
        assemble_overlap (left, right, 0, n - best_overlap, best_overlap, left);
        
        left->data[best_overlap]   = 0;
        left->qscore[best_overlap] = 0;
        *(left->data - 1) = 0;   /* flag that it's one piece */
      }
     else
      {
        return (0);
      }
   }

  /* TODO: Remember to reset *(let->data - 1) */

  /* TODO: Optimize this in assemble_overlap? */

  return (1);
}

double
assemble_overlap (struct read_t * left, struct read_t * right, int base_left, int base_right, int ol_size, struct read_t * ai)
{
  int           i; 
  char          x, y;
  char          qx, qy;
  //double        exp_match  = 0;

  for (i = 0; i < ol_size; ++i)
   {
     x  = left->data[base_left + i]; 
     y  = right->data[base_right + i];
     qx = left->qscore[base_left + i];
     qy = right->qscore[base_right + i];
     if ( (x == 'N' || x == 'n') && (y == 'N' || y == 'n'))
      {
        //exp_match += 0.25; sm_len
        ai->data[base_left + i]          = 'N';
        ai->qscore[base_left + i] = ( qx < qy ) ? qx : qy;
      }
     else if (x == 'N' || x == 'n')
      {
        //exp_match += 0.25; 
        ai->data[base_left + i]          = y;
        ai->qscore[base_left + i] = qy;
      }
     else if (y == 'N' || y == 'n')
      {
        //exp_match += 0.25; 
        ai->data[base_left + i]          = x;
        ai->qscore[base_left + i] = qx;
      }
     else
      {
        if (x == y)
         {
           //exp_match += (sc_eq[(int)qx][(int)qy] / match_score);
           
           ai->data[base_left + i] = x;
           ai->qscore[base_left + i] = (right->qscore[base_right + i] - PHRED_INIT) + (left->qscore[base_left + i] - PHRED_INIT) + PHRED_INIT; //qs_mul[qx][qy];
         }
        else
         {
           //exp_match += (1 - sc_neq[(int)qx][(int)qy] / mismatch_score);
           
           if (qx > qy)
            {
              ai->data[base_left + i]          =  x;
              ai->qscore[base_left + i] = qx;
            }
           else
            {
              ai->data[base_left + i]          =  y;
              ai->qscore[base_left + i] = qy;
            }
         }
      }
   }
  return (0);
}

/*
double
assemble_overlap_ef (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai, struct emp_freq  * ef)
{
  int           i; 
  char          x, y;
  char          qx, qy;
  double        exp_match  = 0;

  for (i = 0; i < ol_size; ++i)
   {
     x  = left->data[base_left + i]; 
     y  = right->data[base_right + i];
     qx = left->qscore[base_left + i];
     qy = right->qscore[base_right + i];
     if ( (x == 'N' || x == 'n') && (y == 'N' || y == 'n'))
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = 'N';
        ai->qscore[base_left + i] = ( qx < qy ) ? qx : qy;
      }
     else if (x == 'N' || x == 'n')
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = y;
        ai->qscore[base_left + i] = qy;
      }
     else if (y == 'N' || y == 'n')
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = x;
        ai->qscore[base_left + i] = qx;
      }
     else
      {
        if (x == y)
         {
           //exp_match += (sc_eq[(int)qx][(int)qy] / match_score);
           switch (x)
           {
             case 'A':
               exp_match += (sc_eqA[(int)qy][(int)qx] / match_score);
               break;
             case 'C':
               exp_match += (sc_eqC[(int)qy][(int)qx] / match_score);
               break;
             case 'G':
               exp_match += (sc_eqG[(int)qy][(int)qx] / match_score);
               break;
             case 'T':
               exp_match += (sc_eqT[(int)qy][(int)qx] / match_score);
               break;
           }
           ai->data[base_left + i] = x;
           ai->qscore[base_left + i] = (right->qscore[base_right + i] - PHRED_INIT) + (left->qscore[base_left + i] - PHRED_INIT) + PHRED_INIT; //qs_mul[qx][qy];
         }
        else
         {
           //exp_match += (1 - sc_neq[(int)qx][(int)qy] / mismatch_score);
           switch  (x)
           {
             case 'A':
               switch (y)
                {
                  case 'C':
                    exp_match += (1 - sc_neqAC[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'G':
                    exp_match += (1 - sc_neqAG[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'T':
                    exp_match += (1 - sc_neqAT[(int)qx][(int)qy]/mismatch_score);
                    break;
                }
               break;
             case 'C':
               switch (y)
                {
                  case 'A':
                    exp_match += (1 - sc_neqCA[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'G':
                    exp_match += (1 - sc_neqCG[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'T':
                    exp_match += (1 - sc_neqCT[(int)qx][(int)qy]/mismatch_score);
                    break;
                }
               break;
             case 'G':
               switch (y)
                {
                  case 'A':
                    exp_match += (1 - sc_neqGA[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'C':
                    exp_match += (1 - sc_neqGC[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'T':
                    exp_match += (1 - sc_neqGT[(int)qx][(int)qy]/mismatch_score);
                    break;
                }
               break;
             case 'T':
               switch (y)
                {
                  case 'A':
                    exp_match += (1 - sc_neqTA[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'C':
                    exp_match += (1 - sc_neqTC[(int)qx][(int)qy]/mismatch_score);
                    break;
                  case 'G':
                    exp_match += (1 - sc_neqTG[(int)qx][(int)qy]/mismatch_score);
                    break;
                }
               break;
           }

           if (qx > qy)
            {
              ai->data[base_left + i]          =  x;
              ai->qscore[base_left + i] = qx;
            }
           else
            {
              ai->data[base_left + i]          =  y;
              ai->qscore[base_left + i] = qy;
            }
         }
      }
   }
  if (ol_size == 0) return (0);
  return (exp_match / (double)ol_size);
}
*/

void mstrrev (char * s)
{
  char * q = s;

  while (q && *q) ++q;

  for (--q; s < q; ++s, --q)
   {
     *s = *s ^ *q;
     *q = *s ^ *q;
     *s = *s ^ *q;
   }
}

void mstrcpl (char * s)
{
  if (!s) return;

  while (*s)
   {
     switch (*s)
      {
        case 'A':
        case 'a':
                  *s = 'T';
                  break;
        case 'C':
        case 'c':
                  *s = 'G';
                  break;
        case 'G':
        case 'g':
                  *s = 'C';
                  break;
        case 'T':
        case 't':
                  *s = 'A';
                  break;
      }
     ++s;
   }
}


int 
validate_input (int nleft, int nright)
{
  if (!nleft || !nright)
   {
     fprintf (stderr, "ERROR: At least one of the input files contains no records.\n");
     return (0);
   }

  if (nleft != nright)
   {
     fprintf (stderr, "ERROR: Number of records in the two input files does not match.\n");
     return (0);
   }

  return (1);
}

char *
makefilename (const char * prefix, const char * suffix)
{
  char * filename;

  filename = (char *) malloc ((strlen(prefix) + strlen(suffix) + 1) * sizeof (char));

  strcpy (filename, prefix);
  strcat (filename, suffix);

  return (filename);
}

void write_data (struct read_t ** fwd, struct read_t ** rev, unsigned int elms, FILE ** fd)
{
  int i;
  char two_piece;
  
  for ( i = 0; i < elms; ++ i)
   {
     two_piece = *(fwd[i]->data - 1);
     *(fwd[i]->data - 1) = 0;

     if (*(fwd[i]->qscore - 1) == 1)   /* assembled */
      {
        *(fwd[i]->qscore - 1) = 0;
        fprintf (fd[0], "%s\n", fwd[i]->header);
        if (!two_piece)
         {
           fprintf (fd[0], "%s\n", fwd[i]->data);
         }
        else
         {
           fprintf (fd[0], "%s",   fwd[i]->data);
           fprintf (fd[0], "%s\n", rev[i]->data);
         }
        fprintf (fd[0], "+\n");

        if (!two_piece)
         {
           fprintf (fd[0], "%s\n", fwd[i]->qscore);
         }
        else
         {
           fprintf (fd[0], "%s",   fwd[i]->qscore);
           fprintf (fd[0], "%s\n", rev[i]->qscore);
         }
      }
     else if (*(fwd[i]->qscore - 1) == 2)                                            /* not assembled */
      {
        *(fwd[i]->qscore - 1) = 0;
           /* discarded reads*/
           /* Maybe consider printing the untrimmed sequences */
           fprintf (fd[3], "%s\n", fwd[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", fwd[i]->data,  fwd[i]->qscore);
           fprintf (fd[3], "%s\n", rev[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", rev[i]->data, rev[i]->qscore); /* printing the reverse compliment of the original sequence */
      }
     else   /* unassembled reads*/
      {
        *(fwd[i]->qscore - 1) = 0;
           fprintf (fd[1], "%s\n", fwd[i]->header);
           fprintf (fd[2], "%s\n", rev[i]->header);
           fprintf (fd[1], "%s\n+\n%s\n", fwd[i]->data,  fwd[i]->qscore);
           fprintf (fd[2], "%s\n+\n%s\n", rev[i]->data, rev[i]->qscore); /* printing the reverse compliment of the original sequence */
      }
   }
}

void flip_list (void)
{
  struct blockinfo_t * elm1;
  struct blockinfo_t * elm2;

  elm1 = thr_global.xblock;
  elm2 = thr_global.yblock;

  thr_global.xblock = elm2;
  thr_global.yblock = elm1;
}

inline int assign_reads (struct blockinfo_t * block, struct thread_local_t * thr_local)
{
  int r;

  if (block->reads == block->processed) return (0);

  r = THREAD_MIN_PACKET_SIZE + (rand() % THREAD_MIN_PACKET_SIZE);
  thr_local->block   = block;
  thr_local->start   = block->processed;

  if (block->reads - block->processed <= r + THREAD_PACKET_SIZE_DELTA)
   {
     thr_local->end    = block->reads;
     block->processed = block->reads;
   }
  else
   {
     thr_local->end    = thr_local->start + r;
     block->processed = thr_local->start + r;
   }
  ++ block->threads ;
  return (1);
}

void * entry_point_ef (void * data)
{
  struct thread_local_t * thr_local;
  int ass, i, sleep, elms;
  double                uncalled_forward, uncalled_reverse;

  /* initialization of values */
  thr_local = (struct thread_local_t *)data;

  
  while (1)
   {
     // TODO: The first condition in the next line is useless?
     if (thr_local->block == thr_global.xblock && thr_global.io_thread == thr_local->id)
      {
        pthread_mutex_lock (&cs_mutex_io);
        write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
        // TODO: read_data ();
        elms = db_get_next_reads (thr_global.xblock->fwd, 
                                  thr_global.xblock->rev, 
                                  thr_global.yblock->fwd, 
                                  thr_global.yblock->rev);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("READ %d elms\n", elms);
           printf ("WAKE_UP_ALL!    (reads: %d processed: %d)\n", thr_global.xblock->reads, thr_global.xblock->processed);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        // wakeup threads
        pthread_cond_broadcast (&cs_mutex_cond);
        pthread_mutex_unlock (&cs_mutex_wnd);
      }

     sleep = 1;

     pthread_mutex_lock (&cs_mutex_wnd);
     if (thr_global.xblock->reads == thr_global.xblock->processed)
      {
        if (thr_global.finish)
         {
           pthread_mutex_unlock (&cs_mutex_wnd);
           if (!thr_global.xblock->threads)
            {
              write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
            }
           break;
         }
        /* is this the last thread using the current buffer? */
        if (thr_global.xblock->threads == 0 && thr_global.io_thread == -1)
         {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED IO THREAD %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
           thr_global.io_thread = thr_local->id;
           thr_local->block = thr_global.xblock;
           sleep = 0;
         }
        else
         {
           if (assign_reads (thr_global.yblock, thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED y READS to %d   (%d - %d)  Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.yblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
         }
      }
     else
      {
        if (assign_reads (thr_global.xblock,thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
      }
     if (sleep)
      {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("WAKING %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
      }
     pthread_mutex_unlock (&cs_mutex_wnd);

     if (!sleep && thr_global.io_thread != thr_local->id)
      {
        // TODO: Process reads and make this an inline function
        for (i = thr_local->start; i < thr_local->end; ++ i)
         {
           mstrrev (thr_local->block->rev->reads[i]->data);    /* reverse the sequence */
           mstrcpl (thr_local->block->rev->reads[i]->data);    /* complement the sequence */
           mstrrev (thr_local->block->rev->reads[i]->qscore);  /* reverse the quality scores */

           //TODO Switch for empirical frequencies
           ass = assembly_ef (thr_local->block->fwd->reads[i], thr_local->block->rev->reads[i], thr_local->match_score, thr_local->mismatch_score, thr_local->ef, thr_local->sw);
           *(thr_local->block->fwd->reads[i]->qscore - 1) = ass;
           if (!ass)
            {
              if (trim     (thr_local->block->fwd->reads[i], thr_local->sw->qual_thres, 
                        &uncalled_forward) < thr_local->sw->min_trim_len ||
                  trim_cpl (thr_local->block->rev->reads[i], thr_local->sw->qual_thres,
                        &uncalled_reverse) < thr_local->sw->min_trim_len ||
                  uncalled_forward >= thr_local->sw->max_uncalled || uncalled_reverse >= thr_local->sw->max_uncalled)
               {
                 *(thr_local->block->fwd->reads[i]->qscore - 1) = 2;
               }
              else
               {
                 *(thr_local->block->fwd->reads[i]->qscore - 1) = 3;
               }
            }
         }

        pthread_mutex_lock (&cs_mutex_wnd);
          -- thr_local->block->threads;
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Finished %d  (Remaining threads: %d)\n", thr_local->id, thr_local->block->threads);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_mutex_unlock (&cs_mutex_wnd);
      }
  //   pthread_mutex_lock (&cs_mutex_out);
  //     printf ("Thread: %d\n", thr_local->id);
  //   pthread_mutex_unlock (&cs_mutex_out);
   }
  return (NULL);
}
void * entry_point (void * data)
{
  struct thread_local_t * thr_local;
  int ass, i, sleep, elms;
  double                uncalled_forward, uncalled_reverse;

  /* initialization of values */
  thr_local = (struct thread_local_t *)data;

  
  while (1)
   {
     // TODO: The first condition in the next line is useless?
     if (thr_local->block == thr_global.xblock && thr_global.io_thread == thr_local->id)
      {
        pthread_mutex_lock (&cs_mutex_io);
        write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
        // TODO: read_data ();
        elms = db_get_next_reads (thr_global.xblock->fwd, 
                                  thr_global.xblock->rev, 
                                  thr_global.yblock->fwd, 
                                  thr_global.yblock->rev);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("READ %d elms\n", elms);
           printf ("WAKE_UP_ALL!    (reads: %d processed: %d)\n", thr_global.xblock->reads, thr_global.xblock->processed);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        // wakeup threads
        pthread_cond_broadcast (&cs_mutex_cond);
        pthread_mutex_unlock (&cs_mutex_wnd);
      }

     sleep = 1;

     pthread_mutex_lock (&cs_mutex_wnd);
     if (thr_global.xblock->reads == thr_global.xblock->processed)
      {
        if (thr_global.finish)
         {
           pthread_mutex_unlock (&cs_mutex_wnd);
           if (!thr_global.xblock->threads)
            {
              write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
            }
           break;
         }
        /* is this the last thread using the current buffer? */
        if (thr_global.xblock->threads == 0 && thr_global.io_thread == -1)
         {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED IO THREAD %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
           thr_global.io_thread = thr_local->id;
           thr_local->block = thr_global.xblock;
           sleep = 0;
         }
        else
         {
           if (assign_reads (thr_global.yblock, thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED y READS to %d   (%d - %d)  Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.yblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
         }
      }
     else
      {
        if (assign_reads (thr_global.xblock,thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
      }
     if (sleep)
      {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("WAKING %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
      }
     pthread_mutex_unlock (&cs_mutex_wnd);

     if (!sleep && thr_global.io_thread != thr_local->id)
      {
        // TODO: Process reads and make this an inline function
        for (i = thr_local->start; i < thr_local->end; ++ i)
         {
           mstrrev (thr_local->block->rev->reads[i]->data);    /* reverse the sequence */
           mstrcpl (thr_local->block->rev->reads[i]->data);    /* complement the sequence */
           mstrrev (thr_local->block->rev->reads[i]->qscore);  /* reverse the quality scores */

           //TODO Switch for empirical frequencies
           ass = assembly (thr_local->block->fwd->reads[i], thr_local->block->rev->reads[i], thr_local->match_score, thr_local->mismatch_score, thr_local->sw);
           *(thr_local->block->fwd->reads[i]->qscore - 1) = ass;
           if (!ass)
            {
              if (trim     (thr_local->block->fwd->reads[i], thr_local->sw->qual_thres, 
                        &uncalled_forward) < thr_local->sw->min_trim_len ||
                  trim_cpl (thr_local->block->rev->reads[i], thr_local->sw->qual_thres,
                        &uncalled_reverse) < thr_local->sw->min_trim_len ||
                  uncalled_forward >= thr_local->sw->max_uncalled || uncalled_reverse >= thr_local->sw->max_uncalled)
               {
                 *(thr_local->block->fwd->reads[i]->qscore - 1) = 2;
               }
              else
               {
                 *(thr_local->block->fwd->reads[i]->qscore - 1) = 3;
               }
            }
         }

        pthread_mutex_lock (&cs_mutex_wnd);
          -- thr_local->block->threads;
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Finished %d  (Remaining threads: %d)\n", thr_local->id, thr_local->block->threads);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_mutex_unlock (&cs_mutex_wnd);
      }
  //   pthread_mutex_lock (&cs_mutex_out);
  //     printf ("Thread: %d\n", thr_local->id);
  //   pthread_mutex_unlock (&cs_mutex_out);
   }
  return (NULL);
}

void init_thr_global (void)
{
  thr_global.xblock = (struct blockinfo_t *) calloc (1,sizeof(struct blockinfo_t));
  thr_global.yblock = (struct blockinfo_t *) calloc (1,sizeof(struct blockinfo_t));
  thr_global.xblock->fwd = (struct block_t *) calloc (1,sizeof(struct block_t));
  thr_global.xblock->rev = (struct block_t *) calloc (1,sizeof(struct block_t));
  thr_global.yblock->fwd = (struct block_t *) calloc (1,sizeof(struct block_t));
  thr_global.yblock->rev = (struct block_t *) calloc (1,sizeof(struct block_t));

  thr_global.xblock->reads     = 0;
  thr_global.xblock->processed = 0;
  thr_global.xblock->threads   = 0;
  thr_global.yblock->reads     = 0;
  thr_global.yblock->processed = 0;
  thr_global.yblock->threads   = 0;

  thr_global.io_thread = -1;
  thr_global.finish = 0;
}

void destroy_thr_global (void)
{
  fclose (thr_global.fd[0]);
  fclose (thr_global.fd[1]);
  fclose (thr_global.fd[2]);
  fclose (thr_global.fd[3]);
   
  free (thr_global.xblock->fwd);
  free (thr_global.xblock->rev);
  free (thr_global.yblock->fwd);
  free (thr_global.yblock->rev);
  free (thr_global.xblock);
  free (thr_global.yblock);
}

/** @brief Initialize output file names
    
    Create output file names based on input parameters, open them for writing
    and store the file pointers to the global thread structures
    
    @param sw
      Parsed command-line parameters given by the user
*/
void
init_files (struct user_args * sw)
{
  char                * out[4];
  /* construct output file names */
  out[0] = makefilename (sw->outfile, ".assembled.fastq");
  out[1] = makefilename (sw->outfile, ".unassembled.forward.fastq");
  out[2] = makefilename (sw->outfile, ".unassembled.reverse.fastq");
  out[3] = makefilename (sw->outfile, ".discarded.fastq");

  thr_global.fd[0] = fopen (out[0], "w");
  thr_global.fd[1] = fopen (out[1], "w");
  thr_global.fd[2] = fopen (out[2], "w");
  thr_global.fd[3] = fopen (out[3], "w");

  free (out[0]);
  free (out[1]);
  free (out[2]);
  free (out[3]);
}

void * emp_entry_point (void * data)
{
  struct thread_local_t * thr_local;
  int i, j, sleep, elms;
  struct emp_freq * ef;
  struct read_t * reads[2];
  char * seq;

  /* initialization of values */
  thr_local = (struct thread_local_t *)data;

  ef = (struct emp_freq *) calloc (1, sizeof (struct emp_freq));
  
  while (1)
   {
     // TODO: The first condition in the next line is useless?
     if (thr_local->block == thr_global.xblock && thr_global.io_thread == thr_local->id)
      {
        pthread_mutex_lock (&cs_mutex_io);
        //write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
        // TODO: read_data ();
        elms = db_get_next_reads (thr_global.xblock->fwd, 
                                  thr_global.xblock->rev, 
                                  thr_global.yblock->fwd, 
                                  thr_global.yblock->rev);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("READ %d elms\n", elms);
           printf ("WAKE_UP_ALL!    (reads: %d processed: %d)\n", thr_global.xblock->reads, thr_global.xblock->processed);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        // wakeup threads
        pthread_cond_broadcast (&cs_mutex_cond);
        pthread_mutex_unlock (&cs_mutex_wnd);
      }

     sleep = 1;

     pthread_mutex_lock (&cs_mutex_wnd);
     if (thr_global.xblock->reads == thr_global.xblock->processed)
      {
        if (thr_global.finish) 
         {
           pthread_mutex_unlock (&cs_mutex_wnd);
           break;
         }
        /*
        if (thr_global.finish)
         {
           pthread_mutex_unlock (&cs_mutex_wnd);
           if (!thr_global.xblock->threads)
            {
              write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
            }
           break;
         }
        */
        /* is this the last thread using the current buffer? */
        if (thr_global.xblock->threads == 0 && thr_global.io_thread == -1)
         {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED IO THREAD %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
           thr_global.io_thread = thr_local->id;
           thr_local->block = thr_global.xblock;
           sleep = 0;
         }
        else
         {
           if (assign_reads (thr_global.yblock, thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED y READS to %d   (%d - %d)  Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.yblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
         }
      }
     else
      {
        if (assign_reads (thr_global.xblock,thr_local)) sleep = 0;
           #ifdef PRINT_DEBUG
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
           #endif
      }
     if (sleep)
      {
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("WAKING %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
      }
     pthread_mutex_unlock (&cs_mutex_wnd);

     if (!sleep && thr_global.io_thread != thr_local->id)
      {
        // TODO: Process reads and make this an inline function
        for (i = thr_local->start; i < thr_local->end; ++ i)
         {
           reads[0] = thr_local->block->fwd->reads[i];
           reads[1] = thr_local->block->rev->reads[i];
           for (j = 0; j < 2; ++ j)
            {
              seq = reads[j]->data;
              while (*seq)
               {
                 switch (*seq)
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
                 ++ seq;
               }
            }
         }

        pthread_mutex_lock (&cs_mutex_wnd);
          -- thr_local->block->threads;
           #ifdef PRINT_DEBUG
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Finished %d  (Remaining threads: %d)\n", thr_local->id, thr_local->block->threads);
           pthread_mutex_unlock (&cs_mutex_out);
           #endif
        pthread_mutex_unlock (&cs_mutex_wnd);
      }
  //   pthread_mutex_lock (&cs_mutex_out);
  //     printf ("Thread: %d\n", thr_local->id);
  //   pthread_mutex_unlock (&cs_mutex_out);
   }
//           pthread_mutex_lock (&cs_mutex_out);
           
//  printf ("thread: %2d   A: %d C: %d G: %d T: %d\n", thr_local->id, ef->freqa, ef->freqc, ef->freqg, ef->freqt);
  //         pthread_mutex_unlock (&cs_mutex_out);
  return (ef);
}


int 
main (int argc, char * argv[])
{
  int                   i;
  struct user_args      sw;
  struct emp_freq * ef;
  struct thread_local_t * thr_data;
  int a, c, g, t;
  pthread_t * tid;
  unsigned int elms;

  /* parse command-line arguments */
  if (!decode_switches (argc, argv, &sw))
   {
     /* TODO: Validate reads names */
     usage ();
     return (EXIT_FAILURE);
   }

  a = c = g = t = 0;
  init_thr_global ();
  thr_data = (struct thread_local_t *) calloc (sw.threads, sizeof (struct thread_local_t));
  tid      = (pthread_t *) malloc (sw.threads * sizeof (pthread_t));
  init_files (&sw);

  /* Initialize read buffers */
  init_fastq_reader_double_buffer (sw.fastq_left, 
                                   sw.fastq_right, 
                                   sw.memory, 
                                   thr_global.xblock->fwd, 
                                   thr_global.xblock->rev, 
                                   thr_global.yblock->fwd, 
                                   thr_global.yblock->rev);

  ef = (struct emp_freq *)malloc (sizeof(struct emp_freq));
  if (sw.emp_freqs)
   {
     
     elms = db_get_next_reads (thr_global.yblock->fwd, 
                               thr_global.yblock->rev,
                               thr_global.xblock->fwd,
                               thr_global.xblock->rev);

     thr_global.yblock->reads     = elms;
     thr_global.yblock->processed = 0;

     for (i = 0; i < sw.threads; ++ i)
      {
        thr_data[i].block          = thr_global.xblock;
        thr_data[i].id             = i;
        thr_data[i].sw             = &sw;
        thr_data[i].match_score    = match_score;
        thr_data[i].mismatch_score = mismatch_score;
        thr_data[i].start          = 0;
        thr_data[i].end            = 0;
        pthread_create (&tid[i], NULL, emp_entry_point, (void *)&thr_data[i]); 
      }
     for (i = 0; i < sw.threads; ++ i)
      {
         pthread_join (tid[i], (void **)&ef);
         a += ef->freqa; c += ef->freqc; g += ef->freqg; t += ef->freqt;
//         printf ("tid: %d ef->a: %d ef->c: %d ef->g: %d ef->t: %d\n", (int)tid[i], ef->freqa, ef->freqc, ef->freqg, ef->freqt);
      }
     ef->freqa = a; ef->freqc = c; ef->freqg = g; ef->freqt = t;
     ef->total = a + c + g + t;
     ef->pa = ef->freqa / ef->total; ef->pc = ef->freqc / ef->total; ef->pg = ef->freqg / ef->total; ef->pt = ef->freqt / ef->total;
     ef->q  = ef->pa * ef->pa + ef->pc * ef->pc + ef->pg * ef->pg + ef->pt * ef->pt;
     printf ("a: %f c: %f g: %f t: %f\n", ef->pa, ef->pc, ef->pg, ef->pt);
  rewind_files ();
  thr_global.xblock->fwd->unread = thr_global.xblock->rev->unread = NULL;
  thr_global.yblock->fwd->unread = thr_global.yblock->rev->unread = NULL;
  thr_global.yblock->reads = thr_global.xblock->reads = 0;
  thr_global.yblock->processed = thr_global.xblock->processed = 0;
  thr_global.finish = 0;
  thr_global.xblock->reads     = 0;
  thr_global.xblock->processed = 0;
  thr_global.xblock->threads   = 0;
  thr_global.yblock->reads     = 0;
  thr_global.yblock->processed = 0;
  thr_global.yblock->threads   = 0;

  thr_global.io_thread = -1;
  thr_global.finish = 0;

   }
  else
   {
     ef->freqa = ef->freqc = ef->freqg = ef->freqt = ef->total = ef->pa = ef->pc = ef->pg = ef->pt = ef->q = 0.25;
//     printf ("Set emp freqs to 0.25\n");
   }
  
//  printf ("!!!! A: %d\nC: %d\nG: %d\nT: %d\n", ef->freqa, ef->freqc, ef->freqg, ef->freqt);
//  exit (1);

  init_scores(match_score, mismatch_score, ef);

  elms = db_get_next_reads (thr_global.yblock->fwd, 
                            thr_global.yblock->rev,
                            thr_global.xblock->fwd,
                            thr_global.xblock->rev);

  thr_global.yblock->reads     = elms;
  thr_global.yblock->processed = 0;

  /* pthreads entry point */
  for (i = 0; i < sw.threads; ++ i)
   {
     thr_data[i].block  = thr_global.xblock;
     thr_data[i].id     = i;
     thr_data[i].sw     = &sw;
     thr_data[i].match_score = match_score;
     thr_data[i].mismatch_score = mismatch_score;
     thr_data[i].start  = 0;
     thr_data[i].end    = 0;
     thr_data[i].ef     = ef;
     if (sw.emp_freqs)
       pthread_create (&tid[i], NULL, entry_point_ef, (void *)&thr_data[i]); 
     else
       pthread_create (&tid[i], NULL, entry_point, (void *)&thr_data[i]); 
   }

  for (i = 0; i < sw.threads; ++ i)
   {
     pthread_join (tid[i], NULL);
   }

  free (ef);

  destroy_reader ();
  destroy_thr_global ();
  
  return (EXIT_SUCCESS);
}
