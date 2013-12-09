#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include "args.h"
#include "emp.h"
#include "reader.h"
#include "async.h"

/** @file pear-pt.c
    @brief Main file containing scoring and assembly related functions (pthreads version)
*/
#define         THREAD_MIN_PACKET_SIZE           500
#define         THREAD_PACKET_SIZE_DELTA         20

#define         PEAR_MATCH_SCORE                 1
#define         PEAR_MISMATCH_SCORE              1

#define         NUM_OF_OUTFILES                  4

#define         PEAR_READ_UNASSEMBLED            3
#define         PEAR_READ_DISCARDED              2
#define         PEAR_READ_ASSEMBLED              1

#define         PEAR_DECODE_OUT_TYPE(x)          (*((x->data) - 1))
#define         PEAR_RESET_OUT_TYPE(x)           *((x->data) - 1) = 0
#define         PEAR_DECODE_ASM_TYPE(x)          (*((x->qscore) - 1))
#define         PEAR_RESET_ASM_TYPE(x)           *((x->qscore) - 1) = 0
#define         PEAR_SET_ASM_TYPE(x,y)           *((x->qscore) - 1) = y

#define         PEAR_READ_OUT_SINGLE             0
#define         PEAR_READ_OUT_BOTH               1
#define         PEAR_SET_OUT_TYPE(x,y)           *((x->data) - 1) = y

#define         PEAR_MERGE_NO_TRIM              1
#define         PEAR_MERGE_TRIM_FORWARD         2
#define         PEAR_MERGE_TRIM_BOTH            3
#define         PEAR_MERGE_TRIM_REVERSE         4


extern void print_number (size_t x);

static int trim_cpl (fastqRead * read, struct user_args * sw, double * uncalled);
static int trim (fastqRead * read, struct user_args * sw, double * uncalled);
static void init_scores (int phred_base, struct emp_freq * ef);
static char * makefilename (const char * prefix, const char * suffix);

static char * outfile_extensions[NUM_OF_OUTFILES] = { ".assembled.fastq", 
                                                      ".unassembled.forward.fastq", 
                                                      ".unassembled.reverse.fastq", 
                                                      ".discarded.fastq" };

static char * sanityCheckMessage[2] = {
 "\n[!!] Forward reads file contains more lines. Merging is done line-by-line on both files and remaining reads in forward file are ignored. You are strongly advised to check that corresponding paired-end reads are located at the same line numbers in your files.\n\n",
 "\n[!!] Reverse reads file contains more lines. Merging is done line-by-line on both files and remaining reads in reverse file are ignored. You are strongly advised to check that corresponding paired-end reads are located at the same line numbers in your files.\n\n" };

static pthread_mutex_t cs_mutex_wnd  = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t cs_mutex_io   = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  cs_mutex_cond = PTHREAD_COND_INITIALIZER;

/* used for locking the screen when debugging */
#ifdef __DEBUG__
static pthread_mutex_t cs_mutex_out  = PTHREAD_MUTEX_INITIALIZER;
#endif

struct thread_global_t thr_global;

int inputFileSanity = 0;        /* Check on whether forward/reverse files have the same number of reads */


static unsigned long g_count_assembled   = 0;
static unsigned long g_count_discarded   = 0;
static unsigned long g_count_unassembled = 0;
static unsigned long g_count_total       = 0;



//int stat_test (double, double, int, double);
int stat_test2 (double, double, int, double);

double
assemble_overlap (fastqRead * left, fastqRead * right, int base_left, int base_right, int ol_size, fastqRead * ai, int phred_base);

/* TODO: 1. Dynamically allocate them 
 *       2. Change 256 to the allowed range and compute only the necessary values */
double      sc_eq[256][256];
double     sc_neq[256][256];
double     sc_eqA[256][256];
double     sc_eqC[256][256];
double     sc_eqG[256][256];
double     sc_eqT[256][256];
double   sc_neqAC[256][256];
double   sc_neqAG[256][256];
double   sc_neqAT[256][256];
double   sc_neqCG[256][256];
double   sc_neqCT[256][256];
double   sc_neqGT[256][256];

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
static int
trim (fastqRead * read, struct user_args * sw, double * uncalled)
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
     if (*(data - 1) - sw->phred_base < sw->qual_thres && *data - sw->phred_base < sw->qual_thres)
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
static int
trim_cpl (fastqRead * read, struct user_args * sw, double * uncalled)
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
     if (qscore[i] - sw->phred_base < sw->qual_thres && qscore[i - 1] - sw->phred_base < sw->qual_thres)
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
static void init_scores (int phred_base, struct emp_freq * ef)
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
        ex = pow (10.0, - (i - phred_base) / 10.0);
        ey = pow (10.0, - (j - phred_base) / 10.0);

        sc_eq[i][j]  =    PEAR_MATCH_SCORE * ((1 - ex) * (1 - ey) + (ex * ey) / 3.0);
        sc_neq[i][j] = PEAR_MISMATCH_SCORE * (1 - (1.0 / 3.0) * (1 - ex) * ey - (1.0 / 3.0) * (1 - ey) * ex - (2.0 / 9.0) * ex * ey);
        //qs_mul[i][j] = ex * ey;

        sc_eqA[i][j] = PEAR_MATCH_SCORE * (1 - ex) * (1 - ey) + (ex * ey) * pcgt2 / p2cgt;
        sc_eqC[i][j] = PEAR_MATCH_SCORE * (1 - ex) * (1 - ey) + (ex * ey) * pagt2 / p2agt;
        sc_eqG[i][j] = PEAR_MATCH_SCORE * (1 - ex) * (1 - ey) + (ex * ey) * pact2 / p2act;
        sc_eqT[i][j] = PEAR_MATCH_SCORE * (1 - ex) * (1 - ey) + (ex * ey) * pacg2 / p2acg;

        sc_neqAC[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pc / pcgt) - (1 - ex) * ey * (ef->pa / pagt) - ex * ey * (pg2 + pt2) / pgt2);
        //sc_neqCA[i][j] PEAR_MISMATCH_SCOREch * (1 - (1 - ey) * ex * (ef->pa / pagt) - (1 - ex) * ey * (ef->pc / pcgt) - ex * ey * (pg2 + pt2) / pgt2);

        sc_neqAG[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pg / pcgt) - (1 - ex) * ey * (ef->pa / pact) - ex * ey * (pc2 + pt2) / pct2);
        //sc_neqGA[i][j] PEAR_MISMATCH_SCOREch * (1 - (1 - ey) * ex * (ef->pa / pact) - (1 - ex) * ey * (ef->pg / pcgt) - ex * ey * (pc2 + pt2) / pct2);

        sc_neqAT[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pt / pcgt) - (1 - ex) * ey * (ef->pa / pacg) - ex * ey * (pc2 + pg2) / pcg2);
        //sc_neqTA[i][j] PEAR_MISMATCH_SCOREch * (1 - (1 - ey) * ex * (ef->pa / pacg) - (1 - ex) * ey * (ef->pt / pcgt) - ex * ey * (pc2 + pg2) / pcg2);

        sc_neqCG[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pg / pagt) - (1 - ex) * ey * (ef->pc / pact) - ex * ey * (pa2 + pt2) / pat2);
        //sc_neqGC[i][j] PEAR_MISMATCH_SCOREch * (1 - (1 - ey) * ex * (ef->pc / pact) - (1 - ex) * ey * (ef->pg / pagt) - ex * ey * (pa2 + pt2) / pat2);

        sc_neqCT[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pt / pagt) - (1 - ex) * ey * (ef->pc / pacg) - ex * ey * (pa2 + pg2) / pag2);
        //sc_neqTC[i][j] PEAR_MISMATCH_SCOREch * (1 - (1 - ey) * ex * (ef->pc / pacg) - (1 - ex) * ey * (ef->pt / pagt) - ex * ey * (pa2 + pg2) / pag2);

        sc_neqGT[i][j] = PEAR_MISMATCH_SCORE * (1 - (1 - ey) * ex * (ef->pt / pact) - (1 - ex) * ey * (ef->pg / pacg) - ex * ey * (pa2 + pc2) / pac2);
        //sc_neqTG[i][j] = mismatch * (1 - (1 - ey) * ex * (ef->pg / pacg) - (1 - ex) * ey * (ef->pt / pact) - ex * ey * (pa2 + pc2) / pac2);

     }
   }
}

#ifdef EXPERIMENTAL
INLINE void
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
#endif

INLINE void
scoring_ef_nm (char dleft, char dright, char qleft, char qright, int score_method, double * score, double * oes, struct emp_freq * ef)
{
  double tmp;

  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += (ef->q - (1 - ef->q));
          *oes    = *score;
          break;
        case 2:
          tmp     = (1 - ef->q);
          *oes   += (ef->q - tmp);
          *score -= tmp; 
          break;
        case 3:
          *score -= 1;
          *oes   += (ef->q - (1 - ef->q));
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
               *score += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft]));
               break;
             case 'C':
               *score += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft]));
               break;
             case 'G':
               *score += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft]));
               break;
             case 'T':
               *score += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft]));
               break;
           }
          *oes    = *score;
          break;
        case 2:
          switch (dleft)
           {
             case 'A':
               *score += sc_eqA[(int)qright][(int)qleft];
               *oes += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft]));
               break;
             case 'C':
               *score += sc_eqC[(int)qright][(int)qleft];
               *oes   += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft]));
               break;
             case 'G':
               *score += sc_eqG[(int)qright][(int)qleft];
               *oes   += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft]));
               break;
             case 'T':
               *score += sc_eqT[(int)qright][(int)qleft];
               *oes   += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft]));
               break;
           }
          break;
        case 3:
          switch (dleft)
           {
             case 'A':
               *oes += (sc_eqA[(int)qright][(int)qleft] - (1 - sc_eqA[(int)qright][(int)qleft]));
               break;
             case 'C':
               *oes += (sc_eqC[(int)qright][(int)qleft] - (1 - sc_eqC[(int)qright][(int)qleft]));
               break;
             case 'G':
               *oes += (sc_eqG[(int)qright][(int)qleft] - (1 - sc_eqG[(int)qright][(int)qleft]));
               break;
             case 'T':
               *oes += (sc_eqT[(int)qright][(int)qleft] - (1 - sc_eqT[(int)qright][(int)qleft]));
               break;
           }
          *score += 1;
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
                    *score = *score - (sc_neqAC[(int)qleft][(int)qright] - (1 - sc_neqAC[(int)qleft][(int)qright]));
                    break;
                  case 'G':
                    *score = *score - (sc_neqAG[(int)qleft][(int)qright] - (1 - sc_neqAG[(int)qleft][(int)qright]));
                    break;
                  case 'T':
                    *score = *score - (sc_neqAT[(int)qleft][(int)qright] - (1 - sc_neqAT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    //*score = *score - (sc_neqCA[(int)qleft][(int)qright] - (1 - sc_neqCA[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqAC[(int)qright][(int)qleft] - (1 - sc_neqAC[(int)qright][(int)qleft]));
                    break;
                  case 'G':
                    *score = *score - (sc_neqCG[(int)qleft][(int)qright] - (1 - sc_neqCG[(int)qleft][(int)qright]));
                    break;
                  case 'T':
                    *score = *score - (sc_neqCT[(int)qleft][(int)qright] - (1 - sc_neqCT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    //*score = *score - (sc_neqGA[(int)qleft][(int)qright] - (1 - sc_neqGA[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqAG[(int)qright][(int)qleft] - (1 - sc_neqAG[(int)qright][(int)qleft]));
                    break;
                  case 'C':
                    //*score = *score - (sc_neqGC[(int)qleft][(int)qright] - (1 - sc_neqGC[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqCG[(int)qright][(int)qleft] - (1 - sc_neqCG[(int)qright][(int)qleft]));
                    break;
                  case 'T':
                    *score = *score - (sc_neqGT[(int)qleft][(int)qright] - (1 - sc_neqGT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    //*score = *score - (sc_neqTA[(int)qleft][(int)qright] - (1 - sc_neqTA[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqAT[(int)qright][(int)qleft] - (1 - sc_neqAT[(int)qright][(int)qleft]));
                    break;
                  case 'C':
                    //*score = *score - (sc_neqTC[(int)qleft][(int)qright] - (1 - sc_neqTC[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqCT[(int)qright][(int)qleft] - (1 - sc_neqCT[(int)qright][(int)qleft]));
                    break;
                  case 'G':
                    //*score = *score - (sc_neqTG[(int)qleft][(int)qright] - (1 - sc_neqTG[(int)qleft][(int)qright]));
                    *score = *score - (sc_neqGT[(int)qright][(int)qleft] - (1 - sc_neqGT[(int)qright][(int)qleft]));
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
                    tmp     = sc_neqAC[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'G':
                    tmp     = sc_neqAG[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'T':
                    tmp     = sc_neqAT[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    tmp     = sc_neqAC[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'G':
                    tmp     = sc_neqCG[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'T':
                    tmp     = sc_neqCT[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    tmp     = sc_neqAG[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'C':
                    tmp     = sc_neqCG[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'T':
                    tmp     = sc_neqGT[(int)qleft][(int)qright];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    tmp     = sc_neqAT[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'C':
                    tmp     = sc_neqCT[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                  case 'G':
                    tmp     = sc_neqGT[(int)qright][(int)qleft];
                    *score -= tmp;
                    *oes    = *oes - (tmp - (1 - tmp));
                    break;
                }
               break;
           }
          break;
        case 3:
          *score -= 1;
           switch  (dleft)
           {
             case 'A':
               switch (dright)
                {
                  case 'C':
                    *oes = *oes - (sc_neqAC[(int)qleft][(int)qright] - (1 - sc_neqAC[(int)qleft][(int)qright]));
                    break;
                  case 'G':
                    *oes = *oes - (sc_neqAG[(int)qleft][(int)qright] - (1 - sc_neqAG[(int)qleft][(int)qright]));
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqAT[(int)qleft][(int)qright] - (1 - sc_neqAT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    //*oes = *oes - (sc_neqCA[(int)qleft][(int)qright] - (1 - sc_neqCA[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqAC[(int)qright][(int)qleft] - (1 - sc_neqAC[(int)qright][(int)qleft]));
                    break;
                  case 'G':
                    *oes = *oes - (sc_neqCG[(int)qleft][(int)qright] - (1 - sc_neqCG[(int)qleft][(int)qright]));
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqCT[(int)qleft][(int)qright] - (1 - sc_neqCT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    //*oes = *oes - (sc_neqGA[(int)qleft][(int)qright] - (1 - sc_neqGA[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqAG[(int)qright][(int)qleft] - (1 - sc_neqAG[(int)qright][(int)qleft]));
                    break;
                  case 'C':
                    //*oes = *oes - (sc_neqGC[(int)qleft][(int)qright] - (1 - sc_neqGC[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqCG[(int)qright][(int)qleft] - (1 - sc_neqCG[(int)qright][(int)qleft]));
                    break;
                  case 'T':
                    *oes = *oes - (sc_neqGT[(int)qleft][(int)qright] - (1 - sc_neqGT[(int)qleft][(int)qright]));
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    //*oes = *oes - (sc_neqTA[(int)qleft][(int)qright] - (1 - sc_neqTA[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqAT[(int)qright][(int)qleft] - (1 - sc_neqAT[(int)qright][(int)qleft]));
                    break;
                  case 'C':
                    //*oes = *oes - (sc_neqTC[(int)qleft][(int)qright] - (1 - sc_neqTC[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqCT[(int)qright][(int)qleft] - (1 - sc_neqCT[(int)qright][(int)qleft]));
                    break;
                  case 'G':
                    //*oes = *oes - (sc_neqTG[(int)qleft][(int)qright] - (1 - sc_neqTG[(int)qleft][(int)qright]));
                    *oes = *oes - (sc_neqGT[(int)qright][(int)qleft] - (1 - sc_neqGT[(int)qright][(int)qleft]));
                    break;
                }
               break;
           }
          break;
      }
   }
}

/* TODO: Remember to speed up this function by doing something with the multiplication and division of match/mismatch */
INLINE void 
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

/* TODO: Remember to speed up this function by doing something with the multiplication and division of match/mismatch */
INLINE void 
scoring_nm (char dleft, char dright, char qleft, char qright, int score_method, double * score, double * oes)
{
  double tmp;

  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += ( - 0.5);
          *oes    = *score;
          break;
        case 2:
          tmp     = 0.75;
          *oes   += ( - 0.5);
          *score -= 0.75; 
          break;
        case 3:
          *score -= 1;
          *oes += ( - 0.5 );
          break;
      }
   }
  else if (dleft == dright)     /* equal */
   {
     switch (score_method)
      {
        case 1:
          *score += (sc_eq[(int)qright][(int)qleft] - (1 - sc_eq[(int)qright][(int)qleft]) );
          *oes    = *score;
          break;
        case 2:
          tmp     = sc_eq[(int)qright][(int)qleft];
          *oes   += (tmp - (1 - sc_eq[(int)qright][(int)qleft] ));
          *score += tmp;
          break;
        case 3:
          *score += 1;
          *oes   += (sc_eq[(int)qright][(int)qleft] - (1 - sc_eq[(int)qright][(int)qleft] ) );
          break;
      }
   }
  else          /* not equal */
   {
     switch (score_method)
      {
        case 1:
          *score = *score - (sc_neq[(int)qright][(int)qleft] - (1 - sc_neq[(int)qright][(int)qleft] ) );
          *oes    = *score;
          break;
        case 2:
          tmp     = sc_neq[(int)qright][(int)qleft];
          *oes    = *oes - (tmp - (1 - sc_neq[(int)qright][(int)qleft] ) );
          *score -= tmp;
          break;
        case 3:
          *score -= 1;
          *oes    = *score - (sc_neq[(int)qright][(int)qleft] - (1 - sc_neq[(int)qright][(int)qleft] ) );
          break;
      }
   }
}

INLINE int
assembly_FORWARD_LONGER (fastqRead * forward, fastqRead * reverse, struct emp_freq * ef, struct user_args  * sw, int nForward, int nReverse)
{
  int 
    i,
    j;
  double
    score = 0,
    oes = 0,
    best_score = 0,
    best_oes = 0,
    uncalled = 0;
  int
    best_overlap = 0,       /* overlap of the best alignment */
    nMatch,
    asm_len = 0,
    st_pass,
    bestScoreCase = 0;

     /*   -------------->
     *                  <------
     *    ....
     *    ------------->
     *           <------
     */
  for (i = 1; i <= nReverse; ++ i)
   {
     nMatch = score = oes = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef_nm (forward->data[nForward - i + j], 
                       reverse->data[j], 
                       forward->qscore[nForward - i + j], 
                       reverse->qscore[j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);

        if (forward->data[nForward - i + j] == reverse->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_NO_TRIM;
        best_overlap  = i;
        best_score    = score;
        best_oes      = oes;
      }
   }
  /* compute score for runthrough case */
  /*   -------------->
  *           <------
  *    ....
  *    ------------->
  *    <------
  */
  for (i = nForward - nReverse  - 1; i >= 0; -- i)
   {
     score = oes = nMatch = 0;
     for (j = 0; j < nReverse; ++ j)
      {
        scoring_ef_nm (forward->data[i + j], 
                       reverse->data[j], 
                       forward->qscore[i + j], 
                       reverse->qscore[j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[i + j] == reverse->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_TRIM_FORWARD;
        best_overlap  = i;
        best_score    = score;
        best_oes      = oes;
      }
   }
  /* compute score for runthrough case */
  /*          -------------->
  *          <------
  *    ....
  *          ------------->
  *    <------
  */
  for (i = nReverse - 1; i > 0; -- i)
   {
     score = oes = nMatch = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef_nm (forward->data[j], 
                       reverse->data[nReverse - i  + j], 
                       forward->qscore[j], 
                       reverse->qscore[nReverse - i + j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[j] == reverse->data[nReverse - i + j]) ++nMatch;
      }
     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_TRIM_BOTH;
        best_overlap  = i;
        best_score    = score;
        best_oes      = oes;
      }
   }

  if (sw->test == 1) 
    st_pass = stat_test2 (sw->p_value, best_oes, sw->min_overlap, ef->q);
  else
    st_pass = stat_test2 (sw->p_value, best_oes, best_overlap, ef->q);

  if (!st_pass) return (0);
  
  /* TODO PART */

  switch (bestScoreCase)
   {
     case PEAR_MERGE_NO_TRIM:
       if (best_overlap == 0)
        {
          assert(0);
          asm_len = nForward + nReverse;

          for (j = 0; j < nForward; ++ j)
            if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;
          for (j = 0; j < nReverse; ++ j)
            if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;
          uncalled /= (nForward + nReverse);

          if (0 >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
           {
             PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);
           }
          else
           {
             return (0);
           }
        }
       else
        {
          asm_len = nForward + nReverse - best_overlap;

          /* count uncalled bases in the non-overlapping high-quality part of the forward read */
          for (j = 0; j < nForward - best_overlap; ++ j)
            if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;

          /* count uncalled bases in the non-overlapping high-quality part of the reverse read */
          for (j = best_overlap; j < nReverse; ++ j)
            if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;

          /* count the uncalled bases in the overlapping part of the two reads */
          for (j = nForward - best_overlap; j < nForward; ++ j)
            if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[j - nForward + best_overlap] == 'N' || reverse->data[j - nForward + best_overlap] == 'n'))  ++uncalled;
          uncalled /= asm_len;

          if (nForward + nReverse - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
           {
             if (asm_len > nForward)
               PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);   /* The merged read will require both memory buffers (forward + reverse) */
             else
               PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE); /* The merged read will fit inside the forward read mem buffer */
               
             
             assemble_overlap (forward, reverse, nForward - best_overlap, 0, best_overlap, forward, sw->phred_base);
             memmove (reverse->data,   reverse->data   + best_overlap,  nReverse - best_overlap);
             memmove (reverse->qscore, reverse->qscore + best_overlap,  nReverse - best_overlap);

             reverse->data[nReverse   - best_overlap] = 0;
             reverse->qscore[nReverse - best_overlap] = 0;
           }
          else
           {
             return (0);
           }
        }
       break;
     case PEAR_MERGE_TRIM_FORWARD:
       /* note that here best_overlap refers to the position of the first byte of the overlapping region
          in the forward read */
       asm_len = best_overlap + nReverse;

       /* count uncalled bases in the non-overlapping high-quality part of the forward read */
       for (j = 0; j < best_overlap; ++ j)
         if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;

       /* count the uncalled bases in the overlapping part of the two reads */
       for (j = best_overlap; j < best_overlap + nReverse; ++ j)
         if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[j - nReverse - best_overlap] == 'N' || reverse->data[j - nReverse - best_overlap] == 'n'))  ++uncalled;
       uncalled /= asm_len;

       if (nReverse >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
        {
          PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE); /* The merged read will fit inside the forward read mem buffer */
          
          assemble_overlap (forward, reverse, best_overlap, 0, nReverse, forward, sw->phred_base);

          forward->data[best_overlap + nReverse]   = 0;
          forward->qscore[best_overlap + nReverse] = 0;
        }
       else
        {
          return (0);
        }
       break;
     case PEAR_MERGE_TRIM_BOTH:
       asm_len = best_overlap;
       
       /* compute uncalled */
       for (j = 0; j < asm_len; ++ j)
         if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[nReverse - best_overlap + j] == 'N' || reverse->data[nReverse - best_overlap + j] == 'n')) ++uncalled;
       uncalled /= asm_len;

       if (asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
        {
          assemble_overlap (forward, reverse, 0, nReverse - best_overlap, best_overlap, forward, sw->phred_base);
          
          forward->data[best_overlap]   = 0;
          forward->qscore[best_overlap] = 0;
          PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE);   /* flag that it's one piece */
        }
       else
        {
          return (0);
        }
       break;

   }

  return (1);

}

INLINE int
assembly_READS_EQUAL (fastqRead * forward, fastqRead * reverse, struct emp_freq * ef, struct user_args  * sw, int n)
{
  int 
    i,
    j;
  double
    score = 0,
    oes = 0,
    best_score = 0,
    best_oes = 0,
    uncalled = 0;
  int
    best_overlap = 0,       /* overlap of the best alignment */
    nMatch,
    asm_len = 0,
    st_pass,
    run_through = 0;

  /* compute score for every overlap */
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {     
     nMatch = score = oes = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef_nm (forward->data[n - i + j], 
                       reverse->data[j], 
                       forward->qscore[n - i + j], 
                       reverse->qscore[j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[n - i + j] == reverse->data[j]) ++nMatch;
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
     score = oes = nMatch = 0;
     for (j = 0; j < i; ++j)
      {
        scoring_ef_nm (forward->data[j],
                       reverse->data[n - i + j], 
                       forward->qscore[j], 
                       reverse->qscore[n - i + j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[n - i + j] == reverse->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        run_through  = 1;
        best_overlap = i;
        best_score   = score;
        best_oes     = oes;
      }
   }

  /** NOW LET'S MERGE */
  if (sw->test == 1)
   {
     st_pass = stat_test2 (sw->p_value, best_oes, sw->min_overlap, ef->q);
   }
  else
   {
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
          if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;
        for (j = 0; j < n; ++ j)
          if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);
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
          if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[j] == 'N' || reverse->data[j] == 'n')) ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE);
           assemble_overlap (forward, reverse, 0, 0, n, forward, sw->phred_base);
           forward->data[n]   = 0;
           forward->qscore[n] = 0;
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
          if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;
        for (j = n - best_overlap; j < n; ++ j)
          if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[j - n + best_overlap] == 'N' || reverse->data[j - n + best_overlap] == 'n'))  ++uncalled;
        for (j = best_overlap; j < n; ++ j)
          if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if (2 * n - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);
           
           assemble_overlap (forward, reverse, n - best_overlap, 0, best_overlap, forward, sw->phred_base);
           memmove (reverse->data,   reverse->data   + best_overlap,  n - best_overlap);
           memmove (reverse->qscore, reverse->qscore + best_overlap,  n - best_overlap);
           /* THIS IS WRONG */
           //memcpy (reverse->data,   reverse->data   + best_overlap,  n - best_overlap);
           //memcpy (reverse->qscore, reverse->qscore + best_overlap,  n - best_overlap);

           reverse->data[n   - best_overlap] = 0;
           reverse->qscore[n - best_overlap] = 0;
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
       if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[n - best_overlap + j] == 'N' || reverse->data[n - best_overlap + j] == 'n')) ++uncalled;
     uncalled /= asm_len;

     if (asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
      {
        assemble_overlap (forward, reverse, 0, n - best_overlap, best_overlap, forward, sw->phred_base);
        
        forward->data[best_overlap]   = 0;
        forward->qscore[best_overlap] = 0;
        PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE);   /* flag that it's one piece */
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

INLINE int
assembly_REVERSE_LONGER (fastqRead * forward, fastqRead * reverse, struct emp_freq * ef, struct user_args  * sw, int nForward, int nReverse)
{
  int 
    i,
    j;
  double
    score = 0,
    oes = 0,
    best_score = 0,
    best_oes = 0,
    uncalled = 0;
  int
    best_overlap = 0,       /* overlap of the best alignment */
    nMatch,
    asm_len = 0,
    st_pass,
    bestScoreCase = 0;

  /*   -------->
  *            <------------
  *    ....
  *    ------->
  *    <------------
  */
  for (i = 0; i <= nForward; ++ i)
   {
     nMatch = score = oes = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef_nm (forward->data[nForward - i + j], 
                       reverse->data[j], 
                       forward->qscore[nForward - i + j], 
                       reverse->qscore[j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[nForward - i + j] == reverse->data[j]) ++nMatch;
      }
     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_NO_TRIM;
        best_overlap  = i;
        best_score    = score;
        best_oes      = oes;
      }
   }
  /*           -------->
  *           <------------
  *    ....
  *         ------->
  *    <------------
  */
  for (i = nReverse - nForward  - 1; i >= 0; -- i)
   {
     score = oes = nMatch = 0;
     for (j = 0; j < nForward; ++ j)
      {
        scoring_ef_nm (forward->data[j], 
                       reverse->data[nReverse - nForward - i + j], 
                       forward->qscore[j], 
                       reverse->qscore[nReverse - nForward - i + j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[j] == reverse->data[nReverse - nForward - i + j]) ++nMatch;
      }
     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_TRIM_REVERSE;
        best_overlap  = i;  /* now this is the length of the remaining part of the reverse read on the right of the overlapping region  */
        best_score    = score;
        best_oes      = oes;
      }
   }
  /*               -------->
  *           <------------
  *    ....
  *                ------->
  *    <------------
  */
  for (i = nForward - 1; i > 0; -- i)
   {
     score = oes = nMatch = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef_nm (forward->data[j],
                       reverse->data[nReverse - i  + j], 
                       forward->qscore[j], 
                       reverse->qscore[nReverse - i + j], 
                       sw->score_method, 
                       &score, 
                       &oes, 
                       ef);
        if (forward->data[j] == reverse->data[nReverse - i + j]) ++nMatch;
      }
     if (score > best_score)
      {
        bestScoreCase = PEAR_MERGE_TRIM_BOTH;
        best_overlap  = i;
        best_score    = score;
        best_oes      = oes;
      }
   }

  /* do a statistical test */
  if (sw->test == 1) 
    st_pass = stat_test2 (sw->p_value, best_oes, sw->min_overlap, ef->q);
  else
    st_pass = stat_test2 (sw->p_value, best_oes, best_overlap, ef->q);

  if (!st_pass) return (0);

  switch (bestScoreCase)
   {
     case PEAR_MERGE_NO_TRIM:
       if (best_overlap == 0)
        {
          assert(0);
          asm_len = nForward + nReverse;

          for (j = 0; j < nForward; ++ j)
            if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;
          for (j = 0; j < nReverse; ++ j)
            if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;
          uncalled /= (nForward + nReverse);

          if (0 >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
           {
             PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);
           }
          else
           {
             return (0);
           }
        }
       else
        {
          asm_len = nForward + nReverse - best_overlap;

          /* count uncalled bases in the non-overlapping high-quality part of the forward read */
          for (j = 0; j < nForward - best_overlap; ++ j)
            if (forward->data[j] == 'N' || forward->data[j] == 'n')  ++uncalled;

          /* count uncalled bases in the non-overlapping high-quality part of the reverse read */
          for (j = best_overlap; j < nReverse; ++ j)
            if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;

          /* count the uncalled bases in the overlapping part of the two reads */
          for (j = nForward - best_overlap; j < nForward; ++ j)
            if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[j - nForward + best_overlap] == 'N' || reverse->data[j - nForward + best_overlap] == 'n'))  ++uncalled;
          uncalled /= asm_len;

          if (nForward + nReverse - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
           {
             if (asm_len > nForward)
               PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);   /* The merged read will require both memory buffers (forward + reverse) */
             else
               PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE); /* The merged read will fit inside the forward read mem buffer */
               
             
             assemble_overlap (forward, reverse, nForward - best_overlap, 0, best_overlap, forward, sw->phred_base);
             memmove (reverse->data,   reverse->data   + best_overlap,  nReverse - best_overlap);
             memmove (reverse->qscore, reverse->qscore + best_overlap,  nReverse - best_overlap);

             reverse->data[nReverse   - best_overlap] = 0;
             reverse->qscore[nReverse - best_overlap] = 0;
           }
          else
           {
             return (0);
           }
        }
       break;
     case PEAR_MERGE_TRIM_REVERSE:
       /* note that here best_overlap refers to the length of the remaining part right of the overlapping regionn
          in the reverse read */
       asm_len = best_overlap + nForward;

       /* count uncalled bases in the non-overlapping high-quality part of the reverse read */
       for (j = nReverse - best_overlap; j < nReverse; ++ j)
         if (reverse->data[j] == 'N' || reverse->data[j] == 'n')  ++uncalled;

       /* count the uncalled bases in the overlapping part of the two reads */
       for (j = 0; j < nForward; ++ j)
         if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[nReverse - best_overlap - nForward + j] == 'N' || reverse->data[nReverse - best_overlap -nForward + j] == 'n'))  ++uncalled;
       uncalled /= asm_len;

       if (nForward >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
        {
          PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_BOTH);
          
          assemble_overlap (forward, reverse, 0, nReverse - nForward - best_overlap, nForward, forward, sw->phred_base);

          memmove (reverse->data,   reverse->data + nReverse - best_overlap,   best_overlap);
          memmove (reverse->qscore, reverse->qscore + nReverse - best_overlap, best_overlap);

          reverse->data[best_overlap]   = 0;
          reverse->qscore[best_overlap] = 0;
        }
       else
        {
          return (0);
        }
       break;
     case PEAR_MERGE_TRIM_BOTH:
       asm_len = best_overlap;
       
       /* compute uncalled */
       for (j = 0; j < asm_len; ++ j)
         if ((forward->data[j] == 'N' || forward->data[j] == 'n') && (reverse->data[nReverse - best_overlap + j] == 'N' || reverse->data[nReverse - best_overlap + j] == 'n')) ++uncalled;
       uncalled /= asm_len;

       if (asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
        {
          assemble_overlap (forward, reverse, 0, nReverse - best_overlap, best_overlap, forward, sw->phred_base);
          
          forward->data[best_overlap]   = 0;
          forward->qscore[best_overlap] = 0;
          PEAR_SET_OUT_TYPE(forward,PEAR_READ_OUT_SINGLE);   /* flag that it's one piece */
        }
       else
        {
          return (0);
        }
       break;
   }

  return (1);
  
}


INLINE int
assembly_ef (fastqRead * forward, fastqRead * reverse, struct emp_freq * ef, struct user_args  * sw)
{
  int
    nForward, 
    nReverse,
    rc;
  
  nForward = strlen (forward->data);
  nReverse = strlen (reverse->data);

  /************************* |FORWARD| > |REVERSE| ********************************/
  if (nForward > nReverse)
   {
     rc = assembly_FORWARD_LONGER (forward, reverse, ef, sw, nForward, nReverse);
   }

  /************************* |FORWARD| == |REVERSE| ********************************/
  else if (nForward == nReverse)
   {
     rc = assembly_READS_EQUAL (forward, reverse, ef, sw, nForward);
   }
  /************************* |FORWARD| < |REVERSE| ********************************/
  else
   {
     rc = assembly_REVERSE_LONGER (forward, reverse, ef, sw, nForward, nReverse);
   }

  return (rc);
}

INLINE int
assembly (fastqRead * left, fastqRead * right, struct user_args  * sw)
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
        //scoring (left->data[n - i + j], right->data[j], left->qscore[n - i + j], right->qscore[j], sw->score_method, &score, &oes, match_score, mismatch_score);
        scoring_nm (left->data[n - i + j], right->data[j], left->qscore[n - i + j], right->qscore[j], sw->score_method, &score, &oes);
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
        //scoring (left->data[j], right->data[n - i + j], left->qscore[j], right->qscore[n - i + j], sw->score_method, &score, &oes, match_score, mismatch_score);
        scoring_nm (left->data[j], right->data[n - i + j], left->qscore[j], right->qscore[n - i + j], sw->score_method, &score, &oes);
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
        asm_len = (n << 1);

        for (j = 0; j < n; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = 0; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if ((n << 1) - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(left,PEAR_READ_OUT_BOTH);
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

        if ((n << 1) - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(left,PEAR_READ_OUT_SINGLE);
           assemble_overlap (left, right, 0, 0, n, left, sw->phred_base);
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
        asm_len = (n << 1) - best_overlap;
        for (j = 0; j < n - best_overlap; ++ j)
          if (left->data[j] == 'N' || left->data[j] == 'n')  ++uncalled;
        for (j = n - best_overlap; j < n; ++ j)
          if ((left->data[j] == 'N' || left->data[j] == 'n') && (right->data[j - n + best_overlap] == 'N' || right->data[j - n + best_overlap] == 'n'))  ++uncalled;
        for (j = best_overlap; j < n; ++ j)
          if (right->data[j] == 'N' || right->data[j] == 'n')  ++uncalled;
        uncalled /= asm_len;

        if ((n << 1) - asm_len >= sw->min_overlap && asm_len >= sw->min_asm_len && asm_len <= sw->max_asm_len && uncalled <= sw->max_uncalled)
         {
           PEAR_SET_OUT_TYPE(left,PEAR_READ_OUT_BOTH);
           assemble_overlap (left, right, n - best_overlap, 0, best_overlap, left, sw->phred_base);
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
        assemble_overlap (left, right, 0, n - best_overlap, best_overlap, left, sw->phred_base);
        
        left->data[best_overlap]   = 0;
        left->qscore[best_overlap] = 0;
        PEAR_SET_OUT_TYPE(left,PEAR_READ_OUT_SINGLE);   /* flag that it's one piece */
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

/** @brief Assemble the overlap

    @param left
    @param right
    
    @param base_left
      Position in the forward read where the overlap starts
    
    @param base_right
      Position in the reverse read where the overlap starts
    
    @param ol_size
      Overlap size

    @param ai
      Buffer where to store the merge read

    @param phred_base
*/

double
assemble_overlap (fastqRead * left, fastqRead * right, int base_left, int base_right, int ol_size, fastqRead * ai, int phred_base)
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
           ai->qscore[base_left + i] = (right->qscore[base_right + i] - phred_base) + (left->qscore[base_left + i] - phred_base) + phred_base; //qs_mul[qx][qy];
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

static char *
makefilename (const char * prefix, const char * suffix)
{
  char * filename;

  filename = (char *) malloc ((strlen(prefix) + strlen(suffix) + 1) * sizeof (char));

  strcpy (filename, prefix);
  strcat (filename, suffix);

  return (filename);
}

/** @brief Write data to file
 *
 *  Write the result of overlapping of a specific read to the corresponding output file
 *
 *  @param fwd
 *     Structure containing the forward reads. Note that this structure is mutable and has
 *     been modified to contain the new read which is the result of the assembly.
 *
 *  @param rev
 *     Structure containing the reverse reads. Note that this structure is mutable and has
 *     been modified to contain the new read which is the result of the assembly.
 *
 *  @param elms
 *     Number of reads in the structure
 *
 *  @param fd
 *     Array of file descriptors of output files
 */
void write_data (fastqRead ** fwd, fastqRead ** rev, unsigned int elms, FILE ** fd)
{
  int i;
  char bothOut;   /* this is set if both fwd[i] and rev[i] contain the resulting assembled read and qscore */
  
  for ( i = 0; i < elms; ++ i)
   {
     bothOut = PEAR_DECODE_OUT_TYPE(fwd[i]);
     PEAR_RESET_OUT_TYPE(fwd[i]);

     if (PEAR_DECODE_ASM_TYPE(fwd[i]) == PEAR_READ_ASSEMBLED)   /* assembled */
      {
        PEAR_RESET_ASM_TYPE(fwd[i]);
        fprintf (fd[0], "%s\n", fwd[i]->header);
        if (!bothOut)
         {
           fprintf (fd[0], "%s\n", fwd[i]->data);
         }
        else
         {
           fprintf (fd[0], "%s",   fwd[i]->data);
           fprintf (fd[0], "%s\n", rev[i]->data);
         }
        fprintf (fd[0], "+\n");

        if (!bothOut)
         {
           fprintf (fd[0], "%s\n", fwd[i]->qscore);
         }
        else
         {
           fprintf (fd[0], "%s",   fwd[i]->qscore);
           fprintf (fd[0], "%s\n", rev[i]->qscore);
         }
        
        ++ g_count_assembled;
      }
     else if (PEAR_DECODE_ASM_TYPE(fwd[i]) == PEAR_READ_DISCARDED)                                            /* discarded */
      {
        PEAR_RESET_ASM_TYPE(fwd[i]);
        /* discarded reads*/
        /* Maybe consider printing the untrimmed sequences */
        fprintf (fd[3], "%s\n", fwd[i]->header);
        fprintf (fd[3], "%s\n+\n%s\n", fwd[i]->data,  fwd[i]->qscore);
        fprintf (fd[3], "%s\n", rev[i]->header);
        fprintf (fd[3], "%s\n+\n%s\n", rev[i]->data, rev[i]->qscore); /* printing the reverse compliment of the original sequence */
        ++ g_count_discarded;
      }
     else   /* unassembled reads*/
      {
        PEAR_RESET_ASM_TYPE(fwd[i]);
        fprintf (fd[1], "%s\n", fwd[i]->header);
        fprintf (fd[2], "%s\n", rev[i]->header);
        fprintf (fd[1], "%s\n+\n%s\n", fwd[i]->data,  fwd[i]->qscore);
        fprintf (fd[2], "%s\n+\n%s\n", rev[i]->data, rev[i]->qscore); /* printing the reverse compliment of the original sequence */
        ++ g_count_unassembled;
      }
   }
}

void flip_list (void)
{
  memBlockInfo * elm1;
  memBlockInfo * elm2;

  elm1 = thr_global.xblock;
  elm2 = thr_global.yblock;

  thr_global.xblock = elm2;
  thr_global.yblock = elm1;
}

INLINE int 
assign_reads (memBlockInfo * block, struct thread_local_t * thr_local)
{
  int r;

  if (block->reads == block->processed) return (0);

  r = THREAD_MIN_PACKET_SIZE + (rand() % THREAD_MIN_PACKET_SIZE);
  thr_local->block   = block;
  thr_local->start   = block->processed;

  if (block->reads - block->processed <= r + THREAD_PACKET_SIZE_DELTA)
   {
     thr_local->end    = block->reads;
     block->processed  = block->reads;
   }
  else
   {
     thr_local->end    = thr_local->start + r;
     block->processed  = thr_local->start + r;
   }
  ++ block->threads;
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
        fprintf (stdout, "."); fflush (stdout);
        write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
        // TODO: read_data ();
        elms = db_get_next_reads (thr_global.xblock->fwd, 
                                  thr_global.xblock->rev, 
                                  thr_global.yblock->fwd, 
                                  thr_global.yblock->rev,
                                  &inputFileSanity);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

#ifdef __DEBUG__
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
           if (!thr_global.xblock->threads)
            {
              thr_global.xblock->threads = 1;
              pthread_mutex_unlock (&cs_mutex_wnd);
//              printf ("!!!!!!!!!! reads: %d processed: %d\n", thr_global.xblock->reads, thr_global.xblock->processed);
//              printf ("!!!!!!!!!! Writing another %d reads\n", thr_global.xblock->reads);
              write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
//              printf ("Finsihed\n");
              break;
            }
           pthread_mutex_unlock (&cs_mutex_wnd);
           break;
         }
        /* is this the last thread using the current buffer? */
        if (thr_global.xblock->threads == 0 && thr_global.io_thread == -1 && thr_global.finish == 0)
         {
#ifdef __DEBUG__
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
#ifdef __DEBUG__
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
#ifdef __DEBUG__
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
#endif
      }
     if (sleep)
      {
#ifdef __DEBUG__
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
#endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
#ifdef __DEBUG__
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
           ass = assembly_ef (thr_local->block->fwd->reads[i], thr_local->block->rev->reads[i], thr_local->ef, thr_local->sw);
           PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],ass);
           if (!ass)
            {
              if (trim     (thr_local->block->fwd->reads[i], thr_local->sw, &uncalled_forward) < thr_local->sw->min_trim_len ||
                  trim_cpl (thr_local->block->rev->reads[i], thr_local->sw, &uncalled_reverse) < thr_local->sw->min_trim_len ||
                  uncalled_forward >= thr_local->sw->max_uncalled || uncalled_reverse >= thr_local->sw->max_uncalled)
               {
                 PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],PEAR_READ_DISCARDED);
               }
              else
               {
                 PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],PEAR_READ_UNASSEMBLED);
               }
            }
         }

        pthread_mutex_lock (&cs_mutex_wnd);
          -- thr_local->block->threads;
#ifdef __DEBUG__
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
{  struct thread_local_t * thr_local;
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
        fprintf (stdout, "."); fflush (stdout);
        write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
        // TODO: read_data ();
        elms = db_get_next_reads (thr_global.xblock->fwd, 
                                  thr_global.xblock->rev, 
                                  thr_global.yblock->fwd, 
                                  thr_global.yblock->rev,
                                  &inputFileSanity);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

#ifdef __DEBUG__
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
           if (!thr_global.xblock->threads)
            {
              thr_global.xblock->threads = 1;
              pthread_mutex_unlock (&cs_mutex_wnd);
              write_data (thr_global.xblock->fwd->reads, thr_global.xblock->rev->reads, thr_global.xblock->reads, thr_global.fd);
            }
           pthread_mutex_unlock (&cs_mutex_wnd);
           break;
         }
        /* is this the last thread using the current buffer? */
        if (thr_global.xblock->threads == 0 && thr_global.io_thread == -1)
         {
#ifdef __DEBUG__
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
#ifdef __DEBUG__
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
#ifdef __DEBUG__
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
#endif
      }
     if (sleep)
      {
#ifdef __DEBUG__
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
#endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
#ifdef __DEBUG__
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
           ass = assembly (thr_local->block->fwd->reads[i], thr_local->block->rev->reads[i], thr_local->sw);
           PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],ass);
           if (!ass)
            {
              if (trim     (thr_local->block->fwd->reads[i], thr_local->sw, &uncalled_forward) < thr_local->sw->min_trim_len ||
                  trim_cpl (thr_local->block->rev->reads[i], thr_local->sw, &uncalled_reverse) < thr_local->sw->min_trim_len ||
                  uncalled_forward >= thr_local->sw->max_uncalled || uncalled_reverse >= thr_local->sw->max_uncalled)
               {
                 PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],PEAR_READ_DISCARDED);
               }
              else
               {
                 PEAR_SET_ASM_TYPE(thr_local->block->fwd->reads[i],PEAR_READ_UNASSEMBLED);
               }
            }
         }

        pthread_mutex_lock (&cs_mutex_wnd);
          -- thr_local->block->threads;
#ifdef __DEBUG__
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
  thr_global.xblock = (memBlockInfo *) calloc (1,sizeof(memBlockInfo));
  thr_global.yblock = (memBlockInfo *) calloc (1,sizeof(memBlockInfo));
  thr_global.xblock->fwd = (memBlock *) calloc (1,sizeof(memBlock));
  thr_global.xblock->rev = (memBlock *) calloc (1,sizeof(memBlock));
  thr_global.yblock->fwd = (memBlock *) calloc (1,sizeof(memBlock));
  thr_global.yblock->rev = (memBlock *) calloc (1,sizeof(memBlock));

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
  int i;
  char * out[4];
  
  /* construct output file names */
  for (i = 0; i < NUM_OF_OUTFILES; ++ i)
   {
     out[i] = makefilename (sw->outfile, outfile_extensions[i]);
     thr_global.fd[i] = fopen (out[i], "w");
     free (out[i]);
   }
}

void * emp_entry_point (void * data)
{
  struct thread_local_t * thr_local;
  int i, j, sleep, elms;
  struct emp_freq * ef;
  fastqRead * reads[2];
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
                                  thr_global.yblock->rev,
                                  &inputFileSanity);
          //read_size = strlen (fwd_block.reads[0]->data);
        pthread_mutex_unlock (&cs_mutex_io);

        pthread_mutex_lock (&cs_mutex_wnd);
        thr_global.xblock->reads      =  elms;
        thr_global.xblock->processed  =  0;
        thr_global.io_thread          = -1;
        thr_global.xblock->threads    =  0;
        if (!elms) thr_global.finish  =  1;
        flip_list ();

#ifdef __DEBUG__
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
#ifdef __DEBUG__
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
#ifdef __DEBUG__
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
#ifdef __DEBUG__
           if (!sleep) {
           pthread_mutex_lock (&cs_mutex_out);
           printf ("ASSIGNED x READS to %d  (%d - %d)   Threads: %d\n", thr_local->id, thr_local->start, thr_local->end, thr_global.xblock->threads);
           pthread_mutex_unlock (&cs_mutex_out);}
#endif
      }
     if (sleep)
      {
#ifdef __DEBUG__
           pthread_mutex_lock (&cs_mutex_out);
           printf ("Sleeping %d\n", thr_local->id);
           pthread_mutex_unlock (&cs_mutex_out);
#endif
        pthread_cond_wait (&cs_mutex_cond, &cs_mutex_wnd);
#ifdef __DEBUG__
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

                    case 'N':
                    case 'n':
                      ++ ef->freqn;
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
#ifdef __DEBUG__
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

void
DisplayInstance (struct user_args * sw)
{
  fprintf (stdout, " ____  _____    _    ____ \n"); 
  fprintf (stdout, "|  _ \\| ____|  / \\  |  _ \\\n");
  fprintf (stdout, "| |_) |  _|   / _ \\ | |_) |\n");
  fprintf (stdout, "|  __/| |___ / ___ \\|  _ <\n");
  fprintf (stdout, "|_|   |_____/_/   \\_\\_| \\_\\\n\n");
  fprintf (stdout, "Forward reads file.................: %s\n", sw->fastq_left);
  fprintf (stdout, "Reverse reads file.................: %s\n", sw->fastq_right);
  fprintf (stdout, "PHRED..............................: %d\n", sw->phred_base);
  fprintf (stdout, "Using empirical frequencies........: %s\n", sw->emp_freqs ? "YES" : "NO");
  fprintf (stdout, "Statistical method.................: %s\n", (sw->test - 1) ? "Acceptance probability" : "OES");
  fprintf (stdout, "Maximum assembly length............: %d\n", sw->max_asm_len);
  fprintf (stdout, "Minimum assembly length............: %d\n", sw->min_asm_len);
  fprintf (stdout, "p-value............................: %f\n", sw->p_value);
  fprintf (stdout, "Quality score threshold (trimming).: %d\n", sw->qual_thres);
  fprintf (stdout, "Minimum read size after trimming...: %d\n", sw->min_trim_len);
  fprintf (stdout, "Maximal ratio of uncalled bases....: %f\n", sw->max_uncalled);
  fprintf (stdout, "Minimum overlap....................: %d\n", sw->min_overlap);
  fprintf (stdout, "Scoring method.....................: ");
  if  (sw->score_method == 1)
   {
     fprintf (stdout, "OES with match = +1 and mismatch = -1\n");
   }
  else if (sw->score_method == 2)
   {
     fprintf (stdout, "Scaled score\n");
   }
  else
   {
     fprintf (stdout, "Ignore quality scores, match = +1 and mismatch = -1\n");
   }
  fprintf (stdout, "Threads............................: %d\n\n", sw->threads);


}

int 
main (int argc, char * argv[])
{
  int                   i;
  struct user_args      sw;
  struct emp_freq * ef;
  struct thread_local_t * thr_data;
  double 
    a = 0, 
    c = 0, 
    g = 0, 
    t = 0, 
    n = 0;
  pthread_t * tid;
  unsigned int 
    blockElements;

  /* parse command-line arguments */
  if (!decode_switches (argc, argv, &sw))
   {
     /* TODO: Validate reads names */
     usage ();
     return (EXIT_FAILURE);
   }

  /* Display PEAR instance settings */
  DisplayInstance (&sw);

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

  if (sw.emp_freqs)
   {
     printf ("Computing empirical frequencies....: ");
     fflush (stdout);
     blockElements = db_get_next_reads (thr_global.yblock->fwd, 
                               thr_global.yblock->rev,
                               thr_global.xblock->fwd,
                               thr_global.xblock->rev,
                               &inputFileSanity);

     thr_global.yblock->reads     = blockElements;
     thr_global.yblock->processed = 0;

     for (i = 0; i < sw.threads; ++ i)
      {
        thr_data[i].block          = thr_global.xblock;
        thr_data[i].id             = i;
        thr_data[i].sw             = &sw;
        thr_data[i].start          = 0;
        thr_data[i].end            = 0;
        pthread_create (&tid[i], NULL, emp_entry_point, (void *)&thr_data[i]); 
      }
     for (i = 0; i < sw.threads; ++ i)
      {
         pthread_join (tid[i], (void **)&ef);
         a += ef->freqa; c += ef->freqc; g += ef->freqg; t += ef->freqt; n += ef->freqn;
         free (ef);
//         printf ("tid: %d ef->a: %d ef->c: %d ef->g: %d ef->t: %d\n", (int)tid[i], ef->freqa, ef->freqc, ef->freqg, ef->freqt);
      }
     ef = (struct emp_freq *)malloc (sizeof(struct emp_freq));
     ef->freqa = a; ef->freqc = c; ef->freqg = g; ef->freqt = t; ef->freqn = n;
     ef->total = a + c + g + t;
     ef->pa = ef->freqa / ef->total; ef->pc = ef->freqc / ef->total; ef->pg = ef->freqg / ef->total; ef->pt = ef->freqt / ef->total;
     ef->q  = ef->pa * ef->pa + ef->pc * ef->pc + ef->pg * ef->pg + ef->pt * ef->pt;
     printf ("DONE\n");
     printf ("  A: %f\n  C: %f\n  G: %f\n  T: %f\n", ef->pa, ef->pc, ef->pg, ef->pt);
     printf ("  %ld uncalled bases\n", ef->freqn);
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

     if (inputFileSanity)
      {
        fprintf (stderr, "%s", sanityCheckMessage[inputFileSanity - 1]);
      }
     inputFileSanity = 0;

   }
  else
   {
     ef = (struct emp_freq *)malloc (sizeof(struct emp_freq));
     ef->freqa = ef->freqc = ef->freqg = ef->freqt = ef->total = ef->pa = ef->pc = ef->pg = ef->pt = ef->q = 0.25;
//     printf ("Set emp freqs to 0.25\n");
   }
  
//  printf ("!!!! A: %d\nC: %d\nG: %d\nT: %d\n", ef->freqa, ef->freqc, ef->freqg, ef->freqt);
//  exit (1);

  init_scores(sw.phred_base, ef);

  fprintf (stdout, "Assemblying reads..................: [");
  fflush (stdout);

  //printf ("END DEBUG TEST\n"); //exit(1);

  blockElements = db_get_next_reads (thr_global.yblock->fwd, 
                            thr_global.yblock->rev,
                            thr_global.xblock->fwd,
                            thr_global.xblock->rev,
                            &inputFileSanity);

  thr_global.yblock->reads     = blockElements;
  thr_global.yblock->processed = 0;

  /* pthreads entry point */
  for (i = 0; i < sw.threads; ++ i)
   {
     thr_data[i].block  = thr_global.xblock;
     thr_data[i].id     = i;
     thr_data[i].sw     = &sw;
//     thr_data[i].match_score = match_score;
//     thr_data[i].mismatch_score = mismatch_score;
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
  printf ("]\n\n");

  g_count_total = g_count_assembled + g_count_discarded + g_count_unassembled;
  printf ("Assembled reads ...................: ");
  print_number (g_count_assembled);
  printf (" / ");
  print_number (g_count_total);
  printf (" (%.3f%%)\n", ((double)g_count_assembled / g_count_total) * 100);

  printf ("Discarded reads ...................: ");
  print_number (g_count_discarded);
  printf (" / ");
  print_number (g_count_total);
  printf (" (%.3f%%)\n", ((double)g_count_discarded / g_count_total) * 100);

  printf ("Not assembled reads ...............: ");
  print_number (g_count_unassembled);
  printf (" / ");
  print_number (g_count_total);
  printf (" (%.3f%%)\n", ((double)g_count_unassembled / g_count_total) * 100);



  printf ("Assembled reads file...............: %s%s\n", sw.outfile, ".assembled.fastq");
  printf ("Discarded reads file...............: %s%s\n", sw.outfile, ".discarded.fastq");
  printf ("Unassembled forward reads file.....: %s%s\n", sw.outfile, ".unassembled.forward.fastq");
  printf ("Unassembled reverse reads file.....: %s%s\n", sw.outfile, ".unassembled.reverse.fastq" );

     if (inputFileSanity)
      {
        fprintf (stderr, "%s", sanityCheckMessage[inputFileSanity - 1]);
      }
  free (ef);
  free (tid);

  destroy_reader ();
  destroy_thr_global ();
  free (thr_data);
  
  return (EXIT_SUCCESS);
}
