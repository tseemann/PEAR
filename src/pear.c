#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fastq.h"
#include "args.h"
#include "emp.h"

/* possibly implement it such that we keep a number of strings for each diff_cnt */

#define         ALPHA           1
#define         BETA            2
#define         MAX_READ_SIZE   300
#define         PHRED_INIT      33

extern long double * precomp;

//int stat_test (double, double, int, double);
int stat_test2 (double, double, int, double);

double
assemble_overlap (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai);

/*
double
assemble_overlap_ef (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai, struct emp_freq  * ef);
*/

struct dp_matrix
 {
   int          cnt_match;
   int          cnt_error;
   double       score;
 };

/* TODO: Dynamically allocate them */
double sc_eq[256][256];
double sc_neq[256][256];

double sc_eqA[256][256];
double sc_eqC[256][256];
double sc_eqG[256][256];
double sc_eqT[256][256];

double sc_neqAC[256][256];
double sc_neqAG[256][256];
double sc_neqAT[256][256];

double sc_neqCA[256][256];
double sc_neqCG[256][256];
double sc_neqCT[256][256];

double sc_neqGA[256][256];
double sc_neqGC[256][256];
double sc_neqGT[256][256];

double sc_neqTA[256][256];
double sc_neqTC[256][256];
double sc_neqTG[256][256];

double qs_mul[256][256];

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
  * @param len
  *   Length of the forward sequence
  *
  * @return
  *   Returns the length of the trimmed sequence
  */
int
trim (struct reads_info * read, int min_quality, int len, double * uncalled)
{
  int                   i;

  *uncalled = 0;

  for (i = 0; i < len - 1; ++ i)
   {
     if (read->data[i] == 'N' || read->data[i] == 'n') ++ (*uncalled);
     if ( read->quality_score[i] - PHRED_INIT < min_quality && read->quality_score[i + 1] - PHRED_INIT < min_quality)
      {
        read->quality_score[i + 1] = 0;
        read->data[i + 1] = 0;
        *uncalled = (*uncalled) / (i + 1);
        return i + 1;
      }
   }

  if (read->data[len - 1] ==  'N' || read->data[len - 1] == 'n') ++ (*uncalled);
  *uncalled = (*uncalled) / len;
  return (len);
}

int
trim_cpl (struct reads_info * read, int min_quality, int len, double * uncalled)
{
  int                   i, j;

  *uncalled = 0;

  for (i = len - 1; i > 0; -- i)
   {
     if (read->quality_score[i] - PHRED_INIT < min_quality && read->quality_score[i - 1] - PHRED_INIT < min_quality)
      {
        read->quality_score[i - 1] = 0;
        read->data[i - 1] = 0;
        memmove (read->data, read->data + i, strlen (read->data + i) + 1);
        memmove (read->quality_score, read->quality_score + i, strlen (read->quality_score + i) + 1);
        for (j = 0; read->data[j]; ++ j)
          if (read->data[j] == 'N' || read->data[j] == 'n') ++ (*uncalled);
        *uncalled = (*uncalled) / j;
        return (len - i);
      }
   }
  for (j = 0; read->data[j]; ++ j)
    if (read->data[j] == 'N' || read->data[j] == 'n') ++ (*uncalled);
  *uncalled = (*uncalled) / j;
  return (len);
}

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

        sc_neqAC[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pc / pcgt) - (1 - ex) * ey * (ef->pa / pagt) - ex * ey * (pg2 + pt2) / pgt2);
        sc_neqCA[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pa / pagt) - (1 - ex) * ey * (ef->pc / pcgt) - ex * ey * (pg2 + pt2) / pgt2);

        sc_neqAG[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pg / pcgt) - (1 - ex) * ey * (ef->pa / pact) - ex * ey * (pc2 + pt2) / pct2);
        sc_neqGA[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pa / pact) - (1 - ex) * ey * (ef->pg / pcgt) - ex * ey * (pc2 + pt2) / pct2);

        sc_neqAT[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pt / pcgt) - (1 - ex) * ey * (ef->pa / pacg) - ex * ey * (pc2 + pg2) / pcg2);
        sc_neqTA[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pa / pacg) - (1 - ex) * ey * (ef->pt / pcgt) - ex * ey * (pc2 + pg2) / pcg2);

        sc_neqCG[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pg / pagt) - (1 - ex) * ey * (ef->pc / pact) - ex * ey * (pa2 + pt2) / pat2);
        sc_neqGC[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pc / pact) - (1 - ex) * ey * (ef->pg / pagt) - ex * ey * (pa2 + pt2) / pat2);

        sc_neqCT[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pt / pagt) - (1 - ex) * ey * (ef->pc / pacg) - ex * ey * (pa2 + pg2) / pag2);
        sc_neqTC[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pc / pacg) - (1 - ex) * ey * (ef->pt / pagt) - ex * ey * (pa2 + pg2) / pag2);

        sc_neqGT[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pt / pact) - (1 - ex) * ey * (ef->pg / pacg) - ex * ey * (pa2 + pc2) / pac2);
        sc_neqTG[i][j] = mismatch * (1- (1 - ey) * ex * (ef->pg / pacg) - (1 - ex) * ey * (ef->pt / pact) - ex * ey * (pa2 + pc2) / pac2);
     }
   }
}

inline void
scoring_ef (char dleft, char dright, char qleft, char qright, int score_method, double * score, int match, int mismatch, struct emp_freq * ef)
{
  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += (ef->q * match - (1 - ef->q) * mismatch);
          break;
        case 2:
          *score -= (1 - ef->q) * mismatch; 
          break;
        case 3:
          *score -= mismatch;
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
          break;
        case 2:
          switch (dleft)
           {
             case 'A':
               *score += sc_eqA[(int)qright][(int)qleft];
               break;
             case 'C':
               *score += sc_eqC[(int)qright][(int)qleft];
               break;
             case 'G':
               *score += sc_eqG[(int)qright][(int)qleft];
               break;
             case 'T':
               *score += sc_eqT[(int)qright][(int)qleft];
               break;
           }
          break;
        case 3:
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
          break;
        case 2:
          switch  (dleft)
           {
             case 'A':
               switch (dright)
                {
                  case 'C':
                    *score -= sc_neqAC[(int)qleft][(int)qright];
                    break;
                  case 'G':
                    *score -= sc_neqAG[(int)qleft][(int)qright];
                    break;
                  case 'T':
                    *score -= sc_neqAT[(int)qleft][(int)qright];
                    break;
                }
               break;
             case 'C':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqCA[(int)qleft][(int)qright];
                    break;
                  case 'G':
                    *score -= sc_neqCG[(int)qleft][(int)qright];
                    break;
                  case 'T':
                    *score -= sc_neqCT[(int)qleft][(int)qright];
                    break;
                }
               break;
             case 'G':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqGA[(int)qleft][(int)qright];
                    break;
                  case 'C':
                    *score -= sc_neqGC[(int)qleft][(int)qright];
                    break;
                  case 'T':
                    *score -= sc_neqGT[(int)qleft][(int)qright];
                    break;
                }
               break;
             case 'T':
               switch (dright)
                {
                  case 'A':
                    *score -= sc_neqTA[(int)qleft][(int)qright];
                    break;
                    break;
                  case 'C':
                    *score -= sc_neqTC[(int)qleft][(int)qright];
                    break;
                  case 'G':
                    *score -= sc_neqTG[(int)qleft][(int)qright];
                    break;
                }
               break;
           }
          break;
        case 3:
          *score -= mismatch;
          break;
      }
   }
}

/* TODO: Remember to speed up this function by doing something with the multiplication and division of match/mismatch */
inline void 
scoring (char dleft, char dright, char qleft, char qright, int score_method, double * score, int match, int mismatch)
{
  if (dleft == 'N' || dright == 'N')       /* one of them is N */
   {
     switch (score_method)
      {
        case 1:
          *score += (0.25 * match - (1 - 0.25) * mismatch);
          break;
        case 2:
          *score -= (1 - 0.25) * mismatch; 
          break;
        case 3:
          *score -= mismatch;
          break;
      }
   }
  else if (dleft == dright)     /* equal */
   {
     switch (score_method)
      {
        case 1:
          *score += (sc_eq[(int)qright][(int)qleft] - (1 - sc_eq[(int)qright][(int)qleft] / match) * mismatch);
          break;
        case 2:
          *score += sc_eq[(int)qright][(int)qleft];
          break;
        case 3:
          *score += match;
          break;
      }
   }
  else          /* not equal */
   {
     switch (score_method)
      {
        case 1:
          *score = *score - (sc_neq[(int)qright][(int)qleft] - (1 - sc_neq[(int)qright][(int)qleft] / mismatch) * match);
          break;
        case 2:
          *score -= sc_neq[(int)qright][(int)qleft];
          break;
        case 3:
          *score -= mismatch;
          break;
      }
   }
}

int
assembly_ef (struct reads_info * left, 
          struct reads_info * right, 
          int               * lps, 
          int               * rps, 
          struct asm_info   * ai,
          int                 score_method,
          double            * uncalled,
          double              p_value,
          int                 min_overlap,
          int               * kassian_result,
          int                 match_score,
          int                 mismatch_score,
          int 				  test_method, 	
          struct emp_freq  * ef)
{
  int                   i,j;
  int                   n;
  double                score;
  double                best_score = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  int                   nMatch;
  //double                exp_match; // This is not used any more
  int asm_len = 0;
  *uncalled = 0;
  
  n = strlen (left->data);
  /* compute score for every overlap */
  score = 0;
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {
     nMatch = 0;
     score = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring_ef (left->data[n - i + j], right->data[j], left->quality_score[n - i + j], right->quality_score[j], score_method, &score, match_score, mismatch_score, ef);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        best_overlap = i;
        best_score   = score;
      }
   }

  /* compute for score for runthrough case */
  for (i = n - 1; i > 0; --i)
   {
     score  = 0;
     nMatch = 0;
     for (j = 0; j < i; ++j)
      {
        scoring_ef (left->data[j], right->data[n - i + j], left->quality_score[j], right->quality_score[n - i + j], score_method, &score, match_score, mismatch_score, ef);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        run_through  = 1;
        best_overlap = i;
        best_score   = score;
      }
   }

  /* do the assembly */
  if (!run_through)
   {
     memcpy (ai->data,          left->data,          n - best_overlap);
     memcpy (ai->quality_score, left->quality_score, n - best_overlap);

     assemble_overlap(left, right, n - best_overlap, 0, best_overlap, ai);
          
     memcpy (ai->data          + n, right->data          + best_overlap,  n - best_overlap);
     memcpy (ai->quality_score + n, right->quality_score + best_overlap,  n - best_overlap);

     ai->data[2 * n - best_overlap]          = 0;
     ai->quality_score[2 * n - best_overlap] = 0;
     asm_len = 2 * n - best_overlap;
   }
  else
   {
     assemble_overlap(left, right, 0, n - best_overlap, best_overlap, ai);
     
     ai->data[best_overlap]          = 0;
     ai->quality_score[best_overlap] = 0;
     asm_len = best_overlap;
    
   }

  for (j = 0; j < asm_len; ++ j)
    if (ai->data[j] == 'N' || ai->data[j] == 'n') ++(*uncalled);
    
  *uncalled = (*uncalled) / asm_len;
  
  
  if (test_method == 1){
		*kassian_result = stat_test2 (p_value, best_score, min_overlap, ef->q);
  }else{
        *kassian_result = stat_test2 (p_value, best_score, best_overlap, ef->q);
  }
  
  return (asm_len);
}

int
assembly (struct reads_info * left, 
          struct reads_info * right, 
          int               * lps, 
          int               * rps, 
          struct asm_info   * ai,
          int                 score_method,
          double            * uncalled,
          double              p_value,
          int                 min_overlap,
          int               * kassian_result,
          int                 match_score,
          int                 mismatch_score,
          int 				  test_method)
{
  int                   i,j;
  int                   n;
  double                score;
  double                best_score = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  int                   nMatch;
  //double                exp_match; // this is not used any more
  int                   asm_len = 0;
  *uncalled = 0;
  
  n = strlen (left->data);
     
  /* compute score for every overlap */
  score = 0;
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {     
     nMatch = 0;
     score = 0;
     for (j = 0; j < i; ++ j)
      {
        scoring (left->data[n - i + j], right->data[j], left->quality_score[n - i + j], right->quality_score[j], score_method, &score, match_score, mismatch_score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }
     if (score > best_score)
      {
        best_overlap = i;
        best_score   = score;
      }
   }

  /* compute for score for runthrough case */
  for (i = n - 1; i > 0; --i)
   {
     score  = 0;
     nMatch = 0;
     for (j = 0; j < i; ++j)
      {
        scoring (left->data[j], right->data[n - i + j], left->quality_score[j], right->quality_score[n - i + j], score_method, &score, match_score, mismatch_score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }

     if (score > best_score)
      {
        run_through  = 1;
        best_overlap = i;
        best_score   = score;
      }
   }

  /* do the assembly */
  if (!run_through)
   {
     memcpy (ai->data,          left->data,          n - best_overlap);
     memcpy (ai->quality_score, left->quality_score, n - best_overlap);

     assemble_overlap (left, right, n - best_overlap, 0, best_overlap, ai);
 
     memcpy (ai->data          + n, right->data          + best_overlap,  n - best_overlap);
     memcpy (ai->quality_score + n, right->quality_score + best_overlap,  n - best_overlap);

     ai->data[2 * n - best_overlap]          = 0;
     ai->quality_score[2 * n - best_overlap] = 0;
     asm_len = 2 * n - best_overlap;
   }
  else
   {
     assemble_overlap (left, right, 0, n - best_overlap, best_overlap, ai);
     
     ai->data[best_overlap]          = 0;
     ai->quality_score[best_overlap] = 0;
     asm_len = best_overlap;
    
   }

  for (j = 0; j < asm_len; ++ j)
    if (ai->data[j] == 'N' || ai->data[j] == 'n') ++(*uncalled);
    
  *uncalled = (*uncalled) / asm_len;


  if (test_method == 1){
		*kassian_result = stat_test2 (p_value, best_score, min_overlap, 0.25);
  }else{
        *kassian_result = stat_test2 (p_value, best_score, best_overlap, 0.25);
  }
  
  return (asm_len);
}

double
assemble_overlap (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai)
{
  int           i; 
  char          x, y;
  char          qx, qy;
  //double        exp_match  = 0;

  for (i = 0; i < ol_size; ++i)
   {
     x  = left->data[base_left + i]; 
     y  = right->data[base_right + i];
     qx = left->quality_score[base_left + i];
     qy = right->quality_score[base_right + i];
     if ( (x == 'N' || x == 'n') && (y == 'N' || y == 'n'))
      {
        //exp_match += 0.25; 
        ai->data[base_left + i]          = 'N';
        ai->quality_score[base_left + i] = ( qx < qy ) ? qx : qy;
      }
     else if (x == 'N' || x == 'n')
      {
        //exp_match += 0.25; 
        ai->data[base_left + i]          = y;
        ai->quality_score[base_left + i] = qy;
      }
     else if (y == 'N' || y == 'n')
      {
        //exp_match += 0.25; 
        ai->data[base_left + i]          = x;
        ai->quality_score[base_left + i] = qx;
      }
     else
      {
        if (x == y)
         {
           //exp_match += (sc_eq[(int)qx][(int)qy] / match_score);
           
           ai->data[base_left + i] = x;
           ai->quality_score[base_left + i] = (right->quality_score[base_right + i] - PHRED_INIT) + (left->quality_score[base_left + i] - PHRED_INIT) + PHRED_INIT; //qs_mul[qx][qy];
         }
        else
         {
           //exp_match += (1 - sc_neq[(int)qx][(int)qy] / mismatch_score);
           
           if (qx > qy)
            {
              ai->data[base_left + i]          =  x;
              ai->quality_score[base_left + i] = qx;
            }
           else
            {
              ai->data[base_left + i]          =  y;
              ai->quality_score[base_left + i] = qy;
            }
         }
      }
   }
  //if (ol_size == 0) return (0);
  //return (exp_match / (double)ol_size);
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
     qx = left->quality_score[base_left + i];
     qy = right->quality_score[base_right + i];
     if ( (x == 'N' || x == 'n') && (y == 'N' || y == 'n'))
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = 'N';
        ai->quality_score[base_left + i] = ( qx < qy ) ? qx : qy;
      }
     else if (x == 'N' || x == 'n')
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = y;
        ai->quality_score[base_left + i] = qy;
      }
     else if (y == 'N' || y == 'n')
      {
        exp_match += ef->q; 
        ai->data[base_left + i]          = x;
        ai->quality_score[base_left + i] = qx;
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
           ai->quality_score[base_left + i] = (right->quality_score[base_right + i] - PHRED_INIT) + (left->quality_score[base_left + i] - PHRED_INIT) + PHRED_INIT; //qs_mul[qx][qy];
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
              ai->quality_score[base_left + i] = qx;
            }
           else
            {
              ai->data[base_left + i]          =  y;
              ai->quality_score[base_left + i] = qy;
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

int 
main (int argc, char * argv[])
{
  int                   i, trim_len_fw, trim_len_rev;
  struct reads_info  ** ri_left;
  struct reads_info  ** ri_right;
  int                   cnt_reads_left;
  int                   cnt_reads_right;
  struct asm_info     * ai;
  int                   read_size;
  int                 * qs_lps;
  int                 * qs_rps;
  char                * out[4];
  FILE                * fd[4];
  int                   asm_len;
  int                 * flags;
  double                uncalled, uncalled_forward, uncalled_reverse;
  struct user_args      sw;
  struct emp_freq * ef;
  int kassian_result;
  //struct block_t fwd_block;
  //struct block_t rev_block;

  if (!decode_switches (argc, argv, &sw))
   {
     /* TODO: Validate reads names */
     usage ();
     return (EXIT_FAILURE);
   }

  /* read the two fastq files containing the left-end and right-end reads */
  ri_left  = read_fastq(sw.fastq_left,  &cnt_reads_left);
  ri_right = read_fastq(sw.fastq_right, &cnt_reads_right);

  if (!validate_input (cnt_reads_left, cnt_reads_right))
   {
     return (EXIT_FAILURE);
   }

  read_size = strlen (ri_left[0]->data);

  ef = get_emp_freq (cnt_reads_right, read_size, ri_left, ri_right);

  /* reverse the right ends */
  for (i = 0; i < cnt_reads_right; ++ i)
   {
     mstrrev (ri_right[i]->data);
     mstrcpl (ri_right[i]->data);
     mstrrev(ri_right[i]->quality_score);
   }

  /* allocate memory for the assembled results */
  ai = (struct asm_info *) malloc (cnt_reads_left * sizeof(struct asm_info));
  for (i = 0; i < cnt_reads_left; ++i)
   {
     ai[i].data          = (char *) malloc ((2 * read_size + 1) * sizeof(char));
     ai[i].quality_score = (char *) malloc ((2 * read_size + 1) * sizeof(char));
   }

  init_scores(match_score, mismatch_score, ef);

  flags = (int *) calloc (cnt_reads_left, sizeof (int));

  #pragma omp parallel shared(ri_left,ri_right,ai) private(qs_lps, qs_rps, i, asm_len, uncalled, kassian_result) 
  {
    /* allocate two auxiliary arrays for storing the prefix sum of quality
       scores of the two reads */
    qs_lps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));
    qs_rps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));

    /* flags[i] = 1 (assembled)  0 (discarded) 2 (unassembled) */
    #pragma omp for
    for (i = 0; i < cnt_reads_left; ++ i)
     {
       if (sw.emp_freqs == 0)
        {
          asm_len = assembly (ri_left[i], ri_right[i], qs_lps, qs_rps, &ai[i], sw.score_method, &uncalled, sw.p_value, sw.min_overlap, &kassian_result, match_score, mismatch_score, sw.test);
        }
       else
        {
          asm_len = assembly_ef (ri_left[i], ri_right[i], qs_lps, qs_rps, &ai[i], sw.score_method, &uncalled, sw.p_value, sw.min_overlap, &kassian_result, match_score, mismatch_score, sw.test, ef);
        }

       if (asm_len < read_size)   /* runthrough case */
        {
          if (asm_len >= sw.min_overlap && asm_len >= sw.min_asm_len && asm_len <= sw.max_asm_len && uncalled <= sw.max_uncalled && kassian_result)
           {
             flags[i] = 1;     /* assembled */
           }
          else
           {
             flags[i] = 2;    /* not assembled */
           }
        }
       else     /* normal case */
        {
          if (2 * read_size - asm_len >= sw.min_overlap && asm_len >= sw.min_asm_len && asm_len <= sw.max_asm_len && uncalled <= sw.max_uncalled && kassian_result)
           {
             flags[i] = 1;    /* assembled */
           }
          else
           {
             flags[i] = 2;   /* not assembled */
           }
        }
     }
    free (qs_lps);
    free (qs_rps);
  }
  
  /* construct output file names */
  out[0] = makefilename (sw.outfile, ".assembled.fastq");
  out[1] = makefilename (sw.outfile, ".unassembled.forward.fastq");
  out[2] = makefilename (sw.outfile, ".unassembled.reverse.fastq");
  out[3] = makefilename (sw.outfile, ".discarded.fastq");

  fd[0] = fopen (out[0], "w");
  fd[1] = fopen (out[1], "w");
  fd[2] = fopen (out[2], "w");
  fd[3] = fopen (out[3], "w");

  for (i = 0; i < cnt_reads_left; ++ i)
   {
      if (flags[i] == 1)        /* assembled reads */
      {
        fprintf (fd[0], "%s\n", ri_left[i]->header);
        fprintf (fd[0], "%s\n", ai[i].data);
        fprintf (fd[0], "+\n");
        fprintf (fd[0], "%s\n", ai[i].quality_score);
      }
     else    /* unassembled part */
      {
        trim_len_fw  = trim (ri_left[i], sw.qual_thres, read_size, &uncalled_forward);
        trim_len_rev = trim_cpl (ri_right[i], sw.qual_thres, read_size, &uncalled_reverse);
        if (trim_len_fw < sw.min_trim_len || trim_len_rev < sw.min_trim_len || uncalled_forward >= sw.max_uncalled || uncalled_reverse >= sw.max_uncalled)
         { /* discarded reads*/
           /* Maybe consider printing the untrimmed sequences */
           fprintf (fd[3], "%s\n", ri_left[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", ri_left[i]->data,  ri_left[i]->quality_score);
           fprintf (fd[3], "%s\n", ri_right[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", ri_right[i]->data, ri_right[i]->quality_score); /* printing the reverse compliment of the original sequence */
         }
        else   /* unassembled reads*/
         {
           fprintf (fd[1], "%s\n", ri_left[i]->header);
           fprintf (fd[2], "%s\n", ri_right[i]->header);
           fprintf (fd[1], "%s\n+\n%s\n", ri_left[i]->data,  ri_left[i]->quality_score);
           fprintf (fd[2], "%s\n+\n%s\n", ri_right[i]->data, ri_right[i]->quality_score); /* printing the reverse compliment of the original sequence */
         }
      }
     free (ri_left[i]->header);
     free (ri_left[i]->data);
     free (ri_left[i]->quality_score);
     free (ri_right[i]->data);
     free (ri_right[i]->header);
     free (ri_right[i]->quality_score);
     free (ai[i].data);
     free (ai[i].quality_score);
   }
  free (ri_right);
  free (ri_left);
  free (ai);
  free (ef);


  fclose (fd[0]);
  fclose (fd[1]);
  fclose (fd[2]);
  fclose (fd[3]);
  
  return (0);
}
