#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fastq.h"
#include "bintest.h"
#include "args.h"
#include "emp.h"

/* possibly implement it such that we keep a number of strings for each diff_cnt */

#define         ALPHA           1
#define         BETA            2
#define         MAX_READ_SIZE   300
#define         PHRED_INIT      33

extern long double * precomp;

struct dp_matrix
 {
   int          cnt_match;
   int          cnt_error;
   double       score;
 };

double sc_eq[256][256];
double sc_neq[256][256];


double qs_mul[256][256];

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

void init_scores (int match, int mismatch)
{
  int           i, j;
  double        ex, ey;

  for (i = 0; i < 256; ++ i)
   {
     for (j = 0; j < 256; ++ j) 
      {
        ex = pow (10.0, - (i - PHRED_INIT) / 10.0);
        ey = pow (10.0, - (j - PHRED_INIT) / 10.0);

        sc_eq[i][j]  =    match * ((1 - ex) * (1 - ey) + (ex * ey) / 3.0);
        sc_neq[i][j] = mismatch * (1 - (1.0 / 3.0) * (1 - ex) * ey - (1.0 / 3.0) * (1 - ey) * ex - (2.0 / 9.0) * ex * ey);
        qs_mul[i][j] = ex * ey;

     }
   }
}

inline void 
scoring (char dleft, char dright, char qleft, char qright, int score_method, double * score)
{
  //double                ex, ey;

  if (dright != 'N' && dleft != 'N' && dright == dleft)
   {
     switch (score_method)
      {
        case 1:
          *score += (qright - PHRED_INIT) + (qleft - PHRED_INIT);
          break;
        case 3:
          //ex      = pow (10.0, - (qright - PHRED_INIT) / 10.0);
          //ey      = pow (10.0, - (qleft  - PHRED_INIT) / 10.0);
          //*score += (1 - ex) * (1 - ey) + (ex * ey) / 3.0;
          *score += sc_eq[(int)qright][(int)qleft];
          break;
        case 4:
          *score += 1;
          break;
      }
   }
  else
   {
     if (dright == 'N' && dleft == 'N')
      {
        switch (score_method)
         {
           case 1:
             *score += abs ((qright - PHRED_INIT) - (qleft - PHRED_INIT));
             break;
           case 3:
             //ex      = pow (10.0, - (qright - PHRED_INIT) / 10.0);
             //ey      = pow (10.0, - (qleft  - PHRED_INIT) / 10.0);
             //*score -= (1 - (1.0 / 3.0) * (1 - ex) * ey - (1.0 / 3.0) * (1 - ey) * ex - (2.0 / 9.0) * ex * ey);
             *score -= sc_neq[(int)qright][(int)qleft];
             break;
           case 4:
             *score -= 1;
             break;
         }
      }
     else
      {
        switch (score_method)
         {
           case 1:
             *score += abs ((qright - PHRED_INIT) - (qleft - PHRED_INIT));
             break;
           case 3:
             *score -= 0.75;
             break;
           case 4:
             *score -= 1;
             break;
         }
      }
   }
}

int
assembly (struct reads_info * left, 
          struct reads_info * right, 
          int               * lps, 
          int               * rps, 
          struct asm_info   * ai,
          int                 score_method,
          double            * uncalled)
{
  int                   i,j;
  int                   n;
  double                score;
  double                best_score = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  int                   nMatch;
  int asm_len = 0;
  *uncalled = 0;
  
  n = strlen (left->data);
  rps[0] = 0;
  lps[0] = 0;

  /* compute the prefix sum of quality scores */
  for (i = 1; i <= n; ++ i)
   {
     lps[i] = lps[i - 1] +  left->quality_score[i - 1] - PHRED_INIT;
     rps[i] = rps[i - 1] + right->quality_score[i - 1] - PHRED_INIT;
   }

  /* compute score for every overlap */
  score = 0;
  for (i = 0; i <= n; ++ i)    /* the size of the overlap */
   {
     /* sum up the quality scores of the nonoverlapping regions */
     switch (score_method)
      {
        case 1:
          score  = lps[n - i];
          score += rps[n] - rps[i];
          break;
        case 3:
        case 4:
          score = 0;
      }
     nMatch = 0;

     /* sum up the quality scores of the overlapping regions */
     for (j = 0; j < i; ++ j)
      {
        scoring (left->data[n - i + j], right->data[j], left->quality_score[n - i + j], right->quality_score[j], score_method, &score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }


     switch (score_method)
      {
        case 1:
          score = score / (2 * n - i);
          break;
        case 3:
          /* binomial test */
          //if (precomp[(i + 1) * i / 2 - 1 + nMatch] < 0.5) score = 0;
          break;
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
        scoring (left->data[j], right->data[n - i + j], left->quality_score[j], right->quality_score[n - i + j], score_method, &score);
        if (left->data[n - i + j] == right->data[j]) ++nMatch;
      }


     switch (score_method)
      {
        case 1:
          score = score / i;
          break;
      }
//     printf ( "runthrough QS: %f overlap: %d\n", score, i );
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
     
     /*
     for (i = 0; i < best_overlap; ++i)
      {
        if (left->quality_score[n - best_overlap + i] > right->quality_score[i])
         {
           ai->data[n - best_overlap + i] = left->data[n - best_overlap + i];
         }
        else
         {
           ai->data[n - best_overlap + i] = right->data[i];
         }

        if (right->data[i] != 'N' && left->data[n - best_overlap + i] != 'N' && right->data[i] == left->data[n - best_overlap + i])
         {
           ai->quality_score[n - best_overlap + i] = (right->quality_score[i] - PHRED_INIT) + (left->quality_score[n - best_overlap + i] - PHRED_INIT) + PHRED_INIT;   // THIS IS OK
         }
        else
         {
           
           ai->quality_score[n - best_overlap + i] = abs ((right->quality_score[i] - PHRED_INIT) - (left->quality_score[n - best_overlap + i] - PHRED_INIT)) + PHRED_INIT;
         }
      }
     */
     
     memcpy (ai->data          + n, right->data          + best_overlap,  n - best_overlap);
     memcpy (ai->quality_score + n, right->quality_score + best_overlap,  n - best_overlap);


     ai->data[2 * n - best_overlap]          = 0;
     ai->quality_score[2 * n - best_overlap] = 0;
     asm_len = 2 * n - best_overlap;
   }
  else
   {
     assemble_overlap (left, right, 0, n - best_overlap, best_overlap, ai);
     
     /*
     for (i = 0; i < best_overlap; ++ i)
      {
        if (left->quality_score[i] > right->quality_score[n - best_overlap + i])
         {
           ai->data[i] = left->data[i];
         }
        else
         {
           ai->data[i] = right->data[n - best_overlap + i];
         }

        if ( right->data[n - best_overlap + i] != 'N' && left->data[i] != 'N' && right->data[n - best_overlap + i] == left->data[i])
         {
           ai->quality_score[i] = (right->quality_score[n - best_overlap + i] - PHRED_INIT) + (left->quality_score[i] - PHRED_INIT) + PHRED_INIT;
         }
        else
         {
           ai->quality_score[i] = abs ((right->quality_score[n - best_overlap + i] - PHRED_INIT) - (left->quality_score[i] - PHRED_INIT)) + PHRED_INIT;
         }
      }
     */

     ai->data[best_overlap]          = 0;
     ai->quality_score[best_overlap] = 0;
     asm_len = best_overlap;
    
   }

  for (j = 0; j < asm_len; ++ j)
    if (ai->data[j] == 'N' || ai->data[j] == 'n') ++(*uncalled);
  *uncalled = (*uncalled) / asm_len;

  return (asm_len);
}

inline void
assemble_overlap (struct reads_info * left, struct reads_info * right, int base_left, int base_right, int ol_size, struct asm_info * ai)
{
  int           i; 
  char          x, y;
  char          qx, qy;

  for (i = 0; i < ol_size; ++i)
   {
     x  = left->data[base_left + i]; 
     y  = right->data[base_right + i];
     //if (x == 0 || y == 0) printf ("ERROR!!!\n");
     qx = left->quality_score[base_left + i];
     qy = right->quality_score[base_right + i];
     if (x == 'N' && y == 'N')
      {
        ai->data[base_left + i]          = 'N';
        ai->quality_score[base_left + i] = ( qx < qy ) ? qx : qy;
      }
     else if (x == 'N')
      {
        ai->data[base_left + i]          = y;
        ai->quality_score[base_left + i] = qy;
      }
     else if (y == 'N')
      {
        ai->data[base_left + i]          = x;
        ai->quality_score[base_left + i] = qx;
      }
     else
      {
        if (x == y)
         {
           ai->data[base_left + i] = x;
           ai->quality_score[base_left + i] = (right->quality_score[base_right + i] - PHRED_INIT) + (left->quality_score[base_left + i] - PHRED_INIT) + PHRED_INIT; //qs_mul[qx][qy];
         }
        else
         {
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
}

char * 
strrev (const char * s)
 {
   char       * rev;
   int          i;
   int          j;

   if (!s) return (NULL);

   rev = strdup (s);
   i   = strlen (s);
   j   = 0;

   while (i)
    {
      rev[j] = s[i - 1];
      ++ j;
      -- i;
    }

   return (rev);
 }

char * 
strcpl (const char * s)
{
  char        * cpl;
  int           i, len;

  if (!s) return (NULL);

  cpl = strdup (s);
  len = strlen (s);

  for (i = 0; i < len; ++ i)
   {
     switch (s[i])
      {
        case 'A':
        case 'a':
                  cpl[i] = 'T';
                  break;
        case 'C':
        case 'c':
                  cpl[i] = 'G';
                  break;
        case 'G':
        case 'g':
                  cpl[i] = 'C';
                  break;
        case 'T':
        case 't':
                  cpl[i] = 'A';
                  break;
        default:
                  cpl[i] = s[i];
      }
   }
  return (cpl);
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
precomp_binom_test (int reads_cnt, int read_len, struct reads_info ** ri_left, struct reads_info ** ri_right, int min_asm_len)
{
  int           i, j, k ;
  int           freqa, freqc, freqg, freqt;
  char          c;

  printf ("=> Count START\n");

  freqa = freqc = freqg = freqt = 0;
  for (i = 0; i < reads_cnt; ++ i)
   {
     for (j = 0; j < read_len; ++ j)
      {
        for (k = 1, c = ri_left[i]->data[j]; k > 0; --k)
         {
           switch (c)
            {
              case 'A':
              case 'a':
                ++freqa;
                break;

              case 'C':
              case 'c':
                ++freqc;
                break;

              case 'G':
              case 'g':
                ++freqg;
                break;

              case 'T':
              case 't':
                ++freqt;
                break;
            }
           c = ri_right[i]->data[j];
         }
      }
   }

  printf ("=> Count END\n");
  return (init_precomp (read_len, freqa, freqc, freqg, freqt, min_asm_len));
}

int 
main (int argc, char * argv[])
{
  int                   i, len;
  struct reads_info  ** ri_left;
  struct reads_info  ** ri_right;
  int                   cnt_reads_left;
  int                   cnt_reads_right;
  char                * rev;
  struct asm_info     * ai;
  int                   read_size;
  int                 * qs_lps;
  int                 * qs_rps;
  char                * out[4];
  FILE                * fd[4];
  int                   asm_len;
  int                 * flags;
  double                p_value, geom_mean, uncalled, uncalled_forward, uncalled_reverse;
  struct user_args      sw;
  struct emp_freq * ef;

  if (!decode_switches (argc, argv, &sw))
   {
     /* TODO: Validate overlap */
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

/*
  printf ("Start precomp %d:\n", sw.min_overlap);
  precomp_binom_test (cnt_reads_right, read_size, ri_left, ri_right, sw.min_overlap);
  printf ("End precomp:\n");
  exit(1);
*/
  ef = get_emp_freq (cnt_reads_right, read_size, ri_left, ri_right);

  /* reverse the right ends */
  for (i = 0; i < cnt_reads_right; ++ i)
   {
     rev = strrev (ri_right[i]->data);
     free (ri_right[i]->data);
     ri_right[i]->data = strcpl (rev);
     free (rev);
     rev = strrev(ri_right[i]->quality_score);
     free (ri_right[i]->quality_score);
     ri_right[i]->quality_score = rev;
   }

  /* allocate memory for the assembled results */
  ai = (struct asm_info *) malloc (cnt_reads_left * sizeof(struct asm_info));
  for (i = 0; i < cnt_reads_left; ++i)
   {
     ai[i].data          = (char *) malloc ((2 * read_size + 1) * sizeof(char));
     ai[i].quality_score = (char *) malloc ((2 * read_size + 1) * sizeof(char));
   }

  init_scores(1, 1);

  flags = (int *) calloc (cnt_reads_left, sizeof (int));

  //#pragma omp parallel shared(ri_left,ri_right,ai) private(qs_lps, qs_rps, i) 
  #pragma omp parallel shared(ri_left,ri_right,ai) private(qs_lps, qs_rps, i, asm_len, uncalled) 
  {
    /* allocate two auxiliary arrays for storing the prefix sum of quality
       scores of the two reads */
    qs_lps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));
    qs_rps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));

    /* flags[i] = 1 (assembled)  0 (discarded) 2 (unassembled) */
    #pragma omp for
    for (i = 0; i < cnt_reads_left; ++ i)
     {
       asm_len = assembly (ri_left[i], ri_right[i], qs_lps, qs_rps, &ai[i], sw.score_method, &uncalled);
       geom_mean = p_value = 100000;
   //    asm_len = strlen (ai[i].data);
       //if ((asm_len >= sw.min_asm_len) && (asm_len <= sw.max_asm_len) && (2 * read_size - asm_len >= sw.min_overlap))
       if (2 * read_size - asm_len >= sw.min_overlap && p_value >= sw.p_value)
        {
          if ( (asm_len >= sw.min_asm_len) && (asm_len <= sw.max_asm_len) && (uncalled <= sw.max_uncalled)) //&& (geom_mean >= sw.geom_mean) )
           {
             flags[i] = 1;    /* read goes to assembled */
           }
          else
           {
             flags[i] = 0;    /* read goes to discarded */
           }
        }
       else
        {
          flags[i] = 2;   /* read is unassembled */
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
     else if (flags[i] == 0)    /* discarded part */
      {
        fprintf (fd[3], "%s\n", ri_left[i]->header);
        fprintf (fd[3], "%s\n+\n%s\n", ri_left[i]->data,  ri_left[i]->quality_score);
        fprintf (fd[3], "%s\n", ri_right[i]->header);
        fprintf (fd[3], "%s\n+\n%s\n", ri_right[i]->data,  ri_right[i]->quality_score);

      }
     else    /* unassembled part */
      {
        len = trim (ri_left[i], sw.qual_thres, read_size, &uncalled_forward);
        len += trim_cpl (ri_right[i], sw.qual_thres, read_size, &uncalled_reverse);
        len = 100;
        if (len <= sw.min_asm_len || (uncalled_forward >= sw.max_uncalled || uncalled_reverse >= sw.max_uncalled)) // && (geom_mean_forward <= sw.geom_mean || geom_mean_reverse <= sw.geom_mean))
         {   /* discarded */
//           printf ("WE HAVE A BUG %f %f %f\n", uncalled_forward, uncalled_reverse, sw.max_uncalled );
           /* Maybe consider printing the untrimmed sequences */
           fprintf (fd[3], "%s\n", ri_left[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", ri_left[i]->data,  ri_left[i]->quality_score);
           fprintf (fd[3], "%s\n", ri_right[i]->header);
           fprintf (fd[3], "%s\n+\n%s\n", ri_right[i]->data, ri_right[i]->quality_score); /* printing the reverse compliment of the original sequence */
         }
        else   /* unassembled */
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
