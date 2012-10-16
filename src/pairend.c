#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fastq.h"
#include "time.h"

/* possibly implement it such that we keep a number of strings for each diff_cnt */

#define         ALPHA           1
#define         BETA            2
#define         MAX_READ_SIZE   300
#define         PHRED_INIT      33


struct dp_matrix
 {
   int          cnt_match;
   int          cnt_error;
   double       score;
 };

double sc_eq[256][256];
double sc_neq[256][256];

void 
trim (struct reads_info * read, int min_quality, int len)
{
  int                   i;

  for (i = 0; i < len - 1; ++ i)
   {
     if ( read->quality_score[i] - PHRED_INIT < min_quality && read->quality_score[i + 1] - PHRED_INIT < min_quality)
      {
        read->quality_score[i + 1] = 0;
        read->data[i + 1] = 0;
        break;
      }
   }
}


void 
trim_cpl (struct reads_info * read, int min_quality, int len)
{
  int                   i;

  for (i = len - 1; i > 0; -- i)
   {
     if (read->quality_score[i] - PHRED_INIT < min_quality && read->quality_score[i - 1] - PHRED_INIT < min_quality)
      {
        read->quality_score[i - 1] = 0;
        read->data[i - 1] = 0;
        memmove (read->data, read->data + i, strlen (read->data + i) + 1);
        memmove (read->quality_score, read->quality_score + i, strlen (read->quality_score + i) + 1);
        break;
      }
   }
}

void init_scores (void)
{
  int           i, j;
  double        ex, ey;

  for (i = 0; i < 256; ++ i)
   {
     for (j = 0; j < 256; ++ j) 
      {
        ex = pow (10.0, - (i - PHRED_INIT) / 10.0);
        ey = pow (10.0, - (j - PHRED_INIT) / 10.0);

        sc_eq[i][j]  = (1 - ex) * (1 - ey) + (ex * ey) / 3.0;
        sc_neq[i][j] = (1 - (1.0 / 3.0) * (1 - ex) * ey - (1.0 / 3.0) * (1 - ey) * ex - (2.0 / 9.0) * ex * ey);

     }
   }
}

void scoring (char dleft, char dright, char qleft, char qright, int score_method, double * score)
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
     if (dright == 'N' || dleft == 'N')
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
     else
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
   }
}

void 
assembly (struct reads_info * left, 
          struct reads_info * right, 
          int               * lps, 
          int               * rps, 
          struct asm_info   * ai,
          int                 score_method) 
{
  int                   i,j;
  int                   n;
  double                score;
  double                best_score = 0;
  int                   best_overlap = 0;       /* overlap of the best alignment */
  int                   run_through = 0;
  
  n = strlen (left->data);
  rps[0] = 0;
  lps[0] = 0;

  /* compute the prefix sum of quality scores */
  for (i = 1; i <= n; ++ i)
   {
     lps[i] = lps[i - 1] +  left->quality_score[i - 1] - PHRED_INIT;
     rps[i] = rps[i - 1] + right->quality_score[i - 1] - PHRED_INIT;
   }


  for (i = 0; i <= n; ++ i)
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

     /* sum up the quality scores of the overlapping regions */
     for (j = 0; j < i; ++ j)
      {
        scoring (left->data[n - i + j], right->data[j], left->quality_score[n - i + j], right->quality_score[j], score_method, &score);
      }

     switch (score_method)
      {
        case 1:
          score = score / (2 * n - i);
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
     score = 0;
     for (j = 0; j < i; ++j)
      {
        scoring (left->data[j], right->data[n - i + j], left->quality_score[j], right->quality_score[n - i + j], score_method, &score);
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
     memcpy (ai->data,          left->data, n - best_overlap);
     memcpy (ai->quality_score, left->quality_score,  n - best_overlap);

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
           ai->quality_score[n - best_overlap + i] = (right->quality_score[i] - PHRED_INIT) + (left->quality_score[n - best_overlap + i] - PHRED_INIT) + PHRED_INIT;
         }
        else
         {
           ai->quality_score[n - best_overlap + i] = abs ((right->quality_score[i] - PHRED_INIT) - (left->quality_score[n - best_overlap + i] - PHRED_INIT)) + PHRED_INIT;
         }
      }

     memcpy (ai->data          + n, right->data          + best_overlap,  n - best_overlap);
     memcpy (ai->quality_score + n, right->quality_score + best_overlap,  n - best_overlap);

     ai->data[2 * n - best_overlap]          = 0;
     ai->quality_score[2 * n - best_overlap] = 0;
   }
  else
   {
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

     ai->data[best_overlap]          = 0;
     ai->quality_score[best_overlap] = 0;
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

void makefilename (char * fn, const char * prefix, time_t * t)
{
  struct tm *tmp;

  tmp=localtime (t);
  strcpy (fn, prefix);
  strftime (fn + strlen(prefix), 50 - strlen(prefix), "%Y%m%d%H%M%S", tmp);
  strcat (fn, ".out");
}



int 
main (int argc, char * argv[])
{
  int                   i;
  struct reads_info  ** ri_left;
  struct reads_info  ** ri_right;
  int                   cnt_reads_left;
  int                   cnt_reads_right;
  char                * rev;
  struct asm_info     * ai;
  FILE                * fd,* fdl, * fdr; 
  int                   read_size;
  int                 * qs_lps;
  int                 * qs_rps;
  int                   min_asm_len, max_asm_len, min_qual_score, score_method, min_overlap;
  char                  fnouta[50];
  char                  fnoutul[50];
  char                  fnoutur[50];
  int                   asm_len;
  time_t t;


  t = time (NULL);

  if (argc != 8)
   {
     fprintf(stderr, "Syntax: %s [LEFT-END-READS-FILE] [RIGHT-END-READS-FILE] [MIN-ASSEMBLY-LEN] [MAX-ASSEMBLY-LEN] [QUALITY-SCORE-THRESHOLD]\n", argv[0]);
     return (1);
   }
  
  min_asm_len    = atoi(argv[3]);
  max_asm_len    = atoi(argv[4]);
  min_qual_score = atoi(argv[5]);
  score_method   = atoi(argv[6]);
  min_overlap    = atoi(argv[7]);

  /* read the two fastq files containing the left-end and right-end reads */
  ri_left  = read_fastq(argv[1], &cnt_reads_left);
  ri_right = read_fastq(argv[2], &cnt_reads_right);

  if (!validate_input (cnt_reads_left, cnt_reads_right))
   {
     return (1);
   }

  read_size = strlen (ri_left[0]->data);

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

  init_scores();


  //#pragma omp parallel shared(ri_left,ri_right,ai) private(qs_lps, qs_rps, i) 
  #pragma omp parallel shared(ri_left,ri_right,ai) private(qs_lps, qs_rps, i) 
  {
    /* allocate two auxiliary arrays for storing the prefix sum of quality
       scores of the two reads */
    qs_lps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));
    qs_rps = (int *) calloc (strlen (ri_left[0]->data) + 1, sizeof(int));

    #pragma omp for
    for (i = 0; i < cnt_reads_left; ++ i)
     {
       assembly (ri_left[i], ri_right[i], qs_lps, qs_rps, &ai[i], score_method);
     }
    free (qs_lps);
    free (qs_rps);
  }

  
  /* construct output file names based on the time the application started */
  makefilename (fnouta, "asm", &t);
  makefilename (fnoutul, "unasm-left", &t);
  makefilename (fnoutur, "unasm-right", &t);
  fd  = fopen (fnouta,  "w");
  fdl = fopen (fnoutul, "w");
  fdr = fopen (fnoutur, "w");

  for (i = 0; i < cnt_reads_left; ++ i)
   {
     asm_len =  strlen(ai[i].data);
     if ((asm_len >= min_asm_len) && (asm_len <= max_asm_len) && (2 * read_size - asm_len >= min_overlap))
      {
        fprintf (fd, "%s\n", ri_left[i]->header);
        fprintf (fd, "%s\n", ai[i].data);
        fprintf (fd, "+\n");
        fprintf (fd, "%s\n", ai[i].quality_score);
      }
     else
      {
        fprintf (fdl, "%s\n", ri_left[i]->header);
        fprintf (fdr, "%s\n", ri_right[i]->header);
        trim (ri_left[i], min_qual_score, read_size);
        trim_cpl (ri_right[i], min_qual_score, read_size);
        fprintf (fdl, "%s\n+\n%s\n", ri_left[i]->data,  ri_left[i]->quality_score);
        fprintf (fdr, "%s\n+\n%s\n", ri_right[i]->data, ri_right[i]->quality_score);
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
  fclose (fd);
  fclose (fdl);
  fclose (fdr);
  
  return (0);
}
