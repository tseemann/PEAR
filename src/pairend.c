#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fastq.h"

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

void 
assembly (struct reads_info * left, 
          struct reads_info * right, 
          int               * lps, 
          int               * rps, 
          struct asm_info   * ai)
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

  for (i = 1; i <= n; ++ i)
   {
     lps[i] = lps[i - 1] +  left->quality_score[i - 1] - PHRED_INIT;
     rps[i] = rps[i - 1] + right->quality_score[i - 1] - PHRED_INIT;
   }


  for (i = 0; i <= n; ++ i)
   {
     /* sum up the quality scores of the nonoverlapping regions */
     score  = lps[n - i];
     score += rps[n] - rps[i];

     /* sum up the quality scores of the overlapping regions */
     for (j = 0; j < i; ++ j)
      {
        if (right->data[j] != 'N' && left->data[n - i + j] != 'N' && right->data[j] == left->data[n - i + j])
         {
           score += (right->quality_score[j] - PHRED_INIT) + (left->quality_score[n - i + j] - PHRED_INIT);
         }
        else
         {
           score += abs ((right->quality_score[j] - PHRED_INIT) - (left->quality_score[n - i + j] - PHRED_INIT));
         }
      }

     score = score / (2 * n - i);
 //    printf ( "QS: %f overlap: %d\n", score, i );

     if (score > best_score)
      {
        best_overlap = i;
        best_score   = score;
      }
   }

  for (i = n - 1; i > 0; --i)
   {
     score = 0;
     for (j = 0; j < i; ++j)
      {
        if (left->data[j] != 'N' && right->data[n - i + j] != 'N' && left->data[j] == right->data[n - i + j])
         {
           score += (left->quality_score[j] - PHRED_INIT) + (right->quality_score[n - i + j] - PHRED_INIT);
         }
        else
         {
           score += abs ((left->quality_score[j] - PHRED_INIT) - (right->quality_score[n - i + j] - PHRED_INIT));
         }
      }


     score = score / i;
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
  FILE                * fd;
  int                   read_size;
  int                 * qs_lps;
  int                 * qs_rps;
  int                   min_asm_len, max_asm_len, min_qual_score;


  if (argc != 6)
   {
     fprintf(stderr, "Syntax: %s [LEFT-END-READS-FILE] [RIGHT-END-READS-FILE] [MIN-ASSEMBLY-LEN] [MAX-ASSEMBLY-LEN] [QUALITY-SCORE-THRESHOLD]\n", argv[0]);
     return (1);
   }
  
  min_asm_len    = atoi(argv[3]);
  max_asm_len    = atoi(argv[4]);
  min_qual_score = atoi(argv[5]);

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
       assembly (ri_left[i], ri_right[i], qs_lps, qs_rps, &ai[i]);
     }
    free (qs_lps);
    free (qs_rps);
  }

  fd = fopen ("output-par.fastq","w");
  for (i = 0; i < cnt_reads_left; ++ i)
   {
     fprintf (fd, "%s\n", ri_left[i]->header);
     if ((strlen (ai[i].data) >= min_asm_len) && (strlen (ai[i].data) <= max_asm_len))
      {
        fprintf (fd, "%s\n", ai[i].data);
        fprintf (fd, "+\n");
        fprintf (fd, "%s\n", ai[i].quality_score);
      }
     else
      {
        trim (ri_left[i], min_qual_score, read_size);
        trim_cpl (ri_right[i], min_qual_score, read_size);
        fprintf (fd, "%s%s\n+\n%s%s\n", ri_left[i]->data, ri_right[i]->data, ri_left[i]->quality_score, ri_right[i]->quality_score);
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
  
  return (0);
}
