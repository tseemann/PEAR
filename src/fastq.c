#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fastq.h"

int 
line_count(const char * s)
{
  int           lc = 0;
  int           not_empty = 0;

  if (*s) not_empty = 1;

  while (*s)
   {
     if (*s == '\n') ++ lc;
     ++s;
   }

  if (!lc && not_empty) lc = 1;

  return (lc);
}

struct reads_info ** 
read_fastq (const char * file, int * cnt_reads)
{
  FILE                * fd;
  char                * data = NULL;
  struct reads_info  ** ri;
  int                   cnt_lines;
  char                * data_seg_start;
  char                * data_seg_end;
  int                   i;
  long                  file_size;

  if (!(fd = fopen (file, "r")))
   {
     fprintf (stderr, "Error: cannot open FASTQ file %s\n", file);
     return (0);
   }

  fseek (fd, 0, SEEK_END);
  file_size = ftell (fd);
  rewind (fd);

  data = (char *) malloc ((file_size + 1) * sizeof(char));
  if (file_size != fread (data, sizeof(char), file_size, fd))
   {
     fprintf (stderr, "Error while reading %s\n", file);
     exit (1);
   }

  data[file_size] = '\0';

  if (!data || file_size < 10)
   {
     fprintf (stderr,"Error: file %s is empty\n", file);
     fclose (fd);
     free (data);
     return (0);
   }
  fclose (fd);

  cnt_lines = line_count (data);

  /* check the validity of the FASTA file */
  if (cnt_lines % 4)
   {
     fprintf (stderr, "Error: invalid FASTQ file\n");
     free (data);
     return (0);
   }

  /* allocate memoryspace for the structure */
  ri = (struct reads_info **) calloc (cnt_lines, sizeof(struct reads_info *));

  /* in case it is of FASTQ format, locate the description header */
  data_seg_start = data;   i = 0;
  while (*data_seg_start)
   {
     if (*data_seg_start == '@')
      {
        ri[i] = (struct reads_info *) malloc (sizeof(struct reads_info));

        /* read header */
        data_seg_end = strchr (data_seg_start, '\n');
        if (data_seg_end)
         {
           *data_seg_end = '\0';
           if (*(data_seg_end - 1) == '\r') *(data_seg_end - 1) = '\0';
           ri[i]->header  = strdup (data_seg_start);
           data_seg_start = data_seg_end + 1;
         }

        /* read sequence */
        if (*data_seg_start)
         {
           data_seg_end = strchr (data_seg_start, '\n');
           if (data_seg_end)
            {
              *data_seg_end = '\0';
              if (*(data_seg_end - 1) == '\r') *(data_seg_end - 1) = '\0';
              ri[i]->data    = strdup (data_seg_start);
              data_seg_start = data_seg_end + 1;
            }
         }

        /* ignore the third line */
        if (*data_seg_start)
         {
           data_seg_end   = strchr (data_seg_start, '\n');
           data_seg_start = data_seg_end + 1;
         }

        /* read quality scores */

        if (*data_seg_start)
         {
           data_seg_end = strchr (data_seg_start, '\n');
           if (data_seg_end)
            {
              *data_seg_end = '\0';
              if (*(data_seg_end - 1) == '\r') *(data_seg_end - 1) = '\0';
              ri[i]->quality_score = strdup (data_seg_start);
              data_seg_start = data_seg_end + 1;
            }
         }
      }
     ++ i;
   }
  
  *cnt_reads = i;

  free (data);

  return (ri);
}
