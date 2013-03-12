#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "reader.h"

#define READ_SIZE 400

static FILE * fp1;
static FILE * fp2;
static char * mempool;

/* READ_SIZE 10 and 3108 was the fastest till now */
//struct block_t fwd_block;
//struct block_t rev_block;

unsigned int rcount = 0;

void print_number (size_t x)
{
  unsigned int digits;
  unsigned int triplets = 0;
  unsigned int i, j, k;
  size_t num;
  unsigned int y[4];

  digits = (x > 0) ? (int) log10 ((double)x) + 1 : 1;
  triplets = ceil (digits / 3.0);
  k = 0;

  for (i = 0; i < triplets; ++ i)
   {
     num = x;
     for (j = 0; j < k; ++ j)
      {
        num = num / 1000;
      }
     y[triplets - i - 1] = num % 1000;
     ++k;
   }
  printf ("%d", y[0]);
  for (k = 1; k < triplets; ++k)
   {
    printf (",%03d", y[k]);
   }
}


void comp_mem (size_t memsize, size_t * reads_count, size_t * rawdata_size)
{
  size_t x,y,z;
  size_t u;

  x = memsize;
  y = READ_SIZE;
  z = sizeof (struct read_t);

  u = (size_t) ((double)x / (sizeof(struct read_t *) + y + z));

  *reads_count = u;
  *rawdata_size = u * y;
}

void rewind_files (void)
{
  rewind (fp1);
  rewind (fp2);
}

void init_fastq_reader_double_buffer (const char * file1, const char * file2, size_t memsize, struct block_t * pri_fwd, struct block_t * pri_rev, 
struct block_t * sec_fwd, struct block_t * sec_rev)
{
  size_t reads_count;
  size_t rawdata_size;
  size_t used_mem;
  size_t unused_mem;
  void * mem_start;
  void * dbmem;
  int i;

  #ifdef PRINT_MEM
  printf ("Mempooling...\n");
  #endif

  comp_mem (memsize / 4, &reads_count, &rawdata_size);

  #ifdef PRINT_MEM
  printf ("tm: %f\n", (double)memsize);
  printf ("Total memory: ");
  print_number(memsize);
  printf ("\n");
  printf ("Total memory: %ld\n", memsize);
  #endif

  #ifdef PRINT_MEM
  printf ("Number of reads: ");
  print_number(reads_count);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of reads structure: ");
  print_number(reads_count * sizeof(struct read_t) + reads_count * sizeof(struct read_t *)); 
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of block: ");
  print_number(rawdata_size);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of used memory: ");
  #endif
  used_mem = 2 * (rawdata_size + reads_count * sizeof(struct read_t) + reads_count * sizeof(struct read_t *));
  #ifdef PRINT_MEM
  print_number(used_mem);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of unused memory: ");
  #endif
  unused_mem = memsize - used_mem;
  #ifdef PRINT_MEM
  print_number(unused_mem);
  printf ("\n");
  #endif

  fp1 = fopen (file1, "r");
  fp2 = fopen (file2, "r");

  if (!fp1 || !fp2)
   {
     fprintf (stderr, "Failed to open files..\n");
     abort();
   }

  pri_fwd->max_reads_count   = pri_rev->max_reads_count   = reads_count;
  pri_fwd->rawdata_size      = pri_rev->rawdata_size      = rawdata_size;

  sec_fwd->max_reads_count   = sec_rev->max_reads_count   = reads_count;
  sec_fwd->rawdata_size      = sec_rev->rawdata_size      = rawdata_size;

  /* allocate memory */
  //mempool = calloc (1,memsize);
  mempool = malloc (memsize);
  if (!mempool)
   {
     fprintf (stderr, "Failed to allocate memory...\n");
     abort ();
   }
  #if defined(__LP64__) || defined(_LP64)
  printf ("Allocating %lu bytes\n", memsize);
  #else
  printf ("Allocating %u bytes\n", memsize);
  #endif



  /* reserve area from mempool for the forwards reads */
  pri_fwd->reads = (struct read_t **) mempool;
  mem_start        = mempool + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     pri_fwd->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  pri_fwd->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  pri_fwd->rawdata_end = pri_fwd->rawdata + pri_fwd->rawdata_size;

  /* reserve area from mempool for the backward reads */
  mem_start  = pri_fwd->rawdata_end;
  pri_rev->reads = (struct read_t **) mem_start;
  mem_start  = mem_start + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     pri_rev->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  pri_rev->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  pri_rev->rawdata_end = pri_rev->rawdata + pri_rev->rawdata_size;

  pri_fwd->unread = pri_rev->unread = NULL;
  

  /* double buffering */

  dbmem = mempool + memsize / 2;
  /* reserve area from mempool for the forwards reads */
  sec_fwd->reads = (struct read_t **) dbmem;
  mem_start        = dbmem + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     sec_fwd->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  sec_fwd->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  sec_fwd->rawdata_end = sec_fwd->rawdata + sec_fwd->rawdata_size;

  /* reserve area from mempool for the backward reads */
  mem_start  = sec_fwd->rawdata_end;
  sec_rev->reads = (struct read_t **) mem_start;
  mem_start  = mem_start + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     sec_rev->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  sec_rev->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  sec_rev->rawdata_end = sec_rev->rawdata + sec_rev->rawdata_size;

  sec_fwd->unread = sec_rev->unread = NULL;
  
  #ifdef PRINT_MEM
  printf ("\n");
  #endif
}


void 
init_fastq_reader (const char * file1, const char * file2, size_t memsize, struct block_t * fwd, struct block_t * rev)
{
  size_t reads_count;
  size_t rawdata_size;
  size_t used_mem;
  size_t unused_mem;
  void * mem_start;
  int i;

  #ifdef PRINT_MEM
  printf ("Mempooling...\n");
  #endif

  comp_mem (memsize / 2, &reads_count, &rawdata_size);

  #ifdef PRINT_MEM
  printf ("Total memory: ");
  print_number(memsize);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Number of reads: ");
  print_number(reads_count);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of reads structure: ");
  print_number(reads_count * sizeof(struct read_t) + reads_count * sizeof(struct read_t *)); 
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of block: ");
  print_number(rawdata_size);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of used memory: ");
  #endif
  used_mem = 2 * (rawdata_size + reads_count * sizeof(struct read_t) + reads_count * sizeof(struct read_t *));
  #ifdef PRINT_MEM
  print_number(used_mem);
  printf ("\n");
  #endif

  #ifdef PRINT_MEM
  printf ("Size of unused memory: ");
  #endif
  unused_mem = memsize - used_mem;
  #ifdef PRINT_MEM
  print_number(unused_mem);
  printf ("\n");
  #endif

  fp1 = fopen (file1, "r");
  fp2 = fopen (file2, "r");

  fwd->max_reads_count   = rev->max_reads_count   = reads_count;
  fwd->rawdata_size      = rev->rawdata_size      = rawdata_size;

  /* allocate memory */
  //mempool = calloc (1,memsize);
  mempool = malloc (memsize);


  /* reserve area from mempool for the forwards reads */
  fwd->reads = (struct read_t **) mempool;
  mem_start        = mempool + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     fwd->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  fwd->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  fwd->rawdata_end = fwd->rawdata + fwd->rawdata_size;

  /* reserve area from mempool for the backward reads */
  mem_start  = fwd->rawdata_end;
  rev->reads = (struct read_t **) mem_start;
  mem_start  = mem_start + reads_count * sizeof (struct read_t *);
  for (i = 0; i < reads_count; ++ i)
   {
     rev->reads[i] = (struct read_t *) (mem_start + i * sizeof (struct read_t));
   }
  rev->rawdata     = (char *)(mem_start + reads_count * sizeof (struct read_t));
  rev->rawdata_end = rev->rawdata + rev->rawdata_size;

  fwd->unread = rev->unread = NULL;



  #ifdef PRINT_MEM
  printf ("\n");
  #endif

}

void
destroy_reader (void)
{
  free (mempool);
  fclose (fp1);
  fclose (fp2);
}

/* Parse a block to a reads struct and return a pointer to the last unprocessed
   read, if such one exists */
inline int
parse_block (struct block_t * block)
{
  int phase;
  char * ptr;
  int elms;
  char * offset;
  char * ignore = NULL;

  phase = elms = 0;

  block->unread = ptr = block->rawdata;
  for (offset = block->rawdata; offset != block->rawdata_end && *offset; ++ offset)
   {
     if (*offset == '\n')
      {
        if (phase == 2 && offset - 1 == ptr && (*(offset - 1) != '+'))
         {
           fprintf (stderr, "Entry is missing\n");
           abort();
         }
        else if (phase != 2 && offset - ptr <=2)
         {
           fprintf (stderr, "Entry is missing\n");
           abort();
         }
        switch (phase)
         {
           case 0:
             block->reads[elms]->header = ptr;
             ++phase;
             break;
           case 1:
             block->reads[elms]->data   = ptr;
             ++phase;
             ignore = offset;
             break;
           case 2:
             ++phase;
             break;
           case 3:
             block->reads[elms]->qscore = ptr;
             /* clear up newlines */
             if (*(block->reads[elms]->data - 1) == '\r' || (*(block->reads[elms]->data - 1) == '\n'))
              {
                *(block->reads[elms]->data - 1) = 0;
              }
             else
              {
                *(block->reads[elms]->data) = 0;
              }
             if (*(ptr - 1) == '\r' || *(ptr - 1) == '\n')
              {
                *(ptr - 1) = 0;
              }
             else
              {
                *ptr = 0;
              }
             if (*(offset - 1) == '\r' || *(offset - 1) == '\n')
              {
                *(offset - 1) = 0;
              }
             else
              {
                *offset = 0;
              }
             phase = 0;
             *ignore = 0;
             ++elms;
             block->unread = offset + 1;
             if (elms == block->max_reads_count)
              {
                return elms;
              }
             break;
         }
        ptr = offset + 1;
      }
   }

  return elms;
}

int db_read_fastq_block (struct block_t * block, FILE * fp, struct block_t * old_block)
{
  int nBytes;
  size_t remainder;

  /* TODO: check if something from the previous block exists */
  /* check also if we do not have a large block that couldnt be read in one read */

  /* Check if we do not have a large block that could not be read in one read */
  if (old_block->unread == old_block->rawdata)
   {
     fprintf (stderr, "Error, too large read? Allocate more mem...");
     abort ();
   }

  remainder    = 0;
  if (old_block->unread && old_block->unread != block->rawdata_end && (*(old_block->unread)) )
   {
     remainder = (size_t) (old_block->rawdata_end - old_block->unread);
     memcpy (block->rawdata, old_block->unread, remainder); 
   }
  
  /* TODO: Check this line again for correctness */
  nBytes = fread (block->rawdata + remainder, sizeof (char), block->rawdata_size - remainder, fp);
  if (!nBytes)
   {
//     fprintf (stderr, "Finished reading (0 bytes)?\n");
     return (0);
   }
  else if (nBytes < block->rawdata_size - remainder)
   {
     block->rawdata[nBytes + remainder] = '\0';
//     fprintf (stderr, "Finished reading (less bytes)?\n");
     return (0);
   }

  return (1);
}

int read_fastq_block (struct block_t * block, FILE * fp)
{
  int nBytes;
  size_t remainder;

  /* TODO: check if something from the previous block exists */
  /* check also if we do not have a large block that couldnt be read in one read */

  /* Check if we do not have a large block that could not be read in one read */
  if (block->unread == block->rawdata)
   {
     fprintf (stderr, "Error, too large read? Allocate more mem...");
     abort ();
   }

  remainder    = 0;
  if (block->unread && block->unread != block->rawdata_end && (*(block->unread)) )
   {
     remainder = (size_t) (block->rawdata_end - block->unread);
     memmove (block->rawdata, block->unread, remainder); 
   }
  
  /* TODO: Check this line again for correctness */
  nBytes = fread (block->rawdata + remainder, sizeof (char), block->rawdata_size - remainder, fp);
  if (!nBytes)
   {
//     fprintf (stderr, "Finished reading (0 bytes)?\n");
     return (0);
   }
  else if (nBytes < block->rawdata_size - remainder)
   {
     block->rawdata[nBytes + remainder] = '\0';
//     fprintf (stderr, "Finished reading (less bytes)?\n");
     return (0);
   }

  return (1);
}

static inline void
do_cpuid(uint32_t selector, uint32_t *data)
{
  asm("cpuid"
      : "=a" (data[0]),
      "=b" (data[1]),
      "=c" (data[2]),
      "=d" (data[3])
      : "a"(selector));
}

void print_reads (struct read_t ** fwd, struct read_t ** rev, int elms)
{
  int i;

  for (i = 0; i < elms; ++ i)
   {
      printf ("%s\n%s\n%s\n\n%s\n%s\n%s\n\n", fwd[i]->header, fwd[i]->data, fwd[i]->qscore, rev[i]->header, rev[i]->data, rev[i]->qscore);
      ++rcount;
   }
}

int get_next_reads (struct block_t * fwd_block, struct block_t * rev_block)
{
  int n1, n2;
  int eof1 = 1, eof2 = 1;

  if (eof1)  eof1 = read_fastq_block (fwd_block, fp1);
  if (eof2)  eof2 = read_fastq_block (rev_block, fp2);
  
  n1 = parse_block (fwd_block);
  n2 = parse_block (rev_block);

//  printf ("Read %d and %d reads\n", n1, n2);
  
  /* align reads if different count selected */
  if (n1 != n2)
   {
     if (!n1 || !n2)
      {
        fprintf (stderr, "Problem, number of reads does not match!\n");
        abort();
      }
     if (n1 > n2)
      {
        fwd_block->unread = fwd_block->reads[n2]->header;
        n1 = n2;
      }
     else
      {
        rev_block->unread = rev_block->reads[n1]->header;
        n2 = n1;
      }
   }
  
  return (n1);
}

int db_get_next_reads (struct block_t * fwd_block, struct block_t * rev_block, struct block_t * old_fwd_block, struct block_t * old_rev_block)
{
  int n1, n2;
  int eof1 = 1, eof2 = 1;

  if (eof1)  eof1 = db_read_fastq_block (fwd_block, fp1, old_fwd_block);
  if (eof2)  eof2 = db_read_fastq_block (rev_block, fp2, old_rev_block);
  
  n1 = parse_block (fwd_block);
  n2 = parse_block (rev_block);
  
//  printf ("Read %d and %d reads\n", n1, n2);
  /* align reads if different count selected */
  if (n1 != n2)
   {
     if (!n1 || !n2)
      {
        fprintf (stderr, "Problem, number of reads does not match!\n");
        abort();
      }
     if (n1 > n2)
      {
        fwd_block->unread = fwd_block->reads[n2]->header;
        n1 = n2;
      }
     else
      {
        rev_block->unread = rev_block->reads[n1]->header;
        n2 = n1;
      }
   }
  
  return (n1);
}

/*
int main (int argc, char * argv[])
{
  uint32_t data[4];
  do_cpuid (0, data);
  int elms;
  struct block_t fwd_block;
  struct block_t rev_block;

  printf("maxcpuid = 0x%x\n", data[0]);
      printf("vendor = %4.4s%4.4s%4.4s\n", (char *) &data[1], (char *)&data[3], (char *)&data[2]);
  if (argc != 4)
   {
     fprintf (stderr, "./%s [MEM-SIZE]\n", argv[0]);
     return (1);
   }
  
  init_fastq_reader (argv[1], argv[2], atoi(argv[3]), &fwd_block, &rev_block);

  while (1)
   {
     elms = get_next_reads(&fwd_block, &rev_block);
     if (!elms) break;
     print_reads (fwd_block.reads, rev_block.reads, elms);
   }

  destroy_reader ();

  printf ("Total reads: %u\n", rcount);

  return (0);
}
*/
