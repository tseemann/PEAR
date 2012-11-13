#ifndef FASTQ_H
#define FASTQ_H

//#define         BUFFER_SIZE             1048576
#define         REALLOC_BLOCK_SIZE      20971520

#define         BUFFER_SIZE             8196
//#define         REALLOC_BLOCK_SIZE      1024

struct reads_info
 {
   char               * header;
   char               * data;
   char               * quality_score;
 };

struct asm_info
 {
   char               * data;
   char               * quality_score;
 };

struct reads_info ** read_fastq (const char * file, int * cnt_reads);

#endif
