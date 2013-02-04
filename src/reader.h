#ifndef READER_H
#define READER_H


struct read_t
 {
   char * header;
   char * data;
   char * qscore;
 };

struct block_t
 {
   struct read_t ** reads;
   char * rawdata;
   char * rawdata_end;
   char * unread;
   int rawdata_size;
   int max_reads_count;
 };


void init_fastq_reader (const char * file1, const char * file2, size_t memsize, struct block_t * fwd, struct block_t * rev);
int get_next_reads (struct block_t * fwd_block, struct block_t * rev_block);
void destroy_reader (void);
void init_fastq_reader_double_buffer (const char * file1, const char * file2, size_t memsize, struct block_t * pri_fwd, struct block_t * pri_rev, 
struct block_t * sec_fwd, struct block_t * sec_rev);
int db_get_next_reads (struct block_t * fwd_block, struct block_t * rev_block, struct block_t * old_fwd_block, struct block_t * old_rev_block);
int db_read_fastq_block (struct block_t * block, FILE * fp, struct block_t * old_block);
#endif
