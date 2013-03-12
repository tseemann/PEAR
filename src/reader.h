#ifndef READER_H
#define READER_H

/** @file reader.h
    @brief Header file for memory pool and reader

    Contains structures for storing the pair-end reads and for managing memory
*/

/** @brief Data structure for storing parsed reads
    
    Consists of three elements. The \a header is a pointer to the header of the read, \a data
    is a pointer to the sequence itself and \a qscore points to the quality score sequence of
    the particular read
*/
struct read_t
 {
   char * header;  /**< @brief Read header */
   char * data;    /**< @brief Read sequence */
   char * qscore;  /**< @brief Quality scores of sequence */
 };

/** @brief A block representing a memory window of the read files
    
    Consists of six elements. \a Reads is the memory space for the
    \a read_t pointers, \a rawdata is the raw data read from the file,
    \a rawdata_end points to the memory cell after the last read byte,
    \a rawdata_size is the number of bytes read and \a max_reads_count
    is the number of \b complete reads parsed after the parsing routine
    parses the block
*/
struct block_t
 {
   struct read_t ** reads;        /**< @brief Array of read_t structures */
   char * rawdata;                /**< @brief Raw data read from file */
   char * rawdata_end;            /**< @brief Pointer to the memory location after the last read byte */
   char * unread;                 /**< @brief Pointer to rawdata, at the start of the incomplete read */
   size_t rawdata_size;           /**< @brief Number of bytes read from file */
   unsigned int max_reads_count;  /**< @brief Number of complete reads parsed from rawdata */
 };


void init_fastq_reader (const char * file1, const char * file2, size_t memsize, struct block_t * fwd, struct block_t * rev);
int get_next_reads (struct block_t * fwd_block, struct block_t * rev_block);
void destroy_reader (void);
void init_fastq_reader_double_buffer (const char * file1, const char * file2, size_t memsize, struct block_t * pri_fwd, struct block_t * pri_rev, 
struct block_t * sec_fwd, struct block_t * sec_rev);
int db_get_next_reads (struct block_t * fwd_block, struct block_t * rev_block, struct block_t * old_fwd_block, struct block_t * old_rev_block);
int db_read_fastq_block (struct block_t * block, FILE * fp, struct block_t * old_block);
void rewind_files (void);
#endif
