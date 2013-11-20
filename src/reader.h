#ifndef READER_H
#define READER_H

#define         PEAR_FORWARD_LARGER              1
#define         PEAR_REVERSE_LARGER              2

#define         PEAR_PARSE_PHASE_HEADER          0
#define         PEAR_PARSE_PHASE_SEQUENCE        1
#define         PEAR_PARSE_PHASE_PLUS_SIGN       2
#define         PEAR_PARSE_PHASE_QUALITY_VALS    3

#ifdef __APPLE__
#define INLINE
#else
#define INLINE inline
#endif


/** @file reader.h
    @brief Header file for memory pool and reader

    Contains structures for storing the pair-end reads and for managing memory
*/

/** @brief Data structure for storing parsed reads
    
    Consists of three elements. The \a header is a pointer to the header of the read, \a data
    is a pointer to the sequence itself and \a qscore points to the quality score sequence of
    the particular read
*/
typedef struct
 {
   char * header;  /**< @brief Read header */
   char * data;    /**< @brief Read sequence */
   char * qscore;  /**< @brief Quality scores of sequence */
 } fastqRead;

/** @brief A block representing a memory window of the read files
    
    Consists of six elements. \a reads is the memory space for the
    \a read_t pointers, \a rawdata is the raw data read from the file,
    \a rawdata_end points to the memory cell after the last read byte,
    \a rawdata_size is the number of bytes read and \a max_reads_count
    is the number of \b complete reads parsed after the parsing routine
    parses the block
*/
typedef struct
 {
   fastqRead ** reads;        /**< @brief Array of read_t structures */
   char * rawdata;                /**< @brief Raw data read from file */
   char * rawdata_end;            /**< @brief Pointer to the memory location after the last read byte */
   char * unread;                 /**< @brief Pointer to rawdata, at the start of the incomplete read */
   int nExtraReads;               /**< @brief Number of complete reads in the unread part */
   size_t rawdata_size;           /**< @brief Number of bytes read from file */
   unsigned int max_reads_count;  /**< @brief Number of complete reads parsed from rawdata */
 } memBlock;


void init_fastq_reader (const char * file1, const char * file2, size_t memsize, memBlock * fwd, memBlock * rev);
int get_next_reads (memBlock * fwd_block, memBlock * rev_block);
void destroy_reader (void);
void init_fastq_reader_double_buffer (const char * file1, const char * file2, size_t memsize, memBlock * pri_fwd, memBlock * pri_rev, 
memBlock * sec_fwd, memBlock * sec_rev);
unsigned int db_get_next_reads (memBlock * fwd_block, memBlock * rev_block, memBlock * old_fwd_block, memBlock * old_rev_block, int *sanity);
int db_read_fastq_block (memBlock * block, FILE * fp, memBlock * old_block);
int read_fastq_block (memBlock * block, FILE * fp);
void rewind_files (void);
#endif
