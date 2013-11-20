#ifndef ASYNC_H
#define ASYNC_H

#include "args.h"
#include "reader.h"
#include "emp.h"

/** @file async.h
    @brief Thread specific data structures

    Header file containing global and local data structures used by threads
*/


/** @brief Wrapper structure for representing reads of processed files

    Consists of 5 elements. \a fwd points to the forward reads block, \a rev
    points to the reverse reads block, \a reads is the number of reads in
    each of the two blocks, \a processed is the number of reads that are
    already processed (assembled, not assembled or discarded) by the threads
    doing the assembly, and \a threads is the number of threads still operating
    with reads from this \a memBlockInfo structure.
*/
typedef struct 
 {
   memBlock * fwd;    /**< @brief Pointer to the forwards reads block */
   memBlock * rev;    /**< @brief Pointer to the reverse reads block */
   unsigned int reads;      /**< @brief Number of reads in each block */
   unsigned int processed;  /**< @brief Number of so far processed reads */
   unsigned int threads;    /**< @brief Number of threads still operating with reads from this structure */
 } memBlockInfo;

/** @brief A global data structure used to control program flow

    In order to speed the implementation the idea of double-buffering, widely used in graphics
    and game programming optimization, is used. The idea is that, since reading from the disk
    is much much slower than accessing RAM, we try to parallelize this process. That means, we
    read two memBlockInfo structures but the threads operate only on one of them, i.e. they
    assemble the reads in only one (\a xblock) datastructure. Once the reads of \a xblock are
    assembled, we flip \a xblock with \a yblock and the threads continue to do the assembly in
    \a xblock (now pointing to the contents of \a yblock). One of the threads though, instead
    of processing the reads, fills \a yblock with the next file contents stored in the disk. Once
    it finishes, it continues to assemble reads with the other threads. The process is repeated
    until all reads are read and processed, therefore not losing any time (except the initial
    lattency) in disk reading.
*/
struct thread_global_t
 {
   memBlockInfo * xblock; /**< @brief Pointer to the first \a memBlockInfo buffer (main processing buffer) */
   memBlockInfo * yblock; /**< @brief Pointer to the second \a memBlockInfo buffer (double-buffer) */
   int io_thread;               /**< @brief Trigger to denote which thread should read the next content of files in \a yblock, otherwise -1 */ 
   int finish;                  /**< @brief Trigger to denote that we read everything from the files */
   FILE * fd[4];                /**< @brief Pointers to the 4 output files */
 };

/** @brief Local data structure for every thread

    Local data structure passed to every thread. It contains the pointer to the \a block
    of reads the current thread is working on, \a sw which is a pointer to the command-line
    arguments, the \a id of the thread, a range (\a start, \a end) of reads that are to
    be processed by this thread, the \a match and mismatch score (\a match_score, \a mismatch_score)
    and the empirical frequencies (\a ef)
*/
struct thread_local_t
 {
   memBlockInfo * block; /**< @brief Block on which the current thread operates */
   struct user_args * sw;      /**< @brief Parsed command-line arguments */
   int id;                     /**< @brief Thread id */
   unsigned int start;         /**< @brief Number of the first read to be processed by current thread */
   unsigned int end;           /**< @brief Number of the last read to be processed by current thread */
//   int match_score;            /**< @brief Match score */
//   int mismatch_score;         /**< @brief Mismatch score */
   struct emp_freq * ef;       /**< @brief Empirical frequencies */
 };

#endif
