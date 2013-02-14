#ifndef ASYNC_H
#define ASYNC_H

#include "args.h"
#include "reader.h"

struct blockinfo_t
 {
   struct block_t * fwd;
   struct block_t * rev;
   unsigned int reads;
   unsigned int processed;
   unsigned int threads;
   struct list_t * next;
 };

struct thread_global_t
 {
   struct blockinfo_t * xblock;
   struct blockinfo_t * yblock;
   int io_thread;
   int finish;
   FILE * fd[4];
 };

struct thread_local_t
 {
   struct blockinfo_t * block;
   //struct block_t * fwd;
   //struct block_t * rev;

   struct user_args * sw;


   int id;

   unsigned int start;
   unsigned int end;
   int match_score;
   int mismatch_score;
 };

#endif
