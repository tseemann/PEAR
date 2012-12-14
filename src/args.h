#ifndef ARGS_H
#define ARGS_H
#include <stdio.h>
#include <getopt.h>

/** @brief User arguments
  * 
  * Contains all values of optional and mandatory arguments entered
  * by the user when executing PAIR. It also contains default values
  * for the parameters that were not specified.
  */
struct user_args 
 {
   char       * fastq_left;   /**< @brief Forward pairend FASTQ filename */
   char       * fastq_right;  /**< @brief Reverse pairend FASTQ filename */
   int          min_asm_len;  /**< @brief Minimum assembly length threshold */
   int          max_asm_len;  /**< @brief Maximum assembly length threshold */
   int          qual_thres;   /**< @brief Quality score threshold */
   int          score_method; /**< @brief Scoring method to use */
   int          min_overlap;  /**< @brief Minimum overlap threshold */
   int          phred_base;   /**< @brief Base Phred quality score, i.e. 33 or 64 */
   double       max_uncalled; /**< @brief Maximum proportion of uncalled bases (N) */
   int          emp_freqs;    /**< @brief Flag whether to compute/use empirical base frequencies */
   double       p_value;      /**< @brief P-value to use */
   double       geom_mean;    /**< @brief Geometric mean */
   char       * outfile;      /**< @brief Output filename to use */
   int          min_trim_len  /**< @brief Minimum trim length */
 };

void usage (void);
int decode_switches (int argc, char * argv[], struct user_args * sw);
#endif
