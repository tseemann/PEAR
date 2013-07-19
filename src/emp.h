#ifndef EMP_H
#define EMP_H

/** @file emp.h
    @brief Header file for empirical frequencies computation

    Header file containing a data structure for storing empirical frequency related data
*/

/** @brief Structure for storing empirical frequency related data

    Stores frequency of each base \a freqa, \a freqc, \a freqg and \a freqt. The total of all
    bases is stored in \a total. \a pa, \a pc, \a pg and \a pt represent the ratio of each
    base to the total number of bases. \a q is the sum of squares of \a pa, \a pc, \a pg and \a pt.
*/
struct emp_freq
 {
   unsigned long freqa;   /**< @brief Frequency of A */
   unsigned long freqc;   /**< @brief Frequency of C */
   unsigned long freqg;   /**< @brief Frequency of G */
   unsigned long freqt;   /**< @brief Frequency of T */
   unsigned long freqn;   /**< @brief Frequency of N */

   double total; /**< @brief \a freqa + \a freqc + \a freqg + \a freqt */
   
   double pa;   /**< @brief \a freqa / total */
   double pc;   /**< @brief \a freqc / total */
   double pg;   /**< @brief \a freqg / total */
   double pt;   /**< @brief \a freqt / total */

   double q;    /**< @brief \a freqa * \a freqa + \a freqc * \a freqc + \a freqg * \a freqg + \a freqt * \a freqt */
 };

#endif
