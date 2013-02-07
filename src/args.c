#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "args.h"

/** @file args.c
    @brief Command-line arguments parsing

    A data-structure and code used for handling the command-line arguments
    parsing phase. 
*/


/** @brief Command-line short and long options
    
    Command-line short and the corresponding long options structure required
    by the \a getopt_long function
*/
static struct option long_options[] =
 {
   { "phred-base",        required_argument, NULL, 'b' },
   { "empirical-freqs",   no_argument,       NULL, 'e' },
   { "left-fastq",        required_argument, NULL, 'f' },
   { "geometric-mean",    required_argument, NULL, 'g' },
   { "help",              no_argument,       NULL, 'h' },
   { "max-asm-length",    required_argument, NULL, 'm' },
   { "min-asm-length",    required_argument, NULL, 'n' },
   { "output",            required_argument, NULL, 'o' },
   { "p-value",           required_argument, NULL, 'p' },
   { "quality-theshold",  required_argument, NULL, 'q' },
   { "right-fastq",       required_argument, NULL, 'r' },
   { "score-method",      required_argument, NULL, 's' },
   { "min-trim-length",   required_argument, NULL, 't' },
   { "max-uncalled-base", required_argument, NULL, 'u' }, 
   { "min-overlap",       required_argument, NULL, 'v' },
   { NULL,                0,                 NULL, 0   }
 };

/** @brief Usage help screen
     
    A help-screen for all command-line options available in PEAR
*/
void usage (void)
{
  fprintf (stdout, " ____  _____    _    ____ \n"); 
  fprintf (stdout, "|  _ \\| ____|  / \\  |  _ \\\n");
  fprintf (stdout, "| |_) |  _|   / _ \\ | |_) |\n");
  fprintf (stdout, "|  __/| |___ / ___ \\|  _ <\n");
  fprintf (stdout, "|_|   |_____/_/   \\_\\_| \\_\\\n");
  fprintf (stdout, "\n.oOo. Paired-End AssembleR .oOo.\n");
  fprintf (stdout, "PEAR v0.8 by Tomas Flouri and Jiajie Zhang\n");
  fprintf (stdout, "Free for academic use, for commercial use or bug report, please contact:\n");
  fprintf (stdout, "flouris@gmail.com and bestzhangjiajie@gmail.com\n");
  fprintf (stdout, "\n\n"); 
  fprintf (stdout, "Usage: pear <options>\n");
  fprintf (stdout, "Standard (mandatory):\n");
  fprintf (stdout, "  -f, --forward-fastq         <str>     Forward paired-end FASTQ file.\n");
  fprintf (stdout, "  -r, --reverse-fastq         <str>     Reverse paired-end FASTQ file.\n");
  fprintf (stdout, "  -o, --output                <str>     Output filename.\n");
  fprintf (stdout, "Optional:\n");
  fprintf (stdout, "  -p, --p-value               <float>   Use a p-value from the set { 1.0, 0.05, 0.01, 0.001, 0.0001 }. If\n"
                   "                                        the p-value of the assembled reads exceeds the specified p-value, the\n"
                   "                                        reads will be output unassembled.\n"
                   "                                        Set 1.0 to disable the test.(default: 0.01)\n");
  fprintf (stdout, "  -v, --min-overlap           <int>     Minimum overlap (default: 10)\n"
                   "                                        If the statistical test is used, the min-overlap can in principal be set to 1,\n"
                   "                                        but in practice setting min-overlap to a proper value will further reduce\n"
                   "                                        false-positive assemlies, since the data will not be perfect.\n"); 	
  fprintf (stdout, "  -m, --max-assembly-length   <int>     Maximum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                        The assembled sequence can be arbitrary long if set\n"
                   "                                        to 0. (default: 0)\n");
  fprintf (stdout, "  -n, --min-assembly-length   <int>     Minimum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                        To disable it set to 0. (default: 50)\n");
  fprintf (stdout, "  -t, --min-trim-length       <int>     Minimum length of reads after trimming the low quality part. If two consecutive\n"
                   "                                        bases quality scores < q, the rest of the reads will be trimed (default: 1)\n");
  fprintf (stdout, "  -q, --quality-threshold     <int>     Quality score threshold used for trimming the low quality\n"
                   "                                        part of the reads (default: 0)\n");
  fprintf (stdout, "  -u, --max-uncalled-base     <float>   Maximal proportion of uncalled bases. A number between 0 and 1.\n"
                   "                                        Set to 0 to discard all reads that contain uncalled bases, or\n"
                   "                                        1 to process all sequences independent on the number of uncalled\n"
                   "                                        bases. (default: 1)\n");
  fprintf (stdout, "  -g, --test-method           <int>     Statistical test method: (default: 1)\n"
                   "                                        1: Test using the highest OES, given the minimum overlap allowed.\n"
                   "                                           Note due to the discret nature of the test, it usually gives a lower p-value\n" 
                   "                                           for assembled sequences than the specified one. For example, set p-value = 0.05,\n" 
                   "                                           using this test, the assembled reads might have an actual p-value of 0.02.\n"
                   "                                        2: Using the acceptance probability. \n"
                   "                                           Test method 2 calculate the same probability as test1, but assumes the minimal overlap\n"
                   "                                           is the observed overlap who has the highest OES, instead of the minimum allowed overlap \n"
                   "                                           predefined as input parameter (-v). Therefore, it is not a valid statistical test, the \n"
                   "                                           'p-value' is really the maximal probability we accpet the assembly. However, we found \n"
                   "                                           in practice, when the actual overlap sizes are small, test 2 can produce more correctly \n"
                   "                                           assembled sequences with only slightly higher false-positive rates.\n");
  fprintf (stdout, "  -e, --empirical-freqs                 Disable empirical base frequencies. (default: use empirical base frequencies)\n");
  fprintf (stdout, "  -s, --score-method          <int>     Scoring method\n"
                   "                                        1: OES with +1 for match and -1 for mismatch.\n"
                   "                                        2: Scaled score, use the probobality of bases been correct or wrong to scale \n"
                   "                                           the score in method 3 (both tests are invalid in use this method).\n"
                   "                                        3: +1 for a match, -1 for a mismatch, ignoring the quality scores.\n"
                   "                                        (both tests are invalid if use this method 2 or 3)(default: 1)\n");				   		
  fprintf (stdout, "  -b, --phred-base            <int>     Base Phred quality score (default: 33)\n");
  fprintf (stdout, "  -h, --help                            This help screen.\n\n");
}

/** @brief Command-line arguments parser
    
    A parser for the command-line options of PEAR. A minimum of the two pair-end
    reads and output filename must be provided.

    @param argc
      The number of command-line parameters given

    @param argv
      The array of command-line parameters

    @param sw
      The structure where the user-defined switches will be stored in
*/
int decode_switches (int argc, char * argv[], struct user_args * sw)
{
  int    opt;
  int    oi;
  char * ep;

  /* initialization */
  sw->fastq_left    = NULL;
  sw->fastq_right   = NULL;
  sw->outfile       = NULL;
  sw->min_asm_len   = 50;
  sw->max_asm_len   = 999999;
  sw->qual_thres    = 0;
  sw->phred_base    = 33;
  sw->max_uncalled  = 1.0;
  sw->min_overlap   = 10;
  sw->emp_freqs     = 1;    /* use empirical base frequencies as default */
  sw->p_value       = 0.01;
  sw->geom_mean     = 0;
  sw->min_trim_len  = 1;
  sw->score_method  = 1;
  sw->test          = 1;

  while ((opt = getopt_long(argc, argv, "b:ef:g:hm:n:o:p:q:r:s:t:u:v:", long_options, &oi)) != -1)
   {
     switch (opt)
      {
        case 'b':
          sw->phred_base  = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0')
           {
             printf ("Please input a correct base Phred quality score.\n");
             return (0);
           }
          break;

        case 'e':
          sw->emp_freqs = 0;
          break;

        case 'f':
          sw->fastq_left  = optarg;
          break;

        case 'g':
          if (!strcmp (optarg, "1") || !strcmp (optarg, "2"))
           {
             sw->test = (int) strtol (optarg, &ep, 10);
           }
          else
           {
             printf ("Invalid testing method.\n");
             return (0);
           }
          break;

        case 'h':
          return (0);

        case 'm':
          sw->max_asm_len = (int)strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid max-assembly-length.\n");
             return (0);
           }
          break;

        case 'n':
          sw->min_asm_len = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid min-assembly-length.\n");
             return (0);
           }
          break;

        case 'o':
          sw->outfile = optarg;
          break;

        case 'p':		
          if (!strcmp (optarg, "1.0") || !strcmp (optarg, "0.05") || !strcmp (optarg, "0.01") || !strcmp (optarg, "0.001") || !strcmp (optarg, "0.0001") )
           {
             sw->p_value = strtod (optarg, &ep);
           }
          else
           {
             printf ("Invalid p-value or minimal OES.\n");
             return (0);
           }
          
          break;

        case 'q':
          sw->qual_thres = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid quality-threshold.\n");
             return (0);
           }
          break;

        case 'r':
          sw->fastq_right = optarg;
          break;

        case 's':
          if (!strcmp (optarg, "1") || !strcmp (optarg, "2") || !strcmp (optarg, "3") )
           {
             sw->score_method = (int) strtol (optarg, &ep, 10);
           }
          else
           {
             printf ("Invalid score-method.\n");
             return (0);
           }
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid score-method.\n");
             return (0);
           }
          break;

        case 't':
          sw->min_trim_len = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid min-trim-length.\n");
             return (0);
           }
          break;

        case 'u':
          sw->max_uncalled = strtod (optarg, &ep);     /* TODO: check this line */
          
          if (ep == optarg || *ep != '\0' || sw->max_uncalled < 0 || sw->max_uncalled > 1)
           {
             printf ("Invalid max-uncalled-base value, must be between 0 and 1.\n");
             return (0);
           }
          break;

        case 'v':
          sw->min_overlap = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Invalid min-overlap length.\n");
             return (0);
           }
          break;
      }
   }
  return (sw->fastq_left && sw->fastq_right && sw->outfile);
}
