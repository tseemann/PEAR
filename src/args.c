#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "pear.h"
#include "args.h"

extern void print_number (size_t x);

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
   { "test-method",       required_argument, NULL, 'g' },
   { "help",              no_argument,       NULL, 'h' },
   { "threads",           required_argument, NULL, 'j' },
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
   { "memory",            required_argument, NULL, 'y' },
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
  fprintf (stdout, "\n.oOo. Paired-End reAd mergeR .oOo.\n");
  fprintf (stdout, "%s v%s released on %s by Tomas Flouri and Jiajie Zhang\n", PROGRAM_NAME, PROGRAM_VERSION, VERSION_DATE);
  fprintf (stdout, "License: %s\n", LICENCE);
  fprintf (stdout, "Bug-reports and requests to: %s\n", CONTACT);
  fprintf (stdout, "\n\n"); 
  fprintf (stdout, "Usage: pear <options>\n");
  fprintf (stdout, "Standard (mandatory):\n");
  fprintf (stdout, "  -f, --forward-fastq         <str>     Forward paired-end FASTQ file.\n");
  fprintf (stdout, "  -r, --reverse-fastq         <str>     Reverse paired-end FASTQ file.\n");
  fprintf (stdout, "  -o, --output                <str>     Output filename.\n");
  fprintf (stdout, "Optional:\n");
  fprintf (stdout, "  -p, --p-value               <float>   Specify  a p-value for the statistical test. If the computed\n"
                   "                                        p-value of a possible assembly exceeds the specified p-value\n"
                   "                                        then  paired-end  read  will not be assembled. Valid options\n"
                   "                                        are: 0.0001, 0.001, 0.01, 0.05 and 1.0. Setting 1.0 disables\n"
                   "                                        the test. (default: 0.01)\n");
  fprintf (stdout, "  -v, --min-overlap           <int>     Specify the minimum overlap size. The minimum overlap may be\n"
                   "                                        set to 1 when the statistical test is used. However, further\n"
                   "                                        restricting  the  minimum overlap size to a proper value may\n"
                   "                                        reduce false-positive assembles. (default: 10)\n"); 	
  fprintf (stdout, "  -m, --max-assembly-length   <int>     Specify   the  maximum  possible  length  of  the  assembled\n"
                   "                                        sequences.  Setting this value to 0 disables the restriction\n"
                   "                                        and assembled sequences may be arbitrary long. (default: 0)\n");
  fprintf (stdout, "  -n, --min-assembly-length   <int>     Specify   the  minimum  possible  length  of  the  assembled\n"
                   "                                        sequences.  Setting this value to 0 disables the restriction\n"
                   "                                        and  assembled  sequences  may be arbitrary short. (default:\n"
                   "                                        50)\n");
  fprintf (stdout, "  -t, --min-trim-length       <int>     Specify  the  minimum length of reads after trimming the low\n"
                   "                                        quality part (see option -q). (default: 1)\n");
  fprintf (stdout, "  -q, --quality-threshold     <int>     Specify  the  quality  score  threshold for trimming the low\n"
                   "                                        quality  part  of  a  read.  If  the  quality  scores of two\n"
                   "                                        consecutive  bases  are  strictly  less  than  the specified\n"
                   "                                        threshold,  the  rest of the read will be trimmed. (default:\n"
                   "                                        0)\n");
  fprintf (stdout, "  -u, --max-uncalled-base     <float>   Specify  the maximal proportion of uncalled bases in a read.\n"
                   "                                        Setting this value to 0 will cause PEAR to discard all reads\n"
                   "                                        containing  uncalled  bases.  The other extreme setting is 1\n"
                   "                                        which  causes  PEAR  to process all reads independent on the\n"
                   "                                        number of uncalled bases. (default: 1)\n");
  fprintf (stdout, "  -g, --test-method           <int>     Specify  the  type  of  statistical  test.  Two  options are\n"
                   "                                        available. (default: 1)\n"
                   "                                        1: Given the minimum allowed overlap, test using the highest\n"
                   "                                        OES. Note that due to its discrete nature, this test usually\n"
                   "                                        yields  a lower p-value for the assembled read than the cut-\n"
                   "                                        off  (specified  by -p). For example, setting the cut-off to\n"
                   "                                        0.05  using  this  test,  the  assembled reads might have an\n"
                   "                                        actual p-value of 0.02.\n\n"
                   "                                        2. Use the acceptance probability (m.a.p). This test methods\n"
                   "                                        computes  the same probability as test method 1. However, it\n"
                   "                                        assumes  that  the  minimal  overlap is the observed overlap\n"
                   "                                        with  the  highest  OES, instead of the one specified by -v.\n"
                   "                                        Therefore,  this  is  not  a  valid statistical test and the\n"
                   "                                        'p-value'  is  in fact the maximal probability for accepting\n"
                   "                                        the assembly. Nevertheless, we observed in practice that for\n"
                   "                                        the case the actual overlap sizes are relatively small, test\n"
                   "                                        2  can  correctly  assemble  more  reads  with only slightly\n"
                   "                                        higher false-positive rate.\n");
  fprintf (stdout, "  -e, --empirical-freqs                 Disable  empirical base frequencies. (default: use empirical\n"
                   "                                        base frequencies)\n");
  fprintf (stdout, "  -s, --score-method          <int>     Specify the scoring method. (default: 2)\n"
                   "                                        1. OES with +1 for match and -1 for mismatch.\n"
                   "                                        2: Assembly score (AS). Use +1 for match and -1 for mismatch\n"
                   "                                        multiplied by base quality scores.\n"
                   "                                        3: Ignore quality scores and use +1 for a match and -1 for a\n"
                   "                                        mismatch.\n");
  fprintf (stdout, "  -b, --phred-base            <int>     Base PHRED quality score. (default: 33)\n");
  fprintf (stdout, "  -y, --memory                <str>     Specify  the  amount of memory to be used. The number may be\n"
                   "                                        followed  by  one  of  the  letters  K,  M,  or  G  denoting\n"
                   "                                        Kilobytes,  Megabytes and Gigabytes, respectively. Bytes are\n"
                   "                                        assumed in case no letter is specified.\n");
  fprintf (stdout, "  -j, --threads               <int>     Number of threads to use\n");
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
  int     oi;
  char *  ep;
  int      n;
  int      x;

  /* initialization */
  sw->fastq_left    =      NULL;
  sw->fastq_right   =      NULL;
  sw->outfile       =      NULL;
  sw->min_asm_len   =        50;
  sw->max_asm_len   =    999999;
  sw->qual_thres    =         0;
  sw->phred_base    =        33;
  sw->max_uncalled  =       1.0;
  sw->min_overlap   =        10;
  sw->emp_freqs     =         1;  /* use empirical base frequencies as default */
  sw->p_value       =      0.01;
  sw->geom_mean     =         0;
  sw->min_trim_len  =         1;
  sw->score_method  =         2;
  sw->test          =         1;
  sw->memory        = 200000000;
  sw->threads       =         1;

  while ((opt = getopt_long(argc, argv, "b:ef:g:hj:m:n:o:p:q:r:s:t:u:v:y:", long_options, &oi)) != -1)
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

        case 'j':
          sw->threads = (int)strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' || !sw->threads)
           {
             printf ("Invalid number of threads.\n");
             return (0);
           }
          break;
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
             fprintf (stderr, "Invalid min-overlap length.\n");
             return (0);
           }
          break;
        case 'y':
          n = strlen (optarg);
          if (optarg[n - 1] != 'k' && optarg[n - 1] != 'K' &&
              optarg[n - 1] != 'M' && optarg[n - 1] != 'm' &&
              optarg[n - 1] != 'G' && optarg[n - 1] != 'g' &&
              optarg[n - 1] < '0' && optarg[n - 1] > '9')
           {
             fprintf (stderr, "Invalid memory size specified\n");
             return (0);
           }
          switch (optarg[n - 1])
           {
             case 'K':
             case 'k':
               x = 1024;
               optarg[n - 1] = 0;
               break;
             case 'M':
             case 'm':
               x = 1048576;
               optarg[n - 1] = 0;
               break;
             case 'G':
             case 'g':
               x = 1073741824;
               optarg[n - 1] = 0;
               break;
             default:
               x = 1;
           }
          sw->memory = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0')
           {
             fprintf (stderr, "Invalid memory size specified\n");
             return (0);
           }
          sw->memory *= x;
          if (!sw->memory) -- sw->memory;
          break;
      }
   }
  return (sw->fastq_left && sw->fastq_right && sw->outfile);
}
