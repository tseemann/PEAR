#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "args.h"

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

void usage (void)
{
  fprintf (stdout, " ____  _____    _    ____ \n"); 
  fprintf (stdout, "|  _ \\| ____|  / \\  |  _ \\\n");
  fprintf (stdout, "| |_) |  _|   / _ \\ | |_) |\n");
  fprintf (stdout, "|  __/| |___ / ___ \\|  _ <\n");
  fprintf (stdout, "|_|   |_____/_/   \\_\\_| \\_\\\n");
  fprintf (stdout, "\n.oOo. Pair-End AssembleR .oOo.\n");
  fprintf (stdout, "PEAR v0.1 by Tomas Flouri and Jiajie Zhang\n");
  fprintf (stdout, "\n\n"); 
  
  fprintf (stdout, "Usage: pear <options>\n");
  fprintf (stdout, "Standard (mandatory):\n");
  fprintf (stdout, "  -f, --forward-fastq         <str>     Forward pairend FASTQ file.\n");
  fprintf (stdout, "  -r, --reverse-fastq         <str>     Reverse pairend FASTQ file.\n");
  fprintf (stdout, "  -o, --output                <str>     Output filename.\n");
  fprintf (stdout, "Optional:\n");
  fprintf (stdout, "  -p, --p-value               <float>   Use a p-value from the set { 0.05, 0.01, 0.001, 0.0001 }. If\n"
                   "                                        the p-value of the assembled reads exceeds the specified p-value, the\n"
                   "                                        reads will be output unassembled. (default: 0.01)\n");
  fprintf (stdout, "  -v, --min-overlap           <int>     Minimum overlap (default: 10)\n");
  fprintf (stdout, "  -m, --max-assembly-length   <int>     Maximum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                        The assembled sequence can be arbitrary long if set\n"
                   "                                        to 0. (default: 0)\n");
  fprintf (stdout, "  -n, --min-assembly-length   <int>     Minimum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                        To disable it set to 0. (default: 0)\n");
  fprintf (stdout, "  -t, --min-trim-length       <int>     Minimum length of reads after trimming the low quality"
                   "                                        part (default: 1)\n");
  fprintf (stdout, "  -q, --quality-threshold     <int>     Quality score threshold used for trimming the low quality"
                   "                                        part of the reads (default: 0)\n");
  fprintf (stdout, "  -u, --max-uncalled-base     <float>   Maximal proportion of uncalled bases. A number between 0 and 1.\n"
                   "                                        Set to 0 to discard all reads that contain uncalled bases, or\n"
                   "                                        1 to process all sequences independent on the number of uncalled\n"
                   "                                        bases. (default: 1)\n");
  fprintf (stdout, "  -g, --geometric-mean        <float>   Minimum value of geometric mean of the quality score. (default: 0)\n");
  fprintf (stdout, "  -e, --empirical-freqs                 Use empirical base frequencies.\n");
  fprintf (stdout, "  -s, --score-method          <int>     Scoring method\n");
  fprintf (stdout, "  -b, --phred-base            <int>     Base Phred quality score (default: 33)\n");
  fprintf (stdout, "  -h, --help                            This help screen.\n\n");
}

int decode_switches (int argc, char * argv[], struct user_args * sw)
{
  int    opt;
  int    oi;
  char * ep;

  /* initialization */
  sw->fastq_left    = NULL;
  sw->fastq_right   = NULL;
  sw->min_asm_len   = 0;
  sw->max_asm_len   = 999999;
  sw->qual_thres    = 0;
  sw->phred_base    = 33;
  sw->max_uncalled  = 1.0;
  sw->min_overlap   = 10;
  sw->emp_freqs     = 0;    /* do not use empirical base frequencies as default */
  sw->p_value       = 0.01;
  sw->geom_mean     = 0;
  sw->min_trim_len  = 1;


  while ((opt = getopt_long(argc, argv, "b:ef:g:hm:n:o:p:q:r:s:t:u:v:", long_options, &oi)) != -1)
   {
     switch (opt)
      {
        case 'b':
          sw->phred_base  = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0')
           {
             printf ("Problem 6\n");
             return (0);
           }
          break;

        case 'e':
          sw->emp_freqs = 1;
          break;

        case 'f':
          sw->fastq_left  = optarg;
          break;

        case 'g':
          sw->geom_mean = strtod (optarg, &ep);     /* TODO: check this line */
          if (ep == optarg || *ep != '\0' || sw->geom_mean < 0 || sw->geom_mean > 1)
           {
             printf ("Problem 7\n");
             return (0);
           }
          break;

        case 'h':
          return (0);

        case 'm':
          sw->max_asm_len = (int)strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 1\n");
             return (0);
           }
          break;

        case 'n':
          sw->min_asm_len = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 2\n");
             return (0);
           }
          break;

        case 'o':
          sw->outfile = optarg;
          break;

        case 'p':
          if (!strcmp (optarg, "0.05") || !strcmp (optarg, "0.01") || !strcmp (optarg, "0.001") || !strcmp (optarg, "0.0001"))
           {
             sw->p_value = strtod (optarg, &ep);
           }
          else
           {
             printf ("Problem 8\n");
             return (0);
           }
          break;

        case 'q':
          sw->qual_thres = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 5\n");
             return (0);
           }
          break;

        case 'r':
          sw->fastq_right = optarg;
          break;

        case 's':
          sw->score_method = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 3\n");
             return (0);
           }
          break;

        case 't':
          sw->min_trim_len = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 10\n");
             return (0);
           }
          break;

        case 'u':
          sw->max_uncalled = strtod (optarg, &ep);     /* TODO: check this line */
          
          if (ep == optarg || *ep != '\0' || sw->max_uncalled < 0 || sw->max_uncalled > 1)
           {
             printf ("Problem 9\n");
             return (0);
           }
          break;

        case 'v':
          sw->min_overlap = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 4\n");
             return (0);
           }
          break;
      }
   }
  return (sw->fastq_left && sw->fastq_right);
}
