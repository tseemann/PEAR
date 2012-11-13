#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "args.h"

static struct option long_options[] =
 {
   { "max-asm-length",   required_argument, NULL, 'm' },
   { "min-asm-length",   required_argument, NULL, 'n' },
   { "quality-theshold", required_argument, NULL, 'q' },
   { "left-fastq",       required_argument, NULL, 'l' },
   { "right-fastq",      required_argument, NULL, 'r' },
   { "score-method",     required_argument, NULL, 's' },
   { "help",             no_argument,       NULL, 'h' },
   { "min-overlap",      required_argument, NULL, 'o' },
   { NULL,               0,                 NULL, 0   }
 };

void usage (void)
{
  fprintf (stdout, " ____  _____    _    ____ \n"); 
  fprintf (stdout, "|  _ \\| ____|  / \\  |  _ \\\n");
  fprintf (stdout, "| |_) |  _|   / _ \\ | |_) |\n");
  fprintf (stdout, "|  __/| |___ / ___ \\|  _ <\n");
  fprintf (stdout, "|_|   |_____/_/   \\_\\_| \\_\\\n");
  fprintf (stdout, "\n.oOo. Pair-End AssembleR .oOo.\n");
  fprintf (stdout, "\n\n"); 
  
  fprintf (stdout, "Usage: pear <options>\n");
  fprintf (stdout, "Standard (mandatory):\n");
  fprintf (stdout, "  -l, --left-fastq          <str>     Left pairend FASTQ file.\n");
  fprintf (stdout, "  -r, --right-fastq         <str>     Right pairend FASTQ file.\n");
  fprintf (stdout, "Optional:\n");
  fprintf (stdout, "  -m, --max-asm-length      <int>     Maximum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                      The assembled sequence can be arbitrary long if set\n"
                   "                                      to 0. (default: 0)\n");
  fprintf (stdout, "  -n, --min-asm-length      <int>     Minimum possible size of the assembled sequence.\n");
  fprintf (stdout, "                                      To disable it set to 0. (default: 0)\n");
  fprintf (stdout, "  -q, --quality-threshold   <int>     Quality score threshold. Dont remember what that\n"
                   "                                      really is. (default: 20)\n");
  fprintf (stdout, "  -s, --score-method        <int>     Scoring method\n");
  fprintf (stdout, "  -o, --min-overlap         <int>     Minimum overlap\n");
  fprintf (stdout, "  -h, --help                          This help screen.\n\n");
}

int decode_switches (int argc, char * argv[], struct user_args * sw)
{
  int    opt;
  int    oi;
  char * ep;

  /* initialization */
  sw->fastq_left  = NULL;
  sw->fastq_right = NULL;
  sw->min_asm_len = 0;
  sw->max_asm_len = 0;
  sw->qual_thres  = 0;

  while ((opt = getopt_long(argc, argv, "m:n:q:l:r:s:o:h", long_options, &oi)) != -1)
   {
     switch (opt)
      {
        case 'l':
          sw->fastq_left  = optarg;
          break;

        case 'r':
          sw->fastq_right = optarg;
          break;

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
        case 's':
          sw->score_method = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 3\n");
             return (0);
           }
          break;

        case 'o':
          sw->min_overlap = (int) strtol (optarg, &ep, 10);
          if (ep == optarg || *ep != '\0' )
           {
             printf ("Problem 4\n");
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

        case 'h':
          return (0);
      }
   }
  return (sw->fastq_left && sw->fastq_right);
}
