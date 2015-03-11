/*

 Copyright (C) 2015 Héctor Condori Alagón.

 This file is part of ALN, the massive Smith-Waterman pairwise aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Foobar is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <unistd.h>
#include <pthread.h>

#include "structs/alignment.h"
#include <avx2_aligners/avx2_sw.h>
#include <avx2_aligners/avx2_swi.h>
#include "common/backtrack.h"
#include "avx_aligners/avx_sw.h"
#include "common/utils.h"
#include "structs/queue.h"

//#include "fasta_parser.h"

/*
 void yyerror(char *s)
 {
 fprintf(stderr, "error: %s\n", s);
 }
 */

#define MAX_BATCH_COUNT 512
#define THRESHOLD_COUNT  (MAX_BATCH_COUNT / 4)

queue* q;

pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t not_empty_cond = PTHREAD_COND_INITIALIZER;
pthread_cond_t need_refill_cond = PTHREAD_COND_INITIALIZER;

int more_batches = 1;
int do_int = 0;
int do_float = 1;

char* seqs1_fname = NULL;
char* seqs2_fname = NULL;
char* out_fname = NULL;
FILE* seqs1_file;
FILE* seqs2_file;
FILE* out_file;

float gap_open = 10.0f,
gap_extend = 0.5f;	//default values
int gap_open_i = 10.0, gap_extend_i = 1;	//default values
float* sm;
int* smi;

//pointer to pairwise aligner
alignment_f32*
(*pw_aligner) (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
	       char* seqs2[8], float* subs_matrix, float gap_open,
	       float gap_extend, int dup_strings);

typedef struct
{
  int pair_count;
  char** seqs1_ids;
  char** seqs2_ids;
  char** seqs1;
  char** seqs2;
} aln_batch;

void
aln_batch_shallow_free (aln_batch* batch)
{
  free (batch->seqs1_ids);
  free (batch->seqs2_ids);
  free (batch->seqs1);
  free (batch->seqs2);
  free (batch);
}

void*
consumer (void* ptr)
{
  alignment_f32* alns;
  alignment_i32* i32_alns;

  aln_batch* batch = NULL;

  while (1)
    {
      pthread_mutex_lock (&queue_mutex);
      while ((batch = dequeue (q)) == NULL)
	{
	  pthread_cond_wait (&not_empty_cond, &queue_mutex);
	}
      //pthread_cond_signal (&need_refill_cond);
      pthread_mutex_unlock (&queue_mutex);

      if (batch->pair_count == 0)	//we are done
	{
	  pthread_mutex_lock (&queue_mutex);
	  enqueue (q, batch);	//reenque it for other processes to see it
	  pthread_mutex_unlock (&queue_mutex);
	  break;
	}

      if (do_float)
	{
	  alns = pw_aligner (batch->seqs1_ids, batch->seqs2_ids, batch->seqs1,
	  batch->seqs2, sm, gap_open, gap_extend, 0);
	  /*alns = avx2_sw_f32_with_match (batch->seqs1_ids, batch->seqs2_ids,
					 batch->seqs1, batch->seqs2, 5, -4,
					 gap_open, gap_extend, 0);*/

	  for (int i = 0; i < batch->pair_count; i++)
	    {
	      pthread_mutex_lock (&print_mutex);
	      print_alignment_f32 (out_file, &alns[i]);
	      pthread_mutex_unlock (&print_mutex);
	    }
	  for (int i = 0; i < 8; i++)
	    alignment_f32_free (&alns[i]);
	  free (alns);
	}
      else
	{
	  /*i32_alns = avx2_sw_i32_with_matrix (batch->seqs1_ids,
	   batch->seqs2_ids, batch->seqs1,
	   batch->seqs2, smi, gap_open_i,
	   gap_extend_i, 0);*/
	  i32_alns = avx2_sw_i32_with_match (batch->seqs1_ids, batch->seqs2_ids,
					     batch->seqs1, batch->seqs2, 5, -4,
					     gap_open_i, gap_extend_i, 0);
	  for (int i = 0; i < batch->pair_count; i++)
	    {
	      pthread_mutex_lock (&print_mutex);
	      print_alignment_i32 (out_file, &i32_alns[i]);
	      pthread_mutex_unlock (&print_mutex);
	    }
	  for (int i = 0; i < 8; i++)
	    alignment_i32_free (&i32_alns[i]);
	  free (i32_alns);
	}
      aln_batch_shallow_free (batch);
      batch = NULL;
    }
  pthread_exit (0);
}

void
print_help ()
{

}

/*
 * FASTA parser
 */
int
get_sequences (FILE* f, char** seqs_ids, char** seqs, int* state, char* ch)
{
  int seq_size = 128;
  char* seq_buffer = (char*) malloc (seq_size * sizeof(char));
  int seqs_n = 0;
  int pos = 0;

  char seqid_buffer[128];
  int bufferid_pos = 0;

  //*state = 0;   //0: start, 1: seq header, 2: seq id, 3: seq

  while (seqs_n < 8 && (*ch = fgetc (f)) != EOF)
    {
      switch (*state)
	{
	case 0:
	  switch (*ch)
	    {
	    case '>':
	      *state = 2;
	      break;
	    default:
	      break;
	    }
	  break;
	case 1:
	  if (*ch == '\n')
	    {
	      *state = 3;
	    }
	  break;
	case 2:
	  switch (*ch)
	    {
	    case '\n':
	      *state = 3;
	      seqid_buffer[bufferid_pos] = '\0';
	      seqs_ids[seqs_n] = strdup (seqid_buffer);
	      bufferid_pos = 0;
	      break;
	    case ' ':
	      *state = 1;
	      seqid_buffer[bufferid_pos] = '\0';
	      seqs_ids[seqs_n] = strdup (seqid_buffer);
	      bufferid_pos = 0;
	      break;
	    default:
	      seqid_buffer[bufferid_pos++] = *ch;
	      break;
	    }
	  break;
	case 3:	//read the sequence
	  if (*ch == '>')	//new record
	    {
	      seq_buffer[pos] = '\0';
	      seqs[seqs_n++] = strdup (seq_buffer);
	      *state = 2;
	      pos = 0;
	    }
	  else if (*ch == '\n')
	    continue;
	  else
	    {
	      if (pos == seq_size - 1)  //we should grow the data buffer
		{
		  seq_size += 128;
		  seq_buffer = (char*) realloc (seq_buffer,
						seq_size * sizeof(char));
		}
	      seq_buffer[pos++] = *ch;
	      //puts ("**************");
	    }
	  break;
	default:
	  break;
	}
    }

  if (*ch == EOF && *state == 3 && pos > 0)
    {
      seq_buffer[pos] = '\0';
      seqs[seqs_n++] = strdup (seq_buffer);
    }

  //init empty chars
  for (int i = seqs_n; i < 8; i++)
    {
      seqs_ids[i] = strdup ("");
      seqs[i] = strdup ("");
    }

  //if (*ch != EOF)
  //fseek (f, -1, SEEK_CUR);

  free (seq_buffer);
  return seqs_n;
}

int verbose_flag = 0;

static struct option long_options[] =
  {
    { "verbose", no_argument, &verbose_flag, 1 },
    { "brief", no_argument, &verbose_flag, 0 },
    { "help", no_argument, 0, 'h' },
    { "jobs", no_argument, 0, 'j' },
    { "gap-open", required_argument, 0, 'O' },
    { "gap_extend", required_argument, 0, 'E' },
    { "matrix", optional_argument, 0, '\0' },
    { "match", optional_argument, 0, '\0' },
    { "mismatch", optional_argument, 0, '\0' },
    { "data-type", required_argument, 0, 'd' },
    { "seqs-a", required_argument, 0, 'a' },
    { "seqs-b", required_argument, 0, 'b' },
    { "query", required_argument, 0, 'Q' },
    { "refs", required_argument, 0, 'R' },
    { "alignment-type", required_argument, 0, 't' },
    { "output", required_argument, 0, 'o' },
    { 0, 0, 0, 0 } };

int
main (int argc, char** argv)
{
  int jobs = 1;
  int jobs_is_set = 0;

  int platform = 0;	//CPU only by default

  int c;

  while (1)
    {
      int option_index = 0;
      c = getopt_long (argc, argv, "a:b:d:hj:O:E:o:", long_options,
		       &option_index);
      if (c == -1)
	break;
      switch (c)
	{
	case 0:
	  if (long_options[option_index].flag != 0)
	    break;
	  printf ("option %s", long_options[option_index].name);
	  if (optarg)
	    printf (" with arg %s", optarg);
	  printf ("\n");
	  break;
	case 'a':
	  seqs1_fname = strdup (optarg);
	  break;
	case 'b':
	  seqs2_fname = strdup (optarg);
	  break;
	case 'd':
	  if (strcmp (optarg, "float") == 0)
	    do_float = 1;
	  else if (strcmp (optarg, "int") == 0)
	    do_float = 0;
	  else
	    {
	      printf ("Error: invalid argument for --dataformat.\n");
	      exit (0);
	    }
	  break;
	case 'j':
	  jobs = atoi (optarg);
	  if (jobs < 1)
	    {
	      printf ("Error: invalid argument for jobs.\n");
	      exit (1);
	    }
	  jobs_is_set = 1;
	  break;
	case 'o':
	  out_fname = strdup (optarg);
	  break;
	case 'h':
	  puts ("Show help");
	  exit (0);
	case 'c':
	  printf ("option -c with value `%s'\n", optarg);
	  break;
	case 'O':
	  gap_open = atof (optarg);
	  gap_open_i = atoi (optarg);
	  break;
	case 'E':
	  gap_extend = atof (optarg);
	  gap_extend_i = atoi (optarg);
	  break;

	case '?':

	  break;

	default:
	  abort ();
	}
    }

  sm = (float*) malloc (sizeof(float) * 128 * 128);

  for (int i = 0; i < 128; ++i)
    {
      for (int j = 0; j < 128; ++j)
	{
	  if (i == j)
	    sm[128 * i + j] = 5;
	  else
	    sm[128 * i + j] = -4;
	}
    }
  *sm = -4;

  smi = (int*) malloc (sizeof(int) * 128 * 128);

  for (int i = 0; i < 128; ++i)
    {
      for (int j = 0; j < 128; ++j)
	{
	  if (i == j)
	    smi[128 * j + i] = 5;
	  else
	    smi[128 * j + i] = -4;
	}
    }
  *smi = -4;

  ssize_t len;
  size_t temp = 0;

  printf ("===============\n");
  printf ("Aln version 0.0\n");
  printf ("===============\n\n");

  int has_sse41 = __builtin_cpu_supports ("sse4.1");
  int has_avx = __builtin_cpu_supports ("avx");
  int has_avx2 = __builtin_cpu_supports ("avx2");

  if (platform == 0)
    {
      if (has_avx2)
	{
	  printf (" * AVX2 support detected.\n");
	  pw_aligner = &avx2_sw_f32_with_matrix;
	}
      else if (has_avx)
	{
	  printf (" * AVX support detected.\n");
	  pw_aligner = &avx_sw_f32_with_matrix;
	}
      else if (has_sse41)
      	{
      	  printf (" * SSE 4.1 support detected.\n");
      	}
      else
	{
	  printf (" * No AVX support found. Using SISD version.");
	}
    }

  int ncpus = sysconf (_SC_NPROCESSORS_ONLN);

  printf (" * %d CPUs detected.\n", ncpus);
  if (jobs_is_set)
    {
      printf (" * Thread number is manually set to %d.\n", jobs);
    }
  else
    {
      jobs = ncpus;
      printf (" * Starting %d threads.\n", jobs);
    }

  seqs1_file = fopen (seqs1_fname, "r");
  if (seqs1_file == NULL)
    {
      printf ("Error at opening file %s.\n", seqs1_fname);
      exit (EXIT_FAILURE);
    }

  seqs2_file = fopen (seqs2_fname, "r");
  if (seqs2_file == NULL)
    {
      printf ("Error at opening file %s.\n", seqs2_fname);
      exit (EXIT_FAILURE);
    }

  if (out_fname == NULL)
    out_file = stdout;
  else
    {
      out_file = fopen (out_fname, "w");
      if (out_file == NULL)
	{
	  printf ("Error at opening file %s.\n", out_fname);
	  exit (EXIT_FAILURE);
	}
    }

  char* seqs1[8];
  char* seqs2[8];

  char* seqs1_ids[8];
  char* seqs2_ids[8];

  int seqs1_n;
  int seqs2_n;

  int pair_count = 0;

  int state1 = 0;
  int state2 = 0;
  char ch1;
  char ch2;

  //has_avx2 = 0;

  q = new_queue ();
  aln_batch* batch;

  pthread_t consumers[jobs];
  for (int i = 0; i < jobs; i++)
    {
      pthread_create (&consumers[i], NULL, consumer, NULL);
    }

  do
    {
      batch = (aln_batch*) malloc (sizeof(aln_batch));
      batch->seqs1_ids = (char**) malloc (8 * sizeof(char*));
      batch->seqs2_ids = (char**) malloc (8 * sizeof(char*));
      batch->seqs1 = (char**) malloc (8 * sizeof(char*));
      batch->seqs2 = (char**) malloc (8 * sizeof(char*));
      seqs1_n = get_sequences (seqs1_file, batch->seqs1_ids, batch->seqs1,
			       &state1, &ch1);
      seqs2_n = get_sequences (seqs2_file, batch->seqs2_ids, batch->seqs2,
			       &state2, &ch2);
      batch->pair_count = MIN(seqs1_n, seqs2_n);

      pthread_mutex_lock (&queue_mutex);
      enqueue (q, batch);
      pthread_cond_signal (&not_empty_cond);
      pthread_mutex_unlock (&queue_mutex);

      pair_count += MIN(seqs1_n, seqs2_n);
    }
  while (seqs1_n > 0 && seqs2_n > 0);

  //no more work to do :)
  pthread_mutex_lock (&queue_mutex);
  //more_batches = 0;
  pthread_cond_broadcast (&not_empty_cond);
  pthread_mutex_unlock (&queue_mutex);

  for (int i = 0; i < jobs; i++)
    {
      pthread_join (consumers[i], NULL);
    }

  free (q);
  fclose (seqs1_file);
  fclose (seqs2_file);

  if (out_file != stdout)
    fclose (out_file);

  if (seqs1_n != seqs2_n)
    {
      printf ("1: %d\n", seqs1_n);
      printf ("2: %d\n", seqs2_n);
      printf ("Warning: sequence count do not match in both files.\n");
    }

  printf ("%d alignments have been performed.\n", pair_count);

  pthread_mutex_destroy (&queue_mutex);
  pthread_mutex_destroy (&print_mutex);
  pthread_cond_destroy (&not_empty_cond);
  pthread_cond_destroy (&need_refill_cond);

  return 0;
}
