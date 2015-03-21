/*

 Copyright (C) 2015 Héctor Condori Alagón.

 This file is part of ALN, the massive Smith-Waterman pairwise sequence aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALN is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALN.  If not, see <http://www.gnu.org/licenses/>.
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <unistd.h>
#include <pthread.h>

#include "structs/alignment.h"
#include "avx_aligners/avx_sw.h"
#include "avx2_aligners/avx2_sw.h"
#include "avx2_aligners/avx2_swi.h"
#include "common/backtrack.h"
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

Queue* q_f32;
Queue* q_i8;
Queue* q_i16;
Queue* q_i32;
Queue* q_ready;

pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t q_ready_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

    pthread_cond_t not_empty_cond = PTHREAD_COND_INITIALIZER;
    pthread_cond_t need_refill_cond = PTHREAD_COND_INITIALIZER;
    pthread_cond_t ready_alignments_cond = PTHREAD_COND_INITIALIZER;

    int more_batches = 1;
    int do_int = 0;
    int do_float = 1;

    int use_match = 0;
    int use_matrix = 1;

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
/*alignment*
 (*pw_aligner) (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
 char* seqs2[8], float* subs_matrix, float gap_open,
 float gap_extend, int dup_strings);*/
void
(*pw_aligner) (alignment** alignments, float* subs_matrix, float gap_open,
	       float gap_extend);

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
worker (void* ptr)
{
  Queue* current_queue = NULL;
  int batch_size;

  if (do_float)
    {
      current_queue = q_f32;
      batch_size = 8;		//we must fetch 8 alignments.
    }
  else
    {
      current_queue = q_i8;
      batch_size = 32;		//we must fetch 32 alignments at the beginning.
    }

  alignment* alns[8];
  alignment* aln = NULL;

  while (1)
    {
      pthread_mutex_lock (&queue_mutex);

      for (int i = 0; i < 8; i++)
	{
	  while ((aln = dequeue (q_f32)) == NULL)
	    pthread_cond_wait (&not_empty_cond, &queue_mutex);
	  alns[i] = aln;
	  //puts(aln->seq2->sequence);
	}
      pthread_mutex_unlock (&queue_mutex);

      if (alignment_is_null (alns[0]))
	{
	  puts ("Done XD");
	  pthread_exit (0);
	}

      if (do_float)
	{
	  pw_aligner (alns, sm, gap_open, gap_extend);
	  pthread_mutex_lock (&q_ready_mutex);
	  for (int i = 0; i < 8; i++)
	    enqueue (q_ready, alns[i]);	//some elements might be null in the last batch
	  pthread_cond_signal (&ready_alignments_cond);
	  pthread_mutex_unlock (&q_ready_mutex);
	}
      else
	{
	  /*alns = avx2_sw_i32_with_matrix (batch->seqs1_ids, batch->seqs2_ids,
	   batch->seqs1, batch->seqs2, smi,
	   gap_open_i, gap_extend_i, 0);*/
	  /*i32_alns = avx2_sw_i32_with_match (batch->seqs1_ids, batch->seqs2_ids,
	   batch->seqs1, batch->seqs2, 5, -4,
	   gap_open_i, gap_extend_i, 0);*/
	}
    }

  pthread_exit (0);
}

/*
 * Dedicated thread for generating output
 */
void*
print_worker (void* ptr)
{
  alignment* aln;

  while (1)
    {
      pthread_mutex_lock (&q_ready_mutex);

      while ((aln = dequeue (q_ready)) == NULL)
	pthread_cond_wait (&ready_alignments_cond, &q_ready_mutex);

      if (alignment_is_null (aln))
	{
	  pthread_mutex_unlock (&q_ready_mutex);
	  pthread_exit (0);
	}

      pthread_mutex_unlock (&q_ready_mutex);

      print_alignment (out_file, aln);
      alignment_free (aln);
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
load_sequence (FILE* f, char** seqs_ids, char** seqs, int* state, char* ch)
{
  int seq_size = 128;
  char* seq_buffer = (char*) malloc (seq_size * sizeof(char));
  int seqs_n = 0;
  int pos = 0;

  char seqid_buffer[128];
  int bufferid_pos = 0;

  //*state = 0;   //0: start, 1: seq header, 2: seq id, 3: seq

  while (seqs_n < 1 && (*ch = fgetc (f)) != EOF)
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
  for (int i = seqs_n; i < 1; i++)
    {
      seqs_ids[i] = NULL;
      seqs[i] = NULL;
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
  //has_avx = 0;

  if (platform == 0)
    {
      if (has_avx2)
	{
	  printf (" * AVX2 support detected.\n");
	  pw_aligner = &avx2_sw_f32_with_matrix_inplace;
	}
      else if (has_avx)
	{
	  printf (" * AVX support detected.\n");
	  pw_aligner = &avx_sw_f32_with_matrix_inplace;
	}
      else if (has_sse41)
	{
	  printf (" * SSE 4.1 support detected.\n");
	}
      else
	{
	  printf (" * No AVX or SSE support found. Using SISD version.");
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

  char* seq1_sequence;
  char* seq2_sequence;

  char* seq1_id;
  char* seq2_id;

  int seq1_n = 0;
  int seq2_n = 0;

  int pair_count = 0;

  int state1 = 0;
  int state2 = 0;
  char ch1;
  char ch2;

  q_f32 = Queue_new ();
  q_i8 = Queue_new ();
  q_i16 = Queue_new ();
  q_i32 = Queue_new ();
  q_ready = Queue_new ();

  pthread_t consumers[jobs];
  pthread_t print_consumer;

  for (int i = 0; i < jobs; i++)
    pthread_create (&consumers[i], NULL, worker, NULL);
  pthread_create (&print_consumer, NULL, print_worker, NULL);

  alignment* aln;
  Sequence* seq1;
  Sequence* seq2;

  do
    {
      aln = (alignment*) malloc (sizeof(alignment));
      seq1_n = load_sequence (seqs1_file, &seq1_id, &seq1_sequence, &state1,
			      &ch1);
      seq2_n = load_sequence (seqs2_file, &seq2_id, &seq2_sequence, &state2,
			      &ch2);
      seq1 = Sequence_new (seq1_id, seq1_sequence, strlen (seq1_sequence));
      seq2 = Sequence_new (seq2_id, seq2_sequence, strlen (seq2_sequence));
      aln->seq1 = seq1;
      aln->seq2 = seq2;
      pthread_mutex_lock (&queue_mutex);
      enqueue (q_f32, aln);
      if ((pair_count + 1) % 8 == 0)
	pthread_cond_signal (&not_empty_cond);
      pthread_mutex_unlock (&queue_mutex);
      pair_count += MIN(seq1_n, seq2_n);
    }
  while (seq1_n && seq2_n);

  //centinel
  alignment* centinel = (alignment*) calloc (1, sizeof(alignment));

  //no more work to do :)
  pthread_mutex_lock (&queue_mutex);
  for (int i = 0; i < (jobs + 1) * 8; i++)
    {
      enqueue (q_f32, centinel);
    }
  //more_batches = 0;
  pthread_cond_broadcast (&not_empty_cond);
  pthread_mutex_unlock (&queue_mutex);

  for (int i = 0; i < jobs; i++)
    pthread_join (consumers[i], NULL);
  pthread_join (print_consumer, NULL);

  free (q_i8);
  free (q_i16);
  free (q_i32);
  free (q_f32);
  free (q_ready);

  if (seqs1_file != NULL)
    fclose (seqs1_file);
  if (seqs2_file != NULL)
    fclose (seqs2_file);

  if (out_file != stdout)
    fclose (out_file);

  if (seq1_n != seq2_n)
    {
      printf ("1: %d\n", seq1_n);
      printf ("2: %d\n", seq2_n);
      printf ("Warning: sequence count do not match in both files.\n");
    }

  printf ("%d alignments have been performed.\n", pair_count);

  pthread_mutex_destroy (&queue_mutex);
  pthread_mutex_destroy (&q_ready_mutex);
  pthread_mutex_destroy (&print_mutex);
  pthread_cond_destroy (&not_empty_cond);
  pthread_cond_destroy (&need_refill_cond);
  pthread_cond_destroy (&ready_alignments_cond);

  return 0;
}
