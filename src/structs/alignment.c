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

#include <stdlib.h>
#include <string.h>

#include "structs/alignment.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/*
 * Frees resources related to this structure.
 */
void
alignment_f32_free (alignment_f32* aln)
{
  free (aln->seq1_id);
  free (aln->seq2_id);
  free (aln->seq1);
  free (aln->seq2);
  free (aln->aln1);
  free (aln->aln2);
}

void
alignment_i32_free (alignment_i32* aln)
{
  free (aln->seq1_id);
  free (aln->seq2_id);
  free (aln->seq1);
  free (aln->seq2);
  free (aln->aln1);
  free (aln->aln2);
}

void
gen_aln (int len, char* aln1, char* aln2, char* matches, int* gaps)
{
  *gaps = 0;
  for (int i = 0; i < len; i++)
    {
      matches[i] =
	  (aln1[i] != '-' && aln2[i] != '-') ?
	      (aln1[i] == aln2[i] ? '|' : '.') : ' ';

      if (aln1[i] == '-')
	(*gaps)++;
      if (aln2[i] == '-')
	(*gaps)++;
    }
}

void
print_alignment_f32 (FILE *f, alignment_f32* aln)
{
  int wrap = 50;
  int gaps;
  char* matches = (char*) malloc ((aln->aln_len + 1) * sizeof(char));
  matches[aln->aln_len] = '\0';
  gen_aln (aln->aln_len, aln->aln1, aln->aln2, matches, &gaps);

  fprintf (f, "#=======================================\n");
  fprintf (f, "# 1: %s\n", aln->seq1_id);
  fprintf (f, "# 2: %s\n", aln->seq2_id);
  fprintf (f, "%d\n", aln->aln_x);
  fprintf (f, "%d\n", aln->aln_y);

  fprintf (f, "# Length: %d\n", aln->aln_len);
  fprintf (f, "# Gaps:     %d/%d (%.1f%%)\n", gaps, aln->aln_len,
	   ((float) gaps / (float) aln->aln_len) * 100);
  fprintf (f, "# Score:    %.1f\n\n\n", aln->score);
  fprintf (f, "#=======================================\n\n\n");

  int len = aln->aln_len;
  int n = (len + wrap - 1) / wrap;

  for (int i = 0; i < n; i++)
    {
      fprintf (f, "1: ");
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (aln->aln1[j], f);
	}
      fputc ('\n', f);
      fprintf (f, "   ", aln->score);
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (matches[j], f);
	}
      fputc ('\n', f);
      fprintf (f, "2: ");
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (aln->aln2[j], f);
	}
      fputc ('\n', f);
      fputc ('\n', f);
    }

  /*
   fprintf (f, "%s\n", aln->aln1);
   fprintf (f, "%s\n", matches);
   fprintf (f, "%s\n", aln->aln2);
   */
  fprintf (f, "\n\n");

  free (matches);
}

void
print_alignment_i32 (FILE *f, alignment_i32* aln)
{
  int wrap = 50;
  int gaps;
  char* matches = (char*) malloc ((aln->aln_len + 1) * sizeof(char));
  matches[aln->aln_len] = '\0';
  gen_aln (aln->aln_len, aln->aln1, aln->aln2, matches, &gaps);

  fprintf (f, "#=======================================\n");

  fprintf (f, "# 1: %s\n", aln->seq1_id);
  fprintf (f, "# 2: %s\n", aln->seq2_id);
  fprintf (f, "%d\n", aln->aln_x);
  fprintf (f, "%d\n", aln->aln_y);

  fprintf (f, "# Length: %d\n", aln->aln_len);
  fprintf (f, "# Gaps:     %d/%d (%.1f%%)\n", gaps, aln->aln_len,
	   ((float) gaps / (float) aln->aln_len) * 100);
  fprintf (f, "# Score:    %d.0\n\n\n", aln->score);
  fprintf (f, "#=======================================\n\n\n");

  int len = aln->aln_len;
  int n = (len + wrap - 1) / wrap;

  for (int i = 0; i < n; i++)
    {
      fprintf (f, "1: ");
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (aln->aln1[j], f);
	}
      fputc ('\n', f);
      fprintf (f, "   ", aln->score);
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (matches[j], f);
	}
      fputc ('\n', f);
      fprintf (f, "2: ");
      for (int j = i * wrap; j < MIN(wrap * (i + 1), len); j++)
	{
	  fputc_unlocked (aln->aln2[j], f);
	}
      fputc ('\n', f);
      fputc ('\n', f);
    }

  fprintf (f, "\n\n");

  free (matches);
}

