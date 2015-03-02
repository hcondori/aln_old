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


#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <stdio.h>

typedef struct
{
  char* seq1_id;
  char* seq2_id;
  char* seq1;
  char* seq2;
  int seq1_len;
  int seq2_len;
  char* aln1;
  char* aln2;
  int aln_len;
  int aln_x0;
  int aln_y0;
  int aln_x;
  int aln_y;
  float score;
} alignment_f32;

typedef struct
{
  char* seq1_id;
  char* seq2_id;
  char* seq1;
  char* seq2;
  int seq1_len;
  int seq2_len;
  char* aln1;
  char* aln2;
  int aln_len;
  int aln_x0;
  int aln_y0;
  int aln_x;
  int aln_y;
  int score;
} alignment_i32;

void
alignment_i32_free (alignment_i32* aln);

void
alignment_f32_free (alignment_f32* aln);

void
print_alignment_f32 (FILE *f, alignment_f32* aln);

void
print_alignment_i32 (FILE *f, alignment_i32* aln);

#endif /* ALIGNMENT_H_ */
