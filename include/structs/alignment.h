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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <stdio.h>
#include <stdint.h>

typedef enum
{
  ALN_SCORE_INT8, ALN_SCORE_INT16, ALN_SCORE_INT32, ALN_SCORE_FLOAT32
} aln_score_type;

typedef union
{
  int8_t i8;
  int16_t i16;
  int32_t i32;
  float f32;
} aln_score;

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
  aln_score score;
  aln_score_type type;
} alignment;

void
alignment_free (alignment* aln);

void
print_alignment (FILE *f, alignment* aln);

#endif /* ALIGNMENT_H_ */
