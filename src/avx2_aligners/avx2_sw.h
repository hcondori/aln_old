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


#ifndef AVX2_SW_H_
#define AVX2_SW_H_

#include <assert.h>
#include <float.h>
#include <immintrin.h>

#include "structs/alignment.h"
#include "common/backtrack.h"
#include "common/utils.h"

alignment_f32*
avx2_sw_f32_with_matrix (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
			 char* seqs2[8], float* subs_matrix, float gap_open,
			 float gap_extend, int dup_strings);

alignment_f32*
avx2_sw_f32_with_match (char* seqs1_id[8], char* seqs2_id[8], char** seqs1,
			char** seqs2, float match, float mismatch,
			float gap_open, float gap_extend, int dup_strings);

#endif /* AVX2_SW_H_ */
