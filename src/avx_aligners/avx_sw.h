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

#ifndef AVX_SW_H_
#define AVX_SW_H_

#include "structs/sequence.h"
#include "structs/alignment.h"

alignment*
sw_1_to_1_f32_mith_matrix (char* seq1, char* seq2, float* subs_matrix,
			   float gap_open, float gap_extend);

void
avx_sw_f32_with_matrix_inplace (alignment** alignments, float* subs_matrix,
				float gap_open, float gap_extend);

alignment*
avx_sw_f32_with_matrix (char** seqs1_id, char** seqs2_id, char** seqs1,
			char** seqs2, float* subs_matrix, float gap_open,
			float gap_extend, int dup_strings);

#endif /* AVX_SW_H_ */
