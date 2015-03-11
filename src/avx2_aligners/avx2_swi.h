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


#ifndef AVX2_SWI_H_
#define AVX2_SWI_H_

alignment_i32*
avx2_sw_i32_with_matrix (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
			 char* seqs2[8], int* subs_matrix, int gap_open,
			 int gap_extend, int dup_strings);

alignment_i32*
avx2_sw_i32_with_match (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
			char* seqs2[8], int match, int mismatch, int gap_open,
			int gap_extend, int dup_strings);

#endif /* AVX2_SWI_H_ */
