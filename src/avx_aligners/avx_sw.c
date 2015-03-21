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

#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "avx_sw.h"
#include "common/backtrack.h"
#include "common/utils.h"
#include "common/avx_common.h"

/*
 *Core function.
 *
 */
int*
avx_fill_table_8_to_8_f32 (int* seqs1, int* seqs2, int x, int y,
			   float* subs_matrix, float gap_open, float gap_extend,
			   float* max_score, int* ipos, int* jpos)
{
  int mask = 0x7FFFFFFF;
  __m256 vmask = _mm256_castsi256_ps (_mm256_set1_epi32 (mask));
  __m256i s1, s2, temp_index, index;
  __m256i v128 = _mm256_set1_epi32 (128);
  __m256 vopen = _mm256_set1_ps (gap_open);
  __m256 vextend = _mm256_set1_ps (gap_extend);
  __m256 vzero = _mm256_setzero_ps ();
  __m256 diff;
  __m256 vepsilon = _mm256_set1_ps (EPSILON);
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256 max = _mm256_setzero_ps ();
  __m256 E, E_sub;
  __m256 F, F_sub;
  __m256 diag, score, H, H_diag, H_left, temp;
  __m256 c_up, c_left, b_up, b_left, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0;

  float buffer[8] __attribute__((aligned(32)));
  int buffer_i[8] __attribute__((aligned(32)));
  int buffer_j[8] __attribute__((aligned(32)));

  float* aF = (float*) aligned_alloc (32, 8 * y * sizeof(float));
  float* aH = (float*) aligned_alloc (32, 8 * y * sizeof(float));
  int* bt_flag = (int*) malloc (x * y * sizeof(int));
  //if we are sure that sequences will be small, we keep them in the stack
  //float aF[8 * y] __attribute__((aligned(32)));
  //float aH[8 * y] __attribute__((aligned(32)));

  memset (aH, 0, 8 * y * sizeof(float));
  memset (bt_flag, 0, y * sizeof(int));

  for (int i = 0; i < 8 * y; i++)
    {
      aF[i] = -FLT_MAX;
    }

  for (int i = 1; i < x; i++)
    {
      E = _mm256_set1_ps (-FLT_MAX);
      H_diag = _mm256_setzero_ps ();
      H = _mm256_setzero_ps ();
      bt_flag[i * y] = 0;
      for (int j = 1; j < y; j++)
	{
	  for (int k = 0; k < 8; k++)
	    {
	      buffer[k] = subs_matrix[128 * (*(seqs1 + 8 * (i - 1) + k))
		  + *(seqs2 + 8 * (j - 1) + k)];
	    }
	  score = _mm256_load_ps (buffer);
	  H_left = _mm256_load_ps (aH + 8 * j);
	  diag = _mm256_add_ps (H_diag, score);
	  H_diag = H_left;

	  E_sub = _mm256_sub_ps (E, vextend);	//for now, E is E_up
	  E = _mm256_sub_ps (H, vopen);		//for now, H is H_up
	  E = _mm256_max_ps (E, E_sub);

	  F_sub = _mm256_load_ps (aF + 8 * j);
	  F_sub = _mm256_sub_ps (F_sub, vextend);
	  F = _mm256_sub_ps (H_left, vopen);
	  F = _mm256_max_ps (F, F_sub);
	  _mm256_store_ps (aF + 8 * j, F);

	  H = _mm256_max_ps (E, F);
	  H = _mm256_max_ps (H, diag);
	  H = _mm256_max_ps (H, vzero);
	  _mm256_store_ps (aH + 8 * j, H);

	  //logic tests

	  diff = _mm256_sub_ps (H, vzero);
	  h_gt_0 = _mm256_movemask_ps (
	      _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ));
	  diff = _mm256_sub_ps (H, E);
	  diff = _mm256_and_ps (diff, vmask);	//absolute value
	  H_eq_E = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
	  h_eq_e = _mm256_movemask_ps (H_eq_E);
	  diff = _mm256_sub_ps (H, F);
	  diff = _mm256_and_ps (diff, vmask);	//absolute value
	  H_eq_F = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
	  h_eq_f = _mm256_movemask_ps (H_eq_F);

	  //********FLAGS********
	  diff = _mm256_sub_ps (E, E_sub);
	  diff = _mm256_and_ps (diff, vmask);	//absolute value
	  c_up = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // E[i,j] == E[i,j-1]-gap_extent ?
	  flag = _mm256_movemask_ps (c_up); // & h_eq_e;		//c_up

	  diff = _mm256_sub_ps (F, F_sub);
	  diff = _mm256_and_ps (diff, vmask);	//absolute value
	  c_left = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // F[i,j] == F[i-1,j]-gap_extent ?
	  flag <<= 8;
	  flag |= _mm256_movemask_ps (c_left); // & h_eq_f;		//c_left

	  diff = _mm256_sub_ps (H, diag);
	  diff = _mm256_and_ps (diff, vmask);	//absolute value
	  H_eq_diag = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
	  h_eq_d = _mm256_movemask_ps (H_eq_diag);
	  flag <<= 8;
	  flag |= (h_eq_e | h_eq_d) & h_gt_0;		//b_up

	  //b_left
	  flag <<= 8;
	  flag = flag | (((h_eq_f & ~h_eq_e) | h_eq_d) & h_gt_0);

	  bt_flag[y * i + j] = flag;

	  diff = _mm256_sub_ps (H, max);
	  temp = _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ);

	  _mm256_store_ps (buffer, temp);
	  for (int k = 0; k < 8; k++)
	    {
	      if (buffer[k])
		{
		  ipos[k] = i;
		  jpos[k] = j;
		}
	    }
	  max = _mm256_max_ps (H, max);
	}
    }
  _mm256_store_ps (max_score, max);
  free (aF);
  free (aH);
  return bt_flag;
}

alignment*
avx_sw_f32_with_matrix (char** seqs1_id, char** seqs2_id, char** seqs1,
			char** seqs2, float* subs_matrix, float gap_open,
			float gap_extend, int dup_strings)
{
  SIMD_SW_8_TO_8(avx_fill_table_8_to_8_f32, float, ALN_FLOAT32)
}

void
avx_sw_f32_with_matrix_inplace (alignment** alignments, float* subs_matrix,
				float gap_open, float gap_extend)
{
  puts("You shouldnt be here!!!!!!!!!!!!");
  SIMD_SW_8_TO_8_INPLACE(avx_fill_table_8_to_8_f32, float, ALN_FLOAT32)
}
