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
#include <limits.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "structs/alignment.h"
#include "common/backtrack.h"
#include "common/utils.h"

/**
 * Core function
 *
 * @param overflow: alignments that need to be rerun with wider registers
 *
 */
int*
avx2_fill_table_8_to_8_i8 (int* seqs1, int* seqs2, int x, int y,
			   int* subs_matrix, int gap_open, int gap_extend,
			   int max_score[32], int ipos[32], int jpos[32],
			   int overflow[32])
{
  __m256i s1, s2, temp_index, index;
  __m256i v128 = _mm256_set1_epi32 (128);
  __m256i vopen = _mm256_set1_epi32 (gap_open);
  __m256i vextend = _mm256_set1_epi32 (gap_extend);
  __m256i vzero = _mm256_setzero_si256 ();
  __m256i vepsilon = _mm256_set1_epi32 (EPSILON);
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256i max = _mm256_setzero_si256 ();
  __m256i E, E_sub;
  __m256i F, F_sub;
  __m256i diag, score, H, H_diag, H_left, temp;
  __m256i c_up, c_left, b_up, b_left, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0;

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
      s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 8 * (i - 1)));
      temp_index = _mm256_mullo_epi32 (s1, v128);
      E = _mm256_set1_epi32 (-10000);
      H_diag = _mm256_setzero_si256 ();
      H = _mm256_setzero_si256 ();
      bt_flag[i * y] = 0;
      for (int j = 1; j < y; j++)
	{
	  s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 8 * (j - 1)));
	  index = _mm256_add_epi32 (temp_index, s2);
	  score = _mm256_i32gather_epi32(subs_matrix, index, 4);
	  H_left = _mm256_load_si256 ((__m256i *) (aH + 8 * j));
	  diag = _mm256_add_epi32 (H_diag, score);

	  E_sub = _mm256_sub_epi32 (E, vextend);	//for now, E is E_up
	  E = _mm256_sub_epi32 (H, vopen);		//for now, H is H_up
	  E = _mm256_max_epi32 (E, E_sub);

	  F_sub = _mm256_load_si256 ((__m256i *) (aF + 8 * j));
	  F_sub = _mm256_sub_epi32 (F_sub, vextend);
	  F = _mm256_sub_epi32 (H_left, vopen);
	  F = _mm256_max_epi32 (F, F_sub);
	  _mm256_store_si256 ((__m256i *) (aF + 8 * j), F);

	  H = _mm256_max_epi32 (E, F);
	  H = _mm256_max_epi32 (H, diag);
	  H = _mm256_max_epi32 (H, vzero);
	  _mm256_store_si256 ((__m256i *) (aH + 8 * j), H);

	  //logic tests
	  h_gt_0 = _mm256_movemask_ps (
	      _mm256_castsi256_ps (_mm256_cmpgt_epi32 (H, vzero)));

	  H_eq_E = _mm256_cmpeq_epi32 (H, E);
	  h_eq_e = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_E));

	  H_eq_F = _mm256_cmpeq_epi32 (H, F);
	  h_eq_f = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_F));

	  //********FLAGS********

	  c_up = _mm256_cmpeq_epi32 (E, E_sub); // E[i,j] == E[i,j-1]-gap_extent ?
	  flag = _mm256_movemask_ps (_mm256_castsi256_ps (c_up)); // & h_eq_e;		//c_up

	  c_left = _mm256_cmpeq_epi32 (F, F_sub); // F[i,j] == F[i-1,j]-gap_extent ?
	  flag <<= 8;
	  flag |= _mm256_movemask_ps (_mm256_castsi256_ps (c_left)); // & h_eq_f;		//c_left

	  H_eq_diag = _mm256_cmpeq_epi32 (H, diag);
	  h_eq_d = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_diag));
	  flag <<= 8;
	  flag |= (h_eq_e | h_eq_d) & h_gt_0;		//b_up

	  //b_left
	  flag <<= 8;
	  flag = flag | (((h_eq_f & ~h_eq_e) | h_eq_d) & h_gt_0);

	  bt_flag[y * i + j] = flag;

	  temp = _mm256_cmpgt_epi32 (H, max);
	  imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi32 (i), temp);
	  jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi32 (j), temp);
	  max = _mm256_max_epi32 (H, max);

	  H_diag = H_left;
	}
    }
  _mm256_store_si256 ((__m256i *) max_score, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
  free (aF);
  free (aH);
  return bt_flag;
}


/**
 * Core function
 *
 * @param overflow: alignments that need to be rerun with wider registers
 *
 */
int*
avx2_fill_table_1_to_32_i8 (int* seqs1, int* seqs2, int x, int y,
			   int* subs_matrix, int gap_open, int gap_extend,
			   int max_score[32], int ipos[32], int jpos[32],
			   int overflow[32])
{
  __m256i s1, s2, temp_index, index;
  __m256i v128 = _mm256_set1_epi32 (128);
  __m256i vopen = _mm256_set1_epi32 (gap_open);
  __m256i vextend = _mm256_set1_epi32 (gap_extend);
  __m256i vzero = _mm256_setzero_si256 ();
  __m256i vepsilon = _mm256_set1_epi32 (EPSILON);
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256i max = _mm256_setzero_si256 ();
  __m256i E, E_sub;
  __m256i F, F_sub;
  __m256i diag, score, H, H_diag, H_left, temp;
  __m256i c_up, c_left, b_up, b_left, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0;

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
      s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 8 * (i - 1)));
      temp_index = _mm256_mullo_epi32 (s1, v128);
      E = _mm256_set1_epi32 (-10000);
      H_diag = _mm256_setzero_si256 ();
      H = _mm256_setzero_si256 ();
      bt_flag[i * y] = 0;
      for (int j = 1; j < y; j++)
	{
	  s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 8 * (j - 1)));
	  index = _mm256_add_epi32 (temp_index, s2);
	  score = _mm256_i32gather_epi32(subs_matrix, index, 4);
	  H_left = _mm256_load_si256 ((__m256i *) (aH + 8 * j));
	  diag = _mm256_add_epi32 (H_diag, score);

	  E_sub = _mm256_sub_epi32 (E, vextend);	//for now, E is E_up
	  E = _mm256_sub_epi32 (H, vopen);		//for now, H is H_up
	  E = _mm256_max_epi32 (E, E_sub);

	  F_sub = _mm256_load_si256 ((__m256i *) (aF + 8 * j));
	  F_sub = _mm256_sub_epi32 (F_sub, vextend);
	  F = _mm256_sub_epi32 (H_left, vopen);
	  F = _mm256_max_epi32 (F, F_sub);
	  _mm256_store_si256 ((__m256i *) (aF + 8 * j), F);

	  H = _mm256_max_epi32 (E, F);
	  H = _mm256_max_epi32 (H, diag);
	  H = _mm256_max_epi32 (H, vzero);
	  _mm256_store_si256 ((__m256i *) (aH + 8 * j), H);

	  //logic tests
	  h_gt_0 = _mm256_movemask_ps (
	      _mm256_castsi256_ps (_mm256_cmpgt_epi32 (H, vzero)));

	  H_eq_E = _mm256_cmpeq_epi32 (H, E);
	  h_eq_e = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_E));

	  H_eq_F = _mm256_cmpeq_epi32 (H, F);
	  h_eq_f = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_F));

	  //********FLAGS********

	  c_up = _mm256_cmpeq_epi32 (E, E_sub); // E[i,j] == E[i,j-1]-gap_extent ?
	  flag = _mm256_movemask_ps (_mm256_castsi256_ps (c_up)); // & h_eq_e;		//c_up

	  c_left = _mm256_cmpeq_epi32 (F, F_sub); // F[i,j] == F[i-1,j]-gap_extent ?
	  flag <<= 8;
	  flag |= _mm256_movemask_ps (_mm256_castsi256_ps (c_left)); // & h_eq_f;		//c_left

	  H_eq_diag = _mm256_cmpeq_epi32 (H, diag);
	  h_eq_d = _mm256_movemask_ps (_mm256_castsi256_ps (H_eq_diag));
	  flag <<= 8;
	  flag |= (h_eq_e | h_eq_d) & h_gt_0;		//b_up

	  //b_left
	  flag <<= 8;
	  flag = flag | (((h_eq_f & ~h_eq_e) | h_eq_d) & h_gt_0);

	  bt_flag[y * i + j] = flag;

	  temp = _mm256_cmpgt_epi32 (H, max);
	  imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi32 (i), temp);
	  jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi32 (j), temp);
	  max = _mm256_max_epi32 (H, max);

	  H_diag = H_left;
	}
    }
  _mm256_store_si256 ((__m256i *) max_score, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
  free (aF);
  free (aH);
  return bt_flag;
}




/*
 *Core function.
 *
 */
alignment_i32*
avx2_sw_i32_with_matrix (char* seqs1_id[8], char* seqs2_id[8], char* seqs1[8],
			 char* seqs2[8], int* subs_matrix, int gap_open,
			 int gap_extend, int dup_strings)
{
  int max_i;
  int max_j;

  int lens_i[8];
  int lens_j[8];

  for (int i = 0; i < 8; i++)
    lens_i[i] = strlen (seqs1[i]);

  for (int j = 0; j < 8; j++)
    lens_j[j] = strlen (seqs2[j]);

  int* packed_seqs1 = avx_pack_8_seqs_with_len (seqs1, lens_i, &max_i);
  int* packed_seqs2 = avx_pack_8_seqs_with_len (seqs2, lens_j, &max_j);

  int max_score[8] __attribute__((aligned(32)));
  int pos_i[8] __attribute__((aligned(32)));
  int pos_j[8] __attribute__((aligned(32)));

  int* flags = avx2_fill_table_8_to_8_i32 (packed_seqs1, packed_seqs2,
					   max_i + 1, max_j + 1, subs_matrix,
					   gap_open, gap_extend, max_score,
					   pos_i, pos_j);

  alignment_i32* alignments = (alignment_i32*) malloc (
      8 * sizeof(alignment_i32));

  int x0, y0;

  for (int i = 0; i < 8; i++)
    {
      //assert(max_i + max_j + 1 > 0);
      char* m = (char*) calloc ((max_i + max_j + 1), sizeof(char)); //worst-case length
      char* n = (char*) calloc ((max_i + max_j + 1), sizeof(char));

      //assert(m != NULL);
      //assert(n != NULL);

      sw_backtrack (i, flags, seqs1[i], seqs2[i], max_i + 1, max_j + 1, m, n,
		    pos_i[i], pos_j[i], &x0, &y0);
      //assert(strlen (m) == strlen (n));

      int new_len = strlen (m) + 1;

      //adjusting to actual sizes
      m = (char*) realloc (m, new_len * sizeof(char));
      n = (char*) realloc (n, new_len * sizeof(char));

      if (dup_strings)
	{
	  alignments[i].seq1_id = strdup (seqs1_id[i]);
	  alignments[i].seq2_id = strdup (seqs2_id[i]);
	  alignments[i].seq1 = strdup (seqs1[i]);
	  alignments[i].seq2 = strdup (seqs2[i]);
	}
      else
	{
	  alignments[i].seq1_id = seqs1_id[i];
	  alignments[i].seq2_id = seqs2_id[i];
	  alignments[i].seq1 = seqs1[i];
	  alignments[i].seq2 = seqs2[i];
	}
      alignments[i].seq1_len = lens_i[i];
      alignments[i].seq2_len = lens_j[i];
      alignments[i].aln1 = m;
      alignments[i].aln2 = n;
      alignments[i].aln_len = strlen (m);
      alignments[i].aln_x0 = x0;
      alignments[i].aln_y0 = y0;
      alignments[i].aln_x = pos_i[i];
      alignments[i].aln_y = pos_j[i];
      alignments[i].score = max_score[i];
    }

  free (packed_seqs1);
  free (packed_seqs2);
  free (flags);
  return alignments;
}
