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

#include "avx_sw.h"

#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "common/backtrack.h"
#include "common/utils.h"

/*
 * Simple Smith-Waterman implementation for compatibility with old processors.
 */
alignment_f32*
sw_1_to_1_f32_mith_matrix (char* seq1, char* seq2, float* subs_matrix,
			   float gap_open, float gap_extend)
{
  int w = strlen (seq1) + 1;
  int h = strlen (seq2) + 1;

  float diagonal, E, E_sub, F, F_sub, H, H_left, H_diag;
  float* aF = (float*) malloc (h * sizeof(float));
  for (int i = 0; i < h; i++)
    {
      aF[i] = -FLT_MAX;
    }
  float* aH = (float*) calloc (h, sizeof(float));

  int flag;
  int *flags = (int*) malloc (h * w * sizeof(int));

  float max = 0;
  int x = 0, y = 0;

  for (int i = 1; i < w; i++)
    {
      E = -FLT_MAX;
      H = 0;
      H_diag = 0;
      for (int j = 1; j < h; j++)
	{
	  H_left = aH[j];
	  F_sub = aF[j];
	  diagonal = H_diag + subs_matrix[seq1[i - 1] * 128 + seq2[j - 1]];

	  E_sub = E - gap_extend;
	  E = H - gap_open;
	  E = MAX(E, E_sub);

	  F_sub = F_sub - gap_extend;
	  F = H_left - gap_open;
	  F = MAX(F, F_sub);
	  aF[j] = F;

	  H = MAX(E, F);
	  H = MAX(H, diagonal);
	  H = MAX(H, 0);
	  aH[j] = H;
	  H_diag = H_left;

	  if (H > max)
	    {
	      max = H;
	      x = i;
	      y = j;
	    }

	  //flags
	  flag = (E == E_sub) ? 1 : 0; //c_up
	  flag <<= 8;
	  flag += (F == F_sub) ? 1 : 0; //c_left
	  flag <<= 8;
	  flag += ((H == E || H == diagonal) && H > 0) ? 1 : 0; //b_up
	  flag <<= 8;
	  flag += (((H == F && H != E) || H == diagonal) && H > 0) ? 1 : 0; //b_left
	  flags[i * w + j] = flag;
	}
    }
  char* aln1 = (char*) calloc (w + h - 1, sizeof(char));
  char* aln2 = (char*) calloc (w + h - 1, sizeof(char));
  int x0, y0;
  sw_backtrack (0, flags, seq1, seq2, w, h, aln1, aln2, x, y, &x0, &y0);
  alignment_f32* alignment = (alignment_f32*) malloc (sizeof(alignment_f32));

  alignment->aln1 = seq1;
  alignment->aln2 = seq2;
  alignment->seq1_len = w - 1;
  alignment->seq2_len = h - 1;
  alignment->aln1 = aln1;
  alignment->aln2 = aln2;
  alignment->aln_len = strlen (aln1);
  alignment->aln_x0 = x0;
  alignment->aln_y0 = y0;
  alignment->aln_x = x;
  alignment->aln_y = y;
  alignment->score = max;
  return alignment;
}

void
avx2_fill_table_1_to_m_f32 (int* bt_flag, char* query, char** refs,
			    int ref_count, int x, int y, float* subs_matrix,
			    float gap_open, float gap_extend, float* max_score,
			    int* ipos, int* jpos)
{

}

/*
 *Core function.
 *
 */
void
avx_fill_table_8_to_8_f32 (int* bt_flag, int* seqs1, int* seqs2, int x, int y,
			   float* subs_matrix, float gap_open, float gap_extend,
			   float* max_score, int* ipos, int* jpos)
{
  int mask = 0x7FFFFFFF;
  __m256 vmask = _mm256_castsi256_ps (_mm256_set1_epi32 (mask));
  __m256i s1, s2, temp_index, index;
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
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0, test;

  float buffer[8] __attribute__((aligned(32)));
  int buffer_i[8] __attribute__((aligned(32)));
  int buffer_j[8] __attribute__((aligned(32)));

  float* aF = (float*) aligned_alloc (32, 8 * y * sizeof(float));
  float* aH = (float*) aligned_alloc (32, 8 * y * sizeof(float));
  //if we are sure that sequences will be small, we keep them in the stack
  //float aF[8 * y] __attribute__((aligned(32)));
  //float aH[8 * y] __attribute__((aligned(32)));

  memset (aH, 0, 8 * y * sizeof(float));
  memset (ipos, 0, 8 * sizeof(int));
  memset (jpos, 0, 8 * sizeof(int));

  for (int i = 0; i < 8 * y; i++)
    {
      aF[i] = -FLT_MAX;
    }

  for (int i = 1; i < x; i++)
    {

      E = _mm256_set1_ps (-FLT_MAX);
      H_diag = _mm256_setzero_ps ();
      H = _mm256_setzero_ps ();
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
	  h_gt_0 = _mm256_movemask_ps (_mm256_cmp_ps(H, vzero, _CMP_GT_OQ));
	  H_eq_E = _mm256_cmp_ps(H, E, _CMP_EQ_OQ);
	  h_eq_e = _mm256_movemask_ps (H_eq_E);
	  H_eq_F = _mm256_cmp_ps(H, F, _CMP_EQ_OQ);
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

	  temp = _mm256_cmp_ps(H, max, _CMP_GT_OQ);
	  //test = _mm256_movemask_ps (temp);
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
	  H_diag = H_left;
	}
    }
  _mm256_store_ps (max_score, max);
  free (aF);
  free (aH);
}

alignment_f32*
avx_sw_f32_with_matrix (char* seqs1_id[8], char* seqs2_id[8], char** seqs1,
			char** seqs2, float* subs_matrix, float gap_open,
			float gap_extend, int dup_strings)
{
  int max_i = 0;
  int max_j = 0;

  int lens_i[8];
  int lens_j[8];

  for (int i = 0; i < 8; i++)
    lens_i[i] = strlen (seqs1[i]);

  for (int j = 0; j < 8; j++)
    lens_j[j] = strlen (seqs2[j]);

  int* packed_seqs1 = avx_pack_8_seqs_with_len (seqs1, lens_i, &max_i);
  int* packed_seqs2 = avx_pack_8_seqs_with_len (seqs2, lens_j, &max_j);

  int* flags = (int*) malloc (sizeof(int) * (max_i + 1) * (max_j + 1));

  float max_score[8] __attribute__((aligned(32)));
  int pos_i[8] __attribute__((aligned(32)));
  int pos_j[8] __attribute__((aligned(32)));

  avx_fill_table_8_to_8_f32 (flags, packed_seqs1, packed_seqs2, max_i + 1,
			     max_j + 1, subs_matrix, gap_open, gap_extend,
			     max_score, pos_i, pos_j);

  alignment_f32* alignments = (alignment_f32*) malloc (
      8 * sizeof(alignment_f32));

  int x0, y0;

  for (int i = 0; i < 8; i++)
    {
      char* m = (char*) calloc ((max_i + max_j + 1), sizeof(char));
      char* n = (char*) calloc ((max_i + max_j + 1), sizeof(char));

      sw_backtrack (i, flags, seqs1[i], seqs2[i], max_i + 1, max_j + 1, m, n,
		    pos_i[i], pos_j[i], &x0, &y0);

      int new_len = strlen (m) + 1;

      //adjusting to actual sizes
      //m = (char*) realloc (m, new_len * sizeof(char));
      //n = (char*) realloc (n, new_len * sizeof(char));

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
