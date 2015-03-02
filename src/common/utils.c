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


#include "utils.h"

inline int*
avx_pack_8_seqs (char** char_seqs, int *max)
{
  *max = 0;
  int temp;
  int lens[8] =
    { 0 };
  for (int i = 0; i < 8; i++)
    {
      lens[i] = strlen (char_seqs[i]);
      if (lens[i] > *max)
	*max = lens[i];
    }
  int* seqs = (int*) aligned_alloc (32, 8 * sizeof(int) * (*max));
  for (int i = 0; i < *max; i++)
    {
      for (int j = 0; j < 8; j++)
	{
	  if (i < lens[j])
	    seqs[8 * i + j] = (int) char_seqs[j][i];
	  else
	    seqs[8 * i + j] = 0;
	}
    }
  return seqs;
}

int*
avx_pack_8_seqs_with_len (char* char_seqs[8], int lens[8], int *max)
{
  *max = 0;
  for (int i = 0; i < 8; i++)
    {
      if (lens[i] > *max)
	*max = lens[i];
    }
  int* seqs = (int*) aligned_alloc (32, 8 * sizeof(int) * (*max));
  for (int i = 0; i < *max; i++)
    {
      for (int j = 0; j < 8; j++)
	{
	  if (i < lens[j])
	    seqs[8 * i + j] = (int) char_seqs[j][i];
	  else
	    seqs[8 * i + j] = 0;
	}
    }
  return seqs;
}

