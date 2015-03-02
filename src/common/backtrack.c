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

#include "backtrack.h"

void
inplace_reverse (char* str)
{
  int len = strlen (str);
  char temp;
  for (int i = 0; i < len / 2; i++)
    {
      temp = str[i];
      str[i] = str[len - i - 1];
      str[len - i - 1] = temp;
    }
}

/*
 *
 */
void
sw_backtrack (int index, int* flags, const char* a, const char* b, int w, int h,
	      char* aln1, char* aln2, int x, int y, int* x0, int* y0)
{
  int d_mask = 0x00000101 << index;	//diagonal
  int bl_mask = 0x00000001 << index;	//begin left
  int bu_mask = 0x00000100 << index;	//begin up
  int cl_mask = 0x00010000 << index;	//continue left
  int cu_mask = 0x01000000 << index;	//continue up

  int c = 0;

  while ((flags[x * h + y] & d_mask) && x > 1 && y > 1)
    {
      if ((flags[x * h + y] & d_mask) == d_mask)
	{
	  aln1[c] = a[--x];
	  aln2[c] = b[--y];
	  c++;
	}
      else if (flags[x * h + y] & bl_mask)
	{
	  do
	    {
	      aln1[c] = a[x - 1];
	      aln2[c] = '-';
	      c++;
	    }
	  while (flags[(x--) * h + y] & cl_mask);
	}
      else
	{
	  do
	    {
	      aln1[c] = '-';
	      aln2[c] = b[y - 1];
	      c++;
	    }
	  while (flags[x * h + (y--)] & cu_mask);
	}
    }
  *x0 = (int) (x + 1);
  *y0 = (int) (y + 1);
  inplace_reverse (aln1);
  inplace_reverse (aln2);
}
