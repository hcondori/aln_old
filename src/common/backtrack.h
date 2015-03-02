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


#ifndef BACKTRACK_H_
#define BACKTRACK_H_

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define U_FEPS 1.192e-6F
#define EPS 0.001f

#define LEFT 1
#define UP 2
#define DIAGONAL 3

void
get_alignment (int* bt_flag, int i, int j, int x, int y, int n, char* path,
	       int* path_len);

void
aln2char (char* seq1, char* seq2, char* path, char* out1, char* out2, int i,
	  int j, int path_len);

void
inplace_reverse (char* str);

void
sw_backtrack (int index, int* flags, const char* a, const char* b, int len_a,
	       int len_b, char* aln1, char* aln2, int x, int y, int* x0,
	       int* y0);

#endif /* BACKTRACK_H_ */
