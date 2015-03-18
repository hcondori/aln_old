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



#ifndef SEQUENCE_H_
#define SEQUENCE_H_

typedef struct
{
  char* id;		//standard null-terminated C string
  char* sequence;	//for performance, this might not be a standard null-terminated C string
  int length;		//and this is its length
} Sequence;

Sequence*
Sequence_new (char* id, char* sequence, int length);

void
Sequence_free (Sequence* sequence);

#endif /* SEQUENCE_H_ */
