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

#ifndef QUEUE_H_
#define QUEUE_H_

typedef struct
{
  void* data;
  void* next;
} Queue_node;

typedef struct
{
  Queue_node* front;
  Queue_node* back;
  int count;
} Queue;

Queue_node*
new_queue_node (void* data);

Queue*
Queue_new ();

void
enqueue (Queue* q, void* data);

void*
dequeue (Queue* q);

int
Queue_is_empty (Queue* q);

#endif /* QUEUE_H_ */
