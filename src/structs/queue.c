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

#include <stdlib.h>

#include "structs/queue.h"

Queue_node*
new_queue_node (void* data)
{
  Queue_node* qn = (Queue_node*) malloc (sizeof(Queue_node));
  qn->data = data;
  qn->next = NULL;
  return qn;
}

Queue*
Queue_new ()
{
  Queue* q = (Queue*) calloc (1, sizeof(Queue));
  return q;
}

void
enqueue (Queue* q, void* data)
{
  Queue_node* qn = new_queue_node (data);
  if (q->back == NULL)	//empty
    {
      q->front = qn;
      q->back = qn;
    }
  else
    {
      q->back->next = qn;
      q->back = qn;
    }
  q->count++;
}

void*
dequeue (Queue* q)
{
  if (q->front == NULL)	//empty
    {
      return NULL;
    }
  else
    {
      void* data = q->front->data;
      Queue_node* qn = q->front;
      q->front = qn->next;
      if (q->back == qn)
	q->back = NULL;
      free (qn);
      q->count--;
      return data;
    }
}

int
Queue_is_empty (Queue* q)
{
  return q->count == 0;
}

void
Queue_free (Queue* q)
{
  if (q != NULL)
    while (q->count > 0)
      dequeue (q);
  free (q);
}
