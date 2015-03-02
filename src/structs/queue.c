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


#include <stdlib.h>

#include "structs/queue.h"

queue_node*
new_queue_node (void* data)
{
  queue_node* qn = (queue_node*) malloc (sizeof(queue_node));
  qn->data = data;
  qn->next = NULL;
  return qn;
}

queue*
new_queue ()
{
  queue* q = (queue*) calloc (1, sizeof(queue));
  return q;
}

void
enqueue (queue* q, void* data)
{
  queue_node* qn = new_queue_node (data);
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
dequeue (queue* q)
{
  if (q->front == NULL)	//empty
    {
      return NULL;
    }
  else
    {
      void* data = q->front->data;
      queue_node* qn = q->front;
      q->front = qn->next;
      if (q->back == qn)
	q->back = NULL;
      free (qn);
      q->count--;
      return data;
    }
}
