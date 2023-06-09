#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../Mpi/multidev.h"

#include "memory_wrapper.h"

struct memory_allocated_t *memory_allocated_base=NULL;

size_t memory_used=0;
size_t max_memory_used=0;

int posix_memalign_wrapper(void **memptr,size_t alignment,size_t size,const char* varname)
{
  int res=posix_memalign(memptr,alignment,size);
  
  // increment the used size, and take the maximum
  memory_used+=size;
  if(memory_used>max_memory_used) max_memory_used=memory_used;
  
  // creates the page and insert on top of the stack
  struct memory_allocated_t *all=(struct memory_allocated_t *)malloc(sizeof(struct memory_allocated_t));
  all->ptr=*memptr;
  all->varname=varname;
  all->size=size;
  all->next=memory_allocated_base;
  memory_allocated_base=all;
  
  return res;
}

void free_wrapper(void *memptr)
{
  // creates the page and insert on top of the stack
  struct memory_allocated_t *all=memory_allocated_base,*prev=NULL;
  
  while(all!=NULL && all->ptr!=memptr)
    {
      prev=all;
      all=all->next;
    };
  
  if(all==NULL)
    {
      fprintf(stderr,"\n  MPI%02d - Failed to find pointer %p in the list while freeing \n\n\n",devinfo.myrank,memptr);
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  
  memory_used-=all->size;
  
  if(prev!=NULL)
    prev->next=all->next;
  
  free(all->ptr);
  free(all);
};
