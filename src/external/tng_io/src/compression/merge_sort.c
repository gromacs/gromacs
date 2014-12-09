/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/merge_sort.h"

static void ms_inner(void *base, const size_t size,
             const size_t start, const size_t end,
             int (*compar)(const void *v1,const void *v2,const void *private),
             const void *private,
             char *workarray)
{
  size_t middle;
  if ((end-start)>1)
    {
      char *cbase=(char *)base;
      middle=start+(end-start)/2;
#if 0
      printf("For start %d end %d obtained new middle: %d\n",start,end,middle);
#endif
      ms_inner(base,size,
           start,middle,
           compar,private,workarray);
      ms_inner(base,size,
           middle,end,
           compar,private,workarray);
#if 0
      printf("For start %d end %d Before merge: Comparing element %d with %d\n",start,end,middle-1,middle);
#endif
      if (compar(cbase+(middle-1)*size,cbase+middle*size,private)>0)
        {
          /* Merge to work array. */
          size_t i, n=end-start;
          size_t ileft=start;
          size_t iright=middle;
          for (i=0; i<n; i++)
            {
              if (ileft==middle)
                {
                  memcpy(workarray+i*size,cbase+iright*size,size);
                  iright++;
                }
              else if (iright==end)
                {
                  memcpy(workarray+i*size,cbase+ileft*size,size);
                  ileft++;
                }
              else
                {
        #if 0
                  printf("For start %d end %d In merge: Comparing element %d with %d\n",start,end,ileft,iright);
        #endif
                  if (compar(cbase+ileft*size,cbase+iright*size,private)>0)
                    {
                      memcpy(workarray+i*size,cbase+iright*size,size);
                      iright++;
                    }
                  else
                    {
                      memcpy(workarray+i*size,cbase+ileft*size,size);
                      ileft++;
                    }
                }
            }
          /* Copy result back. */
          memcpy(cbase+start*size,workarray,(end-start)*size);
        }
    }
}


void Ptngc_merge_sort(void *base, const size_t nmemb, const size_t size,
        int (*compar)(const void *v1,const void *v2,const void *private),
        void *private)
{
  char *warr=warnmalloc(nmemb*size);
  ms_inner(base,size,0,nmemb,compar,private,warr);
  free(warr);
}


#ifdef TEST

static int compint(const void *v1, const void *v2,const void *private)
{
  const int *i1=(const int *)v1;
  const int *i2=(const int *)v2;
  if (*i1<*i2)
    return -1;
  else if (*i1>*i2)
    return 1;
  else
    return 0;
}

static int qcompint(const void *v1, const void *v2)
{
  return compint(v1,v2,NULL);
}

#define N 1000000
int main()
{
  int *arr=warnmalloc(N*sizeof *arr);
  int i;
  for (i=0; i<N; i++)
    scanf("%d",arr+i);
#if 1
  merge_sort(arr,N,sizeof *arr,compint,NULL);
#else
  qsort(arr,N,sizeof *arr,qcompint);
#endif
  for (i=0; i<N; i++)
    printf("%d %d\n",i,arr[i]);
  return 0;
}
#endif
