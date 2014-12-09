/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg and Magnus Lundborg
 * Copyright (c) 2010, 2013-2014 The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/bwt.h"

#if 0
#define SHOWIT
#endif
#if 0
#define SHOWIT2
#endif

static int compare_index(int i1,int i2, const int nvals, unsigned int *vals, unsigned int *nrepeat)
{
  int i,j;
  for (i=0; i<nvals; i++)
    {
      /* If we have repeating patterns, we might be able to start the
         comparison later in the string. */
      /* Do we have a repeating pattern? If so are
         the repeating patterns the same length? */
      int repeat1=(int)(nrepeat[i1]>>8);
      int k1=(int)(nrepeat[i1]&0xFFU);
      int repeat2=(int)(nrepeat[i2]>>8);
      int k2=(int)(nrepeat[i2]&0xFFU);

      if ((repeat1>1) && (repeat2>1) && (k1==k2))
        {
          int maxskip=0;
          /* Yes. Compare the repeating patterns. */
          for (j=0; j<k1; j++)
            {
              unsigned int v1=vals[(i1+j)%nvals];
              unsigned int v2=vals[(i2+j)%nvals];
              if (v1<v2)
                return -1;
              else if (v1>v2)
                return 1;
            }
          /* The repeating patterns are equal. Skip as far as we can
             before continuing. */
          maxskip=repeat1;
          if (repeat2<repeat1)
            maxskip=repeat2;
          i1=(i1+maxskip)%nvals;
          i2=(i2+maxskip)%nvals;
          i+=maxskip-1;
        }
      else
        {
          if (vals[i1]<vals[i2])
            return -1;
          else if (vals[i1]>vals[i2])
            return 1;
          i1++;
          if (i1>=nvals)
            i1=0;
          i2++;
          if (i2>=nvals)
            i2=0;
        }
    }
  return 0;
}

void Ptngc_bwt_merge_sort_inner(int *indices, const int nvals,unsigned int *vals,
                                const int start, const int end,
                                unsigned int *nrepeat,
                                int *workarray)
{
  int middle;
  if ((end-start)>1)
    {
      middle=start+(end-start)/2;
#if 0
      printf("For start %d end %d obtained new middle: %d\n",start,end,middle);
#endif
      Ptngc_bwt_merge_sort_inner(indices,nvals,vals,
                           start,middle,
                           nrepeat,
                           workarray);
      Ptngc_bwt_merge_sort_inner(indices,nvals,vals,
                           middle,end,
                           nrepeat,
                           workarray);
#if 0
      printf("For start %d end %d Before merge: Comparing element %d with %d\n",start,end,middle-1,middle);
#endif
      if (compare_index(indices[middle-1],indices[middle],nvals,vals,nrepeat)>0)
        {
          /* Merge to work array. */
          int i, n=end-start;
          int ileft=start;
          int iright=middle;
          for (i=0; i<n; i++)
            {
              if (ileft==middle)
                {
                  workarray[i]=indices[iright];
                  iright++;
                }
              else if (iright==end)
                {
                  workarray[i]=indices[ileft];
                  ileft++;
                }
              else
                {
#if 0
                  printf("For start %d end %d In merge: Comparing element %d with %d\n",start,end,ileft,iright);
#endif
                  if (compare_index(indices[ileft],indices[iright],nvals,vals,nrepeat)>0)
                    {
                      workarray[i]=indices[iright];
                      iright++;
                    }
                  else
                    {
                      workarray[i]=indices[ileft];
                      ileft++;
                    }
                }
            }
          /* Copy result back. */
          memcpy(indices+start,workarray,(end-start)*sizeof(int));
        }
    }
}

/* Burrows-Wheeler transform. */
void Ptngc_comp_to_bwt(unsigned int *vals, const int nvals,
                       unsigned int *output, int *index)
{
  int i;
  int *indices=warnmalloc(2*nvals*sizeof *indices);
  unsigned int *nrepeat=warnmalloc(nvals*sizeof *nrepeat);
  int *warr=indices+nvals;

  if (nvals>0xFFFFFF)
    {
      fprintf(stderr,"BWT cannot pack more than %d values.\n",0xFFFFFF);
      exit(1);
    }

  /* Also note that repeat pattern k (kmax) cannot be larger than 255. */
#if 0
  printf("Size of arrays is %.2f M\n",4*nvals*sizeof *indices/1024./1024.);
#endif
  for (i=0; i<nvals; i++)
    indices[i]=i;
  /* Find the length of the initial repeating pattern for the strings. */
  /* First mark that the index does not have a found repeating string. */

  memset(nrepeat, 0U, sizeof(unsigned int) * nvals);

#ifdef SHOWIT
  printf("nvals is %d\n",nvals);
#endif
  for (i=0; i<nvals; i++)
    {
      /* If we have not already found a repeating string we must find
         it. */
      if (!nrepeat[i])
        {
          int maxrepeat=nvals*2;
          int j,k,m;
          int good_j=-1, good_k=0;
          int kmax=16;
          /* Track repeating patterns.
             k=1 corresponds to AAAAA...
             k=2 corresponds to ABABAB...
             k=3 corresponds to ABCABCABCABC...
             k=4 corresponds to ABCDABCDABCD...
             etc. */
          for (k=kmax; k>=1; k--)
            {
            try_next_k:
              if (k>=1)
                {
#ifdef SHOWIT
                  printf("Trying k=%d at i=%d\n",k,i);
#endif
                  for (j=k; j<maxrepeat; j+=k)
                    {
                      int is_equal=1;
#ifdef SHOWIT
                      printf("Trying j=%d at i=%d for k %d\n",j,i,k);
#endif
                      for (m=0; m<k; m++)
                        if (vals[(i+m)%nvals]!=vals[(i+j+m)%nvals])
                          {
                            is_equal=0;
                            break;
                          }
                      if (is_equal)
                        {
                          int new_j=j+k;
                          if (new_j>maxrepeat)
                            new_j=j;
                          if ((new_j>good_j) || ((new_j==good_j) && (k<good_k)))
                            {
                              good_j=new_j; /* We have found that
                                             the strings repeat
                                             for this length... */
                              good_k=k;        /* ...and with this
                                                  length of the
                                                  repeating
                                                  pattern. */
#ifdef SHOWIT
                              printf("Best j and k is now %d and %d\n",good_j,good_k);
#endif
                            }
                        }
                      else
                        {
                          /* We know that it is no point in trying
                             with more than m */
                          if (j==0)
                            {
                              k=m;
#ifdef SHOWIT
                              printf("Setting new k to m: %d\n",k);
#endif
                            }
                          else
                            k--;
#ifdef SHOWIT
                          printf("Trying next k\n");
#endif
                          goto try_next_k;
                        }
                    }
                }
            }
          /* From good_j and good_k we know the repeat for a large
             number of strings. The very last repeat length should not
             be assigned, since it can be much longer if a new test is
             done. */
          for (m=0; (m+good_k<good_j) && (i+m<nvals); m+=good_k)
            {
              int repeat=good_j-m;
              if (repeat>nvals)
                repeat=nvals;
              nrepeat[i+m]=((unsigned int) (good_k)) | (((unsigned int) (repeat))<<8);
            }
          /* If no repetition was found for this value signal that here. */
          if (!nrepeat[i])
            nrepeat[i+m]=257U; /* This is 1<<8 | 1 */
        }
    }
#ifdef SHOWIT
  for (i=0; i<nvals; i++)
    if ((nrepeat[i]>>8)!=1)
      printf("String repeats at %d: %d %d\n",i,nrepeat[i]>>8,nrepeat[i]&0xFFU);
#endif

  /* Sort cyclic shift matrix. */
  Ptngc_bwt_merge_sort_inner(indices,nvals,vals,0,nvals,nrepeat,warr);

#if 0
  /* Test that it really is sorted. */
  for (i=0; i<nvals-1; i++)
    if (compare_index(indices[i],indices[i+1],nvals,vals,nrepeat)!=-1)
      fprintf(stderr,"SORTING NOT COMPLETED AT %d. Index %d against %d\n",i,indices[i],indices[i+1]);

#endif


#ifdef SHOWIT2
  for (i=0; i<nvals; i++)
    {
      int j;
      for (j=0; j<nvals; j++)
        printf("%c",vals[(indices[i]+j)%nvals]);
      printf("\n");
    }
#endif
  /* Which one is the original string? */
  for (i=0; i<nvals; i++)
    if (indices[i]==0)
      break;
  *index=i;
  /* Form output. */
  for (i=0; i<nvals; i++)
    {
      int lastchar=indices[i]-1;
      if (lastchar<0)
        lastchar=nvals-1;
      output[i]=vals[lastchar];
    }
  free(nrepeat);
  free(indices);
}

/* Burrows-Wheeler inverse transform. */
void Ptngc_comp_from_bwt(unsigned int *input, const int nvals, int index,
                         unsigned int *vals)
{
  /* Straightforward from the Burrows-Wheeler paper (page 13). */
  int i;
  unsigned int *c=warnmalloc(0x10000*sizeof *c);
  unsigned int *p=warnmalloc(nvals*sizeof *p);
  unsigned int sum=0;

  memset(c, 0, sizeof(unsigned int) * 0x10000);

  for (i=0; i<nvals; i++)
    {
      p[i]=c[input[i]];
      c[input[i]]++;
    }
  for (i=0; i<0x10000; i++)
    {
      sum+=c[i];
      c[i]=sum-c[i];
    }
  for (i=nvals-1; i>=0; i--)
    {
      vals[i]=input[index];
      index=p[index]+c[input[index]];
    }
  free(p);
  free(c);
}

