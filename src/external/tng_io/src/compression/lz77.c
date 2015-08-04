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
#include "../../include/compression/bwt.h"
#include "../../include/compression/lz77.h"

/* This is a simple Lempel-Ziv-77 compressor. It has not been set up
   to remove every possible repeating pattern, but it might be better
   than simple RLE.

   Lempel-Ziv 77 with separate outputs for length, data, and offsets.
 */

#if 0
/* Sort the strings (similar to BWT) to find similar patterns in the
   input data.  The output is an array with two values for each
   input value. The sorted index value is the first in each doublet.
   The second value is the "inverse" of the first value and can be
   used to find the locations of similar strings. */
static void sort_strings(unsigned int *vals, int nvals,
                         unsigned int *output)
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
  for (i=0; i<nvals; i++)
    indices[i]=i;
  /* Find the length of the initial repeating pattern for the strings. */
  /* First mark that the index does not have a found repeating string. */
  for (i=0; i<nvals; i++)
    nrepeat[i]=0U;
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
                  for (j=k; j<maxrepeat; j+=k)
                    {
                      int is_equal=1;
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
                            }
                        }
                      else
                        {
                          /* We know that it is no point in trying
                             with more than m */
                          if (j==0)
                            {
                              k=m;
                            }
                          else
                            k--;
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

  /* Sort cyclic shift matrix. */
  Ptngc_bwt_merge_sort_inner(indices,nvals,vals,0,nvals,nrepeat,warr);

  /* Form output. */
  for (i=0; i<nvals; i++)
    {
      output[i*2]=indices[i];
      output[indices[i]*2+1]=i;
    }
  free(nrepeat);
  free(indices);
}
#endif

#define NUM_PREVIOUS 4
#define MAX_LEN 0xFFFF
#define MAX_OFFSET 0xFFFF
#define MAX_STRING_SEARCH 8

static void add_circular(int *previous, const int v, const int i)
{
  if (previous[(NUM_PREVIOUS+3)*v+2]!=i-1)
    {
      previous[(NUM_PREVIOUS+3)*v]++;
      if (previous[(NUM_PREVIOUS+3)*v]>NUM_PREVIOUS)
        previous[(NUM_PREVIOUS+3)*v]=NUM_PREVIOUS;
      previous[(NUM_PREVIOUS+3)*v+3+previous[(NUM_PREVIOUS+3)*v+1]]=i;
      previous[(NUM_PREVIOUS+3)*v+1]++;
      if (previous[(NUM_PREVIOUS+3)*v+1]>=NUM_PREVIOUS)
        previous[(NUM_PREVIOUS+3)*v+1]=0;
    }
  previous[(NUM_PREVIOUS+3)*v+2]=i;
}

void Ptngc_comp_to_lz77(unsigned int *vals, const int nvals,
                  unsigned int *data, int *ndata,
                  unsigned int *len, int *nlens,
                  unsigned int *offsets, int *noffsets)
{
  int noff=0;
  int ndat=0;
  int nlen=0;
  int i,j;
  int *previous=warnmalloc(0x20000*(NUM_PREVIOUS+3)*sizeof *previous);
#if 0
  unsigned int *info=warnmalloc(2*nvals*sizeof *info);
  sort_strings(vals,nvals,info);
#endif
  for (i=0; i<0x20000; i++)
    {
      previous[(NUM_PREVIOUS+3)*i]=0; /* Number of items in a circular buffer */
      previous[(NUM_PREVIOUS+3)*i+1]=0; /* Pointer to beginning of circular buffer. */
      previous[(NUM_PREVIOUS+3)*i+2]=-2; /* Last offset that had this value. -2 is really never... */
    }
  for (i=0; i<nvals; i++)
    {
      int k;
#if 0
      int kmin,kmax;
#endif
      int firstoffset=i-MAX_OFFSET;
      if (firstoffset<0)
        firstoffset=0;
      if (i!=0)
        {
          int largest_len=0;
          int largest_offset=0;
          int icirc, ncirc;
          /* Is this identical to a previous offset?  Prefer close
             values for offset. Search through circular buffer for the
             possible values for the start of this string. */
          ncirc=previous[(NUM_PREVIOUS+3)*vals[i]];
          for (icirc=0; icirc<ncirc; icirc++)
            {
              int iptr=previous[(NUM_PREVIOUS+3)*vals[i]+1]-icirc-1;
              if (iptr<0)
                iptr+=NUM_PREVIOUS;
              j=previous[(NUM_PREVIOUS+3)*vals[i]+3+iptr];
              if (j<firstoffset)
                break;
#if 0
              fprintf(stderr,"Starting search for %d at %d. Found %d\n",vals[i],j,vals[j]);
#endif
              while ((j<i) && (vals[j]==vals[i]))
                {
                  if (j>=firstoffset)
                    {
                      for (k=0; i+k<nvals; k++)
                        if (vals[j+k]!=vals[i+k])
                          break;
                      if ((k>largest_len) && ((k>=(i-j)+16) || ((k>4) && (i-j==1))))
                        {
                          largest_len=k;
                          largest_offset=j;
                        }
                    }
                  j++;
                }
            }
#if 0
          /* Search in sorted string buffer too. */
          kmin=info[i*2+1]-MAX_STRING_SEARCH;
          kmax=info[i*2+1]+MAX_STRING_SEARCH;
          if (kmin<0)
            kmin=0;
          if (kmax>=nvals)
            kmax=nvals;
          for (k=kmin; k<kmax; k++)
            {
              int m;
              int s=info[k*2];
              if ((s<i) && (s+MAX_OFFSET>=i))
                {
                  for (m=0; i+m<nvals; m++)
                    if (vals[i+m]!=vals[(s+m)%nvals])
                      break;
                  if ((m>largest_len) && (m>4) && (m+2>=(i-s)))
                    {
                      largest_len=m;
                      largest_offset=s;
#if 0
                      fprintf(stderr,"Offset: %d %d\n",m,i-s);
#endif
                    }
                }
            }
#endif
          /* Check how to write this info. */
          if (largest_len>MAX_LEN)
            largest_len=MAX_LEN;
          if (largest_len)
            {
              if (i-largest_offset==1)
                {
                  data[ndat++]=0;
                }
              else
                {
                  data[ndat++]=1;
                  offsets[noff++]=i-largest_offset;
                }
              len[nlen++]=largest_len;
#if 0
              fprintf(stderr,"l:o: %d:%d data=%d i=%d\n",largest_len,i-largest_offset,ndat,i);
              fflush(stderr);
#endif

#if 0
              fprintf(stderr,"Found largest len %d at %d.\n",largest_len,i-largest_offset);
#endif
              /* Add these values to the circular buffer. */
              for (k=0; k<largest_len; k++)
                add_circular(previous,vals[i+k],i+k);
              i+=largest_len-1;
            }
          else
            {
              data[ndat++]=vals[i]+2;
              /* Add this value to circular buffer. */
              add_circular(previous,vals[i],i);
            }
        }
      else
        {
          data[ndat++]=vals[i]+2;
          /* Add this value to circular buffer. */
          add_circular(previous,vals[i],i);
        }
    }
  *noffsets=noff;
  *ndata=ndat;
  *nlens=nlen;
#if 0
  free(info);
#endif
  free(previous);
}

void Ptngc_comp_from_lz77(unsigned int *data, const int ndata,
                    unsigned int *len, const int nlens,
                    unsigned int *offsets, const int noffsets,
                    unsigned int *vals, const int nvals)
{
  int i=0;
  int joff=0;
  int jdat=0;
  int jlen=0;
  (void)ndata;
  (void)nlens;
  (void)noffsets;
  while (i<nvals)
    {
      unsigned int v=data[jdat++];
      if (v<2)
        {
          int offset=1;
          int k;
          int length=(int)len[jlen++];
#if 0
          fprintf(stderr,"len=%d off=%d i=%d\n",length,offset,i);
#endif
          if (v==1)
            offset=offsets[joff++];
          for (k=0; k<length; k++)
            {
              vals[i]=vals[i-offset];
              if (i>=nvals)
                {
                  fprintf(stderr,"too many vals.\n");
                  exit(EXIT_FAILURE);
                }
              i++;
            }
        }
      else
        vals[i++]=v-2;
    }
}
