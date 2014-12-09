/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg and Magnus Lundborg
 * Copyright (c) 2010, 2013-2014 The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

/* This code is heavily influenced by
   http://hpcv100.rc.rug.nl/xdrf.html
   Based on coordinate compression (c) by Frans van Hoesel.
   and GROMACS xtc files (http://www.gromacs.org)
   (c) Copyright (c) Erik Lindahl, David van der Spoel
*/

/* The cost estimates are ripped right out of xtc2.c, so take these
   with a grain (truckload) of salt. */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/widemuldiv.h"
#include "../../include/compression/bwlzh.h"

static const double iflipgaincheck=0.89089871814033927; /*  1./(2**(1./6)) */

#define MAX_LARGE_RLE 1024 /* Maximum number of large atoms for large RLE. */
#define MAX_SMALL_RLE 12 /* Maximum number of small atoms in one group. */

#define TRESHOLD_INTRA_INTER_DIRECT 1.5 /* How much larger can the direct
                                           frame deltas for the small
                                           triplets be and be accepted anyway
                                           as better than the intra/inter frame
                                           deltas. For better instructions/RLEs. */

#define TRESHOLD_INTER_INTRA 5.0 /* How much larger can the intra
                                    frame deltas for the small
                                    triplets be and be accepted anyway
                                    as better than the inter frame
                                    deltas. */

/* Difference in indices used for determining whether to store as
   large or small. A fun detail in this compression algorithm is that
   if everything works fine, large can often be smaller than small, or
   at least not as large as is large in magic.c. This is a key idea of
   xtc3. */
#define QUITE_LARGE 3
#define IS_LARGE 6

#if 0
#define SHOWIT
#endif

#if 0
#define SHOWIT_LIGHT
#endif

#ifdef USE_WINDOWS
#define TNG_INLINE __inline
#else
#define TNG_INLINE inline
#endif

/* These routines are in xtc2.c */
int Ptngc_magic(unsigned int i);
int Ptngc_find_magic_index(const unsigned int maxval);

static TNG_INLINE unsigned int positive_int(const int item)
{
  int s=0;
  if (item>0)
    s=1+(item-1)*2;
  else if (item<0)
    s=2+(-item-1)*2;
  return s;
}

static TNG_INLINE int unpositive_int(const int val)
{
  int s=(val+1)/2;
  if ((val%2)==0)
    s=-s;
  return s;
}


/* Sequence instructions */
#define INSTR_DEFAULT 0U
#define INSTR_SMALL_RUNLENGTH 1U
#define INSTR_ONLY_LARGE 2U
#define INSTR_ONLY_SMALL 3U
#define INSTR_FLIP 4U
#define INSTR_LARGE_RLE 5U
#define INSTR_LARGE_DIRECT 6U
#define INSTR_LARGE_INTRA_DELTA 7U
#define INSTR_LARGE_INTER_DELTA 8U

#define MAXINSTR 9

struct xtc3_context
{
  unsigned int *instructions;
  int ninstr, ninstr_alloc;
  unsigned int *rle;
  int nrle, nrle_alloc;
  unsigned int *large_direct;
  int nlargedir, nlargedir_alloc;
  unsigned int *large_intra_delta;
  int nlargeintra, nlargeintra_alloc;
  unsigned int *large_inter_delta;
  int nlargeinter, nlargeinter_alloc;
  unsigned int *smallintra;
  int nsmallintra, nsmallintra_alloc;
  int minint[3],maxint[3];
  int has_large;
  int has_large_ints[MAX_LARGE_RLE*3]; /* Large cache. */
  int has_large_type[MAX_LARGE_RLE]; /* What kind of type this large
                                        int is. */
  int current_large_type;
};

static void init_xtc3_context(struct xtc3_context *xtc3_context)
{
  xtc3_context->instructions=NULL;
  xtc3_context->ninstr=0;
  xtc3_context->ninstr_alloc=0;
  xtc3_context->rle=NULL;
  xtc3_context->nrle=0;
  xtc3_context->nrle_alloc=0;
  xtc3_context->large_direct=NULL;
  xtc3_context->nlargedir=0;
  xtc3_context->nlargedir_alloc=0;
  xtc3_context->large_intra_delta=NULL;
  xtc3_context->nlargeintra=0;
  xtc3_context->nlargeintra_alloc=0;
  xtc3_context->large_inter_delta=NULL;
  xtc3_context->nlargeinter=0;
  xtc3_context->nlargeinter_alloc=0;
  xtc3_context->smallintra=NULL;
  xtc3_context->nsmallintra=0;
  xtc3_context->nsmallintra_alloc=0;
  xtc3_context->has_large=0;
  xtc3_context->current_large_type=0;
}

static void free_xtc3_context(struct xtc3_context *xtc3_context)
{
  free(xtc3_context->instructions);
  free(xtc3_context->rle);
  free(xtc3_context->large_direct);
  free(xtc3_context->large_intra_delta);
  free(xtc3_context->large_inter_delta);
  free(xtc3_context->smallintra);
}

/* Modifies three integer values for better compression of water */
static void swap_ints(const int *in, int *out)
{
  out[0]=in[0]+in[1];
  out[1]=-in[1];
  out[2]=in[1]+in[2];
}

static void swap_is_better(const int *input, const int *minint, int *sum_normal, int *sum_swapped)
{
  int normal_max=0;
  int swapped_max=0;
  int i,j;
  int normal[3];
  int swapped[3];
  for (i=0; i<3; i++)
    {
      normal[0]=input[i]-minint[i];
      normal[1]=input[3+i]-input[i]; /* minint[i]-minint[i] cancels out */
      normal[2]=input[6+i]-input[3+i]; /* minint[i]-minint[i] cancels out */
      swap_ints(normal,swapped);
      for (j=1; j<3; j++)
        {
          if (positive_int(normal[j])>(unsigned int)normal_max)
            normal_max=positive_int(normal[j]);
          if (positive_int(swapped[j])>(unsigned int)swapped_max)
            swapped_max=positive_int(swapped[j]);
        }
    }
  if (normal_max==0)
    normal_max=1;
  if (swapped_max==0)
    swapped_max=1;
  *sum_normal=normal_max;
  *sum_swapped=swapped_max;
}

static void allocate_enough_memory(unsigned int **ptr, int *nele, int *nele_alloc)
{
  (*nele)++;
  if (*nele>*nele_alloc)
    {
      *nele_alloc=*nele + *nele/2;
      *ptr=warnrealloc(*ptr,*nele_alloc*sizeof **ptr);
    }
}

static void insert_value_in_array(unsigned int **ptr, int *nele, int *nele_alloc,
                                  const unsigned int value,
                                  char *arrayname)
{
#ifndef SHOWIT
  (void)arrayname;
#endif
  allocate_enough_memory(ptr,nele,nele_alloc);
#ifdef SHOWIT
  fprintf(stderr,"Inserting value %u into array %s @ %d\n",value,arrayname,(*nele)-1);
#endif
  (*ptr)[(*nele)-1]=value;
}



static void swapdecide(struct xtc3_context *xtc3_context, int *input,int *swapatoms, int *large_index, int *minint)
{
  int didswap=0;
  int normal,swapped;
  (void)large_index;
  swap_is_better(input,minint,&normal,&swapped);
  /* We have to determine if it is worth to change the behaviour.
     If diff is positive it means that it is worth something to
     swap. But it costs 4 bits to do the change. If we assume that
     we gain 0.17 bit by the swap per value, and the runlength>2
     for four molecules in a row, we gain something. So check if we
     gain at least 0.17 bits to even attempt the swap.
  */
#ifdef SHOWIT
  fprintf(stderr,"Trying Flip: %g %g\n",(double)swapped/normal, (double)normal/swapped);
#endif
  if (((swapped<normal) && (fabs((double)swapped/normal)<iflipgaincheck)) ||
      ((normal<swapped) && (fabs((double)normal/swapped)<iflipgaincheck)))
    {
      if (swapped<normal)
        {
          if (!*swapatoms)
            {
              *swapatoms=1;
              didswap=1;
            }
        }
      else
        {
          if (*swapatoms)
            {
              *swapatoms=0;
              didswap=1;
            }
        }
    }
  if (didswap)
    {
#ifdef SHOWIT
      fprintf(stderr,"Flip. Swapatoms is now %d\n",*swapatoms);
#endif
      insert_value_in_array(&xtc3_context->instructions,
                            &xtc3_context->ninstr,
                            &xtc3_context->ninstr_alloc,
                            INSTR_FLIP,"instr");
    }
}

/* It is "large" if we have to increase the small index quite a
   bit. Not so much to be rejected by the not very large check
   later. */
static int is_quite_large(const int *input, const int small_index, const int max_large_index)
{
  int is=0;
  int i;
  if (small_index+QUITE_LARGE>=max_large_index)
    is=1;
  else
    {
      for (i=0; i<3; i++)
        if (positive_int(input[i])>(unsigned int)Ptngc_magic(small_index+QUITE_LARGE))
          {
            is=1;
            break;
          }
    }
  return is;
}

#ifdef SHOWIT
int nbits_sum;
int nvalues_sum;
#endif

static void insert_batch(const int *input_ptr, const int ntriplets_left, const int *prevcoord, int *encode_ints, const int startenc, int *nenc)
{
  int nencode=startenc*3;
  int tmp_prevcoord[3];

  tmp_prevcoord[0]=prevcoord[0];
  tmp_prevcoord[1]=prevcoord[1];
  tmp_prevcoord[2]=prevcoord[2];

  if (startenc)
    {
      int i;
      for (i=0; i<startenc; i++)
        {
          tmp_prevcoord[0]+=encode_ints[i*3];
          tmp_prevcoord[1]+=encode_ints[i*3+1];
          tmp_prevcoord[2]+=encode_ints[i*3+2];
#ifdef SHOWIT
          fprintf(stderr,"%6d: %6d %6d %6d\t\t%6d %6d %6d\t\t%6d %6d %6d\n",i*3,
                  tmp_prevcoord[0],tmp_prevcoord[1],tmp_prevcoord[2],
                  encode_ints[i*3],
                  encode_ints[i*3+1],
                  encode_ints[i*3+2],
                  positive_int(encode_ints[i*3]),
                  positive_int(encode_ints[i*3+1]),
                  positive_int(encode_ints[i*3+2]));
#endif
        }
    }

#ifdef SHOWIT
  fprintf(stderr,"New batch\n");
#endif
  while ((nencode<3+MAX_SMALL_RLE*3) && (nencode<ntriplets_left*3))
    {
      encode_ints[nencode]=input_ptr[nencode]-tmp_prevcoord[0];
      encode_ints[nencode+1]=input_ptr[nencode+1]-tmp_prevcoord[1];
      encode_ints[nencode+2]=input_ptr[nencode+2]-tmp_prevcoord[2];
#ifdef SHOWIT
      fprintf(stderr,"%6d: %6d %6d %6d\t\t%6d %6d %6d\t\t%6d %6d %6d\n",nencode,
              input_ptr[nencode],
              input_ptr[nencode+1],
              input_ptr[nencode+2],
              encode_ints[nencode],
              encode_ints[nencode+1],
              encode_ints[nencode+2],
              positive_int(encode_ints[nencode]),
              positive_int(encode_ints[nencode+1]),
              positive_int(encode_ints[nencode+2]));
#endif
      tmp_prevcoord[0]=input_ptr[nencode];
      tmp_prevcoord[1]=input_ptr[nencode+1];
      tmp_prevcoord[2]=input_ptr[nencode+2];
      nencode+=3;
    }
  *nenc=nencode;
}

static void large_instruction_change(struct xtc3_context *xtc3_context, const int i)
{
  /* If the first large is of a different kind than the currently used we must
     emit an "instruction" to change the large type. */
  if (xtc3_context->has_large_type[i]!=xtc3_context->current_large_type)
    {
      unsigned int instr;
      xtc3_context->current_large_type=xtc3_context->has_large_type[i];
      if (xtc3_context->current_large_type==0)
        instr=INSTR_LARGE_DIRECT;
      else if (xtc3_context->current_large_type==1)
        instr=INSTR_LARGE_INTRA_DELTA;
      else
        instr=INSTR_LARGE_INTER_DELTA;
      insert_value_in_array(&xtc3_context->instructions,
                            &xtc3_context->ninstr,
                            &xtc3_context->ninstr_alloc,
                            instr,"instr");
    }
}

static void write_three_large(struct xtc3_context *xtc3_context,
                              const int i)
{
  int m;
  if (xtc3_context->current_large_type==0)
    {
      for (m=0; m<3; m++)
        insert_value_in_array(&xtc3_context->large_direct,
                              &xtc3_context->nlargedir,
                              &xtc3_context->nlargedir_alloc,
                              xtc3_context->has_large_ints[i*3+m],"large direct");
    }
  else if (xtc3_context->current_large_type==1)
    {
      for (m=0; m<3; m++)
        insert_value_in_array(&xtc3_context->large_intra_delta,
                              &xtc3_context->nlargeintra,
                              &xtc3_context->nlargeintra_alloc,
                              xtc3_context->has_large_ints[i*3+m],"large intra");
    }
  else
    {
      for (m=0; m<3; m++)
        insert_value_in_array(&xtc3_context->large_inter_delta,
                              &xtc3_context->nlargeinter,
                              &xtc3_context->nlargeinter_alloc,
                              xtc3_context->has_large_ints[i*3+m],"large inter");
    }
}

static void flush_large(struct xtc3_context *xtc3_context,
                        const int n) /* How many to flush. */
{
  int i;
  i=0;
  while (i<n)
    {
      int j,k;
      /* If the first large is of a different kind than the currently used we must
         emit an "instruction" to change the large type. */
      large_instruction_change(xtc3_context,i);
      /* How many large of the same kind in a row? */
      for (j=0;
           (i+j<n) &&
             (xtc3_context->has_large_type[i+j]==xtc3_context->has_large_type[i]);
           j++);
      if (j<3)
        {
          for (k=0; k<j; k++)
            {
              insert_value_in_array(&xtc3_context->instructions,
                                    &xtc3_context->ninstr,
                                    &xtc3_context->ninstr_alloc,
                                    INSTR_ONLY_LARGE,"instr");
              write_three_large(xtc3_context,i+k);
            }
        }
      else
        {
          insert_value_in_array(&xtc3_context->instructions,
                                &xtc3_context->ninstr,
                                &xtc3_context->ninstr_alloc,
                                INSTR_LARGE_RLE,"instr");
          insert_value_in_array(&xtc3_context->rle,
                                &xtc3_context->nrle,
                                &xtc3_context->nrle_alloc,
                                (unsigned int)j,"rle (large)");
          for (k=0; k<j; k++)
            write_three_large(xtc3_context,i+k);
        }
      i+=j;
    }
  if ((xtc3_context->has_large-n)!=0)
    {
      int j;
      for (i=0; i<xtc3_context->has_large-n; i++)
        {
          xtc3_context->has_large_type[i]=xtc3_context->has_large_type[i+n];
          for (j=0; j<3; j++)
            xtc3_context->has_large_ints[i*3+j]=xtc3_context->has_large_ints[(i+n)*3+j];
        }
    }
  xtc3_context->has_large-=n; /* Number of remaining large atoms in buffer */
}

static double compute_intlen(unsigned int *ints)
{
  /* The largest value. */
  unsigned int m=ints[0];
  if (ints[1]>m)
    m=ints[1];
  if (ints[2]>m)
    m=ints[2];
  return (double)m;
}

static void buffer_large(struct xtc3_context *xtc3_context, int *input, const int inpdata,
                         const int natoms, const int intradelta_ok)
{
  unsigned int direct[3], intradelta[3]={0,}, interdelta[3]={0,};
  double minlen;
  int best_type;
  int frame=inpdata/(natoms*3);
  int atomframe=inpdata%(natoms*3);
  /* If it is full we must write them all. */
  if (xtc3_context->has_large==MAX_LARGE_RLE)
    flush_large(xtc3_context,xtc3_context->has_large); /* Flush all. */
  /* Find out which is the best choice for the large integer. Direct coding, or some
     kind of delta coding? */
  /* First create direct coding. */
  direct[0]=(unsigned int)(input[inpdata]-xtc3_context->minint[0]);
  direct[1]=(unsigned int)(input[inpdata+1]-xtc3_context->minint[1]);
  direct[2]=(unsigned int)(input[inpdata+2]-xtc3_context->minint[2]);
  minlen=compute_intlen(direct);
  best_type=0; /* Direct. */
#if 1
  /* Then try intra coding if we can. */
  if ((intradelta_ok) && (atomframe>=3))
    {
      double thislen;
      intradelta[0]=positive_int(input[inpdata]-input[inpdata-3]);
      intradelta[1]=positive_int(input[inpdata+1]-input[inpdata-2]);
      intradelta[2]=positive_int(input[inpdata+2]-input[inpdata-1]);
      thislen=compute_intlen(intradelta);
      if (thislen*TRESHOLD_INTRA_INTER_DIRECT<minlen)
        {
          minlen=thislen;
          best_type=1; /* Intra delta */
        }
    }
#endif
#if 1
  /* Then try inter coding if we can. */
  if (frame>0)
    {
      double thislen;
      interdelta[0]=positive_int(input[inpdata]-input[inpdata-natoms*3]);
      interdelta[1]=positive_int(input[inpdata+1]-input[inpdata-natoms*3+1]);
      interdelta[2]=positive_int(input[inpdata+2]-input[inpdata-natoms*3+2]);
      thislen=compute_intlen(interdelta);
      if (thislen*TRESHOLD_INTRA_INTER_DIRECT<minlen)
        {
          best_type=2; /* Inter delta */
        }
    }
#endif
  xtc3_context->has_large_type[xtc3_context->has_large]=best_type;
  if (best_type==0)
    {
      xtc3_context->has_large_ints[xtc3_context->has_large*3]=direct[0];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+1]=direct[1];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+2]=direct[2];
    }
  else if (best_type==1)
    {
      xtc3_context->has_large_ints[xtc3_context->has_large*3]=intradelta[0];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+1]=intradelta[1];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+2]=intradelta[2];
    }
  else if (best_type==2)
    {
      xtc3_context->has_large_ints[xtc3_context->has_large*3]=interdelta[0];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+1]=interdelta[1];
      xtc3_context->has_large_ints[xtc3_context->has_large*3+2]=interdelta[2];
    }
  xtc3_context->has_large++;
}

static void output_int(unsigned char *output,int *outdata, const unsigned int n)
{
  output[(*outdata)++]=((unsigned int)n)&0xFFU;
  output[(*outdata)++]=(((unsigned int)n)>>8)&0xFFU;
  output[(*outdata)++]=(((unsigned int)n)>>16)&0xFFU;
  output[(*outdata)++]=(((unsigned int)n)>>24)&0xFFU;
}

#if 0
static void printarray(unsigned int *a, int n, char *name)
{
  int i;
  for (i=0; i<n; i++)
    fprintf(stderr,"%u %s\n",a[i],name);
}
#endif

/* The base_compress routine first compresses all x coordinates, then
   y and finally z. The bases used for each can be different. The
   MAXBASEVALS value determines how many coordinates are compressed
   into a single number. Only resulting whole bytes are dealt with for
   simplicity. MAXMAXBASEVALS is the insanely large value to accept
   files written with that value. BASEINTERVAL determines how often a
   new base is actually computed and stored in the output
   file. MAXBASEVALS*BASEINTERVAL values are stored using the same
   base in BASEINTERVAL different integers. Note that the primarily
   the decompression using a large MAXBASEVALS becomes very slow. */
#define MAXMAXBASEVALS 16384U
#define MAXBASEVALS 24U
#define BASEINTERVAL 8

/* How many bytes are needed to store n values in base base */
static int base_bytes(const unsigned int base, const int n)
{
  int i,j;
  unsigned int largeint[MAXMAXBASEVALS+1];
  unsigned int largeint_tmp[MAXMAXBASEVALS+1];
  int numbytes=0;

  memset(largeint, 0U, sizeof(unsigned int) * (n+1));

  for (i=0; i<n; i++)
    {
      if (i!=0)
        {
          Ptngc_largeint_mul(base,largeint,largeint_tmp,n+1);
          memcpy(largeint, largeint_tmp, (n+1)*sizeof *largeint);
        }
      Ptngc_largeint_add(base-1U,largeint,n+1);
    }
  for (i=0; i<n; i++)
    if (largeint[i])
      for (j=0; j<4; j++)
        if ((largeint[i]>>(j*8))&0xFFU)
          numbytes=i*4+j+1;
  return numbytes;
}

static void base_compress(unsigned int *data, const int len, unsigned char *output, int *outlen)
{
  unsigned int largeint[MAXBASEVALS+1];
  unsigned int largeint_tmp[MAXBASEVALS+1];
  int ixyz, i;
  unsigned int j;
  int nwrittenout=0;
  unsigned int numbytes=0;
  /* Store the MAXBASEVALS value in the output. */
  output[nwrittenout++]=(unsigned char)(MAXBASEVALS&0xFFU);
  output[nwrittenout++]=(unsigned char)((MAXBASEVALS>>8)&0xFFU);
  /* Store the BASEINTERVAL value in the output. */
  output[nwrittenout++]=(unsigned char)(BASEINTERVAL&0xFFU);
  for (ixyz=0; ixyz<3; ixyz++)
    {
      unsigned int base=0U;
      int nvals=0;
      int basegiven=0;

      memset(largeint, 0U, sizeof(unsigned int) * (MAXBASEVALS+1));

      for (i=ixyz; i<len; i+=3)
        {
         if (nvals==0)
           {
             int basecheckvals=0;
             int k;
             if (basegiven==0)
               {
                 base=0U;
                 /* Find the largest value for this particular coordinate. */
                 for (k=i; k<len; k+=3)
                   {
                     if (data[k]>base)
                       base=data[k];
                     basecheckvals++;
                     if (basecheckvals==MAXBASEVALS*BASEINTERVAL)
                       break;
                   }
                 /* The base is one larger than the largest values. */
                 base++;
                 if (base<2)
                   base=2;
                 /* Store the base in the output. */
                 output[nwrittenout++]=(unsigned char)(base&0xFFU);
                 output[nwrittenout++]=(unsigned char)((base>>8)&0xFFU);
                 output[nwrittenout++]=(unsigned char)((base>>16)&0xFFU);
                 output[nwrittenout++]=(unsigned char)((base>>24)&0xFFU);
                 basegiven=BASEINTERVAL;
                 /* How many bytes is needed to store MAXBASEVALS values using this base? */
                 numbytes=base_bytes(base,MAXBASEVALS);
               }
             basegiven--;
#ifdef SHOWIT
             fprintf(stderr,"Base for %d is %u. I need %d bytes for %d values.\n",ixyz,base,numbytes,MAXBASEVALS);
#endif
           }
          if (nvals!=0)
            {
              Ptngc_largeint_mul(base,largeint,largeint_tmp,MAXBASEVALS+1);
              for (j=0; j<MAXBASEVALS+1; j++)
                largeint[j]=largeint_tmp[j];
            }
          Ptngc_largeint_add(data[i],largeint,MAXBASEVALS+1);
#ifdef SHOWIT
          fprintf(stderr,"outputting value %u\n",data[i]);
#endif
          nvals++;
          if (nvals==MAXBASEVALS)
            {
#ifdef SHOWIT
              fprintf(stderr,"Writing largeint: ");
#endif
              for (j=0; j<numbytes; j++)
                {
                  int ilarge=j/4;
                  int ibyte=j%4;
                  output[nwrittenout++]=(unsigned char)((largeint[ilarge]>>(ibyte*8))&(0xFFU));
#ifdef SHOWIT
                  fprintf(stderr,"%02x",(unsigned int)output[nwrittenout-1]);
#endif
                }
#ifdef SHOWIT
              fprintf(stderr,"\n");
#endif
              nvals=0;

              memset(largeint, 0U, sizeof(unsigned int) * (MAXBASEVALS+1));
            }
        }
      if (nvals)
        {
          numbytes=base_bytes(base,nvals);
#ifdef SHOWIT
          fprintf(stderr,"Base for %d is %u. I need %d bytes for %d values.\n",ixyz,base,numbytes,nvals);
#endif
          for (j=0; j<numbytes; j++)
            {
              int ilarge=j/4;
              int ibyte=j%4;
              output[nwrittenout++]=(unsigned char)((largeint[ilarge]>>(ibyte*8))&(0xFFU));
            }
        }
    }
  *outlen=nwrittenout;
}

static void base_decompress(unsigned char *input, const int len, unsigned int *output)
{
  unsigned int largeint[MAXMAXBASEVALS+1];
  unsigned int largeint_tmp[MAXMAXBASEVALS+1];
  int ixyz, i, j;
  int maxbasevals=(int)((unsigned int)(input[0])|(((unsigned int)(input[1]))<<8));
  int baseinterval=(int)input[2];
  if (maxbasevals>(int)MAXMAXBASEVALS)
    {
      fprintf(stderr,"Read a larger maxbasevals value from the file than I can handle. Fix"
              " by increasing MAXMAXBASEVALS to at least %d. Although, this is"
              " probably a bug in TRAJNG, since MAXMAXBASEVALS should already be insanely large enough.\n",maxbasevals);
      exit(EXIT_FAILURE);
    }
  input+=3;
  for (ixyz=0; ixyz<3; ixyz++)
    {
      int numbytes=0;
      int nvals_left=len/3;
      int outvals=ixyz;
      int basegiven=0;
      unsigned int base=0U;
#ifdef SHOWIT
      fprintf(stderr,"Base for %d is %u. I need %d bytes for %d values.\n",ixyz,base,numbytes,maxbasevals);
#endif
      while (nvals_left)
        {
          int n;
          if (basegiven==0)
            {
              base=(unsigned int)(input[0])|
                (((unsigned int)(input[1]))<<8)|
                (((unsigned int)(input[2]))<<16)|
                (((unsigned int)(input[3]))<<24);
              input+=4;
              basegiven=baseinterval;
              /* How many bytes is needed to store maxbasevals values using this base? */
              numbytes=base_bytes(base,maxbasevals);
            }
          basegiven--;
          if (nvals_left<maxbasevals)
            {
              numbytes=base_bytes(base,nvals_left);
#ifdef SHOWIT
              fprintf(stderr,"Base for %d is %u. I need %d bytes for %d values.\n",ixyz,base,numbytes,nvals_left);
#endif
            }
          memset(largeint, 0U, sizeof(unsigned int) * (maxbasevals+1));
#ifdef SHOWIT
          fprintf(stderr,"Reading largeint: ");
#endif
          if (numbytes/4 < maxbasevals+1)
            {
              for (j=0; j<numbytes; j++)
                {
                  int ilarge=j/4;
                  int ibyte=j%4;
                  largeint[ilarge]|=((unsigned int)input[j])<<(ibyte*8);
#ifdef SHOWIT
                  fprintf(stderr,"%02x",(unsigned int)input[j]);
#endif
                }
            }
#ifdef SHOWIT
          fprintf(stderr,"\n");
#endif
          input+=numbytes;
          /* Do the long division required to get the output values. */
          n=maxbasevals;
          if (n>nvals_left)
            n=nvals_left;
          for (i=n-1; i>=0; i--)
            {
              output[outvals+i*3]=Ptngc_largeint_div(base,largeint,largeint_tmp,maxbasevals+1);
              for (j=0; j<maxbasevals+1; j++)
                largeint[j]=largeint_tmp[j];
            }
#ifdef SHOWIT
          for (i=0; i<n; i++)
            fprintf(stderr,"outputting value %u\n",output[outvals+i*3]);
#endif
          outvals+=n*3;
          nvals_left-=n;
        }
    }
}

/* If a large proportion of the integers are large (More than 10\% are >14 bits) we return 0, otherwise 1 */
static int heuristic_bwlzh(unsigned int *ints, const int nints)
{
  int i,num;
  num=0;
  for (i=0; i<nints; i++)
    if (ints[i]>=16384)
      num++;
  if (num>nints/10)
    return 0;
  else
    return 1;
}

/* Speed selects how careful to try to find the most efficient compression. The BWLZH algo is expensive!
   Speed <=2 always avoids BWLZH everywhere it is possible.
   Speed 3 and 4 and 5 use heuristics (check proportion of large value). This should mostly be safe.
   Speed 5 enables the LZ77 component of BWLZH.
   Speed 6 always tests if BWLZH is better and if it is uses it. This can be very slow.
 */
unsigned char *Ptngc_pack_array_xtc3(int *input, int *length, const int natoms, int speed)
{
  unsigned char *output=NULL;
  int i,ienc,j;
  int outdata=0;
  /* Pack triplets. */
  int ntriplets=*length/3;
  int intmax;
  int max_small;
  int small_index;
  int max_large_index;
  int large_index[3];
  int prevcoord[3];
  int runlength=0; /* Initial runlength. "Stupidly" set to zero for
                      simplicity and explicity */
  int swapatoms=0; /* Initial guess is that we should not swap the
                      first two atoms in each large+small
                      transition */
  int didswap; /* Whether swapping was actually done. */
  int inpdata=0;
  int encode_ints[3+MAX_SMALL_RLE*3]; /* Up to 3 large + 24 small ints can be encoded at once */
  int nencode;
  int ntriplets_left=ntriplets;
  int refused=0;
  unsigned char *bwlzh_buf=NULL;
  int bwlzh_buf_len;
  unsigned char *base_buf=NULL;
  int base_buf_len;

  struct xtc3_context xtc3_context;
  init_xtc3_context(&xtc3_context);

  memcpy(xtc3_context.maxint, input, 3*sizeof *xtc3_context.maxint);
  memcpy(xtc3_context.minint, input, 3*sizeof *xtc3_context.maxint);

  /* Values of speed should be sane. */
  if (speed<1)
    speed=1;
  if (speed>6)
    speed=6;

#ifdef SHOWIT
  nbits_sum=0;
  nvalues_sum=0;
#endif
  /* Allocate enough memory for output */
  if (*length < 48)
    output=warnmalloc(8*48*sizeof *output);
  else
    output=warnmalloc(8* *length*sizeof *output);


  for (i=1; i<ntriplets; i++)
    for (j=0; j<3; j++)
      {
        if (input[i*3+j]>xtc3_context.maxint[j])
          xtc3_context.maxint[j]=input[i*3+j];
        if (input[i*3+j]<xtc3_context.minint[j])
          xtc3_context.minint[j]=input[i*3+j];
      }

  large_index[0]=Ptngc_find_magic_index(xtc3_context.maxint[0]-xtc3_context.minint[0]+1);
  large_index[1]=Ptngc_find_magic_index(xtc3_context.maxint[1]-xtc3_context.minint[1]+1);
  large_index[2]=Ptngc_find_magic_index(xtc3_context.maxint[2]-xtc3_context.minint[2]+1);
  max_large_index=large_index[0];
  if (large_index[1]>max_large_index)
    max_large_index=large_index[1];
  if (large_index[2]>max_large_index)
    max_large_index=large_index[2];

#ifdef SHOWIT
  for (j=0; j<3; j++)
    fprintf(stderr,"minint[%d]=%d. maxint[%d]=%d large_index[%d]=%d value=%d\n",j,xtc3_context.minint[j],j,xtc3_context.maxint[j],
            j,large_index[j],Ptngc_magic(large_index[j]));
#endif

  /* Guess initial small index */
  small_index=max_large_index/2;

  /* Find the largest value that is not large. Not large is half index of
     large. */
  max_small=Ptngc_magic(small_index);
  intmax=0;
  for (i=0; i<*length; i++)
    {
      int item=input[i];
      int s=positive_int(item);
      if (s>intmax)
        if (s<max_small)
          intmax=s;
    }
  /* This value is not critical, since if I guess wrong, the code will
     just insert instructions to increase this value. */
  small_index=Ptngc_find_magic_index(intmax);
#ifdef SHOWIT
  fprintf(stderr,"initial small_index=%d value=%d\n",small_index,Ptngc_magic(small_index));
#endif

  output_int(output,&outdata,positive_int(xtc3_context.minint[0]));
  output_int(output,&outdata,positive_int(xtc3_context.minint[1]));
  output_int(output,&outdata,positive_int(xtc3_context.minint[2]));

#if 0
#ifdef SHOWIT
  for (i=0; i<ntriplets_left; i++)
    fprintf(stderr,"VALUE:%d %6d %6d %6d\n",
            i,
            input[inpdata+i*3],
            input[inpdata+i*3+1],
            input[inpdata+i*3+2]);
#endif
#endif

  /* Initial prevcoord is the minimum integers. */
  memcpy(prevcoord, xtc3_context.minint, 3*sizeof *prevcoord);
  prevcoord[0]=xtc3_context.minint[0];
  prevcoord[1]=xtc3_context.minint[1];
  prevcoord[2]=xtc3_context.minint[2];

  while (ntriplets_left)
    {
      if (ntriplets_left<0)
        {
          fprintf(stderr,"TRAJNG: BUG! ntriplets_left<0!\n");
          exit(EXIT_FAILURE);
        }
      /* If only less than three atoms left we just write them all as large integers. Here no swapping is done! */
      if (ntriplets_left<3)
        {
          for (ienc=0; ienc<ntriplets_left; ienc++)
            {
              buffer_large(&xtc3_context,input,inpdata,natoms,1);
              inpdata+=3;
              ntriplets_left--;
            }
          flush_large(&xtc3_context,xtc3_context.has_large); /* Flush all */
        }
      else
        {
          int min_runlength=0;
          int largest_required_base;
          int largest_runlength_base;
          int largest_required_index;
          int largest_runlength_index;
          int new_runlength;
          int new_small_index;
          int iter_runlength;
          int iter_small_index;
          int rle_index_dep;
          didswap=0;
          /* Insert the next batch of integers to be encoded into the buffer */
#ifdef SHOWIT
          fprintf(stderr,"Initial batch\n");
#endif
          insert_batch(input+inpdata,ntriplets_left,prevcoord,encode_ints,0,&nencode);

          /* First we must decide if the next value is large (does not reasonably fit in current small encoding)
             Also, if we have not written any values yet, we must begin by writing a large atom. */
          if ((inpdata==0) || (is_quite_large(encode_ints,small_index,max_large_index)) || (refused))
            {
              /* If any of the next two atoms are large we should probably write them as large and not swap them */
              int no_swap=0;
              if ((is_quite_large(encode_ints+3,small_index,max_large_index)) || (is_quite_large(encode_ints+6,small_index,max_large_index)))
                no_swap=1;
#if 1
              if (!no_swap)
                {
                  /* If doing inter-frame coding results in smaller values we should not do any swapping either. */
                  int frame=inpdata/(natoms*3);
                  if (frame>0)
                    {
                      unsigned int delta[3], delta2[3];
                      delta[0]=positive_int(input[inpdata+3]-input[inpdata-natoms*3+3]);
                      delta[1]=positive_int(input[inpdata+4]-input[inpdata-natoms*3+4]);
                      delta[2]=positive_int(input[inpdata+5]-input[inpdata-natoms*3+5]);
                      delta2[0]=positive_int(encode_ints[3]);
                      delta2[1]=positive_int(encode_ints[4]);
                      delta2[2]=positive_int(encode_ints[5]);
#ifdef SHOWIT
                      fprintf(stderr,"A1: inter delta: %u %u %u. intra delta=%u %u %u\n",
                              delta[0],delta[1],delta[2],
                              delta2[0],delta2[1],delta2[2]);
#endif
                      if (compute_intlen(delta)*TRESHOLD_INTER_INTRA<compute_intlen(delta2))
                        {
                          delta[0]=positive_int(input[inpdata+6]-input[inpdata-natoms*3+6]);
                          delta[1]=positive_int(input[inpdata+7]-input[inpdata-natoms*3+7]);
                          delta[2]=positive_int(input[inpdata+8]-input[inpdata-natoms*3+8]);
                          delta2[0]=positive_int(encode_ints[6]);
                          delta2[1]=positive_int(encode_ints[7]);
                          delta2[2]=positive_int(encode_ints[8]);
#ifdef SHOWIT
                          fprintf(stderr,"A2: inter delta: %u %u %u. intra delta=%u %u %u\n",
                                  delta[0],delta[1],delta[2],
                                  delta2[0],delta2[1],delta2[2]);
#endif
                          if (compute_intlen(delta)*TRESHOLD_INTER_INTRA<compute_intlen(delta2))
                            {
                              no_swap=1;
#ifdef SHOWIT
                              fprintf(stderr,"SETTING NO SWAP!\n");
#endif
                            }
                        }
                    }
                }
#endif
              if (!no_swap)
                {
                  /* Next we must decide if we should swap the first
                     two values. */
#if 1
                  swapdecide(&xtc3_context,input+inpdata,&swapatoms,large_index,xtc3_context.minint);
#else
                  swapatoms=0;
#endif
                  /* If we should do the integer swapping manipulation we should do it now */
                  if (swapatoms)
                    {
                      didswap=1;
                      for (i=0; i<3; i++)
                        {
                          int in[3], out[3];
                          in[0]=input[inpdata+i];
                          in[1]=input[inpdata+3+i]-input[inpdata+i];
                          in[2]=input[inpdata+6+i]-input[inpdata+3+i];
                          swap_ints(in,out);
                          encode_ints[i]=out[0];
                          encode_ints[3+i]=out[1];
                          encode_ints[6+i]=out[2];
                        }
                      /* We have swapped atoms, so the minimum run-length is 2 */
#ifdef SHOWIT
                      fprintf(stderr,"Swap atoms results in:\n");
                      for (i=0; i<3; i++)
                        fprintf(stderr,"%d: %6d %6d %6d\t\t%6d %6d %6d\n",i*3,
                                encode_ints[i*3],
                                encode_ints[i*3+1],
                                encode_ints[i*3+2],
                                positive_int(encode_ints[i*3]),
                                positive_int(encode_ints[i*3+1]),
                                positive_int(encode_ints[i*3+2]));

#endif
                      min_runlength=2;
                    }
                }
              /* Cache large value for later possible combination with
                 a sequence of small integers. */
              if ((swapatoms) && (didswap))
                {
                  buffer_large(&xtc3_context,input,inpdata+3,natoms,0); /* This is a swapped integer, so inpdata is one atom later and
                                                                           intra coding is not ok. */
                  for (ienc=0; ienc<3; ienc++)
                    prevcoord[ienc]=input[inpdata+3+ienc];
                }
              else
                {
                  buffer_large(&xtc3_context,input,inpdata,natoms,1);
                  for (ienc=0; ienc<3; ienc++)
                    prevcoord[ienc]=input[inpdata+ienc];
                }


#ifdef SHOWIT
              fprintf(stderr,"Prevcoord after packing of large: %d %d %d\n",
                      prevcoord[0],prevcoord[1],prevcoord[2]);
#endif

              /* We have written a large integer so we have one less atoms to worry about */
              inpdata+=3;
              ntriplets_left--;

              refused=0;

              /* Insert the next batch of integers to be encoded into the buffer */
#ifdef SHOWIT
              fprintf(stderr,"Update batch due to large int.\n");
#endif
              if ((swapatoms) && (didswap))
                {
                  /* Keep swapped values. */
                  for (i=0; i<2; i++)
                    for (ienc=0; ienc<3; ienc++)
                      encode_ints[i*3+ienc]=encode_ints[(i+1)*3+ienc];
                }
              insert_batch(input+inpdata,ntriplets_left,prevcoord,encode_ints,min_runlength,&nencode);
            }
          /* Here we should only have differences for the atom coordinates. */
          /* Convert the ints to positive ints */
          for (ienc=0; ienc<nencode; ienc++)
            {
              int pint=positive_int(encode_ints[ienc]);
              encode_ints[ienc]=pint;
            }
          /* Now we must decide what base and runlength to do. If we have swapped atoms it will be at least 2.
             If even the next atom is large, we will not do anything. */
          largest_required_base=0;
          /* Determine required base */
          for (ienc=0; ienc<min_runlength*3; ienc++)
            if (encode_ints[ienc]>largest_required_base)
              largest_required_base=encode_ints[ienc];
          /* Also compute what the largest base is for the current runlength setting! */
          largest_runlength_base=0;
          for (ienc=0; (ienc<runlength*3) && (ienc<nencode); ienc++)
            if (encode_ints[ienc]>largest_runlength_base)
              largest_runlength_base=encode_ints[ienc];

          largest_required_index=Ptngc_find_magic_index(largest_required_base);
          largest_runlength_index=Ptngc_find_magic_index(largest_runlength_base);

          if (largest_required_index<largest_runlength_index)
            {
              new_runlength=min_runlength;
              new_small_index=largest_required_index;
            }
          else
            {
              new_runlength=runlength;
              new_small_index=largest_runlength_index;
            }

          /* Only allow increase of runlength wrt min_runlength */
          if (new_runlength<min_runlength)
            new_runlength=min_runlength;

          /* If the current runlength is longer than the number of
             triplets left stop it from being so. */
          if (new_runlength>ntriplets_left)
            new_runlength=ntriplets_left;

          /* We must at least try to get some small integers going. */
          if (new_runlength==0)
            {
              new_runlength=1;
              new_small_index=small_index;
            }

          iter_runlength=new_runlength;
          iter_small_index=new_small_index;

          /* Iterate to find optimal encoding and runlength */
#ifdef SHOWIT
          fprintf(stderr,"Entering iterative loop.\n");
          fflush(stderr);
#endif

          do {
            new_runlength=iter_runlength;
            new_small_index=iter_small_index;

#ifdef SHOWIT
            fprintf(stderr,"Test new_small_index=%d Base=%d\n",new_small_index,Ptngc_magic(new_small_index));
#endif
            /* What is the largest runlength
               we can do with the currently
               selected encoding? Also the max supported runlength is MAX_SMALL_RLE triplets! */
            for (ienc=0; ienc<nencode && ienc<MAX_SMALL_RLE*3; ienc++)
              {
                int test_index=Ptngc_find_magic_index(encode_ints[ienc]);
                if (test_index>new_small_index)
                  break;
              }
            if (ienc/3>new_runlength)
              {
                iter_runlength=ienc/3;
#ifdef SHOWIT
                fprintf(stderr,"I found a new possible runlength: %d\n",iter_runlength);
#endif
              }

            /* How large encoding do we have to use? */
            largest_runlength_base=0;
            for (ienc=0; ienc<iter_runlength*3; ienc++)
              if (encode_ints[ienc]>largest_runlength_base)
                largest_runlength_base=encode_ints[ienc];
            largest_runlength_index=Ptngc_find_magic_index(largest_runlength_base);
            if (largest_runlength_index!=new_small_index)
              {
                iter_small_index=largest_runlength_index;
#ifdef SHOWIT
                fprintf(stderr,"I found a new possible small index: %d Base=%d\n",iter_small_index,Ptngc_magic(iter_small_index));
#endif
              }
          } while ((new_runlength!=iter_runlength) ||
                   (new_small_index!=iter_small_index));

#ifdef SHOWIT
          fprintf(stderr,"Exit iterative loop.\n");
          fflush(stderr);
#endif

          /* Verify that we got something good. We may have caught a
             substantially larger atom. If so we should just bail
             out and let the loop get on another lap. We may have a
             minimum runlength though and then we have to fulfill
             the request to write out these atoms! */
          rle_index_dep=0;
          if (new_runlength<3)
            rle_index_dep=IS_LARGE;
          else if (new_runlength<6)
            rle_index_dep=QUITE_LARGE;
          if ((min_runlength)
              || ((new_small_index<small_index+IS_LARGE) && (new_small_index+rle_index_dep<max_large_index))
#if 1
              || (new_small_index+IS_LARGE<max_large_index)
#endif
)
            {
              /* If doing inter-frame coding of large integers results
                 in smaller values than the small value we should not
                 produce a sequence of small values here. */
              int frame=inpdata/(natoms*3);
              int numsmaller=0;
#if 1
              if ((!swapatoms) && (frame>0))
                {
                  for (i=0; i<new_runlength; i++)
                    {
                      unsigned int delta[3];
                      unsigned int delta2[3];
                      delta[0]=positive_int(input[inpdata+i*3]-input[inpdata-natoms*3+i*3]);
                      delta[1]=positive_int(input[inpdata+i*3+1]-input[inpdata-natoms*3+i*3+1]);
                      delta[2]=positive_int(input[inpdata+i*3+2]-input[inpdata-natoms*3+i*3+2]);
                      delta2[0]=positive_int(encode_ints[i*3]);
                      delta2[1]=positive_int(encode_ints[i*3+1]);
                      delta2[2]=positive_int(encode_ints[i*3+2]);
                      if (compute_intlen(delta)*TRESHOLD_INTER_INTRA<compute_intlen(delta2))
                        numsmaller++;
                    }
                }
#endif
              /* Most of the values should become smaller, otherwise
                 we should encode them with intra coding. */
              if ((!swapatoms) && (numsmaller>=2*new_runlength/3))
                {
                  /* Put all the values in large arrays, instead of the small array */
                  if (new_runlength)
                    {
                      for (i=0; i<new_runlength; i++)
                        buffer_large(&xtc3_context,input,inpdata+i*3,natoms,1);
                      for (i=0; i<3; i++)
                        prevcoord[i]=input[inpdata+(new_runlength-1)*3+i];
                      inpdata+=3*new_runlength;
                      ntriplets_left-=new_runlength;
                    }
                }
              else
                {
                  if ((new_runlength!=runlength) || (new_small_index!=small_index))
                    {
                      int change=new_small_index-small_index;

                      if (new_small_index<=0)
                        change=0;

                      if (change<0)
                        {
                          int ixx;
                          for (ixx=0; ixx<new_runlength; ixx++)
                            {
                              int rejected;
                              do {
                                int ixyz;
                                double isum=0.; /* ints can be almost 32 bit so multiplication will overflow. So do doubles. */
                                for (ixyz=0; ixyz<3; ixyz++)
                                  {
                                    /* encode_ints is already positive (and multiplied by 2 versus the original, just as magic ints) */
                                    double id=encode_ints[ixx*3+ixyz];
                                    isum+=id*id;
                                  }
                                rejected=0;
#ifdef SHOWIT
                                fprintf(stderr,"Tested decrease %d of index: %g>=%g?\n",change,isum,(double)Ptngc_magic(small_index+change)*(double)Ptngc_magic(small_index+change));
#endif
                                if (isum>(double)Ptngc_magic(small_index+change)*(double)Ptngc_magic(small_index+change))
                                  {
#ifdef SHOWIT
                                    fprintf(stderr,"Rejected decrease %d of index due to length of vector: %g>=%g\n",change,isum,(double)Ptngc_magic(small_index+change)*(double)Ptngc_magic(small_index+change));
#endif
                                    rejected=1;
                                    change++;
                                  }
                              } while ((change<0) && (rejected));
                              if (change==0)
                                break;
                            }
                        }

                      /* Always accept the new small indices here. */
                      small_index=new_small_index;
                      /* If we have a new runlength emit it */
                      if (runlength!=new_runlength)
                        {
                          runlength=new_runlength;
                          insert_value_in_array(&xtc3_context.instructions,
                                                &xtc3_context.ninstr,
                                                &xtc3_context.ninstr_alloc,
                                                INSTR_SMALL_RUNLENGTH,"instr");
                          insert_value_in_array(&xtc3_context.rle,
                                                &xtc3_context.nrle,
                                                &xtc3_context.nrle_alloc,
                                                (unsigned int)runlength,"rle (small)");
                        }
#ifdef SHOWIT
                      fprintf(stderr,"Current small index: %d Base=%d\n",small_index,Ptngc_magic(small_index));
#endif
                    }
                  /* If we have a large previous integer we can combine it with a sequence of small ints. */
                  if (xtc3_context.has_large)
                    {
                      /* If swapatoms is set to 1 but we did actually not
                         do any swapping, we must first write out the
                         large atom and then the small. If swapatoms is 1
                         and we did swapping we can use the efficient
                         encoding. */
                      if ((swapatoms) && (!didswap))
                        {
#ifdef SHOWIT
                          fprintf(stderr,"Swapatoms was set to 1 but we did not do swapping!\n");
                          fprintf(stderr,"Only one large integer.\n");
#endif
                          /* Flush all large atoms. */
                          flush_large(&xtc3_context,xtc3_context.has_large);
#ifdef SHOWIT
                          fprintf(stderr,"Sequence of only small integers.\n");
#endif
                          insert_value_in_array(&xtc3_context.instructions,
                                                &xtc3_context.ninstr,
                                                &xtc3_context.ninstr_alloc,
                                                INSTR_ONLY_SMALL,"instr");
                        }
                      else
                        {

#ifdef SHOWIT
                          fprintf(stderr,"Sequence of one large and small integers (good compression).\n");
#endif
                          /* Flush all large atoms but one! */
                          if (xtc3_context.has_large>1)
                            flush_large(&xtc3_context,xtc3_context.has_large-1);

                          /* Here we must check if we should emit a large
                             type change instruction. */
                          large_instruction_change(&xtc3_context,0);

                          insert_value_in_array(&xtc3_context.instructions,
                                                &xtc3_context.ninstr,
                                                &xtc3_context.ninstr_alloc,
                                                INSTR_DEFAULT,"instr");

                          write_three_large(&xtc3_context,0);
                          xtc3_context.has_large=0;
                        }
                    }
                  else
                    {
#ifdef SHOWIT
                      fprintf(stderr,"Sequence of only small integers.\n");
#endif
                      insert_value_in_array(&xtc3_context.instructions,
                                            &xtc3_context.ninstr,
                                            &xtc3_context.ninstr_alloc,
                                            INSTR_ONLY_SMALL,"instr");
                    }
                  /* Insert the small integers into the small integer array. */
                  for (ienc=0; ienc<runlength*3; ienc++)
                    insert_value_in_array(&xtc3_context.smallintra,
                                          &xtc3_context.nsmallintra,
                                          &xtc3_context.nsmallintra_alloc,
                                          (unsigned int)encode_ints[ienc],"smallintra");

#ifdef SHOWIT
                  for (ienc=0; ienc<runlength; ienc++)
                    fprintf(stderr,"Small: %d %d %d\n",
                            encode_ints[ienc*3],
                            encode_ints[ienc*3+1],
                            encode_ints[ienc*3+2]);
#endif
                  /* Update prevcoord. */
                  for (ienc=0; ienc<runlength; ienc++)
                    {
#ifdef SHOWIT
                      fprintf(stderr,"Prevcoord in packing: %d %d %d\n",
                              prevcoord[0],prevcoord[1],prevcoord[2]);
#endif
                      prevcoord[0]+=unpositive_int(encode_ints[ienc*3]);
                      prevcoord[1]+=unpositive_int(encode_ints[ienc*3+1]);
                      prevcoord[2]+=unpositive_int(encode_ints[ienc*3+2]);
                    }
#ifdef SHOWIT
                  fprintf(stderr,"Prevcoord in packing: %d %d %d\n",
                          prevcoord[0],prevcoord[1],prevcoord[2]);
#endif

                  inpdata+=3*runlength;
                  ntriplets_left-=runlength;
#if 1
                }
#endif
            }
          else
            {
#ifdef SHOWIT
              fprintf(stderr,"Refused value: %d old is %d max is %d\n",new_small_index,small_index,max_large_index);
              fflush(stderr);
#endif
              refused=1;
            }
        }
#ifdef SHOWIT
      fprintf(stderr,"Number of triplets left is %d\n",ntriplets_left);
#endif
    }

  /* If we have large previous integers we must flush them now. */
  if (xtc3_context.has_large)
    flush_large(&xtc3_context,xtc3_context.has_large);

  /* Now it is time to compress all the data in the buffers with the bwlzh or base algo. */

#if 0
  /* Inspect the data. */
  printarray(xtc3_context.instructions,xtc3_context.ninstr,"A instr");
  printarray(xtc3_context.rle,xtc3_context.nrle,"A rle");
  printarray(xtc3_context.large_direct,xtc3_context.nlargedir,"A largedir");
  printarray(xtc3_context.large_intra_delta,xtc3_context.nlargeintra,"A largeintra");
  printarray(xtc3_context.large_inter_delta,xtc3_context.nlargeinter,"A largeinter");
  printarray(xtc3_context.smallintra,xtc3_context.nsmallintra,"A smallintra");
  exit(0);
#endif

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"instructions: %d\n",xtc3_context.ninstr);
#endif

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
#define bwlzh_compress bwlzh_compress_verbose
#define bwlzh_compress_no_lz77 bwlzh_compress_no_lz77_verbose
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.ninstr);
  if (xtc3_context.ninstr)
    {
      bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.ninstr));
      if (speed>=5)
        bwlzh_compress(xtc3_context.instructions,xtc3_context.ninstr,bwlzh_buf,&bwlzh_buf_len);
      else
        bwlzh_compress_no_lz77(xtc3_context.instructions,xtc3_context.ninstr,bwlzh_buf,&bwlzh_buf_len);
      output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
      memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
      outdata+=bwlzh_buf_len;
      free(bwlzh_buf);
    }

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"rle: %d\n",xtc3_context.nrle);
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.nrle);
  if (xtc3_context.nrle)
    {
      bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.nrle));
      if (speed>=5)
        bwlzh_compress(xtc3_context.rle,xtc3_context.nrle,bwlzh_buf,&bwlzh_buf_len);
      else
        bwlzh_compress_no_lz77(xtc3_context.rle,xtc3_context.nrle,bwlzh_buf,&bwlzh_buf_len);
      output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
      memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
      outdata+=bwlzh_buf_len;
      free(bwlzh_buf);
    }

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"large direct: %d\n",xtc3_context.nlargedir);
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.nlargedir);
  if (xtc3_context.nlargedir)
    {
      if ((speed<=2) || ((speed<=5) && (!heuristic_bwlzh(xtc3_context.large_direct,xtc3_context.nlargedir))))
        {
          bwlzh_buf=NULL;
          bwlzh_buf_len=INT_MAX;
        }
      else
        {
          bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.nlargedir));
          if (speed>=5)
            bwlzh_compress(xtc3_context.large_direct,xtc3_context.nlargedir,bwlzh_buf,&bwlzh_buf_len);
          else
            bwlzh_compress_no_lz77(xtc3_context.large_direct,xtc3_context.nlargedir,bwlzh_buf,&bwlzh_buf_len);
        }
      /* If this can be written smaller using base compression we should do that. */
      base_buf=warnmalloc((xtc3_context.nlargedir+3)*sizeof(int));
      base_compress(xtc3_context.large_direct,xtc3_context.nlargedir,base_buf,&base_buf_len);
#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
      fprintf(stderr,"Large direct: Base len=%d. BWLZH len=%d\n",base_buf_len,bwlzh_buf_len);
#endif
      if (base_buf_len<bwlzh_buf_len)
        {
          output[outdata++]=0U;
          output_int(output,&outdata,(unsigned int)base_buf_len);
          memcpy(output+outdata,base_buf,base_buf_len);
          outdata+=base_buf_len;
        }
      else
        {
          output[outdata++]=1U;
          output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
          memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
          outdata+=bwlzh_buf_len;
        }
      free(bwlzh_buf);
      free(base_buf);
    }

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"large intra: %d\n",xtc3_context.nlargeintra);
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.nlargeintra);
  if (xtc3_context.nlargeintra)
    {
      if ((speed<=2) || ((speed<=5) && (!heuristic_bwlzh(xtc3_context.large_intra_delta,xtc3_context.nlargeintra))))
        {
          bwlzh_buf=NULL;
          bwlzh_buf_len=INT_MAX;
        }
      else
        {
          bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.nlargeintra));
          if (speed>=5)
            bwlzh_compress(xtc3_context.large_intra_delta,xtc3_context.nlargeintra,bwlzh_buf,&bwlzh_buf_len);
          else
            bwlzh_compress_no_lz77(xtc3_context.large_intra_delta,xtc3_context.nlargeintra,bwlzh_buf,&bwlzh_buf_len);
        }
      /* If this can be written smaller using base compression we should do that. */
      base_buf=warnmalloc((xtc3_context.nlargeintra+3)*sizeof(int));
      base_compress(xtc3_context.large_intra_delta,xtc3_context.nlargeintra,base_buf,&base_buf_len);
#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
      fprintf(stderr,"Large intra: Base len=%d. BWLZH len=%d\n",base_buf_len,bwlzh_buf_len);
#endif
      if (base_buf_len<bwlzh_buf_len)
        {
          output[outdata++]=0U;
          output_int(output,&outdata,(unsigned int)base_buf_len);
          memcpy(output+outdata,base_buf,base_buf_len);
          outdata+=base_buf_len;
        }
      else
        {
          output[outdata++]=1U;
          output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
          memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
          outdata+=bwlzh_buf_len;
        }
      free(bwlzh_buf);
      free(base_buf);
    }

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"large inter: %d\n",xtc3_context.nlargeinter);
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.nlargeinter);
  if (xtc3_context.nlargeinter)
    {
      if ((speed<=2) || ((speed<=5) && (!heuristic_bwlzh(xtc3_context.large_inter_delta,xtc3_context.nlargeinter))))
        {
          bwlzh_buf=NULL;
          bwlzh_buf_len=INT_MAX;
        }
      else
        {
          bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.nlargeinter));
          if (speed>=5)
            bwlzh_compress(xtc3_context.large_inter_delta,xtc3_context.nlargeinter,bwlzh_buf,&bwlzh_buf_len);
          else
            bwlzh_compress_no_lz77(xtc3_context.large_inter_delta,xtc3_context.nlargeinter,bwlzh_buf,&bwlzh_buf_len);
        }
      /* If this can be written smaller using base compression we should do that. */
      base_buf=warnmalloc((xtc3_context.nlargeinter+3)*sizeof(int));
      base_compress(xtc3_context.large_inter_delta,xtc3_context.nlargeinter,base_buf,&base_buf_len);
#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
      fprintf(stderr,"Large inter: Base len=%d. BWLZH len=%d\n",base_buf_len,bwlzh_buf_len);
#endif
      if (base_buf_len<bwlzh_buf_len)
        {
          output[outdata++]=0U;
          output_int(output,&outdata,(unsigned int)base_buf_len);
          memcpy(output+outdata,base_buf,base_buf_len);
          outdata+=base_buf_len;
        }
      else
        {
          output[outdata++]=1U;
          output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
          memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
          outdata+=bwlzh_buf_len;
        }
      free(bwlzh_buf);
      free(base_buf);
    }

#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
  fprintf(stderr,"small intra: %d\n",xtc3_context.nsmallintra);
#endif

  output_int(output,&outdata,(unsigned int)xtc3_context.nsmallintra);
  if (xtc3_context.nsmallintra)
    {
      if ((speed<=2) || ((speed<=5) && (!heuristic_bwlzh(xtc3_context.smallintra,xtc3_context.nsmallintra))))
        {
          bwlzh_buf=NULL;
          bwlzh_buf_len=INT_MAX;
        }
      else
        {
          bwlzh_buf=warnmalloc(bwlzh_get_buflen(xtc3_context.nsmallintra));
          if (speed>=5)
            bwlzh_compress(xtc3_context.smallintra,xtc3_context.nsmallintra,bwlzh_buf,&bwlzh_buf_len);
          else
            bwlzh_compress_no_lz77(xtc3_context.smallintra,xtc3_context.nsmallintra,bwlzh_buf,&bwlzh_buf_len);
        }
      /* If this can be written smaller using base compression we should do that. */
      base_buf=warnmalloc((xtc3_context.nsmallintra+3)*sizeof(int));
      base_compress(xtc3_context.smallintra,xtc3_context.nsmallintra,base_buf,&base_buf_len);
#if defined(SHOWIT) || defined(SHOWIT_LIGHT)
      fprintf(stderr,"Small intra: Base len=%d. BWLZH len=%d\n",base_buf_len,bwlzh_buf_len);
#endif
      if (base_buf_len<bwlzh_buf_len)
        {
          output[outdata++]=0U;
          output_int(output,&outdata,(unsigned int)base_buf_len);
          memcpy(output+outdata,base_buf,base_buf_len);
          outdata+=base_buf_len;
        }
      else
        {
          output[outdata++]=1U;
          output_int(output,&outdata,(unsigned int)bwlzh_buf_len);
          memcpy(output+outdata,bwlzh_buf,bwlzh_buf_len);
          outdata+=bwlzh_buf_len;
        }
      free(bwlzh_buf);
      free(base_buf);
    }
  *length=outdata;

  free_xtc3_context(&xtc3_context);
  return output;
}

static void decompress_bwlzh_block(unsigned char **ptr,
                                   const int nvals,
                                   unsigned int **vals)
{
  int bwlzh_buf_len=(int)(((unsigned int)(*ptr)[0]) |
                          (((unsigned int)(*ptr)[1])<<8) |
                          (((unsigned int)(*ptr)[2])<<16) |
                          (((unsigned int)(*ptr)[3])<<24));
  (*ptr)+=4;
  *vals=warnmalloc(nvals*sizeof (**vals));
  bwlzh_decompress(*ptr,nvals,*vals);
  (*ptr)+=bwlzh_buf_len;
}

static void decompress_base_block(unsigned char **ptr,
                                  const int nvals,
                                  unsigned int **vals)
{
  int base_buf_len=(int)(((unsigned int)(*ptr)[0]) |
                         (((unsigned int)(*ptr)[1])<<8) |
                         (((unsigned int)(*ptr)[2])<<16) |
                         (((unsigned int)(*ptr)[3])<<24));
  (*ptr)+=4;
  *vals=warnmalloc(nvals*sizeof (**vals));
  base_decompress(*ptr,nvals,*vals);
  (*ptr)+=base_buf_len;
}

static void unpack_one_large(struct xtc3_context *xtc3_context,
                             int *ilargedir, int *ilargeintra,
                             int *ilargeinter, int *prevcoord,
                             int *minint, int *output,
                             const int outdata, const int didswap,
                             const int natoms, const int current_large_type)
{
  int large_ints[3]={0,0,0};
  if (current_large_type==0 && xtc3_context->large_direct)
    {
      large_ints[0]=(int)xtc3_context->large_direct[(*ilargedir)]+minint[0];
      large_ints[1]=(int)xtc3_context->large_direct[(*ilargedir)+1]+minint[1];
      large_ints[2]=(int)xtc3_context->large_direct[(*ilargedir)+2]+minint[2];
      (*ilargedir)+=3;
    }
  else if (current_large_type==1 && xtc3_context->large_intra_delta)
    {
      large_ints[0]=unpositive_int(xtc3_context->large_intra_delta[(*ilargeintra)])+prevcoord[0];
      large_ints[1]=unpositive_int(xtc3_context->large_intra_delta[(*ilargeintra)+1])+prevcoord[1];
      large_ints[2]=unpositive_int(xtc3_context->large_intra_delta[(*ilargeintra)+2])+prevcoord[2];
      (*ilargeintra)+=3;
    }
  else if (xtc3_context->large_inter_delta)
    {
      large_ints[0]=unpositive_int(xtc3_context->large_inter_delta[(*ilargeinter)])
        +output[outdata-natoms*3+didswap*3];
      large_ints[1]=unpositive_int(xtc3_context->large_inter_delta[(*ilargeinter)+1])
        +output[outdata-natoms*3+1+didswap*3];
      large_ints[2]=unpositive_int(xtc3_context->large_inter_delta[(*ilargeinter)+2])
        +output[outdata-natoms*3+2+didswap*3];
      (*ilargeinter)+=3;
    }
  memcpy(prevcoord, large_ints, 3*sizeof *prevcoord);
  output[outdata]=large_ints[0];
  output[outdata+1]=large_ints[1];
  output[outdata+2]=large_ints[2];
#ifdef SHOWIT
  fprintf(stderr,"Unpack one large: %d %d %d\n",prevcoord[0],prevcoord[1],prevcoord[2]);
#endif
}


int Ptngc_unpack_array_xtc3(unsigned char *packed,int *output, const int length, const int natoms)
{
  int i;
  int minint[3];
  unsigned char *ptr=packed;
  int prevcoord[3];
  int outdata=0;
  int ntriplets_left=length/3;
  int swapatoms=0;
  int runlength=0;
  int current_large_type=0;
  int iinstr=0;
  int irle=0;
  int ilargedir=0;
  int ilargeintra=0;
  int ilargeinter=0;
  int ismallintra=0;

  struct xtc3_context xtc3_context;
  init_xtc3_context(&xtc3_context);

  for (i=0; i<3; i++)
    {
      minint[i]=unpositive_int((int)(((unsigned int)ptr[0]) |
                                     (((unsigned int)ptr[1])<<8) |
                                     (((unsigned int)ptr[2])<<16) |
                                     (((unsigned int)ptr[3])<<24)));
      ptr+=4;
    }

  xtc3_context.ninstr=(int)(((unsigned int)ptr[0]) |
                            (((unsigned int)ptr[1])<<8) |
                            (((unsigned int)ptr[2])<<16) |
                            (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.ninstr)
    decompress_bwlzh_block(&ptr,xtc3_context.ninstr,&xtc3_context.instructions);

  xtc3_context.nrle=(int)(((unsigned int)ptr[0]) |
                          (((unsigned int)ptr[1])<<8) |
                          (((unsigned int)ptr[2])<<16) |
                          (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.nrle)
    decompress_bwlzh_block(&ptr,xtc3_context.nrle,&xtc3_context.rle);

  xtc3_context.nlargedir=(int)(((unsigned int)ptr[0]) |
                               (((unsigned int)ptr[1])<<8) |
                               (((unsigned int)ptr[2])<<16) |
                               (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.nlargedir)
    {
      if (*ptr++==1)
        decompress_bwlzh_block(&ptr,xtc3_context.nlargedir,&xtc3_context.large_direct);
      else
        decompress_base_block(&ptr,xtc3_context.nlargedir,&xtc3_context.large_direct);
    }

  xtc3_context.nlargeintra=(int)(((unsigned int)ptr[0]) |
                                 (((unsigned int)ptr[1])<<8) |
                                 (((unsigned int)ptr[2])<<16) |
                                 (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.nlargeintra)
    {
      if (*ptr++==1)
        decompress_bwlzh_block(&ptr,xtc3_context.nlargeintra,&xtc3_context.large_intra_delta);
      else
        decompress_base_block(&ptr,xtc3_context.nlargeintra,&xtc3_context.large_intra_delta);
    }

  xtc3_context.nlargeinter=(int)(((unsigned int)ptr[0]) |
                                 (((unsigned int)ptr[1])<<8) |
                                 (((unsigned int)ptr[2])<<16) |
                                 (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.nlargeinter)
    {
      if (*ptr++==1)
        decompress_bwlzh_block(&ptr,xtc3_context.nlargeinter,&xtc3_context.large_inter_delta);
      else
        decompress_base_block(&ptr,xtc3_context.nlargeinter,&xtc3_context.large_inter_delta);
    }

  xtc3_context.nsmallintra=(int)(((unsigned int)ptr[0]) |
                                 (((unsigned int)ptr[1])<<8) |
                                 (((unsigned int)ptr[2])<<16) |
                                 (((unsigned int)ptr[3])<<24));
  ptr+=4;
  if (xtc3_context.nsmallintra)
    {
      if (*ptr++==1)
        decompress_bwlzh_block(&ptr,xtc3_context.nsmallintra,&xtc3_context.smallintra);
      else
        decompress_base_block(&ptr,xtc3_context.nsmallintra,&xtc3_context.smallintra);
    }

  /* Initial prevcoord is the minimum integers. */
  memcpy(prevcoord, minint, 3*sizeof *prevcoord);

  while (ntriplets_left>0 && iinstr<xtc3_context.ninstr)
    {
      int instr=xtc3_context.instructions[iinstr++];
#ifdef SHOWIT
      fprintf(stderr,"instr=%d @ %d\n",instr,iinstr-1);
#endif
#ifdef SHOWIT
      fprintf(stderr,"ntriplets left=%d\n",ntriplets_left);
#endif
      if ((instr==INSTR_DEFAULT) /* large+small */
          || (instr==INSTR_ONLY_LARGE) /* only large */
          || (instr==INSTR_ONLY_SMALL)) /* only small */
        {
          if (instr!=INSTR_ONLY_SMALL)
            {
              int didswap=0;
              if ((instr==INSTR_DEFAULT) && (swapatoms))
                didswap=1;
              unpack_one_large(&xtc3_context,&ilargedir, &ilargeintra, &ilargeinter,
                               prevcoord, minint, output, outdata, didswap,
                               natoms, current_large_type);
              ntriplets_left--;
              outdata+=3;
            }
          if (instr!=INSTR_ONLY_LARGE)
            {
              for (i=0; i<runlength; i++)
                {
                  prevcoord[0]+=unpositive_int(xtc3_context.smallintra[ismallintra]);
                  prevcoord[1]+=unpositive_int(xtc3_context.smallintra[ismallintra+1]);
                  prevcoord[2]+=unpositive_int(xtc3_context.smallintra[ismallintra+2]);
                  ismallintra+=3;
                  output[outdata+i*3]=prevcoord[0];
                  output[outdata+i*3+1]=prevcoord[1];
                  output[outdata+i*3+2]=prevcoord[2];
#ifdef SHOWIT
                  fprintf(stderr,"Unpack small: %d %d %d\n",prevcoord[0],prevcoord[1],prevcoord[2]);
#endif
                }
              if ((instr==INSTR_DEFAULT) && (swapatoms))
                {
                  for (i=0; i<3; i++)
                    {
                      int tmp=output[outdata-3+i];
                      output[outdata-3+i]=output[outdata+i];
                      output[outdata+i]=tmp;
                    }
#ifdef SHOWIT
                  fprintf(stderr,"Unswap results in\n");
                  for (i=0; i<3; i++)
                    fprintf(stderr,"%d %d %d\n",output[outdata-3+i*3],output[outdata-2+i*3],output[outdata-1+i*3]);
#endif
                }
              ntriplets_left-=runlength;
              outdata+=runlength*3;
            }
        }
      else if (instr==INSTR_LARGE_RLE && irle<xtc3_context.nrle)
        {
          int large_rle=xtc3_context.rle[irle++];
#ifdef SHOWIT
          fprintf(stderr,"large_rle=%d @ %d\n",large_rle,irle-1);
#endif
          for (i=0; i<large_rle; i++)
            {
              unpack_one_large(&xtc3_context,&ilargedir, &ilargeintra, &ilargeinter,
                               prevcoord, minint, output, outdata, 0,
                               natoms, current_large_type);
              ntriplets_left--;
              outdata+=3;
            }
        }
      else if (instr==INSTR_SMALL_RUNLENGTH && irle<xtc3_context.nrle)
        {
          runlength=xtc3_context.rle[irle++];
#ifdef SHOWIT
          fprintf(stderr,"small_rle=%d @ %d\n",runlength,irle-1);
#endif
        }
      else if (instr==INSTR_FLIP)
        {
          swapatoms=1-swapatoms;
#ifdef SHOWIT
          fprintf(stderr,"new flip=%d\n",swapatoms);
#endif
        }
      else if (instr==INSTR_LARGE_DIRECT)
        {
          current_large_type=0;
#ifdef SHOWIT
          fprintf(stderr,"large direct\n");
#endif
        }
      else if (instr==INSTR_LARGE_INTRA_DELTA)
        {
          current_large_type=1;
#ifdef SHOWIT
          fprintf(stderr,"large intra delta\n");
#endif
        }
      else if (instr==INSTR_LARGE_INTER_DELTA)
        {
          current_large_type=2;
#ifdef SHOWIT
          fprintf(stderr,"large inter delta\n");
#endif
        }
    }
  if (ntriplets_left<0)
    {
      fprintf(stderr,"TRAJNG XTC3: A bug has been found. At end ntriplets_left<0\n");
      exit(EXIT_FAILURE);
    }
  free_xtc3_context(&xtc3_context);
  return 0;
}
