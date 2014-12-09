/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg and Magnus Lundborg
 * Copyright (c) 2010, 2013-2014 The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <stdlib.h>
#include <string.h>
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/mtf.h"

/* "Partial" MTF. Byte based. */
/* Move to front coding.
   Acceptable inputs are max 8 bits (0-0xFF) */
static void comp_conv_to_mtf_byte(unsigned char *vals, const int nvals,
                                  unsigned char *valsmtf)
{
  int i;
  /* Indices into a linked list */
  int list[256];
  int dict[256];
  /* Head of the linked list */
  int head;
  for (i=0; i<256; i++)
    dict[i]=i;
  for (i=0; i<255; i++)
    list[i]=i+1;
  list[255]=-1; /* end. */
  head=0;
  for (i=0; i<nvals; i++)
    {
      int v=(int)vals[i];
      /* Find how early in the dict the value is */
      int ptr=head;
      int oldptr=-1;
      int r=0;
      while (dict[ptr]!=v)
        {
          oldptr=ptr;
          ptr=list[ptr];
          r++;
        }
      valsmtf[i]=(unsigned char)r;
      /* Move it to front in list */
      /* Is it the head? Then it is already at the front. */
      if (oldptr!=-1)
        {
          /* Remove it from inside the list */
          list[oldptr]=list[ptr];
          /* Move it to the front. */
          list[ptr]=head;
          head=ptr;
        }
    }
}

void Ptngc_comp_conv_to_mtf_partial(unsigned int *vals, const int nvals,
                              unsigned int *valsmtf)
{
  unsigned char *tmp=warnmalloc(nvals*2);
  int i, j;

  memset(valsmtf, 0U, sizeof(unsigned int) * nvals);

  for (j=0; j<3; j++)
    {
      for (i=0; i<nvals; i++)
        tmp[i]=(unsigned char)((vals[i]>>(8*j))&0xFF);
      comp_conv_to_mtf_byte(tmp,nvals,tmp+nvals);
      for (i=0; i<nvals; i++)
        valsmtf[i]|=(((unsigned int)(tmp[nvals+i]))<<(8*j));
    }
  free(tmp);
}

void Ptngc_comp_conv_to_mtf_partial3(unsigned int *vals, const int nvals,
                               unsigned char *valsmtf)
{
  unsigned char *tmp=warnmalloc(nvals);
  int i, j;
  for (j=0; j<3; j++)
    {
      for (i=0; i<nvals; i++)
        tmp[i]=(unsigned char)((vals[i]>>(8*j))&0xFF);
      comp_conv_to_mtf_byte(tmp,nvals,valsmtf+j*nvals);
    }
  free(tmp);
}

/* Move to front decoding */
static void comp_conv_from_mtf_byte(unsigned char *valsmtf, const int nvals,
                             unsigned char *vals)
{
  int i;
  /* Indices into a linked list */
  int list[256];
  int dict[256];
  /* Head of the linked list */
  int head;
  for (i=0; i<256; i++)
    dict[i]=i;
  for (i=0; i<255; i++)
    list[i]=i+1;
  list[255]=-1; /* end. */
  head=0;
  for (i=0; i<nvals; i++)
    {
      int r=(int)valsmtf[i];
      /* Find value at position r in the list */
      int ptr=head;
      int oldptr=-1;
      int cnt=0;
      while (cnt<r)
        {
          oldptr=ptr;
          ptr=list[ptr];
          cnt++;
        }
      vals[i]=(unsigned int)dict[ptr];
      /* Move it to front in list */
      /* Is it the head? Then it is already at the front. */
      if (oldptr!=-1)
        {
          /* Remove it from inside the list */
          list[oldptr]=list[ptr];
          /* Move it to the front. */
          list[ptr]=head;
          head=ptr;
        }
    }
}

void Ptngc_comp_conv_from_mtf_partial(unsigned int *valsmtf, const int nvals,
                                unsigned int *vals)
{
  unsigned char *tmp=warnmalloc(nvals*2);
  int i, j;

  memset(vals, 0U, sizeof(unsigned int) * nvals);

  for (j=0; j<3; j++)
    {
      for (i=0; i<nvals; i++)
        tmp[i]=(unsigned char)((valsmtf[i]>>(8*j))&0xFF);
      comp_conv_from_mtf_byte(tmp,nvals,tmp+nvals);
      for (i=0; i<nvals; i++)
        vals[i]|=(((unsigned int)(tmp[nvals+i]))<<(8*j));
    }
  free(tmp);
}

void Ptngc_comp_conv_from_mtf_partial3(unsigned char *valsmtf, const int nvals,
                                 unsigned int *vals)
{
  unsigned char *tmp=warnmalloc(nvals);
  int i, j;

  memset(vals, 0U, sizeof(unsigned int) * nvals);

  for (j=0; j<3; j++)
    {
      comp_conv_from_mtf_byte(valsmtf+j*nvals,nvals,tmp);
      for (i=0; i<nvals; i++)
        vals[i]|=(((unsigned int)(tmp[i]))<<(8*j));
    }
  free(tmp);
}

/* Move to front coding.
   Acceptable inputs are max 24 bits (0-0xFFFFFF) */
void Ptngc_comp_conv_to_mtf(unsigned int *vals, const int nvals,
                      unsigned int *dict, const int ndict,
                      unsigned int *valsmtf)
{
  int i;
  /* Indices into a linked list */
  int *list=warnmalloc(ndict*sizeof *list);
  /* Head of the linked list */
  int head;
  for (i=0; i<ndict-1; i++)
    list[i]=i+1;
  list[ndict-1]=-1; /* end. */
  head=0;
  for (i=0; i<nvals; i++)
    {
      int v=vals[i];
      /* Find how early in the dict the value is */
      int ptr=head;
      int oldptr=-1;
      int r=0;
      while (dict[ptr]!=v)
        {
          oldptr=ptr;
          ptr=list[ptr];
          r++;
        }
      valsmtf[i]=r;
      /* Move it to front in list */
      /* Is it the head? Then it is already at the front. */
      if (oldptr!=-1)
        {
          /* Remove it from inside the list */
          list[oldptr]=list[ptr];
          /* Move it to the front. */
          list[ptr]=head;
          head=ptr;
        }
    }
  free(list);
}

/* Move to front decoding */
void Ptngc_comp_conv_from_mtf(unsigned int *valsmtf, const int nvals,
                        unsigned int *dict, const int ndict,
                        unsigned int *vals)
{
  int i;
  /* Indices into a linked list */
  int *list=warnmalloc(ndict*sizeof *list);
  /* Head of the linked list */
  int head;
  for (i=0; i<ndict-1; i++)
    list[i]=i+1;
  list[ndict-1]=-1; /* end. */
  head=0;
  for (i=0; i<nvals; i++)
    {
      int r=valsmtf[i];
      /* Find value at position r in the list */
      int ptr=head;
      int oldptr=-1;
      int cnt=0;
      while (cnt<r)
        {
          oldptr=ptr;
          ptr=list[ptr];
          cnt++;
        }
      vals[i]=dict[ptr];
      /* Move it to front in list */
      /* Is it the head? Then it is already at the front. */
      if (oldptr!=-1)
        {
          /* Remove it from inside the list */
          list[oldptr]=list[ptr];
          /* Move it to the front. */
          list[ptr]=head;
          head=ptr;
        }
    }
  free(list);
}

