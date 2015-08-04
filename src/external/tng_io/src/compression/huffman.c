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
#include "../../include/compression/huffman.h"

#define MAX_HUFFMAN_LEN 31

enum htree_type { htree_leaf, htree_node };

struct htree_leaf
{
  enum htree_type nodeleaf;
  unsigned int idict; /* Index into input dictionary */
  unsigned int prob;
  unsigned int bit; /* One or zero */
};

struct htree_node
{
  enum htree_type nodeleaf;
  union htree_nodeleaf *n1;
  union htree_nodeleaf *n2;
  unsigned int bit; /* One or zero */
  unsigned int prob;
};

union htree_nodeleaf
{
  enum htree_type nodeleaf;
  struct htree_node node;
  struct htree_leaf leaf;
};

struct codelength
{
  unsigned int code;
  int length;
  unsigned int dict;
  unsigned int prob;
};

static int comp_htree(const void *leafptr1, const void *leafptr2, const void *private)
{
  const union htree_nodeleaf *leaf1=(union htree_nodeleaf *)leafptr1;
  const union htree_nodeleaf *leaf2=(union htree_nodeleaf *)leafptr2;
  int rval=0;
  (void)private;

  if (leaf1->leaf.prob<leaf2->leaf.prob)
    rval=1;
  else if (leaf1->leaf.prob>leaf2->leaf.prob)
    rval=-1;
  return rval;
}

static void assign_codes(union htree_nodeleaf *htree,
                         struct codelength *codelength,
                         unsigned int code,
                         int length,
                         const int top)
{
#if 0
  printf("Assign codes called with code %d length %d\n",code,length);
#endif
  if (htree->nodeleaf==htree_leaf)
    {
      codelength[htree->leaf.idict].length=length+1;
      codelength[htree->leaf.idict].code=(code<<1)|htree->leaf.bit;
#if 0
      printf("I am a leaf: %d %d\n",
             codelength[htree->leaf.idict].length,
             codelength[htree->leaf.idict].code);
#endif
    }
  else
    {
      if (!top)
        {
          code<<=1;
          code|=htree->node.bit;
          length++;
        }
#if 0
      printf("I am a node length: %d\n",length);
      printf("I am a node code: %d\n",code);
#endif
      assign_codes(htree->node.n1,codelength,code,length,0);
      assign_codes(htree->node.n2,codelength,code,length,0);
    }
}

static void free_nodes(union htree_nodeleaf *htree, int top)
{
  if (htree->nodeleaf==htree_leaf)
    {
      if (!top)
        free(htree);
    }
  else
    {
      free_nodes(htree->node.n1,0);
      free_nodes(htree->node.n2,0);
      if (!top)
        free(htree);
    }
}

static void flush_8bits(unsigned int *combine, unsigned char **output, int *bitptr)
{
  while ((*bitptr)>=8)
    {
      unsigned int mask=~(0xFFU<<((*bitptr)-8));
      unsigned char out=(unsigned char)((*combine)>>((*bitptr)-8));
      **output=out;
      (*output)++;
      (*bitptr)-=8;
      (*combine)&=mask;
    }
}

static void writebits(unsigned int value, int length, unsigned char **output, int *bitptr)
{
  unsigned int mask;
  unsigned int combine=(unsigned int)**output;
  if (length>=8)
    mask=0xFFU<<(length-8);
  else
    mask=0xFFU>>(8-length);
  while (length>8)
    {
      /* Make room for the bits. */
      combine<<=8;
      (*bitptr)+=8;
      combine|=(value&mask)>>(length-8);
      flush_8bits(&combine,output,bitptr);
      length-=8;
      mask>>=8;
    }
  if (length)
    {
      /* Make room for the bits. */
      combine<<=length;
      (*bitptr)+=length;
      combine|=value;
      flush_8bits(&combine,output,bitptr);
    }
  **output=(unsigned char)combine;
}

static unsigned int readbits(int length, unsigned char **input, int *bitptr)
{
  unsigned int val=0U;
  unsigned int extract_mask=0x80U>>*bitptr;
  unsigned char thisval=**input;
  while (length--)
    {
      val<<=1;
      val|=((extract_mask & thisval)!=0);
      *bitptr=(*bitptr)+1;
      extract_mask>>=1;
      if (!extract_mask)
        {
          extract_mask=0x80U;
          *input=(*input)+1;
          *bitptr=0;
          thisval=**input;
        }
    }
  return val;
}

static int comp_codes(const void *codeptr1, const void *codeptr2, const void *private)
{
  const struct codelength *code1=(struct codelength *)codeptr1;
  const struct codelength *code2=(struct codelength *)codeptr2;
  int rval=0; /* It shouldn't be possible to get equal here, though. */
  (void)private;
  if (code1->length>code2->length)
    rval=1;
  else if (code1->length<code2->length)
    rval=-1;
  else if (code1->dict>code2->dict)
    rval=1;
  else
    rval=-1;
  return rval;
}

static int comp_codes_value(const void *codeptr1, const void *codeptr2, const void *private)
{
  const struct codelength *code1=(struct codelength *)codeptr1;
  const struct codelength *code2=(struct codelength *)codeptr2;

  int rval=0; /* It shouldn't be possible to get equal here, though. */
  (void)private;
  if (code1->dict>code2->dict)
    rval=1;
  else
    rval=-1;
  return rval;
}

/* The huffman_dict array should be 131077 (0x20005) long. The
huffman_dict_unpacked array should be 131077 long (note five longer than
0x20000) */
void Ptngc_comp_conv_to_huffman(unsigned int *vals, const int nvals,
                          unsigned int *dict, const int ndict,
                          unsigned int *prob,
                          unsigned char *huffman,
                          int *huffman_len,
                          unsigned char *huffman_dict,
                          int *huffman_dictlen,
                          unsigned int *huffman_dict_unpacked,
                          int *huffman_dict_unpackedlen)
{
  int i;
  int nleft;
  union htree_nodeleaf *htree;
  struct codelength *codelength;
  int bitptr;
  unsigned char *huffman_ptr;
  int code;
  int longcodes=1;
  while (longcodes)
    {
      /* Create array of leafs (will be array of nodes/trees during
         buildup of tree. */
      htree=warnmalloc(ndict*sizeof *htree);
      codelength=warnmalloc(ndict*sizeof *codelength);
      bitptr=0;
      huffman_ptr=huffman;
      for (i=0; i<ndict; i++)
        {
          htree[i].nodeleaf=htree_leaf;
          htree[i].leaf.idict=i;
          htree[i].leaf.prob=prob[i];
        }
      /* Sort the leafs wrt probability. */
      Ptngc_merge_sort(htree,ndict,sizeof *htree,comp_htree,NULL);

#if 0
      for (i=0; i<ndict; i++)
        {
          printf("%d %d\n",dict[htree[i].leaf.idict],htree[i].leaf.prob);
        }
#endif

      /* Build tree. */
      if (ndict==1)
        {
          codelength[0].code=1;
          codelength[0].length=1;
        }
      else
        {
          /* Nodes and leafs left. */
          nleft=ndict;

          /* Take the two least probable symbols (which are at the end of the
             array and combine them until there is nothing left. */
          while (nleft>1)
            {
              union htree_nodeleaf *n1=warnmalloc(sizeof *n1);
              union htree_nodeleaf *n2=warnmalloc(sizeof *n2);
              int new_place;
              int p1,p2, new_prob;
              *n1=htree[nleft-1];
              *n2=htree[nleft-2];
              if (n1->nodeleaf==htree_leaf)
                {
                  p1=n1->leaf.prob;
                  n1->leaf.bit=0;
                }
              else
                {
                  p1=n1->node.prob;
                  n1->node.bit=0;
                }
              if (n2->nodeleaf==htree_leaf)
                {
                  p2=n2->leaf.prob;
                  n2->leaf.bit=1;
                }
              else
                {
                  p2=n2->node.prob;
                  n2->node.bit=1;
                }
              nleft--;
              /* Create a new node */
              htree[nleft-1].nodeleaf=htree_node;
              htree[nleft-1].node.n1=n1;
              htree[nleft-1].node.n2=n2;
              new_prob=p1+p2;
              htree[nleft-1].node.prob=new_prob;
              /* Use insertion sort to place this in the correct place in the
                 array. */
              /* Where should it be inserted? */
              new_place=nleft;
              while (new_place>0)
                {
                  int pc;
                  if (htree[new_place-1].nodeleaf==htree_node)
                    pc=htree[new_place-1].node.prob;
                  else
                    pc=htree[new_place-1].leaf.prob;
                  if (new_prob<pc)
                    break;
                  else
                    new_place--;
                }
              if (new_place!=nleft)
                {
                  /* Shift array (overlapping regions!) */
                  union htree_nodeleaf nodecopy=htree[nleft-1];
                  memmove(htree+new_place+1,
                          htree+new_place,
                          (nleft-1-new_place)*sizeof *htree);
                  htree[new_place]=nodecopy;
                }
            }
        }
      /* Create codes from tree */
      assign_codes(htree,codelength,0,0,1);
      /* Canonicalize */
      /* First put values into to the codelength array for sorting. */
      for (i=0; i<ndict; i++)
        {
          codelength[i].dict=dict[i];
          codelength[i].prob=prob[i];
        }
      /* Sort codes wrt length/value */
      Ptngc_merge_sort(codelength,ndict,sizeof *codelength,comp_codes,NULL);
      /* Canonicalize codes. */
      code=0;
      for (i=0; i<ndict; i++)
        {
          codelength[i].code=code;
          if (i<ndict-1)
            code=(code+1)<<(codelength[i+1].length-codelength[i].length);
        }
      /* Free all nodes / leaves. */
      free_nodes(htree,1);
      /* Free array that held nodes/leaves. */
      free(htree);

      longcodes=0;
      /* Check if there is a too long huffman code. */
      for (i=0; i<ndict; i++)
        if (codelength[i].length>MAX_HUFFMAN_LEN)
          longcodes=1;

      /* If the codes are too long alter the probabilities. */
      if (longcodes)
        {
          for (i=0; i<ndict; i++)
            {
              prob[i]>>=1;
              if (prob[i]==0)
                prob[i]=1;
            }

          /* Free codelength. We will compute a new one. */
          free(codelength);
        }
    }

#if 0
  {
    for (i=0; i<ndict; i++)
      {
        printf("%d %d\t\t %d %d ",codelength[i].dict,codelength[i].prob,codelength[i].length,codelength[i].code);
        {
          unsigned int c=codelength[i].code;
          int j;
          unsigned int mask=1<<(codelength[i].length-1);
          for (j=codelength[i].length-1; j>=0; j--)
            {
              int bit=c&mask;
              if (bit)
                printf("1");
              else
                printf("0");
              mask>>=1;
            }
          printf("\n");
        }
      }
  }
#endif

  /* Simply do compression by writing out the bits. */
  for (i=0; i<nvals; i++)
    {
      int r;
      for (r=0; r<ndict; r++)
        if (codelength[r].dict==vals[i])
          break;
      writebits(codelength[r].code,codelength[r].length,&huffman_ptr,&bitptr);
    }
  if (bitptr)
    writebits(0,8-bitptr,&huffman_ptr,&bitptr);
  *huffman_len=(int)(huffman_ptr-huffman);
  /* Output dictionary. */
  /* First the largest symbol value is written in 16 bits. No bits are
     encoded for symbols larger than this.  Then one bit signifies if
     there is a used symbol: 1 If unused entry: 0 If used symbol the 5
     following bits encode the length of the symbol. Worst case is
     thus 6*65538 bits used for the dictionary. That won't happen
     unless there's really that many values in use. If that is so,
     well, either we compress well, or we have many values anyway. */
  /* First sort the dictionary wrt symbol */
  Ptngc_merge_sort(codelength,ndict,sizeof *codelength,comp_codes_value,NULL);
  bitptr=0;
  huffman_ptr=huffman_dict;
  *huffman_ptr++=(unsigned char)(codelength[ndict-1].dict&0xFFU);
  *huffman_ptr++=(unsigned char)((codelength[ndict-1].dict>>8)&0xFFU);
  *huffman_ptr++=(unsigned char)((codelength[ndict-1].dict>>16)&0xFFU);
  huffman_dict_unpacked[0]=(unsigned char)(codelength[ndict-1].dict&0xFFU);
  huffman_dict_unpacked[1]=(unsigned char)((codelength[ndict-1].dict>>8)&0xFFU);
  huffman_dict_unpacked[2]=(unsigned char)((codelength[ndict-1].dict>>16)&0xFFU);
  for (i=0; i<=(int)codelength[ndict-1].dict; i++)
    {
      /* Do I have this value? */
      int ihave=0;
      int j;
      for (j=0; j<ndict; j++)
        if (codelength[j].dict==(unsigned int)i)
          {

            ihave=1;
            writebits(1,1,&huffman_ptr,&bitptr);
            writebits(codelength[j].length,5,&huffman_ptr,&bitptr);
            huffman_dict_unpacked[3+i]=codelength[j].length;
            break;
          }
      if (!ihave)
        {
          writebits(0,1,&huffman_ptr,&bitptr);
          huffman_dict_unpacked[3+i]=0;
        }
    }
  if (bitptr)
    writebits(0,8-bitptr,&huffman_ptr,&bitptr);
  *huffman_dictlen=(int)(huffman_ptr-huffman_dict);
  *huffman_dict_unpackedlen=3+codelength[ndict-1].dict+1;

  /* Free info about codes and length. */
  free(codelength);
}

void Ptngc_comp_conv_from_huffman(unsigned char *huffman,
                            unsigned int *vals, const int nvals,
                            const int ndict,
                            unsigned char *huffman_dict,
                            const int huffman_dictlen,
                            unsigned int *huffman_dict_unpacked,
                            const int huffman_dict_unpackedlen)
{
  struct codelength *codelength=warnmalloc(ndict*sizeof *codelength);
  int i,j;
  int maxdict;
  int code;
  unsigned char *huffman_ptr;
  int bitptr;
  (void)huffman_dictlen;
  (void)huffman_dict_unpackedlen;
  if (huffman_dict_unpacked)
    {
      maxdict=huffman_dict_unpacked[0]|(huffman_dict_unpacked[1]<<8)|(huffman_dict_unpacked[2]<<16);
      j=0;
      for(i=0; i<=maxdict; i++)
        {
          if (huffman_dict_unpacked[3+i]!=0)
            {
              codelength[j].length=huffman_dict_unpacked[3+i];
              codelength[j].dict=i;
#if 0
              printf("%d %d\n",
                     codelength[j].length,
                     codelength[j].dict);
#endif
              j++;
            }
        }
    }
  else
    {
      huffman_ptr=huffman_dict;
      maxdict=((unsigned int)huffman_ptr[0])|(((unsigned int)huffman_ptr[1])<<8)|(((unsigned int)huffman_ptr[2])<<16);
      huffman_ptr+=3;
      bitptr=0;
      j=0;
      for(i=0; i<=maxdict; i++)
        {
          int bit=readbits(1,&huffman_ptr,&bitptr);
          if (bit)
            {
              codelength[j].length=readbits(5,&huffman_ptr,&bitptr);
              codelength[j].dict=i;
#if 0
              printf("%d %d\n",
                     codelength[j].length,
                     codelength[j].dict);
#endif
              j++;
            }
        }
    }
  /* Sort codes wrt length/value. */
  Ptngc_merge_sort(codelength,ndict,sizeof *codelength,comp_codes,NULL);
  /* Canonicalize codes. */
  code=0;
  for (i=0; i<ndict; i++)
    {
      codelength[i].code=code;
      if (i<ndict-1)
        code=(code+1)<<(codelength[i+1].length-codelength[i].length);
    }
#if 0
  {
    for (i=0; i<ndict; i++)
      {
        printf("%d -\t\t %d %d ",codelength[i].dict,codelength[i].length,codelength[i].code);
        {
          unsigned int c=codelength[i].code;
          int j;
          unsigned int mask=1<<(codelength[i].length-1);
          for (j=codelength[i].length-1; j>=0; j--)
            {
              int bit=c&mask;
              if (bit)
                printf("1");
              else
                printf("0");
              mask>>=1;
            }
          printf("\n");
        }
      }
  }
#endif
  /* Decompress data. */
  huffman_ptr=huffman;
  bitptr=0;
  for (i=0; i<nvals; i++)
    {
      unsigned int symbol;
      int len=codelength[0].length;
      symbol=readbits(len,&huffman_ptr,&bitptr);
      j=0;
      while (symbol!=codelength[j].code)
        {
          int newlen;
          j++;
          newlen=codelength[j].length;
          if (newlen!=len)
            {
              symbol<<=(newlen-len);
              symbol|=readbits(newlen-len,&huffman_ptr,&bitptr);
              len=newlen;
            }
        }
      vals[i]=codelength[j].dict;
    }
  /* Free info about codes and length. */
  free(codelength);
}
