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
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/tng_compress.h"
#include "../../include/compression/bwlzh.h"
#include "../../include/compression/huffman.h"
#include "../../include/compression/dict.h"
#include "../../include/compression/rle.h"
#include "../../include/compression/vals16.h"

int Ptngc_comp_huff_buflen(const int nvals)
{
  return 132000+nvals*8;
}

/* the value pointed to by chosen_algo should be sent as -1 for autodetect. */
void Ptngc_comp_huff_compress_verbose(unsigned int *vals, int nvals,
                                unsigned char *huffman, int *huffman_len,
                                int *huffdatalen,
                                int *huffman_lengths,int *chosen_algo,
                                const int isvals16)
{
  unsigned int *dict=warnmalloc(0x20005*sizeof *dict);
  unsigned int *hist=warnmalloc(0x20005*sizeof *hist);
  unsigned int *vals16=NULL;
  unsigned char *huffdict=warnmalloc(0x20005*sizeof *huffdict);
  unsigned int *huffdictunpack=warnmalloc(0x20005*sizeof *huffdictunpack);
  unsigned char *huffman1=warnmalloc(2*0x20005*sizeof *huffman1);
  unsigned char *huffdict1=warnmalloc(0x20005*sizeof *huffdict1);
  unsigned int *huffdictunpack1=warnmalloc(0x20005*sizeof *huffdictunpack1);
  unsigned int *huffdictrle=warnmalloc((3*0x20005+3)*sizeof *huffdictrle);
  unsigned char *huffman2=warnmalloc(6*0x20005*sizeof *huffman2);
  unsigned char *huffdict2=warnmalloc(0x20005*sizeof *huffdict2);
  unsigned int *huffdictunpack2=warnmalloc(0x20005*sizeof *huffdictunpack2);
  int i;
  int ndict,ndict1,ndict2;
  int nhuff,nhuffdict,nhuffdictunpack;
  int nhuff1,nhuffdict1,nhuffdictunpack1;
  int nhuffrle,nhuff2,nhuffdict2,nhuffdictunpack2;
  int nvals16;

  /* Do I need to convert to vals16? */
  if (!isvals16)
  {
    vals16=warnmalloc(nvals*3*sizeof *vals16);
    Ptngc_comp_conv_to_vals16(vals,nvals,vals16,&nvals16);
    nvals=nvals16;
    vals=vals16;
  }
  else
    nvals16=nvals;

  /* Determine probabilities. */
  Ptngc_comp_make_dict_hist(vals,nvals,dict,&ndict,hist);

  /* First compress the data using huffman coding (place it ready for output at 14 (code for algorithm+length etc.). */
  Ptngc_comp_conv_to_huffman(vals,nvals,dict,ndict,hist,
                       huffman+14,&nhuff,
                       huffdict,&nhuffdict,
                       huffdictunpack,&nhuffdictunpack);
  *huffdatalen=nhuff;

  /* Algorithm 0 stores the huffman dictionary directly (+ a code for
     the algorithm) + lengths of the huffman buffer (4) and the huffman dictionary (3). */
  huffman_lengths[0]=nhuff+nhuffdict+1*2+3*4+3+3;
  /* Next we try to compress the huffman dictionary using huffman
     coding ... (algorithm 1) */

  /* Determine probabilities. */
  Ptngc_comp_make_dict_hist(huffdictunpack,nhuffdictunpack,dict,&ndict1,hist);
  /* Pack huffman dictionary */
  Ptngc_comp_conv_to_huffman(huffdictunpack,nhuffdictunpack,
                       dict,ndict1,hist,
                       huffman1,&nhuff1,
                       huffdict1,&nhuffdict1,
                       huffdictunpack1,&nhuffdictunpack1);
  huffman_lengths[1]=nhuff+nhuff1+nhuffdict1+1*2+3*4+3+3+3+3+3;

  /* ... and rle + huffman coding ... (algorithm 2) Pack any repetetitive patterns. */
  Ptngc_comp_conv_to_rle(huffdictunpack,nhuffdictunpack,
                   huffdictrle,&nhuffrle,1);

  /* Determine probabilities. */
  Ptngc_comp_make_dict_hist(huffdictrle,nhuffrle,dict,&ndict2,hist);
  /* Pack huffman dictionary */
  Ptngc_comp_conv_to_huffman(huffdictrle,nhuffrle,
                       dict,ndict2,hist,
                       huffman2,&nhuff2,
                       huffdict2,&nhuffdict2,
                       huffdictunpack2,&nhuffdictunpack2);
  huffman_lengths[2]=nhuff+nhuff2+nhuffdict2+1*2+3*4+3+3+3+3+3+3;

  /* Choose the best algorithm and output the data. */
  if ((*chosen_algo==0) || ((*chosen_algo==-1) &&
                            (((huffman_lengths[0]<huffman_lengths[1]) &&
                              (huffman_lengths[0]<huffman_lengths[2])))))
    {
      *chosen_algo=0;
      *huffman_len=huffman_lengths[0];
      huffman[0]=isvals16;
      huffman[1]=0;
      huffman[2]=((unsigned int)nvals16)&0xFFU;
      huffman[3]=(((unsigned int)nvals16)>>8)&0xFFU;
      huffman[4]=(((unsigned int)nvals16)>>16)&0xFFU;
      huffman[5]=(((unsigned int)nvals16)>>24)&0xFFU;
      huffman[6]=((unsigned int)nvals)&0xFFU;
      huffman[7]=(((unsigned int)nvals)>>8)&0xFFU;
      huffman[8]=(((unsigned int)nvals)>>16)&0xFFU;
      huffman[9]=(((unsigned int)nvals)>>24)&0xFFU;
      huffman[10]=((unsigned int)nhuff)&0xFFU;
      huffman[11]=(((unsigned int)nhuff)>>8)&0xFFU;
      huffman[12]=(((unsigned int)nhuff)>>16)&0xFFU;
      huffman[13]=(((unsigned int)nhuff)>>24)&0xFFU;
      huffman[14+nhuff]=((unsigned int)nhuffdict)&0xFFU;
      huffman[15+nhuff]=(((unsigned int)nhuffdict)>>8)&0xFFU;
      huffman[16+nhuff]=(((unsigned int)nhuffdict)>>16)&0xFFU;
      huffman[17+nhuff]=((unsigned int)ndict)&0xFFU;
      huffman[18+nhuff]=(((unsigned int)ndict)>>8)&0xFFU;
      huffman[19+nhuff]=(((unsigned int)ndict)>>16)&0xFFU;
      for (i=0; i<nhuffdict; i++)
        huffman[20+nhuff+i]=huffdict[i];
    }
  else if ((*chosen_algo==1) || ((*chosen_algo==-1) &&
                                 ((huffman_lengths[1]<huffman_lengths[2]))))
    {
      *chosen_algo=1;
      *huffman_len=huffman_lengths[1];
      huffman[0]=isvals16;
      huffman[1]=1;
      huffman[2]=((unsigned int)nvals16)&0xFFU;
      huffman[3]=(((unsigned int)nvals16)>>8)&0xFFU;
      huffman[4]=(((unsigned int)nvals16)>>16)&0xFFU;
      huffman[5]=(((unsigned int)nvals16)>>24)&0xFFU;
      huffman[6]=((unsigned int)nvals)&0xFFU;
      huffman[7]=(((unsigned int)nvals)>>8)&0xFFU;
      huffman[8]=(((unsigned int)nvals)>>16)&0xFFU;
      huffman[9]=(((unsigned int)nvals)>>24)&0xFFU;
      huffman[10]=((unsigned int)nhuff)&0xFFU;
      huffman[11]=(((unsigned int)nhuff)>>8)&0xFFU;
      huffman[12]=(((unsigned int)nhuff)>>16)&0xFFU;
      huffman[13]=(((unsigned int)nhuff)>>24)&0xFFU;
      huffman[14+nhuff]=((unsigned int)nhuffdictunpack)&0xFFU;
      huffman[15+nhuff]=(((unsigned int)nhuffdictunpack)>>8)&0xFFU;
      huffman[16+nhuff]=(((unsigned int)nhuffdictunpack)>>16)&0xFFU;
      huffman[17+nhuff]=((unsigned int)ndict)&0xFFU;
      huffman[18+nhuff]=(((unsigned int)ndict)>>8)&0xFFU;
      huffman[19+nhuff]=(((unsigned int)ndict)>>16)&0xFFU;
      huffman[20+nhuff]=((unsigned int)nhuff1)&0xFFU;
      huffman[21+nhuff]=(((unsigned int)nhuff1)>>8)&0xFFU;
      huffman[22+nhuff]=(((unsigned int)nhuff1)>>16)&0xFFU;
      huffman[23+nhuff]=((unsigned int)nhuffdict1)&0xFFU;
      huffman[24+nhuff]=(((unsigned int)nhuffdict1)>>8)&0xFFU;
      huffman[25+nhuff]=(((unsigned int)nhuffdict1)>>16)&0xFFU;
      huffman[26+nhuff]=((unsigned int)ndict1)&0xFFU;
      huffman[27+nhuff]=(((unsigned int)ndict1)>>8)&0xFFU;
      huffman[28+nhuff]=(((unsigned int)ndict1)>>16)&0xFFU;
      for (i=0; i<nhuff1; i++)
        huffman[29+nhuff+i]=huffman1[i];
      for (i=0; i<nhuffdict1; i++)
        huffman[29+nhuff+nhuff1+i]=huffdict1[i];
    }
  else
    {
      *chosen_algo=2;
      *huffman_len=huffman_lengths[2];
      huffman[0]=isvals16;
      huffman[1]=2;
      huffman[2]=((unsigned int)nvals16)&0xFFU;
      huffman[3]=(((unsigned int)nvals16)>>8)&0xFFU;
      huffman[4]=(((unsigned int)nvals16)>>16)&0xFFU;
      huffman[5]=(((unsigned int)nvals16)>>24)&0xFFU;
      huffman[6]=((unsigned int)nvals)&0xFFU;
      huffman[7]=(((unsigned int)nvals)>>8)&0xFFU;
      huffman[8]=(((unsigned int)nvals)>>16)&0xFFU;
      huffman[9]=(((unsigned int)nvals)>>24)&0xFFU;
      huffman[10]=((unsigned int)nhuff)&0xFFU;
      huffman[11]=(((unsigned int)nhuff)>>8)&0xFFU;
      huffman[12]=(((unsigned int)nhuff)>>16)&0xFFU;
      huffman[13]=(((unsigned int)nhuff)>>24)&0xFFU;
      huffman[14+nhuff]=((unsigned int)nhuffdictunpack)&0xFFU;
      huffman[15+nhuff]=(((unsigned int)nhuffdictunpack)>>8)&0xFFU;
      huffman[16+nhuff]=(((unsigned int)nhuffdictunpack)>>16)&0xFFU;
      huffman[17+nhuff]=((unsigned int)ndict)&0xFFU;
      huffman[18+nhuff]=(((unsigned int)ndict)>>8)&0xFFU;
      huffman[19+nhuff]=(((unsigned int)ndict)>>16)&0xFFU;
      huffman[20+nhuff]=((unsigned int)nhuffrle)&0xFFU;
      huffman[21+nhuff]=(((unsigned int)nhuffrle)>>8)&0xFFU;
      huffman[22+nhuff]=(((unsigned int)nhuffrle)>>16)&0xFFU;
      huffman[23+nhuff]=((unsigned int)nhuff2)&0xFFU;
      huffman[24+nhuff]=(((unsigned int)nhuff2)>>8)&0xFFU;
      huffman[25+nhuff]=(((unsigned int)nhuff2)>>16)&0xFFU;
      huffman[26+nhuff]=((unsigned int)nhuffdict2)&0xFFU;
      huffman[27+nhuff]=(((unsigned int)nhuffdict2)>>8)&0xFFU;
      huffman[28+nhuff]=(((unsigned int)nhuffdict2)>>16)&0xFFU;
      huffman[29+nhuff]=((unsigned int)ndict2)&0xFFU;
      huffman[30+nhuff]=(((unsigned int)ndict2)>>8)&0xFFU;
      huffman[31+nhuff]=(((unsigned int)ndict2)>>16)&0xFFU;
      for (i=0; i<nhuff2; i++)
        huffman[32+nhuff+i]=huffman2[i];
      for (i=0; i<nhuffdict2; i++)
        huffman[32+nhuff+nhuff2+i]=huffdict2[i];
    }
  if (!isvals16)
    free(vals16);

  free(huffdictunpack2);
  free(huffdict2);
  free(huffman2);
  free(huffdictrle);
  free(huffdictunpack1);
  free(huffdict1);
  free(huffman1);
  free(huffdictunpack);
  free(huffdict);
  free(hist);
  free(dict);
}

void Ptngc_comp_huff_compress(unsigned int *vals, const int nvals,
                        unsigned char *huffman, int *huffman_len)
{
  int huffman_lengths[N_HUFFMAN_ALGO];
  int algo=-1;
  int huffdatalen;
  Ptngc_comp_huff_compress_verbose(vals,nvals,huffman,huffman_len,&huffdatalen,
                             huffman_lengths,&algo,0);
}

void Ptngc_comp_huff_decompress(unsigned char *huffman, const int huffman_len,
                          unsigned int *vals)
{
  int isvals16=(int)huffman[0];
  unsigned int *vals16=NULL;
  int algo=(int)huffman[1];
  int nvals16=(int)((unsigned int)huffman[2]|
                  (((unsigned int)huffman[3])<<8)|
                  (((unsigned int)huffman[4])<<16)|
                  (((unsigned int)huffman[5])<<24));
  int nvals=(int)((unsigned int)huffman[6]|
                  (((unsigned int)huffman[7])<<8)|
                  (((unsigned int)huffman[8])<<16)|
                  (((unsigned int)huffman[9])<<24));
  int nhuff=(int)((unsigned int)huffman[10]|
                  (((unsigned int)huffman[11])<<8)|
                  (((unsigned int)huffman[12])<<16)|
                  (((unsigned int)huffman[13])<<24));
  int ndict=(int)((unsigned int)huffman[17+nhuff]|
                  (((unsigned int)huffman[18+nhuff])<<8)|
                  (((unsigned int)huffman[19+nhuff])<<16));
  (void)huffman_len;
  if (!isvals16)
    vals16=warnmalloc(nvals16*sizeof *vals16);
  else
    {
      vals16=vals;
      nvals16=nvals;
    }
  if (algo==0)
    {
      int nhuffdict=(int)((unsigned int)huffman[14+nhuff]|
                          (((unsigned int)huffman[15+nhuff])<<8)|
                          (((unsigned int)huffman[16+nhuff])<<16));
      Ptngc_comp_conv_from_huffman(huffman+14,vals16,nvals16,ndict,
                             huffman+20+nhuff,nhuffdict,NULL,0);
    }
  else if (algo==1)
    {
      unsigned int *huffdictunpack=warnmalloc(0x20005*sizeof *huffdictunpack);
      /* First the dictionary needs to be uncompressed. */
      int nhuffdictunpack=(int)((unsigned int)huffman[14+nhuff]|
                                (((unsigned int)huffman[15+nhuff])<<8)|
                                (((unsigned int)huffman[16+nhuff])<<16));
      int nhuff1=(int)((unsigned int)huffman[20+nhuff]|
                       (((unsigned int)huffman[21+nhuff])<<8)|
                       (((unsigned int)huffman[22+nhuff])<<16));
      int nhuffdict1=(int)((unsigned int)huffman[23+nhuff]|
                           (((unsigned int)huffman[24+nhuff])<<8)|
                           (((unsigned int)huffman[25+nhuff])<<16));
      int ndict1=(int)((unsigned int)huffman[26+nhuff]|
                       (((unsigned int)huffman[27+nhuff])<<8)|
                       (((unsigned int)huffman[28+nhuff])<<16));
      Ptngc_comp_conv_from_huffman(huffman+29+nhuff,huffdictunpack,
                             nhuffdictunpack,ndict1,
                             huffman+29+nhuff+nhuff1,nhuffdict1,NULL,0);
      /* Then decompress the "real" data. */
      Ptngc_comp_conv_from_huffman(huffman+14,vals16,nvals16,ndict,
                             NULL,0,huffdictunpack,nhuffdictunpack);
      free(huffdictunpack);
    }
  else if (algo==2)
    {
      unsigned int *huffdictunpack=warnmalloc(0x20005*sizeof *huffdictunpack);
      unsigned int *huffdictrle=warnmalloc((3*0x20005+3)*sizeof *huffdictrle);
      /* First the dictionary needs to be uncompressed. */
      int nhuffdictunpack=(int)((unsigned int)huffman[14+nhuff]|
                                (((unsigned int)huffman[15+nhuff])<<8)|
                                (((unsigned int)huffman[16+nhuff])<<16));
      int nhuffrle=(int)((unsigned int)huffman[20+nhuff]|
                         (((unsigned int)huffman[21+nhuff])<<8)|
                         (((unsigned int)huffman[22+nhuff])<<16));
      int nhuff2=(int)((unsigned int)huffman[23+nhuff]|
                       (((unsigned int)huffman[24+nhuff])<<8)|
                       (((unsigned int)huffman[25+nhuff])<<16));
      int nhuffdict2=(int)((unsigned int)huffman[26+nhuff]|
                           (((unsigned int)huffman[27+nhuff])<<8)|
                           (((unsigned int)huffman[28+nhuff])<<16));
      int ndict2=(int)((unsigned int)huffman[29+nhuff]|
                       (((unsigned int)huffman[30+nhuff])<<8)|
                       (((unsigned int)huffman[31+nhuff])<<16));
      Ptngc_comp_conv_from_huffman(huffman+32+nhuff,huffdictrle,
                             nhuffrle,ndict2,
                             huffman+32+nhuff+nhuff2,nhuffdict2,NULL,0);
      /* Then uncompress the rle data */
      Ptngc_comp_conv_from_rle(huffdictrle,huffdictunpack,nhuffdictunpack);
      /* Then decompress the "real" data. */
      Ptngc_comp_conv_from_huffman(huffman+14,vals16,nvals16,ndict,
                             NULL,0,huffdictunpack,nhuffdictunpack);
      free(huffdictrle);
      free(huffdictunpack);
    }

  /* Do I need to convert from vals16? */
  if (!isvals16)
  {
    int nvalsx;
    Ptngc_comp_conv_from_vals16(vals16,nvals16,vals,&nvalsx);
    free(vals16);
  }
}

static char *huff_algo_names[N_HUFFMAN_ALGO]=
  {
    "Huffman (dict=raw)",
    "Huffman (dict=Huffman)",
    "Huffman (dict=RLE+Huffman)"
  };

char *Ptngc_comp_get_huff_algo_name(const int algo)
{
  if (algo<0)
    return NULL;
  else if (algo>=N_HUFFMAN_ALGO)
    return NULL;
  return huff_algo_names[algo];
}
