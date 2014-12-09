/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../include/compression/tng_compress.h"
#include "../../include/compression/coder.h"
#include "../../include/compression/fixpoint.h"

/* Please see tng_compress.h for info on how to call these routines. */

/* This becomes TNGP for positions (little endian) and TNGV for velocities. In ASCII. */
#define MAGIC_INT_POS 0x50474E54
#define MAGIC_INT_VEL 0x56474E54

#define SPEED_DEFAULT 2 /* Default to relatively fast compression. For very good compression it makes sense to
                           choose speed=4 or speed=5 */

#define PRECISION(hi,lo) (Ptngc_i32x2_to_d(hi,lo))

#define MAX_FVAL 2147483647.

static int verify_input_data(double *x, const int natoms, const int nframes, const double precision)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        if (fabs(x[iframe*natoms*3+i*3+j]/precision+0.5)>=MAX_FVAL)
          goto error;
  return 0;
 error:
#if 0
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        if (fabs(x[iframe*natoms*3+i*3+j]/precision+0.5)>=MAX_FVAL)
          printf("ERROR. Too large value: %d %d %d: %g %g %g\n",iframe,i,j,x[iframe*natoms*3+i*3+j],precision,x[iframe*natoms*3+i*3+j]/precision/MAX_FVAL);
#endif
  return 1;
}

static int verify_input_data_float(float *x, const int natoms, const int nframes, const float precision)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        if (fabs(x[iframe*natoms*3+i*3+j]/precision+0.5)>=MAX_FVAL)
          goto error;
  return 0;
 error:
#if 0
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        if (fabs(x[iframe*natoms*3+i*3+j]/precision+0.5)>=MAX_FVAL)
          printf("ERROR. Too large value: %d %d %d: %g %g %g\n",iframe,i,j,x[iframe*natoms*3+i*3+j],precision,x[iframe*natoms*3+i*3+j]/precision/MAX_FVAL);
#endif
  return 1;
}

static int quantize(double *x, const int natoms, const int nframes,
                    const double precision,
                    int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        quant[iframe*natoms*3+i*3+j]=(int)floor((x[iframe*natoms*3+i*3+j]/precision)+0.5);
  return verify_input_data(x,natoms,nframes,precision);
}

static int quantize_float(float *x, int natoms, const int nframes,
                          const float precision,
                          int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        quant[iframe*natoms*3+i*3+j]=(int)floor((x[iframe*natoms*3+i*3+j]/precision)+0.5);
  return verify_input_data_float(x,natoms,nframes,precision);
}

static void quant_inter_differences(int *quant, const int natoms, const int nframes,
                                    int *quant_inter)
{
  int iframe, i, j;
  /* The first frame is used for absolute positions. */
  for (i=0; i<natoms; i++)
    for (j=0; j<3; j++)
      quant_inter[i*3+j]=quant[i*3+j];
  /* For all other frames, the difference to the previous frame is used. */
  for (iframe=1; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        quant_inter[iframe*natoms*3+i*3+j]=quant[iframe*natoms*3+i*3+j]-quant[(iframe-1)*natoms*3+i*3+j];
}

static void quant_intra_differences(int *quant, const int natoms, const int nframes,
                                    int *quant_intra)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    {
      /* The first atom is used with its absolute position. */
      for (j=0; j<3; j++)
        quant_intra[iframe*natoms*3+j]=quant[iframe*natoms*3+j];
      /* For all other atoms the intraframe differences are computed. */
      for (i=1; i<natoms; i++)
        for (j=0; j<3; j++)
          quant_intra[iframe*natoms*3+i*3+j]=quant[iframe*natoms*3+i*3+j]-quant[iframe*natoms*3+(i-1)*3+j];
    }
}

static void unquantize(double *x, int natoms, const int nframes,
                       const double precision,
                       int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        x[iframe*natoms*3+i*3+j]=(double)quant[iframe*natoms*3+i*3+j]*precision;
}

static void unquantize_float(float *x, const int natoms, const int nframes,
                             const float precision,
                             int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (i=0; i<natoms; i++)
      for (j=0; j<3; j++)
        x[iframe*natoms*3+i*3+j]=(float)quant[iframe*natoms*3+i*3+j]*precision;
}

static void unquantize_inter_differences(double *x, int natoms, const int nframes,
                                         const double precision,
                                         int *quant)
{
  int iframe, i, j;
  for (i=0; i<natoms; i++)
    for (j=0; j<3; j++)
      {
        int q=quant[i*3+j]; /* First value. */
        x[i*3+j]=(double)q*precision;
        for (iframe=1; iframe<nframes; iframe++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=(double)q*precision;
          }
      }
}

static void unquantize_inter_differences_float(float *x, const int natoms, const int nframes,
                                               const float precision,
                                               int *quant)
{
  int iframe, i, j;
  for (i=0; i<natoms; i++)
    for (j=0; j<3; j++)
      {
        int q=quant[i*3+j]; /* First value. */
        x[i*3+j]=(float)q*precision;
        for (iframe=1; iframe<nframes; iframe++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=(float)q*precision;
          }
      }
}

static void unquantize_inter_differences_int(int *x, const int natoms, const int nframes,
                                             int *quant)
{
  int iframe, i, j;
  for (i=0; i<natoms; i++)
    for (j=0; j<3; j++)
      {
        int q=quant[i*3+j]; /* First value. */
        x[i*3+j]=q;
        for (iframe=1; iframe<nframes; iframe++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=q;
          }
      }
}

/* In frame update required for the initial frame if intra-frame
   compression was used. */
static void unquant_intra_differences_first_frame(int *quant, const int natoms)
{
  int i,j;
  for (j=0; j<3; j++)
    {
      int q=quant[j];
      for (i=1; i<natoms; i++)
        {
          q+=quant[i*3+j];
          quant[i*3+j]=q;
        }
    }
#if 0
  for (j=0; j<3; j++)
    for (i=0; i<natoms; i++)
      {
        printf("UQ: %d %d %d: %d\n",0,j,i,quant[i*3+j]);
      }
#endif
}

static void unquantize_intra_differences(double *x, const int natoms, const int nframes,
                                         const double precision,
                                         int *quant)
{
  int iframe, i, j;
#if 0
  printf("UQ precision=%g\n",precision);
#endif
  for (iframe=0; iframe<nframes; iframe++)
    for (j=0; j<3; j++)
      {
        int q=quant[iframe*natoms*3+j];
        x[iframe*natoms*3+j]=(double)q*precision;
        for (i=1; i<natoms; i++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=(double)q*precision;
          }
      }
}

static void unquantize_intra_differences_float(float *x, const int natoms, const int nframes,
                                               const float precision,
                                               int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (j=0; j<3; j++)
      {
        int q=quant[iframe*natoms*3+j];
        x[iframe*natoms*3+j]=(float)q*precision;
        for (i=1; i<natoms; i++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=(float)q*precision;
          }
      }
}

static void unquantize_intra_differences_int(int *x, const int natoms, const int nframes,
                                             int *quant)
{
  int iframe, i, j;
  for (iframe=0; iframe<nframes; iframe++)
    for (j=0; j<3; j++)
      {
        int q=quant[iframe*natoms*3+j];
        x[iframe*natoms*3+j]=q;
        for (i=1; i<natoms; i++)
          {
            q+=quant[iframe*natoms*3+i*3+j];
            x[iframe*natoms*3+i*3+j]=q;
          }
      }
}

/* Buffer num 8 bit bytes into buffer location buf */
static void bufferfix(unsigned char *buf, fix_t v, int num)
{
  /* Store in little endian format. */
  unsigned char c; /* at least 8 bits. */
  c=(unsigned char)(v & 0xFFU);
  while (num--)
    {
      *buf++=c;
      v >>= 8;
      c=(unsigned char)(v & 0xFFU);
    }
}

static fix_t readbufferfix(unsigned char *buf, int num)
{
  unsigned char b;
  int shift=0;
  fix_t f=0UL;
  int cnt=0;
  do
    {
      b=buf[cnt++];
      f |= ((fix_t)b & 0xFF)<<shift;
      shift+=8;
    } while (--num);
  return f;
}

/* Perform position compression from the quantized data. */
static void compress_quantized_pos(int *quant, int *quant_inter, int *quant_intra,
                                   const int natoms, const int nframes,
                                   const int speed,
                                   const int initial_coding, const int initial_coding_parameter,
                                   const int coding, const int coding_parameter,
                                   const fix_t prec_hi, const fix_t prec_lo,
                                   int *nitems,
                                   char *data)
{
  int bufloc=0;
  char *datablock=NULL;
  int length=0;
  /* Information needed for decompression. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)MAGIC_INT_POS,4);
  bufloc+=4;
  /* Number of atoms. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)natoms,4);
  bufloc+=4;
  /* Number of frames. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)nframes,4);
  bufloc+=4;
  /* Initial coding. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)initial_coding,4);
  bufloc+=4;
  /* Initial coding parameter. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)initial_coding_parameter,4);
  bufloc+=4;
  /* Coding. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)coding,4);
  bufloc+=4;
  /* Coding parameter. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)coding_parameter,4);
  bufloc+=4;
  /* Precision. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,prec_lo,4);
  bufloc+=4;
  if (data)
    bufferfix((unsigned char*)data+bufloc,prec_hi,4);
  bufloc+=4;
  /* The initial frame */
  if ((initial_coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
      (initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_POS_XTC3))
    {
      struct coder *coder=Ptngc_coder_init();
      length=natoms*3;
      datablock=(char*)Ptngc_pack_array(coder,quant,&length,
                                            initial_coding,initial_coding_parameter,natoms,speed);
      Ptngc_coder_deinit(coder);
    }
  else if ((initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
           (initial_coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
    {
      struct coder *coder=Ptngc_coder_init();
      length=natoms*3;
      datablock=(char*)Ptngc_pack_array(coder,quant_intra,&length,
                                            initial_coding,initial_coding_parameter,natoms,speed);
      Ptngc_coder_deinit(coder);
    }
  /* Block length. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)length,4);
  bufloc+=4;
  /* The actual data block. */
  if (data)
    memcpy(data+bufloc,datablock,length);
  free(datablock);
  bufloc+=length;
  /* The remaining frames */
  if (nframes>1)
    {
      datablock=NULL;
      /* Inter-frame compression? */
      if ((coding==TNG_COMPRESS_ALGO_POS_STOPBIT_INTER) ||
          (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTER) ||
          (coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTER))
        {
          struct coder *coder=Ptngc_coder_init();
          length=natoms*3*(nframes-1);
          datablock=(char*)Ptngc_pack_array(coder,quant_inter+natoms*3,&length,
                                                coding,coding_parameter,natoms,speed);
          Ptngc_coder_deinit(coder);
        }
      /* One-to-one compression? */
      else if ((coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
               (coding==TNG_COMPRESS_ALGO_POS_XTC3) ||
               (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE))
        {
          struct coder *coder=Ptngc_coder_init();
          length=natoms*3*(nframes-1);
          datablock=(char*)Ptngc_pack_array(coder,quant+natoms*3,&length,
                                                coding,coding_parameter,natoms,speed);
          Ptngc_coder_deinit(coder);
        }
      /* Intra-frame compression? */
      else if ((coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
               (coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
        {
          struct coder *coder=Ptngc_coder_init();
          length=natoms*3*(nframes-1);
          datablock=(char*)Ptngc_pack_array(coder,quant_intra+natoms*3,&length,
                                                coding,coding_parameter,natoms,speed);
          Ptngc_coder_deinit(coder);
        }
      /* Block length. */
      if (data)
        bufferfix((unsigned char*)data+bufloc,(fix_t)length,4);
      bufloc+=4;
      if (datablock)
        {
          if (data)
            memcpy(data+bufloc,datablock,length);
          free(datablock);
        }
      bufloc+=length;
    }
  *nitems=bufloc;
}

/* Perform velocity compression from vel into the data block */
static void compress_quantized_vel(int *quant, int *quant_inter,
                                   const int natoms, const int nframes,
                                   const int speed,
                                   const int initial_coding, const int initial_coding_parameter,
                                   const int coding, const int coding_parameter,
                                   const fix_t prec_hi, const fix_t prec_lo,
                                   int *nitems,
                                   char *data)
{
  int bufloc=0;
  char *datablock=NULL;
  int length;
  /* Information needed for decompression. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)MAGIC_INT_VEL,4);
  bufloc+=4;
  /* Number of atoms. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)natoms,4);
  bufloc+=4;
  /* Number of frames. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)nframes,4);
  bufloc+=4;
  /* Initial coding. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)initial_coding,4);
  bufloc+=4;
  /* Initial coding parameter. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)initial_coding_parameter,4);
  bufloc+=4;
  /* Coding. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)coding,4);
  bufloc+=4;
  /* Coding parameter. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)coding_parameter,4);
  bufloc+=4;
  /* Precision. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,prec_lo,4);
  bufloc+=4;
  if (data)
    bufferfix((unsigned char*)data+bufloc,prec_hi,4);
  bufloc+=4;

  length=natoms*3;
  /* The initial frame */
  if ((initial_coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE))
    {
      struct coder *coder=Ptngc_coder_init();
      datablock=(char*)Ptngc_pack_array(coder,quant,&length,
                                            initial_coding,initial_coding_parameter,natoms,speed);
      Ptngc_coder_deinit(coder);
    }
  /* Block length. */
  if (data)
    bufferfix((unsigned char*)data+bufloc,(fix_t)length,4);
  bufloc+=4;
  /* The actual data block. */
  if (data && datablock)
    {
      memcpy(data+bufloc,datablock,length);
      free(datablock);
      bufloc+=length;
    }
  /* The remaining frames */
  if (nframes>1)
    {
      datablock=NULL;
      /* Inter-frame compression? */
      if ((coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_INTER) ||
          (coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER) ||
          (coding==TNG_COMPRESS_ALGO_VEL_BWLZH_INTER))
        {
          struct coder *coder=Ptngc_coder_init();
          length=natoms*3*(nframes-1);
          datablock=(char*)Ptngc_pack_array(coder,quant_inter+natoms*3,&length,
                                                coding,coding_parameter,natoms,speed);
          Ptngc_coder_deinit(coder);
        }
      /* One-to-one compression? */
      else if ((coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE) ||
               (coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE) ||
               (coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE))
        {
          struct coder *coder=Ptngc_coder_init();
          length=natoms*3*(nframes-1);
          datablock=(char*)Ptngc_pack_array(coder,quant+natoms*3,&length,
                                                coding,coding_parameter,natoms,speed);
          Ptngc_coder_deinit(coder);
        }
      /* Block length. */
      if (data)
        bufferfix((unsigned char*)data+bufloc,(fix_t)length,4);
      bufloc+=4;
      if (data)
        memcpy(data+bufloc,datablock,length);
      free(datablock);
      bufloc+=length;
    }
  *nitems=bufloc;
}

static int determine_best_coding_stop_bits(struct coder *coder,int *input, int *length,
                                           int *coding_parameter, const int natoms)
{
  int bits;
  unsigned char *packed;
  int best_length=0;
  int new_parameter=-1;
  int io_length;
  for (bits=1; bits<20; bits++)
    {
      io_length=*length;
      packed=Ptngc_pack_array(coder,input,&io_length,
                                  TNG_COMPRESS_ALGO_STOPBIT,bits,natoms,0);
      if (packed)
        {
          if ((new_parameter==-1) || (io_length<best_length))
            {
              new_parameter=bits;
              best_length=io_length;
            }
          free(packed);
        }
    }
  if (new_parameter==-1)
    return 1;

  *coding_parameter=new_parameter;
  *length=best_length;
  return 0;
}

static int determine_best_coding_triple(struct coder *coder,int *input, int *length,
                                        int *coding_parameter, const int natoms)
{
  int bits;
  unsigned char *packed;
  int best_length=0;
  int new_parameter=-1;
  int io_length;
  for (bits=1; bits<20; bits++)
    {
      io_length=*length;
      packed=Ptngc_pack_array(coder,input,&io_length,
                                  TNG_COMPRESS_ALGO_TRIPLET,bits,natoms,0);
      if (packed)
        {
          if ((new_parameter==-1) || (io_length<best_length))
            {
              new_parameter=bits;
              best_length=io_length;
            }
          free(packed);
        }
    }
  if (new_parameter==-1)
    return 1;

  *coding_parameter=new_parameter;
  *length=best_length;
  return 0;
}

static void determine_best_pos_initial_coding(int *quant, int *quant_intra, const int natoms, const int speed,
                                              const fix_t prec_hi, const fix_t prec_lo,
                                              int *initial_coding, int *initial_coding_parameter)
{
  if (*initial_coding==-1)
    {
      /* Determine all parameters automatically */
      int best_coding;
      int best_coding_parameter;
      int best_code_size;
      int current_coding;
      int current_coding_parameter;
      int current_code_size;
      struct coder *coder;
      /* Start with XTC2, it should always work. */
      current_coding=TNG_COMPRESS_ALGO_POS_XTC2;
      current_coding_parameter=0;
      compress_quantized_pos(quant,NULL,quant_intra,natoms,1,speed,
                             current_coding,current_coding_parameter,
                             0,0,prec_hi,prec_lo,&current_code_size,NULL);
      best_coding=current_coding;
      best_coding_parameter=current_coding_parameter;
      best_code_size=current_code_size;

      /* Determine best parameter for triplet intra. */
      current_coding=TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3;
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant_intra,&current_code_size,&current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Determine best parameter for triplet one-to-one. */
      current_coding=TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3;
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant,&current_code_size,&current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      if (speed>=2)
        {
          current_coding=TNG_COMPRESS_ALGO_POS_XTC3;
          current_coding_parameter=0;
          compress_quantized_pos(quant,NULL,quant_intra,natoms,1,speed,
                                 current_coding,current_coding_parameter,
                                 0,0,prec_hi,prec_lo,&current_code_size,NULL);
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      /* Test BWLZH intra */
      if (speed>=6)
        {
          current_coding=TNG_COMPRESS_ALGO_POS_BWLZH_INTRA;
          current_coding_parameter=0;
          compress_quantized_pos(quant,NULL,quant_intra,natoms,1,speed,
                                 current_coding,current_coding_parameter,
                                 0,0,prec_hi,prec_lo,&current_code_size,NULL);
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
            }
        }
      *initial_coding=best_coding;
      *initial_coding_parameter=best_coding_parameter;
    }
  else
    {
      if (*initial_coding_parameter==-1)
        {
          if ((*initial_coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
              (*initial_coding==TNG_COMPRESS_ALGO_POS_XTC3) ||
              (*initial_coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
            *initial_coding_parameter=0;
          else if (*initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA)
            {
              struct coder *coder=Ptngc_coder_init();
              int current_code_size=natoms*3;
              determine_best_coding_triple(coder,quant_intra,&current_code_size,initial_coding_parameter,natoms);
              Ptngc_coder_deinit(coder);
            }
          else if (*initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE)
            {
              struct coder *coder=Ptngc_coder_init();
              int current_code_size=natoms*3;
              determine_best_coding_triple(coder,quant,&current_code_size,initial_coding_parameter,natoms);
              Ptngc_coder_deinit(coder);
            }
        }
    }
}

static void determine_best_pos_coding(int *quant, int *quant_inter, int *quant_intra, const int natoms, const int nframes, const int speed,
                                      const fix_t prec_hi, const fix_t prec_lo,
                                      int *coding, int *coding_parameter)
{
  if (*coding==-1)
    {
      /* Determine all parameters automatically */
      int best_coding;
      int best_coding_parameter;
      int best_code_size;
      int current_coding;
      int current_coding_parameter;
      int current_code_size;
      int initial_code_size;
      struct coder *coder;
      /* Always use XTC2 for the initial coding. */
      compress_quantized_pos(quant,quant_inter,quant_intra,natoms,1,speed,
                             TNG_COMPRESS_ALGO_POS_XTC2,0,
                             0,0,
                             prec_hi,prec_lo,&initial_code_size,NULL);
      /* Start with XTC2, it should always work. */
      current_coding=TNG_COMPRESS_ALGO_POS_XTC2;
      current_coding_parameter=0;
      compress_quantized_pos(quant,quant_inter,quant_intra,natoms,nframes,speed,
                             TNG_COMPRESS_ALGO_POS_XTC2,0,
                             current_coding,current_coding_parameter,
                             prec_hi,prec_lo,&current_code_size,NULL);
      best_coding=current_coding;
      best_coding_parameter=current_coding_parameter;
      best_code_size=current_code_size-initial_code_size; /* Correct for the use of XTC2 for the first frame. */

      /* Determine best parameter for stopbit interframe coding. */
      current_coding=TNG_COMPRESS_ALGO_POS_STOPBIT_INTER;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      if (!determine_best_coding_stop_bits(coder,quant_inter+natoms*3,&current_code_size,
                                           &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Determine best parameter for triplet interframe coding. */
      current_coding=TNG_COMPRESS_ALGO_POS_TRIPLET_INTER;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant_inter+natoms*3,&current_code_size,
                                        &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Determine best parameter for triplet intraframe coding. */
      current_coding=TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant_intra+natoms*3,&current_code_size,
                                        &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Determine best parameter for triplet one-to-one coding. */
      current_coding=TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant+natoms*3,&current_code_size,
                                        &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Test BWLZH inter */
      if (speed>=4)
        {
          current_coding=TNG_COMPRESS_ALGO_POS_BWLZH_INTER;
          current_coding_parameter=0;
          compress_quantized_pos(quant,quant_inter,quant_intra,natoms,nframes,speed,
                                 TNG_COMPRESS_ALGO_POS_XTC2,0,
                                 current_coding,current_coding_parameter,
                                 prec_hi,prec_lo,&current_code_size,NULL);
          current_code_size-=initial_code_size; /* Correct for the use of XTC2 for the first frame. */
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }

      /* Test BWLZH intra */
      if (speed>=6)
        {
          current_coding=TNG_COMPRESS_ALGO_POS_BWLZH_INTRA;
          current_coding_parameter=0;
          compress_quantized_pos(quant,quant_inter,quant_intra,natoms,nframes,speed,
                                 TNG_COMPRESS_ALGO_POS_XTC2,0,
                                 current_coding,current_coding_parameter,
                                 prec_hi,prec_lo,&current_code_size,NULL);
          current_code_size-=initial_code_size; /* Correct for the use of XTC2 for the first frame. */
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
            }
        }
      *coding=best_coding;
      *coding_parameter=best_coding_parameter;
    }
  else if (*coding_parameter==-1)
    {
      if ((*coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
          (*coding==TNG_COMPRESS_ALGO_POS_XTC3) ||
          (*coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTER) ||
          (*coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
        *coding_parameter=0;
      else if (*coding==TNG_COMPRESS_ALGO_POS_STOPBIT_INTER)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_stop_bits(coder,quant_inter+natoms*3,&current_code_size,
                                          coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTER)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_triple(coder,quant_inter+natoms*3,&current_code_size,
                                       coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_triple(coder,quant_intra+natoms*3,&current_code_size,
                                       coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_triple(coder,quant+natoms*3,&current_code_size,
                                       coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
    }
}

static void determine_best_vel_initial_coding(int *quant, const int natoms, const int speed,
                                              const fix_t prec_hi, const fix_t prec_lo,
                                              int *initial_coding, int *initial_coding_parameter)
{
  if (*initial_coding==-1)
    {
      /* Determine all parameters automatically */
      int best_coding=-1;
      int best_coding_parameter=-1;
      int best_code_size=-1;
      int current_coding;
      int current_coding_parameter;
      int current_code_size;
      struct coder *coder;
      /* Start to determine best parameter for stopbit one-to-one. */
      current_coding=TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE;
      current_code_size=natoms*3;
      current_coding_parameter=0;
      coder=Ptngc_coder_init();
      if (!determine_best_coding_stop_bits(coder,quant,&current_code_size,
                                           &current_coding_parameter,natoms))
        {
          best_coding=current_coding;
          best_coding_parameter=current_coding_parameter;
          best_code_size=current_code_size;
        }
      Ptngc_coder_deinit(coder);

      /* Determine best parameter for triplet one-to-one. */
      current_coding=TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE;
      coder=Ptngc_coder_init();
      current_code_size=natoms*3;
      current_coding_parameter=0;
      if (!determine_best_coding_triple(coder,quant,&current_code_size,&current_coding_parameter,natoms))
        {
          if ((best_coding==-1) || (current_code_size<best_code_size))
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
              best_code_size=current_code_size;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Test BWLZH one-to-one */
      if (speed>=4)
        {
          current_coding=TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE;
          current_coding_parameter=0;
          compress_quantized_vel(quant,NULL,natoms,1,speed,
                                 current_coding,current_coding_parameter,
                                 0,0,prec_hi,prec_lo,&current_code_size,NULL);
          if ((best_coding==-1) || (current_code_size<best_code_size))
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
            }
        }
      *initial_coding=best_coding;
      *initial_coding_parameter=best_coding_parameter;
    }
  else if (*initial_coding_parameter==-1)
    {
      if (*initial_coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE)
        *initial_coding_parameter=0;
      else if (*initial_coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3;
          determine_best_coding_stop_bits(coder,quant,&current_code_size,
                                          initial_coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*initial_coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3;
          determine_best_coding_triple(coder,quant,&current_code_size,initial_coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
    }
}

static void determine_best_vel_coding(int *quant, int *quant_inter, const int natoms, const int nframes, const int speed,
                                      const fix_t prec_hi, const fix_t prec_lo,
                                      int *coding, int *coding_parameter)
{
  if (*coding==-1)
    {
      /* Determine all parameters automatically */
      int best_coding;
      int best_coding_parameter;
      int best_code_size;
      int current_coding;
      int current_coding_parameter;
      int current_code_size;
      int initial_code_size;
      int initial_numbits=5;
      struct coder *coder;
      /* Use stopbits one-to-one coding for the initial coding. */
      compress_quantized_vel(quant,NULL,natoms,1,speed,
                             TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE,initial_numbits,
                             0,0,prec_hi,prec_lo,&initial_code_size,NULL);

      /* Test stopbit one-to-one */
      current_coding=TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE;
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      coder=Ptngc_coder_init();
      determine_best_coding_stop_bits(coder,quant+natoms*3,&current_code_size,
                                      &current_coding_parameter,natoms);
      Ptngc_coder_deinit(coder);
      best_coding=current_coding;
      best_code_size=current_code_size;
      best_coding_parameter=current_coding_parameter;

      /* Test triplet interframe */
      current_coding=TNG_COMPRESS_ALGO_VEL_TRIPLET_INTER;
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      coder=Ptngc_coder_init();
      if (!determine_best_coding_triple(coder,quant_inter+natoms*3,&current_code_size,
                                        &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_code_size=current_code_size;
              best_coding_parameter=current_coding_parameter;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Test triplet one-to-one */
      current_coding=TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE;
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      coder=Ptngc_coder_init();
      if (!determine_best_coding_triple(coder,quant+natoms*3,&current_code_size,
                                        &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_code_size=current_code_size;
              best_coding_parameter=current_coding_parameter;
            }
        }
      Ptngc_coder_deinit(coder);

      /* Test stopbit interframe */
      current_coding=TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER;
      current_code_size=natoms*3*(nframes-1);
      current_coding_parameter=0;
      coder=Ptngc_coder_init();
      if (!determine_best_coding_stop_bits(coder,quant_inter+natoms*3,&current_code_size,
                                           &current_coding_parameter,natoms))
        {
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_code_size=current_code_size;
              best_coding_parameter=current_coding_parameter;
            }
        }
      Ptngc_coder_deinit(coder);

      if (speed>=4)
        {
          /* Test BWLZH inter */
          current_coding=TNG_COMPRESS_ALGO_VEL_BWLZH_INTER;
          current_coding_parameter=0;
          compress_quantized_vel(quant,quant_inter,natoms,nframes,speed,
                                 TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE,initial_numbits,
                                 current_coding,current_coding_parameter,
                                 prec_hi,prec_lo,&current_code_size,NULL);
          current_code_size-=initial_code_size; /* Correct for the initial frame */
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_code_size=current_code_size;
              best_coding_parameter=current_coding_parameter;
            }

          /* Test BWLZH one-to-one */
          current_coding=TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE;
          current_coding_parameter=0;
          compress_quantized_vel(quant,quant_inter,natoms,nframes,speed,
                                 TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE,initial_numbits,
                                 current_coding,current_coding_parameter,
                                 prec_hi,prec_lo,&current_code_size,NULL);
          current_code_size-=initial_code_size; /* Correct for the initial frame */
          if (current_code_size<best_code_size)
            {
              best_coding=current_coding;
              best_coding_parameter=current_coding_parameter;
            }
        }
      *coding=best_coding;
      *coding_parameter=best_coding_parameter;
    }
  else if (*coding_parameter==-1)
    {
      if ((*coding==TNG_COMPRESS_ALGO_VEL_BWLZH_INTER) ||
          (*coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE))
        *coding_parameter=0;
      else if (*coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_stop_bits(coder,quant+natoms*3,&current_code_size,
                                          coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_INTER)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_triple(coder,quant_inter+natoms*3,&current_code_size,
                                       coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_triple(coder,quant+natoms*3,&current_code_size,
                                       coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
      else if (*coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER)
        {
          struct coder *coder=Ptngc_coder_init();
          int current_code_size=natoms*3*(nframes-1);
          determine_best_coding_stop_bits(coder,quant_inter+natoms*3,&current_code_size,
                                          coding_parameter,natoms);
          Ptngc_coder_deinit(coder);
        }
    }
}

char DECLSPECDLLEXPORT *tng_compress_pos_int(int *pos, const int natoms, const int nframes,
                                             const unsigned long prec_hi, const unsigned long prec_lo,
                                             int speed,int *algo,
                                             int *nitems)
{
  char *data=malloc(natoms*nframes*14+11*4); /* 12 bytes are required to store 4 32 bit integers
                                               This is 17% extra. The final 11*4 is to store information
                                               needed for decompression. */
  int *quant=pos; /* Already quantized positions. */
  int *quant_intra=malloc(natoms*nframes*3*sizeof *quant_intra);
  int *quant_inter=malloc(natoms*nframes*3*sizeof *quant_inter);

  int initial_coding, initial_coding_parameter;
  int coding, coding_parameter;
  if (speed==0)
    speed=SPEED_DEFAULT; /* Set the default speed. */
  /* Boundaries of speed. */
  if (speed<1)
    speed=1;
  if (speed>6)
    speed=6;
  initial_coding=algo[0];
  initial_coding_parameter=algo[1];
  coding=algo[2];
  coding_parameter=algo[3];

  quant_inter_differences(quant,natoms,nframes,quant_inter);
  quant_intra_differences(quant,natoms,nframes,quant_intra);

  /* If any of the above codings / coding parameters are == -1, the optimal parameters must be found */
  if (initial_coding==-1)
    {
      initial_coding_parameter=-1;
      determine_best_pos_initial_coding(quant,quant_intra,natoms,speed,prec_hi,prec_lo,
                                        &initial_coding,&initial_coding_parameter);
    }
  else if (initial_coding_parameter==-1)
    {
      determine_best_pos_initial_coding(quant,quant_intra,natoms,speed,prec_hi,prec_lo,
                                        &initial_coding,&initial_coding_parameter);
    }

  if (nframes==1)
    {
      coding=0;
      coding_parameter=0;
    }

  if (nframes>1)
    {
      if (coding==-1)
        {
          coding_parameter=-1;
          determine_best_pos_coding(quant,quant_inter,quant_intra,natoms,nframes,speed,prec_hi,prec_lo,
                                    &coding,&coding_parameter);
        }
      else if (coding_parameter==-1)
        {
          determine_best_pos_coding(quant,quant_inter,quant_intra,natoms,nframes,speed,prec_hi,prec_lo,
                                    &coding,&coding_parameter);
        }
    }

  compress_quantized_pos(quant,quant_inter,quant_intra,natoms,nframes,speed,
                         initial_coding,initial_coding_parameter,
                         coding,coding_parameter,
                         prec_hi,prec_lo,nitems,data);
  free(quant_inter);
  free(quant_intra);
  if (algo[0]==-1)
    algo[0]=initial_coding;
  if (algo[1]==-1)
    algo[1]=initial_coding_parameter;
  if (algo[2]==-1)
    algo[2]=coding;
  if (algo[3]==-1)
    algo[3]=coding_parameter;
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_pos(double *pos, const int natoms, const int nframes,
                                         const double desired_precision,
                                         const int speed,int *algo,
                                         int *nitems)
{
  int *quant=malloc(natoms*nframes*3*sizeof *quant);
  char *data;
  fix_t prec_hi, prec_lo;
  Ptngc_d_to_i32x2(desired_precision,&prec_hi,&prec_lo);

  if (quantize(pos,natoms,nframes,PRECISION(prec_hi,prec_lo),quant))
    data=NULL; /* Error occured. Too large input values. */
  else
    data=tng_compress_pos_int(quant,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
  free(quant);
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_pos_float(float *pos, const int natoms, const int nframes,
                                               const float desired_precision,
                                               const int speed, int *algo,
                                               int *nitems)
{
  int *quant=malloc(natoms*nframes*3*sizeof *quant);
  char *data;
  fix_t prec_hi, prec_lo;
  Ptngc_d_to_i32x2((double)desired_precision,&prec_hi,&prec_lo);

  if (quantize_float(pos,natoms,nframes,(float)PRECISION(prec_hi,prec_lo),quant))
    data=NULL; /* Error occured. Too large input values. */
  else
    data=tng_compress_pos_int(quant,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
  free(quant);
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_pos_find_algo(double *pos, const int natoms, const int nframes,
                                                   const double desired_precision,
                                                   const int speed,
                                                   int *algo,
                                                   int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_pos(pos,natoms,nframes,desired_precision,speed,algo,nitems);
}

char DECLSPECDLLEXPORT *tng_compress_pos_float_find_algo(float *pos, const int natoms, const int nframes,
                                                         const float desired_precision,
                                                         const int speed,
                                                         int *algo,
                                                         int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_pos_float(pos,natoms,nframes,desired_precision,speed,algo,nitems);
}

char DECLSPECDLLEXPORT *tng_compress_pos_int_find_algo(int *pos, const int natoms, const int nframes,
                                                       const unsigned long prec_hi, const unsigned long prec_lo,
                                                       const int speed, int *algo,
                                                       int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_pos_int(pos,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
}



int DECLSPECDLLEXPORT tng_compress_nalgo(void)
{
  return 4; /* There are currently four parameters required:

 1) The compression algorithm for the first frame (initial_coding).
 2) One parameter to the algorithm for the first frame (the initial coding parameter).
 3) The compression algorithm for the remaining frames (coding).
 4) One parameter to the algorithm for the remaining frames (the coding parameter). */
}

char DECLSPECDLLEXPORT *tng_compress_vel_int(int *vel, const int natoms, const int nframes,
                                             const unsigned long prec_hi, const unsigned long prec_lo,
                                             int speed, int *algo,
                                             int *nitems)
{
  char *data=malloc(natoms*nframes*14+11*4); /* 12 bytes are required to store 4 32 bit integers
                                               This is 17% extra. The final 11*4 is to store information
                                               needed for decompression. */
  int *quant=vel;
  int *quant_inter=malloc(natoms*nframes*3*sizeof *quant_inter);

  int initial_coding, initial_coding_parameter;
  int coding, coding_parameter;
  if (speed==0)
    speed=SPEED_DEFAULT; /* Set the default speed. */
  /* Boundaries of speed. */
  if (speed<1)
    speed=1;
  if (speed>6)
    speed=6;
  initial_coding=algo[0];
  initial_coding_parameter=algo[1];
  coding=algo[2];
  coding_parameter=algo[3];

  quant_inter_differences(quant,natoms,nframes,quant_inter);

  /* If any of the above codings / coding parameters are == -1, the optimal parameters must be found */
  if (initial_coding==-1)
    {
      initial_coding_parameter=-1;
      determine_best_vel_initial_coding(quant,natoms,speed,prec_hi,prec_lo,
                                        &initial_coding,&initial_coding_parameter);
    }
  else if (initial_coding_parameter==-1)
    {
      determine_best_vel_initial_coding(quant,natoms,speed,prec_hi,prec_lo,
                                        &initial_coding,&initial_coding_parameter);
    }

  if (nframes==1)
    {
      coding=0;
      coding_parameter=0;
    }

  if (nframes>1)
    {
      if (coding==-1)
        {
          coding_parameter=-1;
          determine_best_vel_coding(quant,quant_inter,natoms,nframes,speed,prec_hi,prec_lo,
                                    &coding,&coding_parameter);
        }
      else if (coding_parameter==-1)
        {
          determine_best_vel_coding(quant,quant_inter,natoms,nframes,speed,prec_hi,prec_lo,
                                    &coding,&coding_parameter);
        }
    }

  compress_quantized_vel(quant,quant_inter,natoms,nframes,speed,
                         initial_coding,initial_coding_parameter,
                         coding,coding_parameter,
                         prec_hi,prec_lo,nitems,data);
  free(quant_inter);
  if (algo[0]==-1)
    algo[0]=initial_coding;
  if (algo[1]==-1)
    algo[1]=initial_coding_parameter;
  if (algo[2]==-1)
    algo[2]=coding;
  if (algo[3]==-1)
    algo[3]=coding_parameter;
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_vel(double *vel, const int natoms, const int nframes,
                                         const double desired_precision,
                                         const int speed, int *algo,
                                         int *nitems)
{
  int *quant=malloc(natoms*nframes*3*sizeof *quant);
  char *data;
  fix_t prec_hi, prec_lo;
  Ptngc_d_to_i32x2(desired_precision,&prec_hi,&prec_lo);
  if (quantize(vel,natoms,nframes,PRECISION(prec_hi,prec_lo),quant))
    data=NULL; /* Error occured. Too large input values. */
  else
    data=tng_compress_vel_int(quant,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
  free(quant);
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_vel_float(float *vel, const int natoms, const int nframes,
                                               const float desired_precision,
                                               const int speed, int *algo,
                                               int *nitems)
{
  int *quant=malloc(natoms*nframes*3*sizeof *quant);
  char *data;
  fix_t prec_hi, prec_lo;
  Ptngc_d_to_i32x2((double)desired_precision,&prec_hi,&prec_lo);
  if (quantize_float(vel,natoms,nframes,(float)PRECISION(prec_hi,prec_lo),quant))
    data=NULL; /* Error occured. Too large input values. */
  else
    data=tng_compress_vel_int(quant,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
  free(quant);
  return data;
}

char DECLSPECDLLEXPORT *tng_compress_vel_find_algo(double *vel, const int natoms, const int nframes,
                                                   const double desired_precision,
                                                   const int speed,
                                                   int *algo,
                                                   int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_vel(vel,natoms,nframes,desired_precision,speed,algo,nitems);
}

char DECLSPECDLLEXPORT *tng_compress_vel_float_find_algo(float *vel, const int natoms, const int nframes,
                                                         const float desired_precision,
                                                         const int speed,
                                                         int *algo,
                                                         int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_vel_float(vel,natoms,nframes,desired_precision,speed,algo,nitems);
}

char DECLSPECDLLEXPORT *tng_compress_vel_int_find_algo(int *vel, const int natoms, const int nframes,
                                                       const unsigned long prec_hi, const unsigned long prec_lo,
                                                       const int speed,
                                                       int *algo,
                                                       int *nitems)
{
  algo[0]=-1;
  algo[1]=-1;
  algo[2]=-1;
  algo[3]=-1;
  return tng_compress_vel_int(vel,natoms,nframes,prec_hi,prec_lo,speed,algo,nitems);
}

int DECLSPECDLLEXPORT tng_compress_inquire(char *data,int *vel, int *natoms,
                                            int *nframes, double *precision,
                                            int *algo)
{
  int bufloc=0;
  fix_t prec_hi, prec_lo;
  int initial_coding, initial_coding_parameter;
  int coding, coding_parameter;
  int magic_int;
  magic_int=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  if (magic_int==MAGIC_INT_POS)
    *vel=0;
  else if (magic_int==MAGIC_INT_VEL)
    *vel=1;
  else
    return 1;
  /* Number of atoms. */
  *natoms=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Number of frames. */
  *nframes=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding. */
  initial_coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding parameter. */
  initial_coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding. */
  coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding parameter. */
  coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Precision. */
  prec_lo=readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  prec_hi=readbufferfix((unsigned char *)data+bufloc,4);
  *precision=PRECISION(prec_hi, prec_lo);
  algo[0]=initial_coding;
  algo[1]=initial_coding_parameter;
  algo[2]=coding;
  algo[3]=coding_parameter;
  return 0;
}

static int tng_compress_uncompress_pos_gen(char *data,double *posd,float *posf,int *posi,unsigned long *prec_hi, unsigned long *prec_lo)
{
  int bufloc=0;
  int length;
  int natoms, nframes;
  int initial_coding, initial_coding_parameter;
  int coding, coding_parameter;
  int *quant=NULL;
  struct coder *coder=NULL;
  int rval=0;
  int magic_int;
  /* Magic integer for positions */
  magic_int=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  if (magic_int!=MAGIC_INT_POS)
    {
      rval=1;
      goto error;
    }
  /* Number of atoms. */
  natoms=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Number of frames. */
  nframes=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding. */
  initial_coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding parameter. */
  initial_coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding. */
  coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding parameter. */
  coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Precision. */
  *prec_lo=readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  *prec_hi=readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Allocate the memory for the quantized positions */
  quant=malloc(natoms*nframes*3*sizeof *quant);
  /* The data block length. */
  length=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* The initial frame */
  coder=Ptngc_coder_init();
  rval=Ptngc_unpack_array(coder,(unsigned char*)data+bufloc,quant,natoms*3,
                         initial_coding,initial_coding_parameter,natoms);
  Ptngc_coder_deinit(coder);
  if (rval)
    goto error;
  /* Skip past the actual data block. */
  bufloc+=length;
  /* Obtain the actual positions for the initial block. */
  if ((initial_coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
      (initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_POS_XTC3))
    {
      if (posd)
        unquantize(posd,natoms,1,PRECISION(*prec_hi,*prec_lo),quant);
      else if (posf)
        unquantize_float(posf,natoms,1,(float)PRECISION(*prec_hi,*prec_lo),quant);
      else if (posi)
        memcpy(posi,quant,natoms*3*sizeof *posi);
    }
  else if ((initial_coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
           (initial_coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
    {
      if (posd)
        unquantize_intra_differences(posd,natoms,1,PRECISION(*prec_hi,*prec_lo),quant);
      else if (posf)
        unquantize_intra_differences_float(posf,natoms,1,(float)PRECISION(*prec_hi,*prec_lo),quant);
      else if (posi)
        unquantize_intra_differences_int(posi,natoms,1,quant);
      unquant_intra_differences_first_frame(quant,natoms);
    }
  /* The remaining frames. */
  if (nframes>1)
    {
      bufloc+=4;
      coder=Ptngc_coder_init();
      rval=Ptngc_unpack_array(coder,(unsigned char *)data+bufloc,quant+natoms*3,(nframes-1)*natoms*3,
                                  coding,coding_parameter,natoms);
      Ptngc_coder_deinit(coder);
      if (rval)
        goto error;
      if ((coding==TNG_COMPRESS_ALGO_POS_STOPBIT_INTER) ||
          (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTER) ||
          (coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTER))
        {
          /* This requires that the first frame is already in one-to-one format, even if intra-frame
             compression was done there. Therefore the unquant_intra_differences_first_frame should be called
             before to convert it correctly. */
          if (posd)
            unquantize_inter_differences(posd,natoms,nframes,PRECISION(*prec_hi,*prec_lo),quant);
          else if (posf)
            unquantize_inter_differences_float(posf,natoms,nframes,(float)PRECISION(*prec_hi,*prec_lo),quant);
          else if (posi)
            unquantize_inter_differences_int(posi,natoms,nframes,quant);
        }
      else if ((coding==TNG_COMPRESS_ALGO_POS_XTC2) ||
               (coding==TNG_COMPRESS_ALGO_POS_XTC3) ||
               (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE))
        {
          if (posd)
            unquantize(posd+natoms*3,natoms,nframes-1,PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (posf)
            unquantize_float(posf+natoms*3,natoms,nframes-1,(float)PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (posi)
            memcpy(posi+natoms*3,quant+natoms*3,natoms*3*(nframes-1)*sizeof *posi);
        }
      else if ((coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
               (coding==TNG_COMPRESS_ALGO_POS_BWLZH_INTRA))
        {
          if (posd)
            unquantize_intra_differences(posd+natoms*3,natoms,nframes-1,PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (posf)
            unquantize_intra_differences_float(posf+natoms*3,natoms,nframes-1,(float)PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (posi)
            unquantize_intra_differences_int(posi+natoms*3,natoms,nframes-1,quant+natoms*3);
        }
    }
 error:
  free(quant);
  return rval;
}

static int tng_compress_uncompress_pos(char *data,double *pos)
{
  unsigned long prec_hi, prec_lo;
  return tng_compress_uncompress_pos_gen(data,pos,NULL,NULL,&prec_hi,&prec_lo);
}

static int tng_compress_uncompress_pos_float(char *data,float *pos)
{
  unsigned long prec_hi, prec_lo;
  return tng_compress_uncompress_pos_gen(data,NULL,pos,NULL,&prec_hi,&prec_lo);
}

static int tng_compress_uncompress_pos_int(char *data,int *pos, unsigned long *prec_hi, unsigned long *prec_lo)
{
  return tng_compress_uncompress_pos_gen(data,NULL,NULL,pos,prec_hi,prec_lo);
}

static int tng_compress_uncompress_vel_gen(char *data,double *veld,float *velf,int *veli,unsigned long *prec_hi, unsigned long *prec_lo)
{
  int bufloc=0;
  int length;
  int natoms, nframes;
  int initial_coding, initial_coding_parameter;
  int coding, coding_parameter;
  int *quant=NULL;
  struct coder *coder=NULL;
  int rval=0;
  int magic_int;
  /* Magic integer for velocities */
  magic_int=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  if (magic_int!=MAGIC_INT_VEL)
    {
      rval=1;
      goto error;
    }
  /* Number of atoms. */
  natoms=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Number of frames. */
  nframes=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding. */
  initial_coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Initial coding parameter. */
  initial_coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding. */
  coding=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Coding parameter. */
  coding_parameter=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Precision. */
  *prec_lo=readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  *prec_hi=readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* Allocate the memory for the quantized positions */
  quant=malloc(natoms*nframes*3*sizeof *quant);
  /* The data block length. */
  length=(int)readbufferfix((unsigned char *)data+bufloc,4);
  bufloc+=4;
  /* The initial frame */
  coder=Ptngc_coder_init();
  rval=Ptngc_unpack_array(coder,(unsigned char*)data+bufloc,quant,natoms*3,
                         initial_coding,initial_coding_parameter,natoms);
  Ptngc_coder_deinit(coder);
  if (rval)
    goto error;
  /* Skip past the actual data block. */
  bufloc+=length;
  /* Obtain the actual positions for the initial block. */
  if ((initial_coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE) ||
      (initial_coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE))
    {
      if (veld)
        unquantize(veld,natoms,1,PRECISION(*prec_hi,*prec_lo),quant);
      else if (velf)
        unquantize_float(velf,natoms,1,(float)PRECISION(*prec_hi,*prec_lo),quant);
      else if (veli)
        memcpy(veli,quant,natoms*3*sizeof *veli);
    }
  /* The remaining frames. */
  if (nframes>1)
    {
      bufloc+=4;
      coder=Ptngc_coder_init();
      rval=Ptngc_unpack_array(coder,(unsigned char *)data+bufloc,quant+natoms*3,(nframes-1)*natoms*3,
                                  coding,coding_parameter,natoms);
      Ptngc_coder_deinit(coder);
      if (rval)
        goto error;
      /* Inter-frame compression? */
      if ((coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_INTER) ||
          (coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER) ||
          (coding==TNG_COMPRESS_ALGO_VEL_BWLZH_INTER))
        {
          /* This requires that the first frame is already in one-to-one format. */
          if (veld)
            unquantize_inter_differences(veld,natoms,nframes,PRECISION(*prec_hi,*prec_lo),quant);
          else if (velf)
            unquantize_inter_differences_float(velf,natoms,nframes,(float)PRECISION(*prec_hi,*prec_lo),quant);
          else if (veli)
            unquantize_inter_differences_int(veli,natoms,nframes,quant);
        }
      /* One-to-one compression? */
      else if ((coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_ONETOONE) ||
               (coding==TNG_COMPRESS_ALGO_VEL_TRIPLET_ONETOONE) ||
               (coding==TNG_COMPRESS_ALGO_VEL_BWLZH_ONETOONE))
        {
          if (veld)
            unquantize(veld+natoms*3,natoms,nframes-1,PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (velf)
            unquantize_float(velf+natoms*3,natoms,nframes-1,(float)PRECISION(*prec_hi,*prec_lo),quant+natoms*3);
          else if (veli)
            memcpy(veli+natoms*3,quant+natoms*3,natoms*3*(nframes-1)*sizeof *veli);
        }
    }
 error:
  free(quant);
  return rval;
}

static int tng_compress_uncompress_vel(char *data,double *vel)
{
  unsigned long prec_hi, prec_lo;
  return tng_compress_uncompress_vel_gen(data,vel,NULL,NULL,&prec_hi,&prec_lo);
}

static int tng_compress_uncompress_vel_float(char *data,float *vel)
{
  unsigned long prec_hi, prec_lo;
  return tng_compress_uncompress_vel_gen(data,NULL,vel,NULL,&prec_hi,&prec_lo);
}

static int tng_compress_uncompress_vel_int(char *data,int *vel, unsigned long *prec_hi, unsigned long *prec_lo)
{
  return tng_compress_uncompress_vel_gen(data,NULL,NULL,vel,prec_hi,prec_lo);
}

/* Uncompresses any tng compress block, positions or velocities. It determines whether it is positions or velocities from the data buffer. The return value is 0 if ok, and 1 if not.
*/
int DECLSPECDLLEXPORT tng_compress_uncompress(char *data,double *posvel)
{
  int magic_int;
  magic_int=(int)readbufferfix((unsigned char *)data,4);
  if (magic_int==MAGIC_INT_POS)
    return tng_compress_uncompress_pos(data,posvel);
  else if (magic_int==MAGIC_INT_VEL)
    return tng_compress_uncompress_vel(data,posvel);
  else
    return 1;
}

int DECLSPECDLLEXPORT tng_compress_uncompress_float(char *data,float *posvel)
{
  int magic_int;
  magic_int=(int)readbufferfix((unsigned char *)data,4);
  if (magic_int==MAGIC_INT_POS)
    return tng_compress_uncompress_pos_float(data,posvel);
  else if (magic_int==MAGIC_INT_VEL)
    return tng_compress_uncompress_vel_float(data,posvel);
  else
    return 1;
}

int DECLSPECDLLEXPORT tng_compress_uncompress_int(char *data,int *posvel, unsigned long *prec_hi, unsigned long *prec_lo)
{
  int magic_int;
  magic_int=(int)readbufferfix((unsigned char *)data,4);
  if (magic_int==MAGIC_INT_POS)
    return tng_compress_uncompress_pos_int(data,posvel,prec_hi,prec_lo);
  else if (magic_int==MAGIC_INT_VEL)
    return tng_compress_uncompress_vel_int(data,posvel,prec_hi,prec_lo);
  else
    return 1;
}

void DECLSPECDLLEXPORT tng_compress_int_to_double(int *posvel_int, const unsigned long prec_hi, const unsigned long prec_lo,
                                                  const int natoms, const int nframes,
                                                  double *posvel_double)
{
  unquantize(posvel_double,natoms,nframes,PRECISION(prec_hi,prec_lo),posvel_int);
}

void DECLSPECDLLEXPORT tng_compress_int_to_float(int *posvel_int, const unsigned long prec_hi, const unsigned long prec_lo,
                                                 const int natoms, const int nframes,
                                                 float *posvel_float)
{
  unquantize_float(posvel_float,natoms,nframes,(float)PRECISION(prec_hi,prec_lo),posvel_int);
}

static char *compress_algo_pos[TNG_COMPRESS_ALGO_MAX]={
  "Positions invalid algorithm",
  "Positions stopbits interframe",
  "Positions triplet interframe",
  "Positions triplet intraframe",
  "Positions invalid algorithm",
  "Positions XTC2",
  "Positions invalid algorithm",
  "Positions triplet one to one",
  "Positions BWLZH interframe",
  "Positions BWLZH intraframe",
  "Positions XTC3"
};

static char *compress_algo_vel[TNG_COMPRESS_ALGO_MAX]={
  "Velocities invalid algorithm",
  "Velocities stopbits one to one",
  "Velocities triplet interframe",
  "Velocities triplet one to one",
  "Velocities invalid algorithm",
  "Velocities invalid algorithm",
  "Velocities stopbits interframe",
  "Velocities invalid algorithm",
  "Velocities BWLZH interframe",
  "Velocities BWLZH one to one",
  "Velocities invalid algorithm"
};

char DECLSPECDLLEXPORT *tng_compress_initial_pos_algo(int *algo)
{
  int i=algo[0];
  if (i<0)
    i=0;
  if (i>=TNG_COMPRESS_ALGO_MAX)
    i=0;
  return compress_algo_pos[i];
}

char DECLSPECDLLEXPORT *tng_compress_pos_algo(int *algo)
{
  int i=algo[2];
  if (i<0)
    i=0;
  if (i>=TNG_COMPRESS_ALGO_MAX)
    i=0;
  return compress_algo_pos[i];
}

char DECLSPECDLLEXPORT *tng_compress_initial_vel_algo(int *algo)
{
  int i=algo[0];
  if (i<0)
    i=0;
  if (i>=TNG_COMPRESS_ALGO_MAX)
    i=0;
  return compress_algo_vel[i];
}

char DECLSPECDLLEXPORT *tng_compress_vel_algo(int *algo)
{
  int i=algo[2];
  if (i<0)
    i=0;
  if (i>=TNG_COMPRESS_ALGO_MAX)
    i=0;
  return compress_algo_vel[i];
}
