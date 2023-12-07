/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2010,2013, The GROMACS development team.
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
 */

/* This code is part of the tng compression routines
 * Written by Daniel Spangberg
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/compression/warnmalloc.h"
#include "../../include/compression/tng_compress.h"
#include "../../include/compression/bwlzh.h"
#include "../../include/compression/dict.h"
#include "../../include/compression/vals16.h"
#include "../../include/compression/rle.h"
#include "../../include/compression/mtf.h"
#include "../../include/compression/bwt.h"
#include "../../include/compression/lz77.h"

#if 0
#    define SHOWIT
#endif

#if 0
#    define SHOWTEST
#endif

#define MAX_VALS_PER_BLOCK 200000

#if 1
#    define PARTIAL_MTF3
#endif

#if 0
#    define PARTIAL_MTF
#endif

int bwlzh_get_buflen(const int nvals)
{
    return 132000 + nvals * 8 + 12 * ((nvals + MAX_VALS_PER_BLOCK) / MAX_VALS_PER_BLOCK);
}

#ifdef SHOWIT
static void printvals(const char* name, unsigned int* vals, const int nvals)
{
    int i;
    int nvalsmax = nvals;
    if (nvalsmax > 99)
        nvalsmax = 99;
#    if 0
  for (i=0; i<nvalsmax; i++)
    fprintf(stderr,"%s %d %u\n",name,i,vals[i]);
#    else
    fprintf(stderr, "%s\n", name);
    {
        unsigned char* x = (unsigned char*)vals;
        for (i = 0; i < nvalsmax * 4; i++)
            fprintf(stderr, "%02x", (unsigned int)x[i]);
        fprintf(stderr, "\n");
    }
#    endif
}
#endif


static void bwlzh_compress_gen(unsigned int*  vals,
                               const int      nvals,
                               unsigned char* output,
                               int*           output_len,
                               const int      enable_lz77,
                               const int      verbose)
{
    unsigned int* vals16;
    int           nvals16;
    int           huffdatalen;
    int           nhufflen[N_HUFFMAN_ALGO];
    int           huffalgo;
    int           bwt_index;
    unsigned int* bwt = NULL;
#ifdef PARTIAL_MTF3
    unsigned char* mtf3 = NULL;
    int            imtfinner;
#endif
    unsigned int*  mtf     = NULL;
    unsigned int*  rle     = NULL;
    unsigned int*  offsets = NULL;
    unsigned int*  lens    = NULL;
    unsigned int*  dict    = warnmalloc(0x20004 * sizeof *dict);
    unsigned int*  hist    = warnmalloc(0x20004 * sizeof *hist);
    int            nrle;
    int            noffsets;
    int            nlens;
    unsigned char* bwlzhhuff = NULL;
    int            bwlzhhufflen;
    int            max_vals_per_block = MAX_VALS_PER_BLOCK;
    int            valsleft;
    int            thisvals;
    int            valstart;
    int            outdata = 0;

    unsigned int* tmpmem = warnmalloc(max_vals_per_block * 18 * sizeof *tmpmem);

#if 0
  verbose=1;
#endif

    bwlzhhuff = warnmalloc(Ptngc_comp_huff_buflen(3 * nvals));
    vals16    = tmpmem;
    bwt       = tmpmem + max_vals_per_block * 3;
    mtf       = tmpmem + max_vals_per_block * 6;
    rle       = tmpmem + max_vals_per_block * 9;
    offsets   = tmpmem + max_vals_per_block * 12;
    lens      = tmpmem + max_vals_per_block * 15;
#ifdef PARTIAL_MTF3
    mtf3 = warnmalloc(max_vals_per_block * 3 * 3
                      * sizeof *mtf3); /* 3 due to expansion of 32 bit to 16 bit, 3 due to up to 3
                                          bytes per 16 value. */
#endif
    if (verbose)
    {
        fprintf(stderr, "Number of input values: %d\n", nvals);
    }

    /* Store the number of real values in the whole block. */
    output[outdata++] = ((unsigned int)nvals) & 0xFFU;
    output[outdata++] = (((unsigned int)nvals) >> 8) & 0xFFU;
    output[outdata++] = (((unsigned int)nvals) >> 16) & 0xFFU;
    output[outdata++] = (((unsigned int)nvals) >> 24) & 0xFFU;

    valsleft = nvals;
    valstart = 0;
    while (valsleft)
    {
        int reducealgo = 1; /* Reduce algo is LZ77. */
        if (!enable_lz77)
        {
            reducealgo = 0;
        }
        thisvals = valsleft;
        if (thisvals > max_vals_per_block)
        {
            thisvals = max_vals_per_block;
        }
        valsleft -= thisvals;
        if (verbose)
        {
            fprintf(stderr, "Creating vals16 block from %d values.\n", thisvals);
        }

#ifdef SHOWIT
        printvals("vals", vals + valstart, thisvals);
#endif

        Ptngc_comp_conv_to_vals16(vals + valstart, thisvals, vals16, &nvals16);
        valstart += thisvals;

#ifdef SHOWTEST
        nvals16 = 99;
#endif

#ifdef SHOWIT
        printvals("vals16", vals16, nvals16);
#endif

        if (verbose)
        {
            fprintf(stderr, "Resulting vals16 values: %d\n", nvals16);
        }
        if (verbose)
        {
            fprintf(stderr, "BWT\n");
        }
        Ptngc_comp_to_bwt(vals16, nvals16, bwt, &bwt_index);

#ifdef SHOWIT
        printvals("bwt", bwt, nvals16);
        fprintf(stderr, "BWT INDEX is %d\n", bwt_index);
#endif

        /* Store the number of real values in this block. */
        output[outdata++] = ((unsigned int)thisvals) & 0xFFU;
        output[outdata++] = (((unsigned int)thisvals) >> 8) & 0xFFU;
        output[outdata++] = (((unsigned int)thisvals) >> 16) & 0xFFU;
        output[outdata++] = (((unsigned int)thisvals) >> 24) & 0xFFU;

        /* Store the number of nvals16 values in this block. */
        output[outdata++] = ((unsigned int)nvals16) & 0xFFU;
        output[outdata++] = (((unsigned int)nvals16) >> 8) & 0xFFU;
        output[outdata++] = (((unsigned int)nvals16) >> 16) & 0xFFU;
        output[outdata++] = (((unsigned int)nvals16) >> 24) & 0xFFU;

        /* Store the BWT index. */
        output[outdata++] = ((unsigned int)bwt_index) & 0xFFU;
        output[outdata++] = (((unsigned int)bwt_index) >> 8) & 0xFFU;
        output[outdata++] = (((unsigned int)bwt_index) >> 16) & 0xFFU;
        output[outdata++] = (((unsigned int)bwt_index) >> 24) & 0xFFU;

        if (verbose)
        {
            fprintf(stderr, "MTF\n");
        }
#ifdef PARTIAL_MTF3
        Ptngc_comp_conv_to_mtf_partial3(bwt, nvals16, mtf3);
        for (imtfinner = 0; imtfinner < 3; imtfinner++)
        {
            if (verbose)
            {
                fprintf(stderr, "Doing partial MTF: %d\n", imtfinner);
            }
            for (int j = 0; j < nvals16; j++)
            {
                mtf[j] = (unsigned int)mtf3[imtfinner * nvals16 + j];
            }
#else
#    ifdef PARTIAL_MTF
        Ptngc_comp_conv_to_mtf_partial(bwt, nvals16, mtf);
#    else
        int ndict;
        Ptngc_comp_canonical_dict(dict, &ndict);
        Ptngc_comp_conv_to_mtf(bwt, nvals16, dict, ndict, mtf);
#    endif

#    ifdef SHOWIT
        printvals("mtf", mtf, nvals16);
#    endif
#endif


            if (reducealgo == 1)
            {
                if (verbose)
                {
                    fprintf(stderr, "LZ77\n");
                }
                reducealgo = 1;
                Ptngc_comp_to_lz77(mtf, nvals16, rle, &nrle, lens, &nlens, offsets, &noffsets);

                if (verbose)
                {
                    fprintf(stderr, "Resulting LZ77 values: %d\n", nrle);
                    fprintf(stderr, "Resulting LZ77 lens: %d\n", nlens);
                    fprintf(stderr, "Resulting LZ77 offsets: %d\n", noffsets);
                }
#ifdef SHOWIT
                printvals("lz77 table", rle, nrle);
                printvals("lz77 lengths", lens, nlens);
                printvals("lz77 offsets", offsets, noffsets);
#endif

#if 0
              if (noffsets)
                {
                  unsigned int thist[0x20004];
                  unsigned int coarse[17]={0,};
                  int jj;
                  Ptngc_comp_make_dict_hist(lens,nlens,dict,&ndict,thist);
                  for (jj=0; jj<ndict; jj++)
                    fprintf(stderr,"%d %u %u L\n",jj,dict[jj],thist[jj]);

                  Ptngc_comp_make_dict_hist(offsets,noffsets,dict,&ndict,thist);
                  for (jj=0; jj<ndict; jj++)
                    {
                      unsigned int v=dict[jj];
                      int numbits=0;
                      while (v)
                        {
                          numbits++;
                          v>>=1;
                        }
                      coarse[numbits-1]+=thist[jj];
                    }
#    if 1
                  for (jj=0; jj<ndict; jj++)
                    fprintf(stderr,"%d %u %u O\n",jj,dict[jj],thist[jj]);
#    else
                  for (jj=0; jj<17; jj++)
                    fprintf(stderr,"%d %u\n",jj+1,coarse[jj]);
#    endif

                }
              exit(0);
#endif

                if (nlens < 2)
                {
                    reducealgo = 0;
                }

#ifdef SHOWTEST
                reducealgo = 1;
#endif
            }
            if (reducealgo == 0)
            {
                if (verbose)
                {
                    fprintf(stderr, "RLE\n");
                }
                /* Do RLE. For any repetetitive characters. */
                Ptngc_comp_conv_to_rle(mtf, nvals16, rle, &nrle, 1);

#ifdef SHOWIT
                printvals("rle", rle, nrle);
#endif
                if (verbose)
                {
                    fprintf(stderr, "Resulting RLE values: %d\n", nrle);
                }
            }

            /* reducealgo: RLE == 0, LZ77 == 1 */
            output[outdata++] = reducealgo;

            if (verbose)
            {
                fprintf(stderr, "Huffman\n");
            }

            huffalgo = -1;
            Ptngc_comp_huff_compress_verbose(rle, nrle, bwlzhhuff, &bwlzhhufflen, &huffdatalen,
                                             nhufflen, &huffalgo, 1);
#ifdef SHOWTEST
            {
                int i;
                fprintf(stderr, "Huffman\n");
                for (i = 0; i < bwlzhhufflen; i++)
                    fprintf(stderr, "%02x", (unsigned int)bwlzhhuff[i]);
                fprintf(stderr, "\n");
                exit(0);
            }
#endif
            if (verbose)
            {
                int i;
                fprintf(stderr, "Huffman data length is %d B.\n", huffdatalen);
                for (i = 0; i < N_HUFFMAN_ALGO; i++)
                {
                    fprintf(stderr, "Huffman dictionary for algorithm %s is %d B.\n",
                            Ptngc_comp_get_huff_algo_name(i), nhufflen[i] - huffdatalen);
                }
                fprintf(stderr, "Resulting algorithm: %s. Size=%d B\n",
                        Ptngc_comp_get_huff_algo_name(huffalgo), bwlzhhufflen);
            }

            /* Store the number of huffman values in this block. */
            output[outdata++] = ((unsigned int)nrle) & 0xFFU;
            output[outdata++] = (((unsigned int)nrle) >> 8) & 0xFFU;
            output[outdata++] = (((unsigned int)nrle) >> 16) & 0xFFU;
            output[outdata++] = (((unsigned int)nrle) >> 24) & 0xFFU;

            /* Store the size of the huffman block. */
            output[outdata++] = ((unsigned int)bwlzhhufflen) & 0xFFU;
            output[outdata++] = (((unsigned int)bwlzhhufflen) >> 8) & 0xFFU;
            output[outdata++] = (((unsigned int)bwlzhhufflen) >> 16) & 0xFFU;
            output[outdata++] = (((unsigned int)bwlzhhufflen) >> 24) & 0xFFU;

            /* Store the huffman block. */
            memcpy(output + outdata, bwlzhhuff, bwlzhhufflen);
            outdata += bwlzhhufflen;

            if (reducealgo == 1)
            {
                /* Store the number of values in this block. */
                output[outdata++] = ((unsigned int)noffsets) & 0xFFU;
                output[outdata++] = (((unsigned int)noffsets) >> 8) & 0xFFU;
                output[outdata++] = (((unsigned int)noffsets) >> 16) & 0xFFU;
                output[outdata++] = (((unsigned int)noffsets) >> 24) & 0xFFU;

                if (noffsets > 0)
                {
                    if (verbose)
                    {
                        fprintf(stderr, "Huffman for offsets\n");
                    }

                    huffalgo = -1;
                    Ptngc_comp_huff_compress_verbose(offsets, noffsets, bwlzhhuff, &bwlzhhufflen,
                                                     &huffdatalen, nhufflen, &huffalgo, 1);
                    if (verbose)
                    {
                        int i;
                        fprintf(stderr, "Huffman data length is %d B.\n", huffdatalen);
                        for (i = 0; i < N_HUFFMAN_ALGO; i++)
                        {
                            fprintf(stderr, "Huffman dictionary for algorithm %s is %d B.\n",
                                    Ptngc_comp_get_huff_algo_name(i), nhufflen[i] - huffdatalen);
                        }
                        fprintf(stderr, "Resulting algorithm: %s. Size=%d B\n",
                                Ptngc_comp_get_huff_algo_name(huffalgo), bwlzhhufflen);
                    }

                    /* If huffman was bad for these offsets, just store the offsets as pairs of bytes. */
                    if (bwlzhhufflen < noffsets * 2)
                    {
                        output[outdata++] = 0;

                        /* Store the size of the huffman block. */
                        output[outdata++] = ((unsigned int)bwlzhhufflen) & 0xFFU;
                        output[outdata++] = (((unsigned int)bwlzhhufflen) >> 8) & 0xFFU;
                        output[outdata++] = (((unsigned int)bwlzhhufflen) >> 16) & 0xFFU;
                        output[outdata++] = (((unsigned int)bwlzhhufflen) >> 24) & 0xFFU;

                        /* Store the huffman block. */
                        memcpy(output + outdata, bwlzhhuff, bwlzhhufflen);
                        outdata += bwlzhhufflen;
                    }
                    else
                    {
                        int i;
                        output[outdata++] = 1;
                        for (i = 0; i < noffsets; i++)
                        {
                            output[outdata++] = (offsets[i]) & 0xFFU;
                            output[outdata++] = ((offsets[i]) >> 8) & 0xFFU;
                        }
                        if (verbose)
                        {
                            fprintf(stderr, "Store raw offsets: %d B\n", noffsets * 2);
                        }
                    }
                }

#if 0
              {
                int i,ndict;
                FILE *f=fopen("len.dict","w");
                Ptngc_comp_make_dict_hist(lens,nlens,dict,&ndict,hist);
                for (i=0; i<ndict; i++)
                  fprintf(f,"%d %d %d\n",i,dict[i],hist[i]);
                fclose(f);
                f=fopen("off.dict","w");
                Ptngc_comp_make_dict_hist(offsets,noffsets,dict,&ndict,hist);
                for (i=0; i<ndict; i++)
                  fprintf(f,"%d %d %d\n",i,dict[i],hist[i]);
                fclose(f);
                f=fopen("len.time","w");
                for (i=0; i<ndict; i++)
                  fprintf(f,"%d\n",lens[i]);
                fclose(f);
                f=fopen("off.time","w");
                for (i=0; i<ndict; i++)
                  fprintf(f,"%d\n",offsets[i]);
                fclose(f);
              }
#endif

                if (verbose)
                {
                    fprintf(stderr, "Huffman for lengths\n");
                }

                huffalgo = -1;
                Ptngc_comp_huff_compress_verbose(lens, nlens, bwlzhhuff, &bwlzhhufflen,
                                                 &huffdatalen, nhufflen, &huffalgo, 1);
                if (verbose)
                {
                    int i;
                    fprintf(stderr, "Huffman data length is %d B.\n", huffdatalen);
                    for (i = 0; i < N_HUFFMAN_ALGO; i++)
                    {
                        fprintf(stderr, "Huffman dictionary for algorithm %s is %d B.\n",
                                Ptngc_comp_get_huff_algo_name(i), nhufflen[i] - huffdatalen);
                    }
                    fprintf(stderr, "Resulting algorithm: %s. Size=%d B\n",
                            Ptngc_comp_get_huff_algo_name(huffalgo), bwlzhhufflen);
                }

                /* Store the number of values in this block. */
                output[outdata++] = ((unsigned int)nlens) & 0xFFU;
                output[outdata++] = (((unsigned int)nlens) >> 8) & 0xFFU;
                output[outdata++] = (((unsigned int)nlens) >> 16) & 0xFFU;
                output[outdata++] = (((unsigned int)nlens) >> 24) & 0xFFU;

                /* Store the size of the huffman block. */
                output[outdata++] = ((unsigned int)bwlzhhufflen) & 0xFFU;
                output[outdata++] = (((unsigned int)bwlzhhufflen) >> 8) & 0xFFU;
                output[outdata++] = (((unsigned int)bwlzhhufflen) >> 16) & 0xFFU;
                output[outdata++] = (((unsigned int)bwlzhhufflen) >> 24) & 0xFFU;

                /* Store the huffman block. */
                memcpy(output + outdata, bwlzhhuff, bwlzhhufflen);
                outdata += bwlzhhufflen;
            }
#ifdef PARTIAL_MTF3
        }
#endif
    }

    *output_len = outdata;
    free(hist);
    free(dict);
    free(bwlzhhuff);
#ifdef PARTIAL_MTF3
    free(mtf3);
#endif
    free(tmpmem);
}


void DECLSPECDLLEXPORT bwlzh_compress(unsigned int* vals, const int nvals, unsigned char* output, int* output_len)
{
    bwlzh_compress_gen(vals, nvals, output, output_len, 1, 0);
}

void DECLSPECDLLEXPORT bwlzh_compress_verbose(unsigned int*  vals,
                                              const int      nvals,
                                              unsigned char* output,
                                              int*           output_len)
{
    bwlzh_compress_gen(vals, nvals, output, output_len, 1, 1);
}


void DECLSPECDLLEXPORT bwlzh_compress_no_lz77(unsigned int*  vals,
                                              const int      nvals,
                                              unsigned char* output,
                                              int*           output_len)
{
    bwlzh_compress_gen(vals, nvals, output, output_len, 0, 0);
}

void DECLSPECDLLEXPORT bwlzh_compress_no_lz77_verbose(unsigned int*  vals,
                                                      const int      nvals,
                                                      unsigned char* output,
                                                      int*           output_len)
{
    bwlzh_compress_gen(vals, nvals, output, output_len, 0, 1);
}


static void bwlzh_decompress_gen(unsigned char* input, const int nvals, unsigned int* vals, const int verbose)
{
    unsigned int* vals16;
    int           nvals16;
    int           bwt_index;
    unsigned int* bwt = NULL;
    unsigned int* mtf = NULL;
#ifdef PARTIAL_MTF3
    unsigned char* mtf3 = NULL;
    int            imtfinner;
#endif
    unsigned int*  rle     = NULL;
    unsigned int*  offsets = NULL;
    unsigned int*  lens    = NULL;
    unsigned int*  dict    = warnmalloc(0x20004 * sizeof *dict);
    unsigned int*  hist    = warnmalloc(0x20004 * sizeof *hist);
    int            nrle, noffsets, nlens;
    unsigned char* bwlzhhuff = NULL;
    int            bwlzhhufflen;
    int            max_vals_per_block = MAX_VALS_PER_BLOCK;
    int            valsleft;
    int            thisvals;
    int            valstart;
    int            inpdata = 0;
    int            nvalsfile;

    unsigned int* tmpmem = warnmalloc(max_vals_per_block * 18 * sizeof *tmpmem);

#if 0
  verbose=1;
#endif


    bwlzhhuff = warnmalloc(Ptngc_comp_huff_buflen(3 * nvals));
    vals16    = tmpmem;
    bwt       = tmpmem + max_vals_per_block * 3;
    mtf       = tmpmem + max_vals_per_block * 6;
    rle       = tmpmem + max_vals_per_block * 9;
    offsets   = tmpmem + max_vals_per_block * 12;
    lens      = tmpmem + max_vals_per_block * 15;
#ifdef PARTIAL_MTF3
    mtf3 = warnmalloc(max_vals_per_block * 3 * 3
                      * sizeof *mtf3); /* 3 due to expansion of 32 bit to 16 bit, 3 due to up to 3
                                          bytes per 16 value. */
#endif

    if (verbose)
    {
        fprintf(stderr, "Number of input values: %d\n", nvals);
    }

    /* Read the number of real values in the whole block. */
    nvalsfile = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                      | (((unsigned int)input[inpdata + 2]) << 16)
                      | (((unsigned int)input[inpdata + 3]) << 24));
    inpdata += 4;

    if (nvalsfile != nvals)
    {
        fprintf(stderr,
                "BWLZH: The number of values found in the file is different from the number of "
                "values expected.\n");
        exit(EXIT_FAILURE);
    }

    valsleft = nvals;
    valstart = 0;
    while (valsleft)
    {
        int valsnew;
        int reducealgo;
        /* Read the number of real values in this block. */
        thisvals = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                         | (((unsigned int)input[inpdata + 2]) << 16)
                         | (((unsigned int)input[inpdata + 3]) << 24));
        inpdata += 4;

        valsleft -= thisvals;

        /* Read the number of nvals16 values in this block. */
        nvals16 = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                        | (((unsigned int)input[inpdata + 2]) << 16)
                        | (((unsigned int)input[inpdata + 3]) << 24));
        inpdata += 4;

        /* Read the BWT index. */
        bwt_index = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                          | (((unsigned int)input[inpdata + 2]) << 16)
                          | (((unsigned int)input[inpdata + 3]) << 24));
        inpdata += 4;

        if (thisvals > max_vals_per_block)
        {
            /* More memory must be allocated for decompression. */
            max_vals_per_block = thisvals;
            if (verbose)
            {
                fprintf(stderr, "Allocating more memory: %d B\n",
                        (int)(max_vals_per_block * 15 * sizeof *tmpmem));
            }
            tmpmem  = warnrealloc(tmpmem, max_vals_per_block * 18 * sizeof *tmpmem);
            vals16  = tmpmem;
            bwt     = tmpmem + max_vals_per_block * 3;
            mtf     = tmpmem + max_vals_per_block * 6;
            rle     = tmpmem + max_vals_per_block * 9;
            offsets = tmpmem + max_vals_per_block * 12;
            lens    = tmpmem + max_vals_per_block * 15;
#ifdef PARTIAL_MTF3
            mtf3 = warnrealloc(mtf3, max_vals_per_block * 3 * 3
                                             * sizeof *mtf3); /* 3 due to expansion of 32 bit to 16 bit, 3
                                                                 due to up to 3 bytes per 16 value. */
#endif
        }

#ifdef PARTIAL_MTF3
        for (imtfinner = 0; imtfinner < 3; imtfinner++)
        {
            if (verbose)
            {
                fprintf(stderr, "Doing partial MTF: %d\n", imtfinner);
            }
#endif

            reducealgo = (int)input[inpdata];
            inpdata++;

            /* Read the number of huffman values in this block. */
            nrle = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                         | (((unsigned int)input[inpdata + 2]) << 16)
                         | (((unsigned int)input[inpdata + 3]) << 24));
            inpdata += 4;

            /* Read the size of the huffman block. */
            bwlzhhufflen = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                                 | (((unsigned int)input[inpdata + 2]) << 16)
                                 | (((unsigned int)input[inpdata + 3]) << 24));
            inpdata += 4;

            if (verbose)
            {
                fprintf(stderr, "Decompressing huffman block of length %d.\n", bwlzhhufflen);
            }
            /* Decompress the huffman block. */
            Ptngc_comp_huff_decompress(input + inpdata, bwlzhhufflen, rle);
            inpdata += bwlzhhufflen;

            if (reducealgo == 1) /* LZ77 */
            {
                int offstore;
                /* Read the number of huffman values in this block. */
                noffsets = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                                 | (((unsigned int)input[inpdata + 2]) << 16)
                                 | (((unsigned int)input[inpdata + 3]) << 24));
                inpdata += 4;

                if (noffsets > 0)
                {
                    /* How are the offsets stored? */
                    offstore = (int)input[inpdata++];
                    if (offstore == 0)
                    {
                        /* Read the size of the huffman block. */
                        bwlzhhufflen = (int)(((unsigned int)input[inpdata])
                                             | (((unsigned int)input[inpdata + 1]) << 8)
                                             | (((unsigned int)input[inpdata + 2]) << 16)
                                             | (((unsigned int)input[inpdata + 3]) << 24));
                        inpdata += 4;

                        if (verbose)
                        {
                            fprintf(stderr, "Decompressing offset huffman block.\n");
                        }

                        /* Decompress the huffman block. */
                        Ptngc_comp_huff_decompress(input + inpdata, bwlzhhufflen, offsets);
                        inpdata += bwlzhhufflen;
                    }
                    else
                    {
                        if (verbose)
                        {
                            fprintf(stderr, "Reading offset block.\n");
                        }
                        for (int i = 0; i < noffsets; i++)
                        {
                            offsets[i] = (int)(((unsigned int)input[inpdata])
                                               | (((unsigned int)input[inpdata + 1]) << 8));
                            inpdata += 2;
                        }
                    }
                }
#if 0
              {
                int i;
                for (i=0; i<nrle; i++)
                  fprintf(stderr,"RLE %d: %d\n",i,rle[i]);
                for (i=0; i<noffsets; i++)
                  fprintf(stderr,"OFFSET %d: %d\n",i,offsets[i]);
              }
#endif


                /* Read the number of huffman values in this block. */
                nlens = (int)(((unsigned int)input[inpdata]) | (((unsigned int)input[inpdata + 1]) << 8)
                              | (((unsigned int)input[inpdata + 2]) << 16)
                              | (((unsigned int)input[inpdata + 3]) << 24));
                inpdata += 4;

                /* Read the size of the huffman block. */
                bwlzhhufflen = (int)(((unsigned int)input[inpdata])
                                     | (((unsigned int)input[inpdata + 1]) << 8)
                                     | (((unsigned int)input[inpdata + 2]) << 16)
                                     | (((unsigned int)input[inpdata + 3]) << 24));
                inpdata += 4;

                if (verbose)
                {
                    fprintf(stderr, "Decompressing length huffman block.\n");
                }

                /* Decompress the huffman block. */
                Ptngc_comp_huff_decompress(input + inpdata, bwlzhhufflen, lens);
                inpdata += bwlzhhufflen;

                if (verbose)
                {
                    fprintf(stderr, "Decompressing LZ77.\n");
                }

                Ptngc_comp_from_lz77(rle, nrle, lens, nlens, offsets, noffsets, mtf, nvals16);
            }
            else if (reducealgo == 0) /* RLE */
            {
#ifdef SHOWIT
                printvals("rle", rle, nrle);
#endif

                if (verbose)
                {
                    fprintf(stderr, "Decompressing rle block.\n");
                }
                Ptngc_comp_conv_from_rle(rle, mtf, nvals16);
            }

#ifdef PARTIAL_MTF3
            for (int i = 0; i < nvals16; i++)
            {
                mtf3[imtfinner * nvals16 + i] = (unsigned char)mtf[i];
            }
        }
#else
#    ifdef SHOWIT
        printvals("mtf", mtf, nvals16);
#    endif

#endif


        if (verbose)
        {
            fprintf(stderr, "Inverse MTF.\n");
        }
#ifdef PARTIAL_MTF3
        Ptngc_comp_conv_from_mtf_partial3(mtf3, nvals16, bwt);
#else
#    ifdef PARTIAL_MTF
        Ptngc_comp_conv_from_mtf_partial(mtf, nvals16, bwt);
#    else
        int ndict;
        Ptngc_comp_canonical_dict(dict, &ndict);
        Ptngc_comp_conv_from_mtf(mtf, nvals16, dict, ndict, bwt);
#    endif
#endif

#ifdef SHOWIT
        printvals("bwt", bwt, nvals16);
        fprintf(stderr, "BWT INDEX is %d\n", bwt_index);
#endif

        if (verbose)
        {
            fprintf(stderr, "Inverse BWT.\n");
        }
        Ptngc_comp_from_bwt(bwt, nvals16, bwt_index, vals16);

#ifdef SHOWIT
        printvals("vals16", vals16, nvals16);
#endif

        if (verbose)
        {
            fprintf(stderr, "Decompressing vals16 block.\n");
        }
        Ptngc_comp_conv_from_vals16(vals16, nvals16, vals + valstart, &valsnew);

#ifdef SHOWIT
        printvals("vals", vals + valstart, thisvals);
#endif

        if (valsnew != thisvals)
        {
            fprintf(stderr, "BWLZH: Block contained different number of values than expected.\n");
            exit(EXIT_FAILURE);
        }
        valstart += thisvals;
    }
    free(hist);
    free(dict);
    free(bwlzhhuff);
#ifdef PARTIAL_MTF3
    free(mtf3);
#endif
    free(tmpmem);
}


void DECLSPECDLLEXPORT bwlzh_decompress(unsigned char* input, const int nvals, unsigned int* vals)
{
    bwlzh_decompress_gen(input, nvals, vals, 0);
}

void DECLSPECDLLEXPORT bwlzh_decompress_verbose(unsigned char* input, const int nvals, unsigned int* vals)
{
    bwlzh_decompress_gen(input, nvals, vals, 1);
}
