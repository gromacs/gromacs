/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/utility/futil.h"

/* This is just for clarity - it can never be anything but 4! */
#define XDR_INT_SIZE 4

/* same order as the definition of xdr_datatype */
const char *xdr_datatype_names[] =
{
    "int",
    "float",
    "double",
    "large int",
    "char",
    "string"
};


/*___________________________________________________________________________
 |
 | what follows are the C routine to read/write compressed coordinates together
 | with some routines to assist in this task (those are marked
 | static and cannot be called from user programs)
 */
#define MAXABS INT_MAX-2

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
static const int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216
};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


/*____________________________________________________________________________
 |
 | sendbits - encode num into buf using the specified number of bits
 |
 | This routines appends the value of num to the bits already present in
 | the array buf. You need to give it the number of bits to use and you
 | better make sure that this number of bits is enough to hold the value
 | Also num must be positive.
 |
 */

static void sendbits(int buf[], int num_of_bits, int num)
{

    unsigned int    cnt, lastbyte;
    int             lastbits;
    unsigned char * cbuf;

    cbuf     = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt      = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte = (unsigned int) buf[2];
    while (num_of_bits >= 8)
    {
        lastbyte     = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
        cbuf[cnt++]  = lastbyte >> lastbits;
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        lastbyte  = (lastbyte << num_of_bits) | num;
        lastbits += num_of_bits;
        if (lastbits >= 8)
        {
            lastbits   -= 8;
            cbuf[cnt++] = lastbyte >> lastbits;
        }
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits > 0)
    {
        cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

/*_________________________________________________________________________
 |
 | sizeofint - calculate bitsize of an integer
 |
 | return the number of bits needed to store an integer with given max size
 |
 */

static int sizeofint(const int size)
{
    int num         = 1;
    int num_of_bits = 0;

    while (size >= num && num_of_bits < 32)
    {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

/*___________________________________________________________________________
 |
 | sizeofints - calculate 'bitsize' of compressed ints
 |
 | given the number of small unsigned integers and the maximum value
 | return the number of bits needed to read or write them with the
 | routines receiveints and sendints. You need this parameter when
 | calling these routines. Note that for many calls I can use
 | the variable 'smallidx' which is exactly the number of bits, and
 | So I don't need to call 'sizeofints for those calls.
 */

static int sizeofints( const int num_of_ints, unsigned int sizes[])
{
    int          i, num;
    int          bytes[32];
    unsigned int num_of_bytes, num_of_bits, bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0]     = 1;
    num_of_bits  = 0;
    for (i = 0; i < num_of_ints; i++)
    {
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp            = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp          >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp            >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num)
    {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}

/*____________________________________________________________________________
 |
 | sendints - send a small set of small integers in compressed format
 |
 | this routine is used internally by xdr3dfcoord, to send a set of
 | small integers to the buffer.
 | Multiplication with fixed (specified maximum ) sizes is used to get
 | to one big, multibyte integer. Allthough the routine could be
 | modified to handle sizes bigger than 16777216, or more than just
 | a few integers, this is not done, because the gain in compression
 | isn't worth the effort. Note that overflowing the multiplication
 | or the byte buffer (32 bytes) is unchecked and causes bad results.
 |
 */

static void sendints(int buf[], const int num_of_ints, const int num_of_bits,
                     unsigned int sizes[], unsigned int nums[])
{

    int          i, num_of_bytes, bytecnt;
    unsigned int bytes[32], tmp;

    tmp          = nums[0];
    num_of_bytes = 0;
    do
    {
        bytes[num_of_bytes++] = tmp & 0xff;
        tmp                 >>= 8;
    }
    while (tmp != 0);

    for (i = 1; i < num_of_ints; i++)
    {
        if (nums[i] >= sizes[i])
        {
            fprintf(stderr, "major breakdown in sendints num %u doesn't "
                    "match size %u\n", nums[i], sizes[i]);
            exit(1);
        }
        /* use one step multiply */
        tmp = nums[i];
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp            = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp          >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp            >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8)
    {
        for (i = 0; i < num_of_bytes; i++)
        {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    }
    else
    {
        for (i = 0; i < num_of_bytes-1; i++)
        {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}


/*___________________________________________________________________________
 |
 | receivebits - decode number from buf using specified number of bits
 |
 | extract the number of bits from the array buf and construct an integer
 | from it. Return that value.
 |
 */

static int receivebits(int buf[], int num_of_bits)
{

    int             cnt, num, lastbits;
    unsigned int    lastbyte;
    unsigned char * cbuf;
    int             mask = (1 << num_of_bits) -1;

    cbuf     = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt      = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];

    num = 0;
    while (num_of_bits >= 8)
    {
        lastbyte     = ( lastbyte << 8 ) | cbuf[cnt++];
        num         |=  (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        if (lastbits < num_of_bits)
        {
            lastbits += 8;
            lastbyte  = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num      |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num   &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}

/*____________________________________________________________________________
 |
 | receiveints - decode 'small' integers from the buf array
 |
 | this routine is the inverse from sendints() and decodes the small integers
 | written to buf by calculating the remainder and doing divisions with
 | the given sizes[]. You need to specify the total number of bits to be
 | used from buf in num_of_bits.
 |
 */

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
                        unsigned int sizes[], int nums[])
{
    int bytes[32];
    int i, j, num_of_bytes, p, num;

    bytes[0]     = bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8)
    {
        bytes[num_of_bytes++] = receivebits(buf, 8);
        num_of_bits          -= 8;
    }
    if (num_of_bits > 0)
    {
        bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--)
    {
        num = 0;
        for (j = num_of_bytes-1; j >= 0; j--)
        {
            num      = (num << 8) | bytes[j];
            p        = num / sizes[i];
            bytes[j] = p;
            num      = num - p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

/*____________________________________________________________________________
 |
 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
 |
 | this routine reads or writes (depending on how you opened the file with
 | xdropen() ) a large number of 3d coordinates (stored in *fp).
 | The number of coordinates triplets to write is given by *size. On
 | read this number may be zero, in which case it reads as many as were written
 | or it may specify the number if triplets to read (which should match the
 | number written).
 | Compression is achieved by first converting all floating numbers to integer
 | using multiplication by *precision and rounding to the nearest integer.
 | Then the minimum and maximum value are calculated to determine the range.
 | The limited range of integers so found, is used to compress the coordinates.
 | In addition the differences between succesive coordinates is calculated.
 | If the difference happens to be 'small' then only the difference is saved,
 | compressing the data even more. The notion of 'small' is changed dynamically
 | and is enlarged or reduced whenever needed or possible.
 | Extra compression is achieved in the case of GROMOS and coordinates of
 | water molecules. GROMOS first writes out the Oxygen position, followed by
 | the two hydrogens. In order to make the differences smaller (and thereby
 | compression the data better) the order is changed into first one hydrogen
 | then the oxygen, followed by the other hydrogen. This is rather special, but
 | it shouldn't harm in the general case.
 |
 */

int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision)
{
    int     *ip  = NULL;
    int     *buf = NULL;
    gmx_bool bRead;

    /* preallocate a small buffer and ip on the stack - if we need more
       we can always malloc(). This is faster for small values of size: */
    unsigned     prealloc_size = 3*16;
    int          prealloc_ip[3*16], prealloc_buf[3*20];
    int          we_should_free = 0;

    int          minint[3], maxint[3], mindiff, *lip, diff;
    int          lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int          minidx, maxidx;
    unsigned     sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int          flag, k;
    int          smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float       *lfp, lf;
    int          tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int          bufsize, xdrid, lsize;
    unsigned int bitsize;
    float        inv_precision;
    int          errval = 1;
    int          rc;

    bRead         = (xdrs->x_op == XDR_DECODE);
    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;

    if (!bRead)
    {
        /* xdrs is open for writing */

        if (xdr_int(xdrs, size) == 0)
        {
            return 0;
        }
        size3 = *size * 3;
        /* when the number of coordinates is small, don't try to compress; just
         * write them as floats using xdr_vector
         */
        if (*size <= 9)
        {
            return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3,
                               (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
        }

        if (xdr_float(xdrs, precision) == 0)
        {
            return 0;
        }

        if (size3 <= prealloc_size)
        {
            ip  = prealloc_ip;
            buf = prealloc_buf;
        }
        else
        {
            we_should_free = 1;
            bufsize        = size3 * 1.2;
            ip             = (int *)malloc((size_t)(size3 * sizeof(*ip)));
            buf            = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
            if (ip == NULL || buf == NULL)
            {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }
        /* buf[0-2] are special and do not contain actual data */
        buf[0]    = buf[1] = buf[2] = 0;
        minint[0] = minint[1] = minint[2] = INT_MAX;
        maxint[0] = maxint[1] = maxint[2] = INT_MIN;
        prevrun   = -1;
        lfp       = fp;
        lip       = ip;
        mindiff   = INT_MAX;
        oldlint1  = oldlint2 = oldlint3 = 0;
        while (lfp < fp + size3)
        {
            /* find nearest integer */
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint1 = lf;
            if (lint1 < minint[0])
            {
                minint[0] = lint1;
            }
            if (lint1 > maxint[0])
            {
                maxint[0] = lint1;
            }
            *lip++ = lint1;
            lfp++;
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint2 = lf;
            if (lint2 < minint[1])
            {
                minint[1] = lint2;
            }
            if (lint2 > maxint[1])
            {
                maxint[1] = lint2;
            }
            *lip++ = lint2;
            lfp++;
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint3 = lf;
            if (lint3 < minint[2])
            {
                minint[2] = lint3;
            }
            if (lint3 > maxint[2])
            {
                maxint[2] = lint3;
            }
            *lip++ = lint3;
            lfp++;
            diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
            if (diff < mindiff && lfp > fp + 3)
            {
                mindiff = diff;
            }
            oldlint1 = lint1;
            oldlint2 = lint2;
            oldlint3 = lint3;
        }
        if ( (xdr_int(xdrs, &(minint[0])) == 0) ||
             (xdr_int(xdrs, &(minint[1])) == 0) ||
             (xdr_int(xdrs, &(minint[2])) == 0) ||
             (xdr_int(xdrs, &(maxint[0])) == 0) ||
             (xdr_int(xdrs, &(maxint[1])) == 0) ||
             (xdr_int(xdrs, &(maxint[2])) == 0))
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        if ((float)maxint[0] - (float)minint[0] >= MAXABS ||
            (float)maxint[1] - (float)minint[1] >= MAXABS ||
            (float)maxint[2] - (float)minint[2] >= MAXABS)
        {
            /* turning value in unsigned by subtracting minint
             * would cause overflow
             */
            errval = 0;
        }
        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = sizeofint(sizeint[0]);
            bitsizeint[1] = sizeofint(sizeint[1]);
            bitsizeint[2] = sizeofint(sizeint[2]);
            bitsize       = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = sizeofints(3, sizeint);
        }
        lip      = ip;
        luip     = (unsigned int *) ip;
        smallidx = FIRSTIDX;
        while (smallidx < LASTIDX && magicints[smallidx] < mindiff)
        {
            smallidx++;
        }
        if (xdr_int(xdrs, &smallidx) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        maxidx       = MIN(LASTIDX, smallidx + 8);
        minidx       = maxidx - 8; /* often this equal smallidx */
        smaller      = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
        smallnum     = magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        larger       = magicints[maxidx] / 2;
        i            = 0;
        while (i < *size)
        {
            is_small  = 0;
            thiscoord = (int *)(luip) + i * 3;
            if (smallidx < maxidx && i >= 1 &&
                abs(thiscoord[0] - prevcoord[0]) < larger &&
                abs(thiscoord[1] - prevcoord[1]) < larger &&
                abs(thiscoord[2] - prevcoord[2]) < larger)
            {
                is_smaller = 1;
            }
            else if (smallidx > minidx)
            {
                is_smaller = -1;
            }
            else
            {
                is_smaller = 0;
            }
            if (i + 1 < *size)
            {
                if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
                    abs(thiscoord[1] - thiscoord[4]) < smallnum &&
                    abs(thiscoord[2] - thiscoord[5]) < smallnum)
                {
                    /* interchange first with second atom for better
                     * compression of water molecules
                     */
                    tmp          = thiscoord[0]; thiscoord[0] = thiscoord[3];
                    thiscoord[3] = tmp;
                    tmp          = thiscoord[1]; thiscoord[1] = thiscoord[4];
                    thiscoord[4] = tmp;
                    tmp          = thiscoord[2]; thiscoord[2] = thiscoord[5];
                    thiscoord[5] = tmp;
                    is_small     = 1;
                }

            }
            tmpcoord[0] = thiscoord[0] - minint[0];
            tmpcoord[1] = thiscoord[1] - minint[1];
            tmpcoord[2] = thiscoord[2] - minint[2];
            if (bitsize == 0)
            {
                sendbits(buf, bitsizeint[0], tmpcoord[0]);
                sendbits(buf, bitsizeint[1], tmpcoord[1]);
                sendbits(buf, bitsizeint[2], tmpcoord[2]);
            }
            else
            {
                sendints(buf, 3, bitsize, sizeint, tmpcoord);
            }
            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];
            thiscoord    = thiscoord + 3;
            i++;

            run = 0;
            if (is_small == 0 && is_smaller == -1)
            {
                is_smaller = 0;
            }
            while (is_small && run < 8*3)
            {
                if (is_smaller == -1 && (
                        SQR(thiscoord[0] - prevcoord[0]) +
                        SQR(thiscoord[1] - prevcoord[1]) +
                        SQR(thiscoord[2] - prevcoord[2]) >= smaller * smaller))
                {
                    is_smaller = 0;
                }

                tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
                tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
                tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;

                prevcoord[0] = thiscoord[0];
                prevcoord[1] = thiscoord[1];
                prevcoord[2] = thiscoord[2];

                i++;
                thiscoord = thiscoord + 3;
                is_small  = 0;
                if (i < *size &&
                    abs(thiscoord[0] - prevcoord[0]) < smallnum &&
                    abs(thiscoord[1] - prevcoord[1]) < smallnum &&
                    abs(thiscoord[2] - prevcoord[2]) < smallnum)
                {
                    is_small = 1;
                }
            }
            if (run != prevrun || is_smaller != 0)
            {
                prevrun = run;
                sendbits(buf, 1, 1); /* flag the change in run-length */
                sendbits(buf, 5, run+is_smaller+1);
            }
            else
            {
                sendbits(buf, 1, 0); /* flag the fact that runlength did not change */
            }
            for (k = 0; k < run; k += 3)
            {
                sendints(buf, 3, smallidx, sizesmall, &tmpcoord[k]);
            }
            if (is_smaller != 0)
            {
                smallidx += is_smaller;
                if (is_smaller < 0)
                {
                    smallnum = smaller;
                    smaller  = magicints[smallidx-1] / 2;
                }
                else
                {
                    smaller  = smallnum;
                    smallnum = magicints[smallidx] / 2;
                }
                sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
            }
        }
        if (buf[1] != 0)
        {
            buf[0]++;
        }
        /* buf[0] holds the length in bytes */
        if (xdr_int(xdrs, &(buf[0])) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }


        rc = errval * (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]));
        if (we_should_free)
        {
            free(ip);
            free(buf);
        }
        return rc;

    }
    else
    {

        /* xdrs is open for reading */

        if (xdr_int(xdrs, &lsize) == 0)
        {
            return 0;
        }
        if (*size != 0 && lsize != *size)
        {
            fprintf(stderr, "wrong number of coordinates in xdr3dfcoord; "
                    "%d arg vs %d in file", *size, lsize);
        }
        *size = lsize;
        size3 = *size * 3;
        if (*size <= 9)
        {
            *precision = -1;
            return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3,
                               (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
        }
        if (xdr_float(xdrs, precision) == 0)
        {
            return 0;
        }

        if (size3 <= prealloc_size)
        {
            ip  = prealloc_ip;
            buf = prealloc_buf;
        }
        else
        {
            we_should_free = 1;
            bufsize        = size3 * 1.2;
            ip             = (int *)malloc((size_t)(size3 * sizeof(*ip)));
            buf            = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
            if (ip == NULL || buf == NULL)
            {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }

        buf[0] = buf[1] = buf[2] = 0;

        if ( (xdr_int(xdrs, &(minint[0])) == 0) ||
             (xdr_int(xdrs, &(minint[1])) == 0) ||
             (xdr_int(xdrs, &(minint[2])) == 0) ||
             (xdr_int(xdrs, &(maxint[0])) == 0) ||
             (xdr_int(xdrs, &(maxint[1])) == 0) ||
             (xdr_int(xdrs, &(maxint[2])) == 0))
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = sizeofint(sizeint[0]);
            bitsizeint[1] = sizeofint(sizeint[1]);
            bitsizeint[2] = sizeofint(sizeint[2]);
            bitsize       = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = sizeofints(3, sizeint);
        }

        if (xdr_int(xdrs, &smallidx) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        maxidx       = MIN(LASTIDX, smallidx + 8);
        minidx       = maxidx - 8; /* often this equal smallidx */
        smaller      = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
        smallnum     = magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        larger       = magicints[maxidx];

        /* buf[0] holds the length in bytes */

        if (xdr_int(xdrs, &(buf[0])) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }


        if (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }



        buf[0] = buf[1] = buf[2] = 0;

        lfp           = fp;
        inv_precision = 1.0 / *precision;
        run           = 0;
        i             = 0;
        lip           = ip;
        while (i < lsize)
        {
            thiscoord = (int *)(lip) + i * 3;

            if (bitsize == 0)
            {
                thiscoord[0] = receivebits(buf, bitsizeint[0]);
                thiscoord[1] = receivebits(buf, bitsizeint[1]);
                thiscoord[2] = receivebits(buf, bitsizeint[2]);
            }
            else
            {
                receiveints(buf, 3, bitsize, sizeint, thiscoord);
            }

            i++;
            thiscoord[0] += minint[0];
            thiscoord[1] += minint[1];
            thiscoord[2] += minint[2];

            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];


            flag       = receivebits(buf, 1);
            is_smaller = 0;
            if (flag == 1)
            {
                run        = receivebits(buf, 5);
                is_smaller = run % 3;
                run       -= is_smaller;
                is_smaller--;
            }
            if (run > 0)
            {
                thiscoord += 3;
                for (k = 0; k < run; k += 3)
                {
                    receiveints(buf, 3, smallidx, sizesmall, thiscoord);
                    i++;
                    thiscoord[0] += prevcoord[0] - smallnum;
                    thiscoord[1] += prevcoord[1] - smallnum;
                    thiscoord[2] += prevcoord[2] - smallnum;
                    if (k == 0)
                    {
                        /* interchange first with second atom for better
                         * compression of water molecules
                         */
                        tmp          = thiscoord[0]; thiscoord[0] = prevcoord[0];
                        prevcoord[0] = tmp;
                        tmp          = thiscoord[1]; thiscoord[1] = prevcoord[1];
                        prevcoord[1] = tmp;
                        tmp          = thiscoord[2]; thiscoord[2] = prevcoord[2];
                        prevcoord[2] = tmp;
                        *lfp++       = prevcoord[0] * inv_precision;
                        *lfp++       = prevcoord[1] * inv_precision;
                        *lfp++       = prevcoord[2] * inv_precision;
                    }
                    else
                    {
                        prevcoord[0] = thiscoord[0];
                        prevcoord[1] = thiscoord[1];
                        prevcoord[2] = thiscoord[2];
                    }
                    *lfp++ = thiscoord[0] * inv_precision;
                    *lfp++ = thiscoord[1] * inv_precision;
                    *lfp++ = thiscoord[2] * inv_precision;
                }
            }
            else
            {
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
            smallidx += is_smaller;
            if (is_smaller < 0)
            {
                smallnum = smaller;
                if (smallidx > FIRSTIDX)
                {
                    smaller = magicints[smallidx - 1] /2;
                }
                else
                {
                    smaller = 0;
                }
            }
            else if (is_smaller > 0)
            {
                smaller  = smallnum;
                smallnum = magicints[smallidx] / 2;
            }
            sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        }
    }
    if (we_should_free)
    {
        free(ip);
        free(buf);
    }
    return 1;
}



/******************************************************************

   XTC files have a relatively simple structure.
   They have a header of 16 bytes and the rest are
   the compressed coordinates of the files. Due to the
   compression 00 is not present in the coordinates.
   The first 4 bytes of the header are the magic number
   1995 (0x000007CB). If we find this number we are guaranteed
   to be in the header, due to the presence of so many zeros.
   The second 4 bytes are the number of atoms in the frame, and is
   assumed to be constant. The third 4 bytes are the frame number.
   The last 4 bytes are a floating point representation of the time.

 ********************************************************************/

/* Must match definition in xtcio.c */
#ifndef XTC_MAGIC
#define XTC_MAGIC 1995
#endif

static const int header_size = 16;

/* Check if we are at the header start.
   At the same time it will also read 1 int */
static int xtc_at_header_start(FILE *fp, XDR *xdrs,
                               int natoms, int * timestep, float * time)
{
    int       i_inp[3];
    float     f_inp[10];
    int       i;
    gmx_off_t off;


    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }
    /* read magic natoms and timestep */
    for (i = 0; i < 3; i++)
    {
        if (!xdr_int(xdrs, &(i_inp[i])))
        {
            gmx_fseek(fp, off+XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* quick return */
    if (i_inp[0] != XTC_MAGIC)
    {
        if (gmx_fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        return 0;
    }
    /* read time and box */
    for (i = 0; i < 10; i++)
    {
        if (!xdr_float(xdrs, &(f_inp[i])))
        {
            gmx_fseek(fp, off+XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* Make a rigourous check to see if we are in the beggining of a header
       Hopefully there are no ambiguous cases */
    /* This check makes use of the fact that the box matrix has 3 zeroes on the upper
       right triangle and that the first element must be nonzero unless the entire matrix is zero
     */
    if (i_inp[1] == natoms &&
        ((f_inp[1] != 0 && f_inp[6] == 0) ||
         (f_inp[1] == 0 && f_inp[5] == 0 && f_inp[9] == 0)))
    {
        if (gmx_fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        *time     = f_inp[0];
        *timestep = i_inp[2];
        return 1;
    }
    if (gmx_fseek(fp, off+XDR_INT_SIZE, SEEK_SET))
    {
        return -1;
    }
    return 0;
}

static int
xtc_get_next_frame_number(FILE *fp, XDR *xdrs, int natoms)
{
    gmx_off_t off;
    int       step;
    float     time;
    int       ret;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
        }
    }
}


static float xtc_get_next_frame_time(FILE *fp, XDR *xdrs, int natoms,
                                     gmx_bool * bOK)
{
    gmx_off_t off;
    float     time;
    int       step;
    int       ret;
    *bOK = 0;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = 1;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = 0;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
    }
}


static float
xtc_get_current_frame_time(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK)
{
    gmx_off_t off;
    int       step;
    float     time;
    int       ret;
    *bOK = 0;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = 1;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = 0;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (gmx_fseek(fp, -2*XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}

/* Currently not used, just for completeness */
static int
xtc_get_current_frame_number(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK)
{
    gmx_off_t off;
    int       ret;
    int       step;
    float     time;
    *bOK = 0;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }


    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = 1;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = 0;
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;

        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (gmx_fseek(fp, -2*XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}



static gmx_off_t xtc_get_next_frame_start(FILE *fp, XDR *xdrs, int natoms)
{
    int       inp;
    gmx_off_t res;
    int       ret;
    int       step;
    float     time;
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (1)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            if ((res = gmx_ftell(fp)) >= 0)
            {
                return res - XDR_INT_SIZE;
            }
            else
            {
                return res;
            }
        }
        else if (ret == -1)
        {
            return -1;
        }
    }
}


static
float
xdr_xtc_estimate_dt(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK)
{
    float     res;
    float     tinit;
    gmx_off_t off;

    *bOK = 0;
    if ((off   = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    tinit = xtc_get_current_frame_time(fp, xdrs, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res = xtc_get_next_frame_time(fp, xdrs, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res -= tinit;
    if (0 != gmx_fseek(fp, off, SEEK_SET))
    {
        *bOK = 0;
        return -1;
    }
    return res;
}

int
xdr_xtc_seek_frame(int frame, FILE *fp, XDR *xdrs, int natoms)
{
    gmx_off_t low = 0;
    gmx_off_t high, pos;


    /* round to 4 bytes */
    int        fr;
    gmx_off_t  offset;
    if (gmx_fseek(fp, 0, SEEK_END))
    {
        return -1;
    }

    if ((high = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    /* round to 4 bytes  */
    high  /= XDR_INT_SIZE;
    high  *= XDR_INT_SIZE;
    offset = ((high/2)/XDR_INT_SIZE)*XDR_INT_SIZE;

    if (gmx_fseek(fp, offset, SEEK_SET))
    {
        return -1;
    }

    while (1)
    {
        fr = xtc_get_next_frame_number(fp, xdrs, natoms);
        if (fr < 0)
        {
            return -1;
        }
        if (fr != frame && llabs(low-high) > header_size)
        {
            if (fr < frame)
            {
                low = offset;
            }
            else
            {
                high = offset;
            }
            /* round to 4 bytes */
            offset = (((high+low)/2)/4)*4;

            if (gmx_fseek(fp, offset, SEEK_SET))
            {
                return -1;
            }
        }
        else
        {
            break;
        }
    }
    if (offset <= header_size)
    {
        offset = low;
    }

    if (gmx_fseek(fp, offset, SEEK_SET))
    {
        return -1;
    }

    if ((pos = xtc_get_next_frame_start(fp, xdrs, natoms)) < 0)
    {
        /* we probably hit an end of file */
        return -1;
    }

    if (gmx_fseek(fp, pos, SEEK_SET))
    {
        return -1;
    }

    return 0;
}



int xdr_xtc_seek_time(real time, FILE *fp, XDR *xdrs, int natoms, gmx_bool bSeekForwardOnly)
{
    float     t;
    float     dt;
    gmx_bool  bOK = FALSE;
    gmx_off_t low = 0;
    gmx_off_t high, offset, pos;
    int       res;
    int       dt_sign = 0;

    if (bSeekForwardOnly)
    {
        low = gmx_ftell(fp)-header_size;
    }
    if (gmx_fseek(fp, 0, SEEK_END))
    {
        return -1;
    }

    if ((high = gmx_ftell(fp)) < 0)
    {
        return -1;
    }
    /* round to int  */
    high  /= XDR_INT_SIZE;
    high  *= XDR_INT_SIZE;
    offset = (((high-low) / 2) / XDR_INT_SIZE) * XDR_INT_SIZE;

    if (gmx_fseek(fp, offset, SEEK_SET))
    {
        return -1;
    }


    /*
     * No need to call xdr_xtc_estimate_dt here - since xdr_xtc_estimate_dt is called first thing in the loop
       dt = xdr_xtc_estimate_dt(fp, xdrs, natoms, &bOK);

       if (!bOK)
       {
        return -1;
       }
     */

    while (1)
    {
        dt = xdr_xtc_estimate_dt(fp, xdrs, natoms, &bOK);
        if (!bOK)
        {
            return -1;
        }
        else
        {
            if (dt > 0)
            {
                if (dt_sign == -1)
                {
                    /* Found a place in the trajectory that has positive time step while
                       other has negative time step */
                    return -2;
                }
                dt_sign = 1;
            }
            else if (dt < 0)
            {
                if (dt_sign == 1)
                {
                    /* Found a place in the trajectory that has positive time step while
                       other has negative time step */
                    return -2;
                }
                dt_sign = -1;
            }
        }
        t = xtc_get_next_frame_time(fp, xdrs, natoms, &bOK);
        if (!bOK)
        {
            return -1;
        }

        /* If we are before the target time and the time step is positive or 0, or we have
           after the target time and the time step is negative, or the difference between
           the current time and the target time is bigger than dt and above all the distance between high
           and low is bigger than 1 frame, then do another step of binary search. Otherwise stop and check
           if we reached the solution */
        if ((((t < time && dt_sign >= 0) || (t > time && dt_sign == -1)) ||
             ((t - time) >= dt && dt_sign >= 0) || ((time - t) >= -dt && dt_sign < 0)) &&
            (llabs(low - high) > header_size))
        {
            if (dt >= 0 && dt_sign != -1)
            {
                if (t < time)
                {
                    low = offset;
                }
                else
                {
                    high = offset;
                }
            }
            else if (dt <= 0 && dt_sign == -1)
            {
                if (t >= time)
                {
                    low = offset;
                }
                else
                {
                    high = offset;
                }
            }
            else
            {
                /* We should never reach here */
                return -1;
            }
            /* round to 4 bytes and subtract header*/
            offset = (((high + low) / 2) / XDR_INT_SIZE) * XDR_INT_SIZE;
            if (gmx_fseek(fp, offset, SEEK_SET))
            {
                return -1;
            }
        }
        else
        {
            if (llabs(low - high) <= header_size)
            {
                break;
            }
            /* re-estimate dt */
            if (xdr_xtc_estimate_dt(fp, xdrs, natoms, &bOK) != dt)
            {
                if (bOK)
                {
                    dt = xdr_xtc_estimate_dt(fp, xdrs, natoms, &bOK);
                }
            }
            if (t >= time && t - time < dt)
            {
                break;
            }
        }
    }

    if (offset <= header_size)
    {
        offset = low;
    }

    gmx_fseek(fp, offset, SEEK_SET);

    if ((pos = xtc_get_next_frame_start(fp, xdrs, natoms)) < 0)
    {
        return -1;
    }

    if (gmx_fseek(fp, pos, SEEK_SET))
    {
        return -1;
    }
    return 0;
}

float
xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK)
{
    float      time;
    gmx_off_t  off;
    int        res;
    *bOK = 1;
    off  = gmx_ftell(fp);
    if (off < 0)
    {
        *bOK = 0;
        return -1;
    }

    if ( (res = gmx_fseek(fp, -3*XDR_INT_SIZE, SEEK_END)) != 0)
    {
        *bOK = 0;
        return -1;
    }

    time = xtc_get_current_frame_time(fp, xdrs, natoms, bOK);
    if (!(*bOK))
    {
        return -1;
    }

    if ( (res = gmx_fseek(fp, off, SEEK_SET)) != 0)
    {
        *bOK = 0;
        return -1;
    }
    return time;
}


int
xdr_xtc_get_last_frame_number(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK)
{
    int        frame;
    gmx_off_t  off;
    int        res;
    *bOK = 1;

    if ((off = gmx_ftell(fp)) < 0)
    {
        *bOK = 0;
        return -1;
    }


    if (gmx_fseek(fp, -3*XDR_INT_SIZE, SEEK_END))
    {
        *bOK = 0;
        return -1;
    }

    frame = xtc_get_current_frame_number(fp, xdrs, natoms, bOK);
    if (!bOK)
    {
        return -1;
    }


    if (gmx_fseek(fp, off, SEEK_SET))
    {
        *bOK = 0;
        return -1;
    }

    return frame;
}
