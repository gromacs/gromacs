/*
 * xdr_stdio.c, XDR implementation on standard i/o file.
 *
 * Copyright (c) 2010, Oracle America, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the "Oracle America, Inc." nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This set of routines implements a XDR on a stdio stream.
 * XDR_ENCODE serializes onto the stream, XDR_DECODE de-serializes
 * from the stream.
 */

/* This file has been modified in the GROMACS distribution by, e.g.:
 * - removing macros such as:
 * "#ifdef EXPORT_RPC_SYMBOLS
 *  libc_hidden_def (xdr_free)
 *  #else
 *  libc_hidden_nolink_sunrpc (xdr_free, GLIBC_2_0)
 *  #endif"
 * - removing macros for fflush, fread, ftell and fwrite.
 * - removing an explicit cast to (struct xdr_ops *).
 * - changing u_int to unsigned int, u_long to unsigned long,
 *   u_short to unsigned short and u_char to unsigned char.
 * - changing caddr_t to char*.
 * - Changing
 * "static const struct xdr_ops xdrstdio_ops = ..."
 * to
 * "static struct XDR::xdr_ops xdrstdio_ops = ..."
 * - Adding xdr_swapbytes, xdr_htonl and xdr_ntohl functions
 *   from GROMACS.
 * - removing headers that are no longer necessary.
 */

#include "types.h"
#include <stdio.h>
#include "xdr.h"

static bool_t xdrstdio_getlong (XDR *, long *);
static bool_t xdrstdio_putlong (XDR *, const long *);
static bool_t xdrstdio_getbytes (XDR *, char *, unsigned int);
static bool_t xdrstdio_putbytes (XDR *, const char *, unsigned int);
static unsigned int xdrstdio_getpos (const XDR *);
static bool_t xdrstdio_setpos (XDR *, unsigned int);
static int32_t *xdrstdio_inline (XDR *, unsigned int);
static void xdrstdio_destroy (XDR *);
static bool_t xdrstdio_getint32 (XDR *, int32_t *);
static bool_t xdrstdio_putint32 (XDR *, const int32_t *);

/*
 * Ops vector for stdio type XDR
 */
static struct XDR::xdr_ops xdrstdio_ops =
{
  xdrstdio_getlong,		/* deserialize a long int */
  xdrstdio_putlong,		/* serialize a long int */
  xdrstdio_getbytes,		/* deserialize counted bytes */
  xdrstdio_putbytes,		/* serialize counted bytes */
  xdrstdio_getpos,		/* get offset in the stream */
  xdrstdio_setpos,		/* set offset in the stream */
  xdrstdio_inline,		/* prime stream for inline macros */
  xdrstdio_destroy,		/* destroy stream */
  xdrstdio_getint32,		/* deserialize a int */
  xdrstdio_putint32		/* serialize a int */
};

/* Copyright The GROMACS Authors */
static uint32_t xdr_swapbytes(uint32_t x)
{
    uint32_t y;
    int          i;
    char*        px = reinterpret_cast<char*>(&x);
    char*        py = reinterpret_cast<char*>(&y);

    for (i = 0; i < 4; i++)
    {
        py[i] = px[3 - i];
    }

    return y;
}

/* Copyright The GROMACS Authors */
static uint32_t xdr_htonl(uint32_t x)
{
    short s = 0x0F00;
    if (*(reinterpret_cast<char*>(&s)) == static_cast<char>(0x0F))
    {
        /* bigendian, do nothing */
        return x;
    }
    else
    {
        /* smallendian,swap bytes */
        return xdr_swapbytes(x);
    }
}

/* Copyright The GROMACS Authors */
static uint32_t xdr_ntohl(uint32_t x)
{
    short s = 0x0F00;
    if (*(reinterpret_cast<char*>(&s)) == static_cast<char>(0x0F))
    {
        /* bigendian, do nothing */
        return x;
    }
    else
    {
        /* smallendian, swap bytes */
        return xdr_swapbytes(x);
    }
}

/*
 * Initialize a stdio xdr stream.
 * Sets the xdr stream handle xdrs for use on the stream file.
 * Operation flag is set to op.
 */
void
xdrstdio_create (XDR *xdrs, FILE *file, enum xdr_op op)
{
  xdrs->x_op = op;
  /* We have to add the const since the `struct xdr_ops' in `struct XDR'
     is not `const'.  */
  xdrs->x_ops = &xdrstdio_ops;
  xdrs->x_private = (char *) file;
  xdrs->x_handy = 0;
  xdrs->x_base = 0;
}

/*
 * Destroy a stdio xdr stream.
 * Cleans up the xdr stream handle xdrs previously set up by xdrstdio_create.
 */
static void
xdrstdio_destroy (XDR *xdrs)
{
  (void) fflush ((FILE *) xdrs->x_private);
  /* xx should we close the file ?? */
};

static bool_t
xdrstdio_getlong (XDR *xdrs, long *lp)
{
  uint32_t mycopy;

  if (fread ((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  *lp = (long) xdr_ntohl (mycopy);
  return TRUE;
}

static bool_t
xdrstdio_putlong (XDR *xdrs, const long *lp)
{
  int32_t mycopy = xdr_htonl ((uint32_t) *lp);

  if (fwrite ((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  return TRUE;
}

static bool_t
xdrstdio_getbytes (XDR *xdrs, char * addr, unsigned int len)
{
  if ((len != 0) && (fread (addr, (int) len, 1,
			    (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static bool_t
xdrstdio_putbytes (XDR *xdrs, const char *addr, unsigned int len)
{
  if ((len != 0) && (fwrite (addr, (int) len, 1,
			     (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static unsigned int
xdrstdio_getpos (const XDR *xdrs)
{
  return (unsigned int) ftell ((FILE *) xdrs->x_private);
}

static bool_t
xdrstdio_setpos (XDR *xdrs, unsigned int pos)
{
  return fseek ((FILE *) xdrs->x_private, (long) pos, 0) < 0 ? FALSE : TRUE;
}

static int32_t *
xdrstdio_inline (XDR *xdrs, unsigned int len)
{
  /*
   * Must do some work to implement this: must insure
   * enough data in the underlying stdio buffer,
   * that the buffer is aligned so that we can indirect through a
   * long *, and stuff this pointer in xdrs->x_buf.  Doing
   * a fread or fwrite to a scratch buffer would defeat
   * most of the gains to be had here and require storage
   * management on this buffer, so we don't do this.
   */
  return NULL;
}

static bool_t
xdrstdio_getint32 (XDR *xdrs, int32_t *ip)
{
  int32_t mycopy;

  if (fread ((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  *ip = xdr_ntohl (mycopy);
  return TRUE;
}

static bool_t
xdrstdio_putint32 (XDR *xdrs, const int32_t *ip)
{
  int32_t mycopy = xdr_htonl (*ip);

  ip = &mycopy;
  if (fwrite ((char *) ip, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  return TRUE;
}
