/*
 * xdr.c, Generic XDR routines implementation.
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
 * These are the "generic" xdr routines used to serialize and de-serialize
 * most common data items.  See xdr.h for more info on the interface to
 * xdr.
 */

/* This file has been modified in the GROMACS distribution by, e.g.:
 * - removing macros such as:
 * "#ifdef EXPORT_RPC_SYMBOLS
 *  libc_hidden_def (xdr_free)
 *  #else
 *  libc_hidden_nolink_sunrpc (xdr_free, GLIBC_2_0)
 *  #endif"
 * - using fputs() instead of __fxprintf()
 * - changing u_int to unsigned int, u_long to unsigned long,
 *   u_short to unsigned short and u_char to unsigned char.
 * - changing caddr_t to char*.
 * - removing functions xdr_hyper(), xdr_uhyper(),
 *   xdr_longlong_t() and xdr_u_longlong_t().
 * - removing headers that are no longer necessary.
 */

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>

#include "types.h"
#include "xdr.h"


/*
 * constants specific to the xdr "protocol"
 */
#define XDR_FALSE	((long) 0)
#define XDR_TRUE	((long) 1)
#define LASTUNSIGNED	((unsigned int) 0-1)

/*
 * for unit alignment
 */
static const char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};

/*
 * Free a data structure using XDR
 * Not a filter, but a convenient utility nonetheless
 */
void
xdr_free (xdrproc_t proc, char *objp)
{
  XDR x;

  x.x_op = XDR_FREE;
  (*proc) (&x, objp);
}

/*
 * XDR nothing
 */
bool_t
xdr_void (void)
{
  return TRUE;
}

/*
 * XDR integers
 */
bool_t
xdr_int (XDR *xdrs, int *ip)
{

#if INT_MAX < LONG_MAX
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (long) *ip;
      return XDR_PUTLONG (xdrs, &l);

    case XDR_DECODE:
      if (!XDR_GETLONG (xdrs, &l))
	{
	  return FALSE;
	}
      *ip = (int) l;
      /* Fall through.  */
    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#elif INT_MAX == LONG_MAX
  return xdr_long (xdrs, (long *) ip);
#elif INT_MAX == SHRT_MAX
  return xdr_short (xdrs, (short *) ip);
#else
#error unexpected integer sizes in_xdr_int()
#endif
}

/*
 * XDR unsigned integers
 */
bool_t
xdr_u_int (XDR *xdrs, unsigned int *up)
{
#if UINT_MAX < ULONG_MAX
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (unsigned long) * up;
      return XDR_PUTLONG (xdrs, &l);

    case XDR_DECODE:
      if (!XDR_GETLONG (xdrs, &l))
	{
	  return FALSE;
	}
      *up = (unsigned int) (unsigned long) l;
      /* Fall through.  */
    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#elif UINT_MAX == ULONG_MAX
  return xdr_u_long (xdrs, (unsigned long *) up);
#elif UINT_MAX == USHRT_MAX
  return xdr_short (xdrs, (short *) up);
#else
#error unexpected integer sizes in_xdr_u_int()
#endif
}

/*
 * XDR long integers
 * The definition of xdr_long() is kept for backward
 * compatibility. Instead xdr_int() should be used.
 */
bool_t
xdr_long (XDR *xdrs, long *lp)
{

  if (xdrs->x_op == XDR_ENCODE
      && (sizeof (int32_t) == sizeof (long)
	  || (int32_t) *lp == *lp))
    return XDR_PUTLONG (xdrs, lp);

  if (xdrs->x_op == XDR_DECODE)
    return XDR_GETLONG (xdrs, lp);

  if (xdrs->x_op == XDR_FREE)
    return TRUE;

  return FALSE;
}

/*
 * XDR unsigned long integers
 * The definition of xdr_u_long() is kept for backward
 * compatibility. Instead xdr_u_int() should be used.
 */
bool_t
xdr_u_long (XDR *xdrs, unsigned long *ulp)
{
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      {
	long int tmp;

	if (XDR_GETLONG (xdrs, &tmp) == FALSE)
	  return FALSE;

	*ulp = (uint32_t) tmp;
	return TRUE;
      }

    case XDR_ENCODE:
      if (sizeof (uint32_t) != sizeof (unsigned long)
	  && (uint32_t) *ulp != *ulp)
	return FALSE;

      return XDR_PUTLONG (xdrs, (long *) ulp);

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}

/*
 * XDR short integers
 */
bool_t
xdr_short (XDR *xdrs, short *sp)
{
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (long) *sp;
      return XDR_PUTLONG (xdrs, &l);

    case XDR_DECODE:
      if (!XDR_GETLONG (xdrs, &l))
	{
	  return FALSE;
	}
      *sp = (short) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}

/*
 * XDR unsigned short integers
 */
bool_t
xdr_u_short (XDR *xdrs, unsigned short *usp)
{
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (unsigned long) * usp;
      return XDR_PUTLONG (xdrs, &l);

    case XDR_DECODE:
      if (!XDR_GETLONG (xdrs, &l))
	{
	  return FALSE;
	}
      *usp = (unsigned short) (unsigned long) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR a char
 */
bool_t
xdr_char (XDR *xdrs, char *cp)
{
  int i;

  i = (*cp);
  if (!xdr_int (xdrs, &i))
    {
      return FALSE;
    }
  *cp = i;
  return TRUE;
}

/*
 * XDR an unsigned char
 */
bool_t
xdr_u_char (XDR *xdrs, unsigned char *cp)
{
  unsigned int u;

  u = (*cp);
  if (!xdr_u_int (xdrs, &u))
    {
      return FALSE;
    }
  *cp = u;
  return TRUE;
}

/*
 * XDR booleans
 */
bool_t
xdr_bool (XDR *xdrs, bool_t *bp)
{
  long lb;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      lb = *bp ? XDR_TRUE : XDR_FALSE;
      return XDR_PUTLONG (xdrs, &lb);

    case XDR_DECODE:
      if (!XDR_GETLONG (xdrs, &lb))
	{
	  return FALSE;
	}
      *bp = (lb == XDR_FALSE) ? FALSE : TRUE;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}

/*
 * XDR enumerations
 */
bool_t
xdr_enum (XDR *xdrs, enum_t *ep)
{
  enum sizecheck
    {
      SIZEVAL
    };				/* used to find the size of an enum */

  /*
   * enums are treated as ints
   */
  if (sizeof (enum sizecheck) == 4)
    {
#if INT_MAX < LONG_MAX
      long l;

      switch (xdrs->x_op)
	{
	case XDR_ENCODE:
	  l = *ep;
	  return XDR_PUTLONG (xdrs, &l);

	case XDR_DECODE:
	  if (!XDR_GETLONG (xdrs, &l))
	    {
	      return FALSE;
	    }
	  *ep = l;
	  /* Fall through.  */
	case XDR_FREE:
	  return TRUE;

	}
      return FALSE;
#else
      return xdr_long (xdrs, (long *) ep);
#endif
    }
  else if (sizeof (enum sizecheck) == sizeof (short))
    {
      return xdr_short (xdrs, (short *) ep);
    }
  else
    {
      return FALSE;
    }
}

/*
 * XDR opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
bool_t
xdr_opaque (XDR *xdrs, char* cp, unsigned int cnt)
{
  unsigned int rndup;
  static char crud[BYTES_PER_XDR_UNIT];

  /*
   * if no data we are done
   */
  if (cnt == 0)
    return TRUE;

  /*
   * round byte count to full xdr units
   */
  rndup = cnt % BYTES_PER_XDR_UNIT;
  if (rndup > 0)
    rndup = BYTES_PER_XDR_UNIT - rndup;

  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (!XDR_GETBYTES (xdrs, cp, cnt))
	{
	  return FALSE;
	}
      if (rndup == 0)
	return TRUE;
      return XDR_GETBYTES (xdrs, (char *)crud, rndup);

    case XDR_ENCODE:
      if (!XDR_PUTBYTES (xdrs, cp, cnt))
	{
	  return FALSE;
	}
      if (rndup == 0)
	return TRUE;
      return XDR_PUTBYTES (xdrs, xdr_zero, rndup);

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}

/*
 * XDR counted bytes
 * *cpp is a pointer to the bytes, *sizep is the count.
 * If *cpp is NULL maxsize bytes are allocated
 */
bool_t
xdr_bytes (XDR *xdrs, char **cpp, unsigned int *sizep, unsigned int maxsize)
{
  char *sp = *cpp;	/* sp is the actual string pointer */
  unsigned int nodesize;

  /*
   * first deal with the length since xdr bytes are counted
   */
  if (!xdr_u_int (xdrs, sizep))
    {
      return FALSE;
    }
  nodesize = *sizep;
  if ((nodesize > maxsize) && (xdrs->x_op != XDR_FREE))
    {
      return FALSE;
    }

  /*
   * now deal with the actual bytes
   */
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (nodesize == 0)
	{
	  return TRUE;
	}
      if (sp == NULL)
	{
	  *cpp = sp = (char *) mem_alloc (nodesize);
	}
      if (sp == NULL)
	{
          (void)fputs("xdr_bytes: out of memory\n", stderr);
	  return FALSE;
	}
      /* Fall through.  */

    case XDR_ENCODE:
      return xdr_opaque (xdrs, sp, nodesize);

    case XDR_FREE:
      if (sp != NULL)
	{
	  mem_free (sp, nodesize);
	  *cpp = NULL;
	}
      return TRUE;
    }
  return FALSE;
}

/*
 * Implemented here due to commonality of the object.
 */
bool_t
xdr_netobj (XDR *xdrs, struct netobj *np)
{

  return xdr_bytes (xdrs, &np->n_bytes, &np->n_len, MAX_NETOBJ_SZ);
}

/*
 * XDR a discriminated union
 * Support routine for discriminated unions.
 * You create an array of xdrdiscrim structures, terminated with
 * an entry with a null procedure pointer.  The routine gets
 * the discriminant value and then searches the array of xdrdiscrims
 * looking for that value.  It calls the procedure given in the xdrdiscrim
 * to handle the discriminant.  If there is no specific routine a default
 * routine may be called.
 * If there is no specific or default routine an error is returned.
 */
bool_t
xdr_union (XDR *xdrs,
	   /* enum to decide which arm to work on */
	   enum_t *dscmp,
	   /* the union itself */
	   char *unp,
	   /* [value, xdr proc] for each arm */
	   const struct xdr_discrim *choices,
	   /* default xdr routine */
	   xdrproc_t dfault)
{
  enum_t dscm;

  /*
   * we deal with the discriminator;  it's an enum
   */
  if (!xdr_enum (xdrs, dscmp))
    {
      return FALSE;
    }
  dscm = *dscmp;

  /*
   * search choices for a value that matches the discriminator.
   * if we find one, execute the xdr routine for that value.
   */
  for (; choices->proc != NULL_xdrproc_t; choices++)
    {
      if (choices->value == dscm)
	return (*(choices->proc)) (xdrs, unp, LASTUNSIGNED);
    }

  /*
   * no match - execute the default xdr routine if there is one
   */
  return ((dfault == NULL_xdrproc_t) ? FALSE :
	  (*dfault) (xdrs, unp, LASTUNSIGNED));
}


/*
 * Non-portable xdr primitives.
 * Care should be taken when moving these routines to new architectures.
 */


/*
 * XDR null terminated ASCII strings
 * xdr_string deals with "C strings" - arrays of bytes that are
 * terminated by a NUL character.  The parameter cpp references a
 * pointer to storage; If the pointer is null, then the necessary
 * storage is allocated.  The last parameter is the max allowed length
 * of the string as specified by a protocol.
 */
bool_t
xdr_string (XDR *xdrs, char **cpp, unsigned int maxsize)
{
  char *sp = *cpp;	/* sp is the actual string pointer */
  /* Initialize to silence the compiler.  It is not really needed because SIZE
     never actually gets used without being initialized.  */
  unsigned int size = 0;
  unsigned int nodesize;

  /*
   * first deal with the length since xdr strings are counted-strings
   */
  switch (xdrs->x_op)
    {
    case XDR_FREE:
      if (sp == NULL)
	{
	  return TRUE;		/* already free */
	}
      /* fall through... */
    case XDR_ENCODE:
      if (sp == NULL)
	return FALSE;
      size = strlen (sp);
      break;
    case XDR_DECODE:
      break;
    }
  if (!xdr_u_int (xdrs, &size))
    {
      return FALSE;
    }
  if (size > maxsize)
    {
      return FALSE;
    }
  nodesize = size + 1;
  if (nodesize == 0)
    {
      /* This means an overflow.  It a bug in the caller which
	 provided a too large maxsize but nevertheless catch it
	 here.  */
      return FALSE;
    }

  /*
   * now deal with the actual bytes
   */
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (sp == NULL)
	*cpp = sp = (char *) mem_alloc (nodesize);
      if (sp == NULL)
	{
          (void)fputs("xdr_string: out of memory\n", stderr);
	  return FALSE;
	}
      sp[size] = 0;
      /* Fall through.  */

    case XDR_ENCODE:
      return xdr_opaque (xdrs, sp, size);

    case XDR_FREE:
      mem_free (sp, nodesize);
      *cpp = NULL;
      return TRUE;
    }
  return FALSE;
}

/*
 * Wrapper for xdr_string that can be called directly from
 * routines like clnt_call
 */
bool_t
xdr_wrapstring (XDR *xdrs, char **cpp)
{
  if (xdr_string (xdrs, cpp, LASTUNSIGNED))
    {
      return TRUE;
    }
  return FALSE;
}
