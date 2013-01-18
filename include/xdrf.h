/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _xdrf_h
#define _xdrf_h


#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"

#ifdef __PGI    /*Portland group compiler*/
#define int64_t long long
#endif

#include "gmx_header_config.h"
#if (defined GMX_NATIVE_WINDOWS || defined GMX_CYGWIN || defined GMX_INTERNAL_XDR)
#include "gmx_system_xdr.h"
#else
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* THESE 3 FUNCTIONS (xdropen, xdrclose and xdr_get_fp) ARE NOW OBSOLETE
   AND ONLY PROVIDED FOR BACKWARD COMPATIBILITY OF 3D PARTY TOOLS.
   THEY SHOULD NOT BE USED ANYWHERE IN GROMACS ITSELF.
   int xdropen(XDR *xdrs, const char *filename, const char *type);
   int xdrclose(XDR *xdrs);
 */

/* the xdr data types; note that there is no data type 'real' because
   here we deal with the types as they are actually written to disk.  */
typedef enum
{
    xdr_datatype_int,
    xdr_datatype_float,
    xdr_datatype_double,
    xdr_datatype_large_int,
    xdr_datatype_char,
    xdr_datatype_string
} xdr_datatype;

/* names corresponding to the xdr_datatype enum */
GMX_LIBGMX_EXPORT
extern const char *xdr_datatype_names[];

/* Read or write reduced precision *float* coordinates */
int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision);


/* Read or write a *real* value (stored as float) */
int xdr_real(XDR *xdrs, real *r);


/* Read or write reduced precision *real* coordinates */
int xdr3drcoord(XDR *xdrs, real *fp, int *size, real *precision);


int xdr_gmx_large_int(XDR *xdrs, gmx_large_int_t *i, const char *warn);
/* Read or write a gmx_large_int_t value.
 * 32bit code reading a 64bit gmx_large_int_t value from xdrs could
 * lead to values out of int range.
 * When warn!=NULL a warning will be written to stderr
 * when a value does not fit,
 * the first line is:
 * "WARNING during %s:", where warn is printed in %s.
 */

int xdr_xtc_seek_time(real time, FILE *fp, XDR *xdrs, int natoms, gmx_bool bSeekForwardOnly);


int xdr_xtc_seek_frame(int frame, FILE *fp, XDR *xdrs, int natoms);


GMX_LIBGMX_EXPORT
float xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK);


int xdr_xtc_get_last_frame_number(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK);

#ifdef __cplusplus
}
#endif

#endif
