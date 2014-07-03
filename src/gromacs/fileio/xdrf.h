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
#ifndef GMX_FILEIO_XDRF_H
#define GMX_FILEIO_XDRF_H

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __PGI    /*Portland group compiler*/
#define int64_t long long
#endif

#include "config.h"

#ifdef GMX_INTERNAL_XDR
#include "gromacs/fileio/gmx_system_xdr.h"
#else
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct t_fileio;

/* Read or write reduced precision *float* coordinates */
int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision);


/* Read or write a *real* value (stored as float) */
int xdr_real(XDR *xdrs, real *r);


/* Read or write reduced precision *real* coordinates */
int xdr3drcoord(XDR *xdrs, real *fp, int *size, real *precision);


int xdr_int64(XDR *xdrs, gmx_int64_t *i);
/* Read or write a gmx_int64_t value.
 * When warn!=NULL a warning will be written to stderr
 * when a value does not fit,
 * the first line is:
 * "WARNING during %s:", where warn is printed in %s.
 */

int xdr_xtc_seek_time(real time, FILE *fp, XDR *xdrs, int natoms, gmx_bool bSeekForwardOnly);


int xdr_xtc_seek_frame(int frame, FILE *fp, XDR *xdrs, int natoms);


float xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK);


int xdr_xtc_get_last_frame_number(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK);


/* Defined in gmxfio.c.
 * TODO: It would be nice to decouple this header from t_fileio completely,
 * and not need the XDR struct in gmxfio.h, but that would require some
 * extra code that is not warranted for this single function.
 * Can be reconsidered if the file I/O gets refactored in the future.
 */
XDR *gmx_fio_getxdr(struct t_fileio *fio);
/* Return the file pointer itself */

#ifdef __cplusplus
}
#endif

#endif
