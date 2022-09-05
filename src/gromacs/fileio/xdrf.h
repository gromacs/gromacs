/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_FILEIO_XDRF_H
#define GMX_FILEIO_XDRF_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __PGI /*Portland group compiler*/
#    define int64_t long long
#endif

#include "config.h"

#if GMX_INTERNAL_XDR
#    include "gromacs/fileio/gmx_internal_xdr.h"
#else
#    include <rpc/rpc.h>
#    include <rpc/xdr.h>
#endif

/* Read or write reduced precision *float* coordinates */
int xdr3dfcoord(XDR* xdrs, float* fp, int* size, float* precision);


/* Read or write a *real* value (stored as float) */
int xdr_real(XDR* xdrs, real* r);


/* Read or write reduced precision *real* coordinates */
int xdr3drcoord(XDR* xdrs, real* fp, int* size, real* precision);


//! Read or write a int32_t value.
int xdr_int32(XDR* xdrs, int32_t* i);

//! Read or write a int64_t value.
int xdr_int64(XDR* xdrs, int64_t* i);

int xdr_xtc_seek_time(real time, FILE* fp, XDR* xdrs, int natoms, gmx_bool bSeekForwardOnly);


int xdr_xtc_seek_frame(int frame, FILE* fp, XDR* xdrs, int natoms);


float xdr_xtc_get_last_frame_time(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK);


int xdr_xtc_get_last_frame_number(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK);

#endif
