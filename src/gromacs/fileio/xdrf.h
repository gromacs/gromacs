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

#include "config.h"

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#if GMX_INTERNAL_XDR
#    include "rpc_xdr/xdr.h"
#else
#    include <rpc/rpc.h>
#    include <rpc/xdr.h>
#endif

// Magic number used for traditional XTC files with 32-bit data buffer size
#define XTC_MAGIC 1995
// New magic number used for (very) large XTC files with 64-bit data buffer size
#define XTC_NEW_MAGIC 2023

/* Until june 2023, the old XDR format could only store up to ~300M atoms.
 * To handle larger systems, we use a newer magic number (2023 instead of 1995).
 * By default we only need this for large systems, but to enable us to check
 * the code without using gigantic test files, we select the version to
 * read/write by providing either XDR_MAGIC XDR_NEW_MAGIC as the last argument.
 * This should always match what you read/write in the xdr file header.
 *
 * The reason for the format update is that XTC logic allocates a buffer 1.2x
 * larger than the natoms*3 and the size - in bytes - of the USED
 * part of this 32-bit-integer-containing buffer in bytes needs to be stored as
 * an integer. This means we MIGHT need 64-bit sizing if natoms is larger than
 * (2^32)/(3*4*1.2)=298261617 atoms, and since the exact size will vary from
 * frame-to-frame, we need 64-bit indexing whenever it might happen,
 * so in that case the first entry will use 8 bytes and be incompatible with
 * older 32-bit reading code. It is up to the calling code to ensure you select
 * the new magic number for systems larger than this - but we will double-check
 * it internally and exit if you forgot.
 */
#define XTC_1995_MAX_NATOMS 298261617

/* Read or write reduced precision *float* coordinates */
int xdr3dfcoord(XDR* xdrs, float* fp, int* size, float* precision, int magic_number);


/* Read or write a *real* value (stored as float) */
int xdr_real(XDR* xdrs, real* r);

/* Read or write reduced precision *real* coordinates.
 */
int xdr3drcoord(XDR* xdrs, real* fp, int* size, real* precision, int magic_number);


//! Read or write a int32_t value.
int xdr_int32(XDR* xdrs, int32_t* i);

//! Read or write a int64_t value.
int xdr_int64(XDR* xdrs, int64_t* i);

int xdr_xtc_seek_time(real time, FILE* fp, XDR* xdrs, int natoms, gmx_bool bSeekForwardOnly);


int xdr_xtc_seek_frame(int frame, FILE* fp, XDR* xdrs, int natoms);


float xdr_xtc_get_last_frame_time(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK);


int xdr_xtc_get_last_frame_number(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK);

#endif
