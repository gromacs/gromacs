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
/*! \internal \file
 * \brief
 * Internal definitions shared by gmxfio*.c files.
 */
#ifndef GMX_FILEIO_GMXFIO_IMPL_H
#define GMX_FILEIO_GMXFIO_IMPL_H

/* This is the new improved and thread safe version of gmxfio.  */


/* WARNING WARNING WARNING WARNING
   The data types used here are PRIVATE to gmxfio routines. DO NOT use them
   directly in your own code, but use the external functions provided in
   include/gmxfio.h

   If you don't heed this warning, your code will suddenly stop working
   at some point in the not-so-distant future.

   WARNING WARNING WARNING WARNING */

#include <filesystem>

#include "thread_mpi/lock.h"

#include "gromacs/fileio/xdrf.h"

struct t_fileio
{
    FILE*    fp;                   /* the file pointer */
    gmx_bool bRead,                /* the file is open for reading */
            bDouble,               /* write doubles instead of floats */
            bReadWrite;            /* the file is open for reading and writing */
    std::filesystem::path fn;      /* the file name */
    XDR*                  xdr;     /* the xdr data pointer */
    enum xdr_op           xdrmode; /* the xdr mode */
    int                   iFTP;    /* the file type identifier */

    t_fileio *next, *prev; /* next and previous file pointers in the
                              linked list */
    tMPI_Lock_t mtx;       /* content locking mutex. This is a fast lock
                              for performance reasons: in some cases every
                              single byte that gets read/written requires
                              a lock */
};

/** lock the mutex associated with a fio  */
void gmx_fio_lock(t_fileio* fio);
/** unlock the mutex associated with a fio  */
void gmx_fio_unlock(t_fileio* fio);

#endif
