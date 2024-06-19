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
#include "gmxpre.h"

#include "mtxio.h"

#include <cstdio>

/* This module provides routines to read/write sparse or full storage
 * matrices from/to files. It is normally used for the Hessian matrix
 * in normal mode analysis.
 */

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* Just a number to identify our file type */
#define GMX_MTXIO_MAGIC_NUMBER 0x34ce8fd2

#define GMX_MTXIO_FULL_MATRIX 0
#define GMX_MTXIO_SPARSE_MATRIX 1


/* Matrix file format definition:
 *
 * All entries are stored in XDR format.
 *
 * 1. Magic number integer, should be GMX_MTXIO_MAGIC_NUMBER
 * 2. An XDR string specifying the Gromacs version used to generate the file.
 * 3. Integer to denote precision. 1 if double, 0 if single precision.
 * 4. Two integers specifying number of rows and columns.
 * 5. Integer to denote storage type:
 *    GMX_MTXIO_FULL_MATRIX or GMX_MTXIO_SPARSE_MATRIX
 *
 * 6. Matrix data.
 *    a) In case of full matrix, this is nrow*ncol floating-point values.
 *    b) In case of sparse matrix the data is:
 *       - Integer specifying compressed_symmetric format (1=yes, 0=no)
 *       - Integer specifying number of rows (again)
 *       - nrow integers specifying the number of data entries on each row ("ndata")
 *       - All the actual entries for each row, stored contiguous.
 *         Each entry consists of an integer column index and floating-point data value.
 */

void gmx_mtxio_write(const std::filesystem::path& filename,
                     int                          nrow,
                     int                          ncol,
                     real*                        full_matrix,
                     gmx_sparsematrix_t*          sparse_matrix)
{
    t_fileio* fio;
    int       i, j, prec;
    size_t    sz;

    if (full_matrix != nullptr && sparse_matrix != nullptr)
    {
        gmx_fatal(FARGS, "Both full AND sparse matrix specified to gmx_mtxio_write().\n");
    }

    fio = gmx_fio_open(filename, "w");

    /* Write magic number */
    i = GMX_MTXIO_MAGIC_NUMBER;
    gmx_fio_do_int(fio, i);

    /* Write generating Gromacs version */
    gmx_fio_write_string(fio, gmx_version());

    /* Write 1 for double, 0 for single precision */
    if (sizeof(real) == sizeof(double))
    {
        prec = 1;
    }
    else
    {
        prec = 0;
    }
    gmx_fio_do_int(fio, prec);

    gmx_fio_do_int(fio, nrow);
    gmx_fio_do_int(fio, ncol);

    if (full_matrix != nullptr)
    {
        /* Full matrix storage format */
        i = GMX_MTXIO_FULL_MATRIX;
        gmx_fio_do_int(fio, i);
        sz = nrow * ncol;
        gmx_fio_ndo_real(fio, full_matrix, sz);
    }
    else
    {
        /* Sparse storage */
        i = GMX_MTXIO_SPARSE_MATRIX;
        gmx_fio_do_int(fio, i);

        gmx_fio_do_gmx_bool(fio, sparse_matrix->compressed_symmetric);
        gmx_fio_do_int(fio, sparse_matrix->nrow);
        if (sparse_matrix->nrow != nrow)
        {
            gmx_fatal(FARGS, "Internal inconsistency in sparse matrix.\n");
        }
        gmx_fio_ndo_int(fio, sparse_matrix->ndata, sparse_matrix->nrow);
        for (i = 0; i < sparse_matrix->nrow; i++)
        {
            for (j = 0; j < sparse_matrix->ndata[i]; j++)
            {
                gmx_fio_do_int(fio, sparse_matrix->data[i][j].col);
                gmx_fio_do_real(fio, sparse_matrix->data[i][j].value);
            }
        }
    }
    gmx_fio_close(fio);
}


void gmx_mtxio_read(const std::filesystem::path& filename,
                    int*                         nrow,
                    int*                         ncol,
                    real**                       full_matrix,
                    gmx_sparsematrix_t**         sparse_matrix)
{
    t_fileio* fio;
    int       i, j, prec;
    char      gmxver[256];
    size_t    sz;

    fio = gmx_fio_open(filename, "r");

    /* Read and check magic number */
    i = GMX_MTXIO_MAGIC_NUMBER;
    gmx_fio_do_int(fio, i);

    if (i != GMX_MTXIO_MAGIC_NUMBER)
    {
        gmx_fatal(FARGS,
                  "No matrix data found in file. Note that the Hessian matrix format changed\n"
                  "in GROMACS 3.3 to enable portable files and sparse matrix storage.\n");
    }

    /* Read generating Gromacs version */
    gmx_fio_do_string(fio, gmxver);

    /* Write 1 for double, 0 for single precision */
    if (sizeof(real) == sizeof(double))
    {
        prec = 1;
    }
    else
    {
        prec = 0;
    }
    gmx_fio_do_int(fio, prec);

    fprintf(stderr,
            "Reading %s precision matrix generated by GROMACS %s\n",
            (prec == 1) ? "double" : "single",
            gmxver);

    gmx_fio_do_int(fio, i);
    *nrow = i;
    gmx_fio_do_int(fio, i);
    *ncol = i;

    gmx_fio_do_int(fio, i);

    if (i == GMX_MTXIO_FULL_MATRIX && nullptr != full_matrix)
    {
        printf("Full matrix storage format, nrow=%d, ncols=%d\n", *nrow, *ncol);

        sz = (*nrow) * (*ncol);
        snew((*full_matrix), sz);
        gmx_fio_ndo_real(fio, (*full_matrix), sz);
    }
    else if (nullptr != sparse_matrix)
    {
        /* Sparse storage */
        printf("Sparse matrix storage format, nrow=%d, ncols=%d\n", *nrow, *ncol);

        snew((*sparse_matrix), 1);
        gmx_fio_do_gmx_bool(fio, (*sparse_matrix)->compressed_symmetric);
        gmx_fio_do_int(fio, (*sparse_matrix)->nrow);
        if ((*sparse_matrix)->nrow != *nrow)
        {
            gmx_fatal(FARGS, "Internal inconsistency in sparse matrix.\n");
        }
        snew((*sparse_matrix)->ndata, (*sparse_matrix)->nrow);
        snew((*sparse_matrix)->nalloc, (*sparse_matrix)->nrow);
        snew((*sparse_matrix)->data, (*sparse_matrix)->nrow);
        gmx_fio_ndo_int(fio, (*sparse_matrix)->ndata, (*sparse_matrix)->nrow);

        for (i = 0; i < (*sparse_matrix)->nrow; i++)
        {
            (*sparse_matrix)->nalloc[i] = (*sparse_matrix)->ndata[i] + 10;
            snew(((*sparse_matrix)->data[i]), (*sparse_matrix)->nalloc[i]);

            for (j = 0; j < (*sparse_matrix)->ndata[i]; j++)
            {
                gmx_fio_do_int(fio, (*sparse_matrix)->data[i][j].col);
                gmx_fio_do_real(fio, (*sparse_matrix)->data[i][j].value);
            }
        }
    }
    gmx_fio_close(fio);
}
