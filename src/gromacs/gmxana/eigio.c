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

#include "eigio.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

void read_eigenvectors(const char *file, int *natoms, gmx_bool *bFit,
                       rvec **xref, gmx_bool *bDMR,
                       rvec **xav, gmx_bool *bDMA,
                       int *nvec, int **eignr,
                       rvec ***eigvec, real **eigval)
{
    t_trnheader head;
    int         i, snew_size;
    t_fileio   *status;
    rvec       *x;
    matrix      box;
    gmx_bool    bOK;

    *bDMR = FALSE;

    /* read (reference (t=-1) and) average (t=0) structure */
    status = open_trn(file, "r");
    fread_trnheader(status, &head, &bOK);
    *natoms = head.natoms;
    snew(*xav, *natoms);
    fread_htrn(status, &head, box, *xav, NULL, NULL);

    if ((head.t >= -1.1) && (head.t <= -0.9))
    {
        snew(*xref, *natoms);
        for (i = 0; i < *natoms; i++)
        {
            copy_rvec((*xav)[i], (*xref)[i]);
        }
        *bDMR = (head.lambda > 0.5);
        *bFit = (head.lambda > -0.5);
        if (*bFit)
        {
            fprintf(stderr, "Read %smass weighted reference structure with %d atoms from %s\n", *bDMR ? "" : "non ", *natoms, file);
        }
        else
        {
            fprintf(stderr, "Eigenvectors in %s were determined without fitting\n", file);
            sfree(*xref);
            *xref = NULL;
        }
        fread_trnheader(status, &head, &bOK);
        fread_htrn(status, &head, box, *xav, NULL, NULL);
    }
    else
    {
        *bFit = TRUE;
        *xref = NULL;
    }
    *bDMA = (head.lambda > 0.5);
    if ((head.t <= -0.01) || (head.t >= 0.01))
    {
        fprintf(stderr, "WARNING: %s does not start with t=0, which should be the "
                "average structure. This might not be a eigenvector file. "
                "Some things might go wrong.\n",
                file);
    }
    else
    {
        fprintf(stderr,
                "Read %smass weighted average/minimum structure with %d atoms from %s\n",
                *bDMA ? "" : "non ", *natoms, file);
    }

    snew(x, *natoms);
    snew_size = 10;
    snew(*eignr, snew_size);
    snew(*eigval, snew_size);
    snew(*eigvec, snew_size);

    *nvec = 0;
    while (fread_trnheader(status, &head, &bOK))
    {
        fread_htrn(status, &head, box, x, NULL, NULL);
        if (*nvec >= snew_size)
        {
            snew_size += 10;
            srenew(*eignr, snew_size);
            srenew(*eigval, snew_size);
            srenew(*eigvec, snew_size);
        }
        i                = head.step;
        (*eigval)[*nvec] = head.t;
        (*eignr)[*nvec]  = i-1;
        snew((*eigvec)[*nvec], *natoms);
        for (i = 0; i < *natoms; i++)
        {
            copy_rvec(x[i], (*eigvec)[*nvec][i]);
        }
        (*nvec)++;
    }
    sfree(x);
    fprintf(stderr, "Read %d eigenvectors (for %d atoms)\n\n", *nvec, *natoms);
}


void write_eigenvectors(const char *trnname, int natoms, real mat[],
                        gmx_bool bReverse, int begin, int end,
                        int WriteXref, rvec *xref, gmx_bool bDMR,
                        rvec xav[], gmx_bool bDMA, real eigval[])
{
    t_fileio *trnout;
    int       ndim, i, j, d, vec;
    matrix    zerobox;
    rvec     *x;

    ndim = natoms*DIM;
    clear_mat(zerobox);
    snew(x, natoms);

    fprintf (stderr,
             "\nWriting %saverage structure & eigenvectors %d--%d to %s\n",
             (WriteXref == eWXR_YES) ? "reference, " : "",
             begin, end, trnname);

    trnout = open_tpx(trnname, "w");
    if (WriteXref == eWXR_YES)
    {
        /* misuse lambda: 0/1 mass weighted fit no/yes */
        fwrite_trn(trnout, -1, -1, bDMR ? 1.0 : 0.0, zerobox, natoms, xref, NULL, NULL);
    }
    else if (WriteXref == eWXR_NOFIT)
    {
        /* misuse lambda: -1 no fit */
        fwrite_trn(trnout, -1, -1, -1.0, zerobox, natoms, x, NULL, NULL);
    }

    /* misuse lambda: 0/1 mass weighted analysis no/yes */
    fwrite_trn(trnout, 0, 0, bDMA ? 1.0 : 0.0, zerobox, natoms, xav, NULL, NULL);

    for (i = 0; i <= (end-begin); i++)
    {

        if (!bReverse)
        {
            vec = i;
        }
        else
        {
            vec = ndim-i-1;
        }

        for (j = 0; j < natoms; j++)
        {
            for (d = 0; d < DIM; d++)
            {
                x[j][d] = mat[vec*ndim+DIM*j+d];
            }
        }

        /* Store the eigenvalue in the time field */
        fwrite_trn(trnout, begin+i, eigval[vec], 0, zerobox, natoms, x, NULL, NULL);
    }
    close_trn(trnout);

    sfree(x);
}
