/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file defines functions for implementing dispersion
 * corrections
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "dispersioncorrection.h"

#include <cstdio>

#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

t_forcetable *makeDispersionCorrectionTable(FILE *fp,
                                            t_forcerec *fr, real rtab,
                                            const char *tabfn)
{
    t_forcetable *dispersionCorrectionTable = NULL;

    if (tabfn == NULL)
    {
        if (debug)
        {
            fprintf(debug, "No table file name passed, can not read table, can not do non-bonded interactions\n");
        }
        return dispersionCorrectionTable;
    }

    t_forcetable *fullTable = make_tables(fp, fr, tabfn, rtab, 0);
    /* Copy the contents of the table to one that has just dispersion
     * and repulsion, to improve cache performance. We want the table
     * data to be aligned to 32-byte boundaries. The pointers could be
     * freed but currently aren't. */
    snew(dispersionCorrectionTable, 1);
    dispersionCorrectionTable->interaction   = GMX_TABLE_INTERACTION_VDWREP_VDWDISP;
    dispersionCorrectionTable->format        = fullTable->format;
    dispersionCorrectionTable->r             = fullTable->r;
    dispersionCorrectionTable->n             = fullTable->n;
    dispersionCorrectionTable->scale         = fullTable->scale;
    dispersionCorrectionTable->formatsize    = fullTable->formatsize;
    dispersionCorrectionTable->ninteractions = 2;
    dispersionCorrectionTable->stride        = dispersionCorrectionTable->formatsize * dispersionCorrectionTable->ninteractions;
    snew_aligned(dispersionCorrectionTable->data, dispersionCorrectionTable->stride*(dispersionCorrectionTable->n+1), 32);

    for (int i = 0; i <= fullTable->n; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            dispersionCorrectionTable->data[8*i+j] = fullTable->data[12*i+4+j];
        }
    }
    sfree_aligned(fullTable->data);
    sfree(fullTable);

    return dispersionCorrectionTable;
}
