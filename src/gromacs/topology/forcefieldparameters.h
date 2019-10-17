/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_FORCEFIELDPARAMETERS_H
#define GMX_TOPOLOGY_FORCEFIELDPARAMETERS_H

#include <cstdio>

#include <vector>

#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

/*! \brief Struct that holds all force field parameters for the simulated system */
struct gmx_ffparams_t
{
    /*! \brief Returns the number of function types, which matches the number of elements in iparams */
    int numTypes() const
    {
        GMX_ASSERT(iparams.size() == functype.size(), "Parameters and function types go together");

        return gmx::ssize(functype);
    }

    /* TODO: Consider merging functype and iparams, either by storing
     *       the functype in t_iparams or by putting both in a single class.
     */
    int                     atnr = 0;    /**< The number of non-bonded atom types */
    std::vector<t_functype> functype;    /**< The function type per type */
    std::vector<t_iparams>  iparams;     /**< Force field parameters per type */
    double                  reppow  = 0; /**< The repulsion power for VdW: C12*r^-reppow   */
    real                    fudgeQQ = 0; /**< The scaling factor for Coulomb 1-4: f*q1*q2  */
    gmx_cmap_t              cmap_grid;   /**< The dihedral correction maps                 */
};

void pr_ffparams(FILE* fp, int indent, const char* title, const gmx_ffparams_t* ffparams, gmx_bool bShowNumbers);

#endif
