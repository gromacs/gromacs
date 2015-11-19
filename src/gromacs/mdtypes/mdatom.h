/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \brief
 * Declares mdatom data structure.
 *
 * \inpublicapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_MDATOM_H
#define GMX_MDTYPES_MDATOM_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct t_mdatoms {
    //! \brief Total mass in state A
    real                   tmassA;
    //! \brief Total mass in state B
    real                   tmassB;
    //! \brief Total mass
    real                   tmass;
    //! \brief Number of atoms in arrays
    int                    nr;
    //! \brief Number of elements in arrays
    int                    nalloc;
    //! \brief Number of energy groups
    int                    nenergrp;
    //! \brief Do we have multiple center of mass motion removal groups
    gmx_bool               bVCMgrps;
    //! \brief Number of perturbed atoms
    int                    nPerturbed;
    //! \brief Number of atoms for which the mass is perturbed
    int                    nMassPerturbed;
    //! \brief Number of atoms for which the charge is perturbed
    int                    nChargePerturbed;
    //! \brief Number of atoms for which the type is perturbed
    int                    nTypePerturbed;
    //! \brief Do we have orientation restraints
    gmx_bool               bOrires;
    //! \brief Atomic mass in A state
    real                  *massA;
    //! \brief Atomic mass in B state
    real                  *massB;
    //! \brief Atomic mass in present state
    real                  *massT;
    //! \brief Inverse atomic mass
    real                  *invmass;
    //! \brief Atomic charge in A state
    real                  *chargeA;
    //! \brief Atomic charge in B state
    real                  *chargeB;
    //! \brief Dispersion constant C6 in A state
    real                  *sqrt_c6A;
    //! \brief Dispersion constant C6 in A state
    real                  *sqrt_c6B;
    //! \brief Van der Waals radius sigma in the A state
    real                  *sigmaA;
    //! \brief Van der Waals radius sigma in the B state
    real                  *sigmaB;
    //! \brief Van der Waals radius sigma^3 in the A state
    real                  *sigma3A;
    //! \brief Van der Waals radius sigma^3 in the B state
    real                  *sigma3B;
    //! \brief Is this atom perturbed
    gmx_bool              *bPerturbed;
    //! \brief Type of atom in the A state
    int                   *typeA;
    //! \brief Type of atom in the B state
    int                   *typeB;
    //! \brief Particle type
    unsigned short        *ptype;
    //! \brief Group index for temperature coupling
    unsigned short        *cTC;
    //! \brief Group index for energy matrix
    unsigned short        *cENER;
    //! \brief Group index for acceleration
    unsigned short        *cACC;
    //! \brief Group index for freezing
    unsigned short        *cFREEZE;
    //! \brief Group index for center of mass motion removal
    unsigned short        *cVCM;
    //! \brief Group index for user 1
    unsigned short        *cU1;
    //! \brief Group index for user 2
    unsigned short        *cU2;
    //! \brief Group index for orientation restraints
    unsigned short        *cORF;
    //! \brief QMMM atoms
    gmx_bool              *bQM;
    //! \brief Number of atoms on this processor. TODO is this still used?
    int                    homenr;
    //! \brief The lambda value used to create the contents of the struct
    real                   lambda;
} t_mdatoms;

#ifdef __cplusplus
}
#endif

#endif
