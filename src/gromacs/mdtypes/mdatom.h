/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

typedef struct t_mdatoms {
    //! Total mass in state A
    real                   tmassA;
    //! Total mass in state B
    real                   tmassB;
    //! Total mass
    real                   tmass;
    //! Number of atoms in arrays
    int                    nr;
    //! Number of elements in arrays
    int                    nalloc;
    //! Number of energy groups
    int                    nenergrp;
    //! Do we have multiple center of mass motion removal groups
    gmx_bool               bVCMgrps;
    //! Do we have any virtual sites?
    gmx_bool               haveVsites;
    //! Do we have atoms that are frozen along 1 or 2 (not 3) dimensions?
    gmx_bool               havePartiallyFrozenAtoms;
    //! Number of perturbed atoms
    int                    nPerturbed;
    //! Number of atoms for which the mass is perturbed
    int                    nMassPerturbed;
    //! Number of atoms for which the charge is perturbed
    int                    nChargePerturbed;
    //! Number of atoms for which the type is perturbed
    int                    nTypePerturbed;
    //! Do we have orientation restraints
    gmx_bool               bOrires;
    //! Atomic mass in A state
    real                  *massA;
    //! Atomic mass in B state
    real                  *massB;
    //! Atomic mass in present state
    real                  *massT;
    //! Inverse atomic mass per atom, 0 for vsites and shells
    real                  *invmass;
    //! Inverse atomic mass per atom and dimension, 0 for vsites, shells and frozen dimensions
    rvec                  *invMassPerDim;
    //! Atomic charge in A state
    real                  *chargeA;
    //! Atomic charge in B state
    real                  *chargeB;
    //! Dispersion constant C6 in A state
    real                  *sqrt_c6A;
    //! Dispersion constant C6 in A state
    real                  *sqrt_c6B;
    //! Van der Waals radius sigma in the A state
    real                  *sigmaA;
    //! Van der Waals radius sigma in the B state
    real                  *sigmaB;
    //! Van der Waals radius sigma^3 in the A state
    real                  *sigma3A;
    //! Van der Waals radius sigma^3 in the B state
    real                  *sigma3B;
    //! Is this atom perturbed
    gmx_bool              *bPerturbed;
    //! Type of atom in the A state
    int                   *typeA;
    //! Type of atom in the B state
    int                   *typeB;
    //! Particle type
    unsigned short        *ptype;
    //! Group index for temperature coupling
    unsigned short        *cTC;
    //! Group index for energy matrix
    unsigned short        *cENER;
    //! Group index for acceleration
    unsigned short        *cACC;
    //! Group index for freezing
    unsigned short        *cFREEZE;
    //! Group index for center of mass motion removal
    unsigned short        *cVCM;
    //! Group index for user 1
    unsigned short        *cU1;
    //! Group index for user 2
    unsigned short        *cU2;
    //! Group index for orientation restraints
    unsigned short        *cORF;
    //! QMMM atoms
    gmx_bool              *bQM;
    //! Number of atoms on this processor. TODO is this still used?
    int                    homenr;
    //! The lambda value used to create the contents of the struct
    real                   lambda;
} t_mdatoms;

#endif
