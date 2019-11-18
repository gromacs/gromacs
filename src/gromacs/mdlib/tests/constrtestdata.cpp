/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief SHAKE and LINCS tests.
 *
 * \todo Better tests for virial are needed.
 * \todo Tests for bigger systems to test threads synchronization,
 *       reduction, etc. on the GPU.
 * \todo Tests for algorithms for derivatives.
 * \todo Free-energy perturbation tests
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "constrtestdata.h"

#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{

ConstraintsTestData::ConstraintsTestData(const std::string&       title,
                                         int                      numAtoms,
                                         std::vector<real>        masses,
                                         std::vector<int>         constraints,
                                         std::vector<real>        constraintsR0,
                                         bool                     computeVirial,
                                         tensor                   virialScaledRef,
                                         bool                     compute_dHdLambda,
                                         float                    dHdLambdaRef,
                                         real                     initialTime,
                                         real                     timestep,
                                         const std::vector<RVec>& x,
                                         const std::vector<RVec>& xPrime,
                                         const std::vector<RVec>& v,
                                         real                     shakeTolerance,
                                         gmx_bool                 shakeUseSOR,
                                         int                      lincsNumIterations,
                                         int                      lincsExpansionOrder,
                                         real                     lincsWarnAngle)
{
    title_    = title;    // Human-friendly name of the system
    numAtoms_ = numAtoms; // Number of atoms

    // Masses of atoms
    masses_ = masses;
    invmass_.resize(numAtoms); // Vector of inverse masses

    for (int i = 0; i < numAtoms; i++)
    {
        invmass_[i] = 1.0 / masses.at(i);
    }

    // Saving constraints to check if they are satisfied after algorithm was applied
    constraints_   = constraints;   // Constraints indices (in type-i-j format)
    constraintsR0_ = constraintsR0; // Equilibrium distances for each type of constraint

    invdt_ = 1.0 / timestep; // Inverse timestep

    // Communication record
    cr_.nnodes = 1;
    cr_.dd     = nullptr;

    // Multisim data
    ms_.sim  = 0;
    ms_.nsim = 1;

    // Input record - data that usually comes from configuration file (.mdp)
    ir_.efep    = 0;
    ir_.init_t  = initialTime;
    ir_.delta_t = timestep;
    ir_.eI      = 0;

    // MD atoms data
    md_.nMassPerturbed = 0;
    md_.lambda         = 0.0;
    md_.invmass        = invmass_.data();
    md_.nr             = numAtoms;
    md_.homenr         = numAtoms;

    // Virial evaluation
    computeVirial_ = computeVirial;
    if (computeVirial)
    {
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                virialScaled_[i][j]    = 0;
                virialScaledRef_[i][j] = virialScaledRef[i][j];
            }
        }
    }


    // Free energy evaluation
    compute_dHdLambda_ = compute_dHdLambda;
    dHdLambda_         = 0;
    if (compute_dHdLambda_)
    {
        ir_.efep      = efepYES;
        dHdLambdaRef_ = dHdLambdaRef;
    }
    else
    {
        ir_.efep      = efepNO;
        dHdLambdaRef_ = 0;
    }

    // Constraints and their parameters (local topology)
    for (int i = 0; i < F_NRE; i++)
    {
        idef_.il[i].nr = 0;
    }
    idef_.il[F_CONSTR].nr = constraints.size();

    snew(idef_.il[F_CONSTR].iatoms, constraints.size());
    int maxType = 0;
    for (index i = 0; i < ssize(constraints); i++)
    {
        if (i % 3 == 0)
        {
            if (maxType < constraints.at(i))
            {
                maxType = constraints.at(i);
            }
        }
        idef_.il[F_CONSTR].iatoms[i] = constraints.at(i);
    }
    snew(idef_.iparams, maxType + 1);
    for (index i = 0; i < ssize(constraints) / 3; i++)
    {
        idef_.iparams[constraints.at(3 * i)].constr.dA = constraintsR0.at(constraints.at(3 * i));
        idef_.iparams[constraints.at(3 * i)].constr.dB = constraintsR0.at(constraints.at(3 * i));
    }

    // Constraints and their parameters (global topology)
    InteractionList interactionList;
    interactionList.iatoms.resize(constraints.size());
    std::copy(constraints.begin(), constraints.end(), interactionList.iatoms.begin());
    InteractionList interactionListEmpty;
    interactionListEmpty.iatoms.resize(0);

    gmx_moltype_t molType;
    molType.atoms.nr             = numAtoms;
    molType.ilist.at(F_CONSTR)   = interactionList;
    molType.ilist.at(F_CONSTRNC) = interactionListEmpty;
    mtop_.moltype.push_back(molType);

    gmx_molblock_t molBlock;
    molBlock.type = 0;
    molBlock.nmol = 1;
    mtop_.molblock.push_back(molBlock);

    mtop_.natoms = numAtoms;
    mtop_.ffparams.iparams.resize(maxType + 1);
    for (int i = 0; i <= maxType; i++)
    {
        mtop_.ffparams.iparams.at(i) = idef_.iparams[i];
    }
    mtop_.bIntermolecularInteractions = false;

    // Coordinates and velocities
    x_.resizeWithPadding(numAtoms);
    xPrime_.resizeWithPadding(numAtoms);
    xPrime0_.resizeWithPadding(numAtoms);
    xPrime2_.resizeWithPadding(numAtoms);

    v_.resizeWithPadding(numAtoms);
    v0_.resizeWithPadding(numAtoms);

    std::copy(x.begin(), x.end(), x_.begin());
    std::copy(xPrime.begin(), xPrime.end(), xPrime_.begin());
    std::copy(xPrime.begin(), xPrime.end(), xPrime0_.begin());
    std::copy(xPrime.begin(), xPrime.end(), xPrime2_.begin());

    std::copy(v.begin(), v.end(), v_.begin());
    std::copy(v.begin(), v.end(), v0_.begin());

    // SHAKE-specific parameters
    ir_.shake_tol = shakeTolerance;
    ir_.bShakeSOR = shakeUseSOR;

    // LINCS-specific parameters
    ir_.nLincsIter     = lincsNumIterations;
    ir_.nProjOrder     = lincsExpansionOrder;
    ir_.LincsWarnAngle = lincsWarnAngle;
}

/*! \brief
 * Reset the data structure so it can be reused.
 *
 * Set the coordinates and velocities back to their values before
 * constraining. The scaled virial tensor and dHdLambda are zeroed.
 *
 */
void ConstraintsTestData::reset()
{
    xPrime_  = xPrime0_;
    xPrime2_ = xPrime0_;
    v_       = v0_;

    if (computeVirial_)
    {
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                virialScaled_[i][j] = 0;
            }
        }
    }
    dHdLambda_ = 0;
}

/*! \brief
 * Cleaning up the memory.
 */
ConstraintsTestData::~ConstraintsTestData()
{
    sfree(idef_.il[F_CONSTR].iatoms);
    sfree(idef_.iparams);
}

} // namespace test
} // namespace gmx
