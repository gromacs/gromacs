/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*!\internal
 * \file
 * \brief
 * Implements outputselector class.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "outputselector.h"

#include <cstddef>

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

/*! \brief
 * Modify atoms information in coordinate frame to fit output selection.
 *
 * Changes the information contained in the coordinate frame t_atoms struct to match
 * the selection provided to the module.
 *
 * \param[in] atoms Pointer to original t_atoms.
 * \param[in] selectionAtoms Pointer to local atoms.
 * \param[in]     sel   Reference to selection.
 */
static void adjustAtomInformation(t_atoms* atoms, t_atoms* selectionAtoms, const Selection& sel)
{
    int natoms = sel.atomCount();

    selectionAtoms->nr   = natoms;
    selectionAtoms->nres = atoms->nres;

    srenew(selectionAtoms->resinfo, atoms->nres);
    srenew(selectionAtoms->atom, natoms);
    srenew(selectionAtoms->atomname, natoms);

    selectionAtoms->haveType    = atoms->haveType;
    selectionAtoms->haveBState  = atoms->haveBState;
    selectionAtoms->haveMass    = atoms->haveMass;
    selectionAtoms->haveCharge  = atoms->haveCharge;
    selectionAtoms->havePdbInfo = atoms->havePdbInfo;

    if (atoms->haveType)
    {
        srenew(selectionAtoms->atomtype, natoms);
    }
    if (atoms->haveBState)
    {
        srenew(selectionAtoms->atomtypeB, natoms);
    }
    if (atoms->havePdbInfo)
    {
        srenew(selectionAtoms->pdbinfo, natoms);
    }

    for (int i = 0; i < natoms; i++)
    {
        int pos                     = sel.position(i).refId();
        selectionAtoms->atom[i]     = atoms->atom[pos];
        selectionAtoms->atomname[i] = atoms->atomname[pos];
        if (selectionAtoms->haveType)
        {
            selectionAtoms->atomtype[i] = atoms->atomtype[pos];
        }
        if (selectionAtoms->haveBState)
        {
            selectionAtoms->atomtypeB[i] = atoms->atomtypeB[pos];
        }
        if (selectionAtoms->havePdbInfo)
        {
            selectionAtoms->pdbinfo[i] = atoms->pdbinfo[pos];
        }
    }
    // Copy residue information unconditionally.
    for (int i = 0; i < atoms->nres; i++)
    {
        selectionAtoms->resinfo[i] = atoms->resinfo[i];
    }
}

void OutputSelector::processFrame(const int /*framenumber*/, t_trxframe* input)
{
    size_t natoms = sel_.atomCount();

    input->natoms = natoms;

    localX_.resize(natoms);
    if (input->bV)
    {
        localV_.resize(natoms);
    }
    if (input->bF)
    {
        localF_.resize(natoms);
    }
    if (input->index)
    {
        localIndex_.resize(natoms);
    }

    for (size_t i = 0; i < natoms; i++)
    {
        int pos = sel_.position(i).refId();

        copy_rvec(input->x[pos], localX_[i]);
        if (input->bV)
        {
            copy_rvec(input->v[pos], localV_[i]);
        }
        if (input->bF)
        {
            copy_rvec(input->f[pos], localF_[i]);
        }
        if (input->index)
        {
            localIndex_[i] = input->index[pos];
        }
    }
    input->x = as_rvec_array(localX_.data());

    input->index = localIndex_.data();

    if (input->bV)
    {
        input->v = as_rvec_array(localV_.data());
    }
    if (input->bF)
    {
        input->f = as_rvec_array(localF_.data());
    }

    if (input->bAtoms)
    {
        if (selectionAtoms_ == nullptr)
        {
            selectionAtoms_.reset(new t_atoms);
            init_t_atoms(selectionAtoms_.get(), natoms, false);
        }
        adjustAtomInformation(input->atoms, selectionAtoms_.get(), sel_);
        input->atoms = selectionAtoms_.get();
    }
}

} // namespace gmx
