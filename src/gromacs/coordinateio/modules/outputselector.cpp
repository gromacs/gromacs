/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Implements outputselector class.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "outputselector.h"

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/*! \brief
 * Modify atoms information in coordinate frame to fit output selection.
 *
 * Changes the information contained in the coordinate frame t_atoms struct to match
 * the selection provided to the module.
 *
 * \param[in,out] atoms Pointer to t_atoms to modify.
 * \param[in]     sel   Reference to selection.
 */
static t_atoms *adjustAtomInformation(t_atoms *atoms, const Selection &sel)
{
    int      natoms = sel.atomCount();

    t_atoms *newAtoms;
    snew(newAtoms, 1);
    init_t_atoms(newAtoms, natoms, atoms->havePdbInfo);
    srenew(newAtoms->resinfo, atoms->nres);
    newAtoms->haveMass   = atoms->haveMass;
    newAtoms->haveCharge = atoms->haveCharge;
    newAtoms->haveType   = atoms->haveType;
    newAtoms->haveBState = atoms->haveBState;

    for (int i = 0; i < natoms; i++)
    {
        int pos = sel.position(i).refId();
        newAtoms->atom[i]     = atoms->atom[pos];
        newAtoms->atomname[i] = atoms->atomname[pos];
        if (newAtoms->haveType)
        {
            newAtoms->atomtype[i] = atoms->atomtype[pos];
        }
        else
        {
            sfree(newAtoms->atomtype);
        }
        if (newAtoms->haveBState)
        {
            newAtoms->atomtypeB[i] = atoms->atomtypeB[pos];
        }
        else
        {
            sfree(newAtoms->atomtypeB);
        }
        if (newAtoms->havePdbInfo)
        {
            newAtoms->pdbinfo[i] = atoms->pdbinfo[pos];
        }
    }
    newAtoms->nr = natoms;
    // Copy residue information unconditionally.
    for (int i = 0; i < atoms->nres; i++)
    {
        newAtoms->resinfo[i] = atoms->resinfo[i];
    }

    sfree(atoms);

    return newAtoms;
}

void
OutputSelector::processFrame(const int /*framenumber*/, t_trxframe *input)
{
    int              natoms = sel_.atomCount();

    input->natoms = natoms;

    rvec *xmem = nullptr;
    rvec *vmem = nullptr;
    rvec *fmem = nullptr;
    snew(xmem, natoms);
    // Free memory allocated before
    if (input->bV)
    {
        snew(vmem, natoms);
    }
    if (input->bF)
    {
        snew(fmem, natoms);
    }

    for (int i = 0; i < natoms; i++)
    {
        int pos = sel_.position(i).refId();
        copy_rvec(input->x[pos], xmem[i]);
        if (input->bV)
        {
            copy_rvec(input->v[pos], vmem[i]);
        }
        if (input->bF)
        {
            copy_rvec(input->f[pos], fmem[i]);
        }
    }
    sfree(input->x);
    input->x = xmem;
    if (input->bV)
    {
        sfree(input->v);
        input->v = vmem;
    }
    if (input->bF)
    {
        sfree(input->f);
        input->f = fmem;
    }

    if (input->bAtoms)
    {
        input->atoms = adjustAtomInformation(input->atoms, sel_);
    }

}

} // namespace gmx
