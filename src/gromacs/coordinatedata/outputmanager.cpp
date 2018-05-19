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
 * Helper classe for handling output of trajectory data to files.
 *
 * \author
 * \ingroup module_coordinatedata
 */


#include "gmxpre.h"

#include "outputmanager.h"

#include <algorithm>

#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

void
OutputManager::checkOptions()
{
}

void
OutputManager::setFiletype()
{
    if (!name_.empty())
    {
        filetype_ = fn2ftp(name_.c_str());
    }
    else
    {
        GMX_THROW(InvalidInputError("Can not open file with an empty name"));
    }
}

void
OutputManager::closeFile()
{
    if (trr_ != nullptr)
    {
        close_trx(trr_);
    }
    trr_ = nullptr;
}

t_trxstatus *
OutputManager::trjOpenTng()
{
    ArrayRef<const int> index     = sel_->atomIndices();
    int                 natoms    = sel_->atomCount();
    const char         *indexname = sel_->name();
    if (!haveMtop() && !haveAtoms())
    {
        GMX_THROW(InvalidInputError("Need either atoms or topology information to write TNG file"));
    }
    gmx_mtop_t  mtop;
    gmx_mtop_t *pMtop = &mtop;
    init_mtop(pMtop);
    if (haveMtop())
    {
        pMtop      = const_cast<gmx_mtop_t*>(mtop_);
    }
    else if (haveAtoms())
    {

        convertAtomsToMtop(nullptr, nullptr, &atoms_, pMtop);
    }


    int                 size = index.size();
    int                *localindex;
    snew(localindex, size);
    int                 runner = 0;

    for (const int *ai = index.begin(); ai != index.end(); ++ai)
    {
        localindex[runner++] = *ai;
    }

    t_trxstatus *out =
        trjtools_gmx_prepare_tng_writing(name_.c_str(),
                                         filemode_[0],
                                         nullptr, //infile_, //how to get the input file here?
                                         nullptr,
                                         natoms,
                                         pMtop,
                                         localindex,
                                         indexname);

    sfree(localindex);
    return out;
}

t_trxstatus *
OutputManager::trjOpenTrr()
{
    return open_trx(name_.c_str(), filemode_.c_str());
}

t_trxstatus *
OutputManager::trjOpenPdb()
{
    if (!haveAtoms())
    {
        GMX_THROW(InvalidInputError("Need trajectory frame with atom names for pdb file writing"));
    }
    return open_trx(name_.c_str(), filemode_.c_str());
}

t_trxstatus *
OutputManager::trjOpenGro()
{
    if (!haveAtoms())
    {
        GMX_THROW(InvalidInputError("Need trajectory frame with atom names for pdb file writing"));
    }
    return open_trx(name_.c_str(), filemode_.c_str());
}

void
OutputManager::clearCoordinateFrame(t_trxframe *frame) const
{
    sfree(frame->x);
    sfree(frame->v);
    sfree(frame->f);
}

bool
OutputManager::haveAtoms() const
{
    return atoms_.nr > 0 ? true : false;
}

bool
OutputManager::haveMtop() const
{
    return mtop_ != nullptr ? true : false;
}

void
OutputManager::initOutput()
{
    if (trr_ == nullptr)
    {
        setFiletype();
        switch (filetype_)
        {
            case (efTNG):
                trr_ = trjOpenTng();
                break;
            case (efPDB):
                trr_ = trjOpenPdb();
                break;
            case (efGRO):
                trr_ = trjOpenGro();
                break;
            case (efTRR):
            case (efXTC):
            case (efG96):
                trr_ = trjOpenTrr();
                break;
            default:
                GMX_THROW(InvalidInputError("Invalid file type"));
        }
    }
}

void
OutputManager::writeFrame(const t_trxframe write) const
{
    t_trxframe        local  = (*const_cast<t_trxframe*>(&write));

    write_trxframe(trr_, &local, nullptr);
    clearCoordinateFrame(&local);
}





} // namespace gmx
