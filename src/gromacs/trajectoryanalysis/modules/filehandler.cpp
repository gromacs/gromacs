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
/*! \internal \file
 * \brief
 * Helper classes for file handling of trajectory files.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */

#include <algorithm>

#include "filehandler.h"
#include "writesettings.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/unique_cptr.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"

namespace gmx
{

void
Filehandler::initFileOptions(IOptionsContainer *options)
{
}



void
Filehandler::modifyFrame(const t_trxframe *oldFrame)
{
    const Selection *sel    = &sel_;
    int              natoms = sel->atomCount();

    newFrame_ = *oldFrame;

    newFrame_.time   = writeDoubles_.time;
    newFrame_.bV     = (oldFrame->bV && writeBool_.bV);
    newFrame_.bF     = (oldFrame->bF && writeBool_.bF);
    newFrame_.natoms = natoms;
    newFrame_.bPrec  = (oldFrame->bPrec && writeBool_.bP);
    newFrame_.prec   = writeDoubles_.prec;
    newFrame_.atoms  = &atoms_;
    newFrame_.bAtoms = writeBool_.bA;

    rvec *xmem = nullptr;
    rvec *vmem = nullptr;
    rvec *fmem = nullptr;
    snew(xmem, natoms);
    if (newFrame_.bV)
    {
        snew(vmem, natoms);
    }
    if (newFrame_.bF)
    {
        snew(fmem, natoms);
    }
    newFrame_.x = xmem;
    newFrame_.v = vmem;
    newFrame_.f = fmem;


    for (int i = 0; i < natoms; i++)
    {
        int pos = sel->position(i).refId();
        copy_rvec(oldFrame->x[pos], newFrame_.x[i]);
        if (newFrame_.bV)
        {
            copy_rvec(oldFrame->v[pos], newFrame_.v[i]);
        }
        if (newFrame_.bF)
        {
            copy_rvec(oldFrame->f[pos], newFrame_.f[i]);
        }
    }
}

void Filehandler::setBox(const matrix box)
{
    copy_mat(box, newFrame_.box);
}

void
Filehandler::closeFile()
{
    if (trr_ != nullptr)
    {
        close_trx(trr_);
    }
    trr_ = nullptr;
}

void
Filehandler::setLegacyInformation(t_atoms *local)
{
    t_atoms inputAtoms = gmx_mtop_global_atoms(mtop_);
    init_t_atoms(local, inputAtoms.nr, inputAtoms.havePdbInfo);
    sfree(local->resinfo);
    local->resinfo = inputAtoms.resinfo;
    int natoms = sel_.atomCount();
    for (int i = 0; (i < natoms); i++)
    {
        local->atomname[i] = inputAtoms.atomname[sel_.position(i).refId()];
        local->atom[i]     = inputAtoms.atom[sel_.position(i).refId()];
        if (inputAtoms.havePdbInfo)
        {
            local->pdbinfo[i] = inputAtoms.pdbinfo[sel_.position(i).refId()];
        }
        local->nres = std::max(local->nres, local->atom[i].resind+1);
    }
    local->nr     = natoms;
    writeBool_.bA = true;
}

t_trxstatus *
Filehandler::trjOpenTng(const char *filename, const char *mode) const
{
    ArrayRef<const int> index     = sel_.atomIndices();
    int                 natoms    = sel_.atomCount();
    const char         *indexname = sel_.name();

    int                 size = index.size();
    int                *localindex;
    snew(localindex, size);
    int                 runner = 0;

    for (const int *ai = index.begin(); ai != index.end(); ++ai)
    {
        localindex[runner++] = *ai;
    }

    gmx_mtop_t mtop = (*const_cast<gmx_mtop_t*>(mtop_));

    t_trxstatus *out =
    trjtools_gmx_prepare_tng_writing(filename,
                                     mode[0],
                                     nullptr, //infile_, //how to get the input file here?
                                     nullptr,
                                     natoms,
                                     &mtop,
                                     localindex,
                                     indexname);

    sfree(localindex);
    return out;
}

t_trxstatus *
Filehandler::trjOpenTrr(const char *filename, const char *mode) const
{
    return open_trx(filename, mode);
}

void
Filehandler::setConnections()
{
    gmx_mtop_t *mtop = const_cast<gmx_mtop_t*>(mtop_);
    t_topology  top  = gmx_mtop_t_to_t_topology(mtop, false);
    connections_ = gmx_conect_generate(&top);
}


void
Filehandler::initOutput(TrajectoryWriteSettings *settings)
{
    settings_     = settings;
    name_         = settings_->getName();
    sel_          = *(settings_->getSel());
    mtop_         = settings_->getMtop();
    writeBool_    = settings_->getWriteBool();
    writeDoubles_ = settings_->getWriteDoubles();
    if (!name_.empty())
    {
        filetype_ = fn2ftp(name_.c_str());
        switch (filetype_)
        {
            case (efTNG):
                trr_ = trjOpenTng(name_.c_str(), filemode_);
                break;
            case (efPDB):
            case (efGRO):
		writeBool_.bGC = true;
            case (efTRR):
            case (efXTC):
            case (efG96):
                trr_ = trjOpenTrr(name_.c_str(), filemode_);
                setLegacyInformation(&atoms_);
                break;
            default:
                gmx_incons("Invalid file type");
        }
    }
    if (writeBool_.bGC)
    {
        setConnections();
    }
}

void
Filehandler::writeFrame() const
{
    const t_trxframe *plocal = &newFrame_;
    t_trxframe        local  = (*const_cast<t_trxframe*>(plocal));

    switch (filetype_)
    {
        case (efTNG):
        case (efTRR):
        case (efXTC):
        case (efPDB):
        case (efGRO):
        case (efG96):
            write_trxframe(trr_, &local, connections_);
            break;
        default:
            gmx_incons("Illegal output file format");
    }
}

t_trxframe *
Filehandler::getFrameForModification()
{
	return &newFrame_;
}

} // namespace gmx

