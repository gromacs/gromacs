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
 * Implements classes in filehandler.h.
 *
 * \ingroup module_analysisdata
 * \author
 */
#include "gmxpre.h"

#include "filehandler.h"

#include <cstdio>
#include <cstring>

#include <string>
#include <vector>

#include "gromacs/analysisdata/modules/settings.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"


namespace gmx
{
Filehandler::Filehandler(TrajectoryDataWriteSettings *settings) : settings_(settings)
{
    strcpy(filemode_, "w");
}
Filehandler::Filehandler()
{
    strcpy(filemode_, "w");
}
Filehandler::~Filehandler()
{
    closeFile();
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
Filehandler::setLegacyInformation()
{
    const Selection  *sel  = sel_;
    const gmx_mtop_t *mtop =  settings_->getMtop();
    // get pointer to storage of atoms data in settings
    t_atoms          *local      = settings_->getAtoms();
    t_atoms           inputAtoms = gmx_mtop_global_atoms(mtop);
    init_t_atoms(local, inputAtoms.nr, inputAtoms.havePdbInfo);
    sfree(local->resinfo);
    int natoms = sel->atomCount();
    snew(local->resinfo, natoms);
    for (int i = 0; (i < natoms); i++)
    {
        int pos = sel->position(i).refId();
        local->atomname[i] = inputAtoms.atomname[pos];
        local->atom[i]     = inputAtoms.atom[pos];
        int localResind = local->atom[i].resind;
        int inputResind = inputAtoms.atom[pos].resind;
        local->resinfo[localResind] = inputAtoms.resinfo[inputResind];
        if (inputAtoms.havePdbInfo)
        {
            local->pdbinfo[i] = inputAtoms.pdbinfo[pos];
        }
        local->nres = std::max(local->nres, local->atom[i].resind+1);
    }
    local->nr = sel->atomCount();
    settings_->setbAtoms(true);
    settings_->setAtoms(local);
    done_atom(&inputAtoms);
}

void
Filehandler::setSettings(TrajectoryDataWriteSettings *settings)
{
    settings_ = settings;
    sel_      = settings_->getInputSel();
}

t_trxstatus *
Filehandler::trjOpenTng(const char *filename, const char *mode) const
{
    ArrayRef<const int> index     = sel_->atomIndices();
    int                 natoms    = sel_->atomCount();
    const char         *indexname = sel_->name();

    int                 size = index.size();
    int                *localindex;
    snew(localindex, size);
    int                 runner = 0;

    for (const int *ai = index.begin(); ai != index.end(); ++ai)
    {
        localindex[runner++] = *ai;
    }

    gmx_mtop_t mtop = (*const_cast<gmx_mtop_t*>(settings_->getMtop()));

//    t_trxstatus &localInput = impl_->infile_;
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
    settings_->setConnections();
}


void
Filehandler::openFile()
{
    std::string name = settings_->getName();
    if (!name.empty())
    {
        int filetype = settings_->getFiletype();

        switch (filetype)
        {
            case (efTNG):
                trr_ = trjOpenTng(name.c_str(), filemode_);
                break;
            case (efPDB):
            case (efGRO):
                settings_->setbGC(true);
            case (efTRR):
            case (efXTC):
            case (efG96):
                trr_ = trjOpenTrr(name.c_str(), filemode_);
                setLegacyInformation();
                break;
            default:
                gmx_incons("Invalid file type");
                // handle this error
        }
    }
    if (settings_->getbGC())
    {
        setConnections();
    }
}

/*! \cond libapi */
bool
Filehandler::isFileOpen() const
{
    return trr_ != nullptr;
}


void
Filehandler::writeValue(t_trxframe &coord) const
{
    t_trxframe *plocal   = &coord;
    int         filetype = settings_->getFiletype();
    switch (filetype)
    {
        case (efTNG):
        case (efTRR):
        case (efXTC):
        case (efPDB):
        case (efGRO):
        case (efG96):
            write_trxframe(trr_, plocal, settings_->getConnections());
            break;
        default:
            gmx_incons("Illegal output file format");
    }
    done_coord(plocal);
}

} // namespace gmx
