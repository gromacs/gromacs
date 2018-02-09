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
 * Helper classe for handling output of trajectory data to files.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "outputmanager.h"

#include <algorithm>

#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

void
OutputManager::initFileOptions(IOptionsContainer *options)
{
    options->addOption(FileNameOption("o").filetype(eftTrajectory).outputFile()
                           .store(&name_).defaultBasename("trajout")
                           .required()
                           .description("Output trajectory after conversion"));
}

void
OutputManager::checkOptions(const TopologyInformation *top)
{
    setFiletype();
    if (metadata_->hasFlag(efRequireAtoms) && !metadata_->hasFlag(efRequireMtop))
    {
        if (!top->topology())
        {
            GMX_THROW(InvalidInputError("Need topology available to write atom information to text files"));
        }
        metadata_->setAtoms(top->topology()->atoms);
        // check if atoms is accesible from frame or from mtop
    }
    if (metadata_->hasFlag(efRequireMtop))
    {
        if (!top->mtop())
        {
            GMX_THROW(InvalidInputError("Full topology needs to be accesible"));
        }
        metadata_->setMtop(top->mtop());
        // check if mtop is accesible
    }
    if (metadata_->hasFlag(efRequireConnectionCheck))
    {
        if (!metadata_->haveAtoms() && !metadata_->haveMtop())
        {
            GMX_THROW(InvalidInputError("Need to have some topology information to write connections"));
        }
        // check if we can write connections
    }

}

void
OutputManager::setFiletype()
{
    if (!name_.empty())
    {
        filetype_ = fn2ftp(name_.c_str());
        switch (filetype_)
        {
            case (efTNG):
                // effectively needs mtop but not really
                metadata_->setFlag(TrajectoryOutputMetadata::efRequireMtop, true);
                break;
            case (efPDB):
                metadata_->setFlag(TrajectoryOutputMetadata::efRequireAtoms, true);
                metadata_->setFlag(TrajectoryOutputMetadata::efRequireConnectionCheck, true);
                break;
            // needs atoms from topology, so have to add mtop
            case (efGRO):
                metadata_->setFlag(TrajectoryOutputMetadata::efRequireAtoms, true);
                // needs atoms from topology, so have to add mtop
                break;
            case (efTRR):
            case (efXTC):
            case (efG96):
                break;
            default:
                gmx_incons("Invalid file type");
        }
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
OutputManager::trjOpenTng() const
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

    t_trxstatus *out =
        trjtools_gmx_prepare_tng_writing(name_.c_str(),
                                         filemode_[0],
                                         nullptr, //infile_, //how to get the input file here?
                                         nullptr,
                                         natoms,
                                         nullptr,
                                         localindex,
                                         indexname);

    sfree(localindex);
    return out;
}

t_trxstatus *
OutputManager::trjOpenTrr() const
{
    return open_trx(name_.c_str(), filemode_);
}

t_trxstatus *
OutputManager::trjOpenPdb() const
{
    return open_trx(name_.c_str(), filemode_);
}

t_trxstatus *
OutputManager::trjOpenGro() const
{
    return open_trx(name_.c_str(), filemode_);
}

void
OutputManager::initOutput()
{
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
            gmx_incons("Invalid file type");
    }
}

void
OutputManager::writeFrame(const t_trxframe *write) const
{
    t_trxframe        local  = (*const_cast<t_trxframe*>(write));

    write_trxframe(trr_, &local, nullptr);
}





} // namespace gmx
