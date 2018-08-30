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
 * Implements methods from outputmanager.h
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "outputmanager.h"

#include <algorithm>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

void
OutputManager::addFlagModule(CoordinateOutputPointer module)
{
    // first check if the module dependencies are satisfied
    // by the OutputManager object that has been constructed.
    if (module->checkModuleFlag(moduleFlags_))
    {
        flagModules_.emplace_back(module);
    }
    else
    {
        GMX_THROW(InconsistentInputError("Module does not support this output"));
    }
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
    moduleFlags_ |= efCoordinateOutput;
    switch (filetype_)
    {
        case (efTNG):
            moduleFlags_ |= (efForceOutput | efVelocityOutput | efAtomOutput | efCustomPrecision);
            break;
        case (efPDB):
            moduleFlags_ |= (efConnectionOutput | efAtomOutput);
            break;
        case (efGRO):
            moduleFlags_ |= (efAtomOutput | efVelocityOutput);
            break;
        case (efTRR):
            moduleFlags_ |= (efForceOutput | efVelocityOutput | efCustomPrecision);
            break;
        case (efXTC):
            moduleFlags_ |= (efCustomPrecision);
            break;
        case (efG96):
            break;
        default:
            GMX_THROW(InvalidInputError("Invalid file type"));
    }
}

void
OutputManager::closeFile()
{
    if (outputFile_ != nullptr)
    {
        close_trx(outputFile_);
    }
    outputFile_ = nullptr;
}

t_trxstatus *
OutputManager::trjOpenTng()
{
    ArrayRef<const int> index     = sel_.atomIndices();
    int                 natoms    = sel_.atomCount();
    const char         *indexname = sel_.name();
    gmx_mtop_t          mtop;
    gmx_mtop_t         *pMtop = &mtop;
    init_mtop(pMtop);
    if (haveMtop())
    {
        pMtop      = const_cast<gmx_mtop_t*>(getMtop());
    }
    else
    {
        GMX_THROW(InvalidInputError("Need topology information to write TNG files"));
    }


    std::vector<int>    localindex;
    for (const int *ai = index.begin(); ai != index.end(); ++ai)
    {
        localindex.emplace_back(*ai);
    }

    t_trxstatus *out =
        trjtools_gmx_prepare_tng_writing(name_.c_str(),
                                         filemode_[0],
                                         nullptr, //infile_, //how to get the input file here?
                                         nullptr,
                                         natoms,
                                         pMtop,
                                         localindex.data(),
                                         indexname);

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
    if (!haveMtop())
    {
        GMX_THROW(InvalidInputError("Need trajectory frame with atom names for pdb file writing"));
    }
    return open_trx(name_.c_str(), filemode_.c_str());
}

t_trxstatus *
OutputManager::trjOpenGro()
{
    if (!haveMtop())
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

void
OutputManager::initOutput()
{
    if (outputFile_ == nullptr)
    {
        switch (filetype_)
        {
            case (efTNG):
                outputFile_          = trjOpenTng();
                break;
            case (efPDB):
                outputFile_          = trjOpenPdb();
                break;
            case (efGRO):
                outputFile_          = trjOpenGro();
                break;
            case (efTRR):
            case (efXTC):
            case (efG96):
                outputFile_ = trjOpenTrr();
                break;
            default:
                GMX_THROW(InvalidInputError("Invalid file type"));
        }
    }
}

void
OutputManager::prepareFrame(const int framenumber, const t_trxframe &input)
{
    if (outputFile_ == nullptr)
    {
        initOutput();
    }
    if (!flagModules_.empty())
    {
        t_trxframe local = input;
        for (const auto &method : flagModules_)
        {
            method.module->processFrame(framenumber, &local);
        }
        writeFrame(local);
    }
    else
    {
        writeFrame(input);
    }
}

void
OutputManager::writeFrame(const t_trxframe write)
{
    t_trxframe        local  = write;
    if (outputFile_ == nullptr)
    {
        initOutput();
    }
    write_trxframe(outputFile_, &local, nullptr);
}

} // namespace gmx
