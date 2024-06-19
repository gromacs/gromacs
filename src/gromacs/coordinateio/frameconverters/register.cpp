/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
/*!\file
 * \internal
 * \brief
 * Implements gmx::ProcessFrameConversion.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "register.h"

#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

ProcessFrameConversion::ProcessFrameConversion() = default;

ProcessFrameConversion::~ProcessFrameConversion() = default;

void ProcessFrameConversion::addAndCheckGuarantee(const unsigned long flag)
{
    listOfGuarantees_ |= flag;
    if (((flag & convertFlag(FrameConverterFlags::NewSystemCenter)) != 0U)
        || ((flag & convertFlag(FrameConverterFlags::UnitCellIsCompact)) != 0U)
        || ((flag & convertFlag(FrameConverterFlags::UnitCellIsRectangular)) != 0U)
        || ((flag & convertFlag(FrameConverterFlags::UnitCellIsTriclinic)) != 0U))
    {
        listOfGuarantees_ &= ~convertFlag(FrameConverterFlags::AtomsInBox);
        listOfGuarantees_ &= ~convertFlag(FrameConverterFlags::MoleculeCOMInBox);
        listOfGuarantees_ &= ~convertFlag(FrameConverterFlags::ResidueCOMInBox);
        listOfGuarantees_ &= ~convertFlag(FrameConverterFlags::SystemIsCenteredInBox);
    }
}

void ProcessFrameConversion::addFrameConverter(FrameConverterPointer module)
{
    moduleChain_.emplace_back(std::move(module));
}

t_trxframe* ProcessFrameConversion::prepareAndTransformCoordinates(const t_trxframe* inputFrame)
{
    if (!frame_)
    {
        frame_ = std::make_unique<t_trxframe>();
        clear_trxframe(frame_.get(), true);
    }
    prepareNewCoordinates(inputFrame);
    convertFrame(frame_.get());
    GMX_ASSERT(inputFrame->natoms == frame_->natoms,
               "Frame conversion methods need to conserve the number of atoms");

    return frame_.get();
}


void ProcessFrameConversion::prepareNewCoordinates(const t_trxframe* inputFrame)
{
    *frame_ = *inputFrame;
    localX_.resize(inputFrame->natoms);
    frame_->x = as_rvec_array(localX_.data());
    if (inputFrame->bV)
    {
        localV_.resize(inputFrame->natoms);
        frame_->v = as_rvec_array(localV_.data());
    }
    if (inputFrame->bF)
    {
        localF_.resize(inputFrame->natoms);
        frame_->f = as_rvec_array(localF_.data());
    }
    for (int i = 0; i < frame_->natoms; i++)
    {
        copy_rvec(inputFrame->x[i], frame_->x[i]);
        if (frame_->bV)
        {
            copy_rvec(inputFrame->v[i], frame_->v[i]);
        }
        if (frame_->bF)
        {
            copy_rvec(inputFrame->f[i], frame_->f[i]);
        }
    }
}

void ProcessFrameConversion::convertFrame(t_trxframe* input)
{
    addAndCheckGuarantee(guarantee());
    for (auto& method : moduleChain_)
    {
        method.module_->convertFrame(input);
        addAndCheckGuarantee(method.module_->guarantee());
    }
}
} // namespace gmx
