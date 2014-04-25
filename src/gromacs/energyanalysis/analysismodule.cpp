/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Implements classes in analysismodule.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "analysismodule.h"

#include <utility>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * EnergyAnalysisModule::Impl
 */

/*! \internal \brief
 * Private implementation class for EnergyAnalysisModule.
 *
 * \ingroup module_trajectoryanalysis
 */
class EnergyAnalysisModule::Impl
{
    public:
        //! Initializes analysis module data with given name and description.
        Impl(const char *name, const char *description)
            : name_(name), description_(description)
        {
        }

        //! Name of the module.
        std::string                     name_;
        //! Description of the module.
        std::string                     description_;
};

/********************************************************************
 * EnergyAnalysisModule
 */

EnergyAnalysisModule::EnergyAnalysisModule(const char *name,
                                           const char *description)
    : impl_(new Impl(name, description))
{
}


EnergyAnalysisModule::~EnergyAnalysisModule()
{
}


const char *EnergyAnalysisModule::name() const
{
    return impl_->name_.c_str();
}


const char *EnergyAnalysisModule::description() const
{
    return impl_->description_.c_str();
}


} // namespace gmx
