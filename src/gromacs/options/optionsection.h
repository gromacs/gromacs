/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::OptionSection and gmx::OptionSectionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSECTION_H
#define GMX_OPTIONS_OPTIONSECTION_H

#include "gromacs/utility/classhelpers.h"

#include "abstractsection.h"

namespace gmx
{

class OptionSectionHandle;

/*! \brief
 * Declares a simple option section.
 *
 * This class declares a simple section that only provides structure for
 * grouping the options, but does not otherwise influence the behavior of the
 * contained options.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class OptionSection : public AbstractOptionSection
{
public:
    //! AbstractOptionSectionHandle corresponding to this option type.
    typedef OptionSectionHandle HandleType;

    //! Creates a section with the given name.
    explicit OptionSection(const char* name) : AbstractOptionSection(name) {}

private:
    std::unique_ptr<IOptionSectionStorage> createStorage() const override;
};

/*! \brief
 * Allows adding options to an OptionSection.
 *
 * An instance of this class is returned from
 * IOptionsContainerWithSections::addSection(), and supports adding options and
 * subsections to a section created with OptionSection.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class OptionSectionHandle : public AbstractOptionSectionHandle
{
public:
    //! Wraps a given section storage object.
    explicit OptionSectionHandle(internal::OptionSectionImpl* section) :
        AbstractOptionSectionHandle(section)
    {
    }
};

class OptionSectionInfo : public AbstractOptionSectionInfo
{
public:
    //! Wraps a given section storage object.
    explicit OptionSectionInfo(internal::OptionSectionImpl* section) :
        AbstractOptionSectionInfo(section)
    {
    }
};

} // namespace gmx

#endif
