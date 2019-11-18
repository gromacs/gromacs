/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2019, by the GROMACS development team, led by
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
 * Declares gmx::IOptionSectionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ISECTIONSTORAGE_H
#define GMX_OPTIONS_ISECTIONSTORAGE_H

namespace gmx
{

/*! \internal
 * \brief
 * Provides behavior specific to a certain option section type.
 *
 * \ingroup module_options
 */
class IOptionSectionStorage
{
public:
    virtual ~IOptionSectionStorage();

    /*! \brief
     * Called once before the first call to startSection().
     *
     * This is called once all options have been added to the section.
     * The current implementation does not call this if startSection() is
     * never called.
     */
    virtual void initStorage() = 0;
    /*! \brief
     * Called when option assignment enters this section.
     */
    virtual void startSection() = 0;
    /*! \brief
     * Called when option assignment leaves this section.
     */
    virtual void finishSection() = 0;
};

} // namespace gmx

#endif
