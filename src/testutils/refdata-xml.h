/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Declares functions for reference data XML persistence.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_XML_H
#define GMX_TESTUTILS_REFDATA_XML_H

#include <string>

#include "testutils/refdata-impl.h"

namespace gmx
{
namespace test
{

class ReferenceDataEntry;

//! \cond internal
/*! \internal
 * \brief
 * Loads reference data from an XML file.
 *
 * \param[in] path  Path to the file from which the data is loaded.
 * \returns   Root entry for the reference data parsed from the file.
 * \throws    TestException if there is a problem reading the file.
 *
 * \ingroup module_testutils
 */
ReferenceDataEntry::EntryPointer
readReferenceDataFile(const std::string &path);
/*! \internal
 * \brief
 * Saves reference data to an XML file.
 *
 * \param[in] path  Path to the file to which the data is saved.
 * \param[in] root  Root entry for the reference data to write.
 * \throws    TestException if there is a problem writing the file.
 *
 * \ingroup module_testutils
 */
void writeReferenceDataFile(const std::string        &path,
                            const ReferenceDataEntry &root);
//! \endcond

} // namespace test
} // namespace gmx

#endif
