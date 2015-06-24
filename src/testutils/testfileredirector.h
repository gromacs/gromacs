/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares generic mock implementations for interfaces in fileredirector.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTFILEREDIRECTOR_H
#define GMX_TESTUTILS_TESTFILEREDIRECTOR_H

#include <set>
#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/fileredirector.h"

namespace gmx
{
namespace test
{

/*! \libinternal \brief
 * In-memory implementation for FileInputRedirectorInterface for tests.
 *
 * By default, this implementation will return `false` for all file existence
 * checks.  To return `true` for a specific path, use addExistingFile().
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class TestFileInputRedirector : public FileInputRedirectorInterface
{
    public:
        TestFileInputRedirector();
        virtual ~TestFileInputRedirector();

        /*! \brief
         * Marks the provided path as an existing file.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * Further checks for existence of the given path will return `true`.
         */
        void addExistingFile(const char *filename);

        // From FileInputRedirectorInterface
        virtual bool fileExists(const char *filename) const;

    private:
        std::set<std::string> existingFiles_;

        GMX_DISALLOW_COPY_AND_ASSIGN(TestFileInputRedirector);
};

} // namespace test
} // namespace gmx

#endif
