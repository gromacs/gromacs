/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares gmx::test::PosixMemstream
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_POSIXMEMSTREAM_H
#define GMX_TESTUTILS_POSIXMEMSTREAM_H

#include <cstdio>

#include <string>

namespace gmx
{

namespace test
{

/*! \brief Open an in-memory string stream that works like a FILE*
 * handle.
 *
 * This permits writing tests based on the contents of that string
 * after the stream has closed.
 *
 * When open_memstream() is not supported, provides a minimal
 * functional implementation that writes to stdout.
 * This lets us have test coverage of legacy code on POSIX-compliant
 * systems, which is much better than nothing.
 */
class PosixMemstream
{
public:
    PosixMemstream();
    ~PosixMemstream();
    //! Get the stream
    FILE* stream();
    //! Close the string stream (if supported)
    void closeStream();
    //! Return whether checking the buffer contents is supported
    static bool canCheckBufferContents();
    /*! \brief Close the stream and return a string
     *
     * \return A string containing the stream contents (if supported),
     * else an emptry string. */
    std::string toString();

private:
    char*  buffer_     = nullptr;
    size_t bufferSize_ = 0;
    FILE*  stream_     = nullptr;
    bool   isOpen_     = true;
};

} // namespace test
} // namespace gmx

#endif
