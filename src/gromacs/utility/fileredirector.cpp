/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements classes and functions from fileredirector.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/fileredirector.h"

#include <filesystem>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

IFileInputRedirector::~IFileInputRedirector() {}

IFileOutputRedirector::~IFileOutputRedirector() {}

namespace
{

/*! \internal
 * \brief
 * Implements the redirector returned by defaultFileInputRedirector().
 *
 * Does not redirect anything, but uses the file system as requested.
 *
 * \ingroup module_utility
 */
class DefaultInputRedirector : public IFileInputRedirector
{
public:
    bool fileExists(const std::filesystem::path& filename, const File::NotFoundHandler& onNotFound) const override
    {
        return File::exists(filename, onNotFound);
    }
};

/*! \internal
 * \brief
 * Implements the redirector returned by defaultFileOutputRedirector().
 *
 * Does not redirect anything, but instead opens the files exactly as
 * requested.
 *
 * \ingroup module_utility
 */
class DefaultOutputRedirector : public IFileOutputRedirector
{
public:
    TextOutputStream&       standardOutput() override { return TextOutputFile::standardOutput(); }
    TextOutputStreamPointer openTextOutputFile(const std::filesystem::path& filename) override
    {
        return TextOutputStreamPointer(new TextOutputFile(filename));
    }
};

} // namespace

//! \cond libapi
IFileInputRedirector& defaultFileInputRedirector()
{
    static DefaultInputRedirector instance;
    return instance;
}

IFileOutputRedirector& defaultFileOutputRedirector()
{
    static DefaultOutputRedirector instance;
    return instance;
}
//! \endcond

} // namespace gmx
