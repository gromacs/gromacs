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
 * Implements classes from testfileredirector.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testfileredirector.h"

#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textstream.h"

#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{
class TestReferenceChecker;

/********************************************************************
 * TestFileInputRedirector
 */

TestFileInputRedirector::TestFileInputRedirector() {}

TestFileInputRedirector::~TestFileInputRedirector() {}

void TestFileInputRedirector::addExistingFile(const std::filesystem::path& filename)
{
    existingFiles_.insert(filename);
}

bool TestFileInputRedirector::fileExists(const std::filesystem::path& filename,
                                         const File::NotFoundHandler& onNotFound) const
{
    if (existingFiles_.count(filename) == 0)
    {
        File::NotFoundInfo info(filename, "File not present in test", nullptr, false, 0);
        onNotFound(info);
        return false;
    }
    return true;
}

/********************************************************************
 * TestFileOutputRedirector::Impl
 */

class TestFileOutputRedirector::Impl
{
public:
    typedef std::shared_ptr<StringOutputStream>         StringStreamPointer;
    typedef std::pair<std::string, StringStreamPointer> FileListEntry;

    StringStreamPointer        stdoutStream_;
    std::vector<FileListEntry> fileList_;
};

/********************************************************************
 * TestFileOutputRedirector
 */

TestFileOutputRedirector::TestFileOutputRedirector() : impl_(new Impl) {}

TestFileOutputRedirector::~TestFileOutputRedirector() {}

TextOutputStream& TestFileOutputRedirector::standardOutput()
{
    if (!impl_->stdoutStream_)
    {
        impl_->stdoutStream_.reset(new StringOutputStream);
        impl_->fileList_.emplace_back("<stdout>", impl_->stdoutStream_);
    }
    return *impl_->stdoutStream_;
}

TextOutputStreamPointer TestFileOutputRedirector::openTextOutputFile(const std::filesystem::path& filename)
{
    Impl::StringStreamPointer stream(new StringOutputStream);
    impl_->fileList_.emplace_back(filename.string(), stream);
    return stream;
}

void TestFileOutputRedirector::checkRedirectedFiles(TestReferenceChecker* checker)
{
    std::vector<Impl::FileListEntry>::const_iterator i;
    for (i = impl_->fileList_.begin(); i != impl_->fileList_.end(); ++i)
    {
        StringTestBase::checkText(checker, i->second->toString(), i->first.c_str());
    }
}

} // namespace test
} // namespace gmx
