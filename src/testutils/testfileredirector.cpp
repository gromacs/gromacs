/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Implements classes from testfileredirector.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testfileredirector.h"

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/stringstream.h"

#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * TestFileInputRedirector
 */

TestFileInputRedirector::TestFileInputRedirector()
{
}

TestFileInputRedirector::~TestFileInputRedirector()
{
}

void TestFileInputRedirector::addExistingFile(const char *filename)
{
    existingFiles_.insert(filename);
}

bool TestFileInputRedirector::fileExists(const char            *filename,
                                         File::NotFoundHandler  onNotFound) const
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
        typedef std::shared_ptr<StringOutputStream> StringStreamPointer;
        typedef std::pair<std::string, StringStreamPointer> FileListEntry;

        StringStreamPointer         stdoutStream_;
        std::vector<FileListEntry>  fileList_;
};

/********************************************************************
 * TestFileOutputRedirector
 */

TestFileOutputRedirector::TestFileOutputRedirector()
    : impl_(new Impl)
{
}

TestFileOutputRedirector::~TestFileOutputRedirector()
{
}

TextOutputStream &TestFileOutputRedirector::standardOutput()
{
    if (!impl_->stdoutStream_)
    {
        impl_->stdoutStream_.reset(new StringOutputStream);
        impl_->fileList_.emplace_back("<stdout>", impl_->stdoutStream_);
    }
    return *impl_->stdoutStream_;
}

TextOutputStreamPointer
TestFileOutputRedirector::openTextOutputFile(const char *filename)
{
    Impl::StringStreamPointer stream(new StringOutputStream);
    impl_->fileList_.emplace_back(filename, stream);
    return stream;
}

void TestFileOutputRedirector::checkRedirectedFiles(TestReferenceChecker *checker)
{
    std::vector<Impl::FileListEntry>::const_iterator i;
    for (i = impl_->fileList_.begin(); i != impl_->fileList_.end(); ++i)
    {
        StringTestBase::checkText(checker, i->second->toString(), i->first.c_str());
    }
}

} // namespace test
} // namespace gmx
