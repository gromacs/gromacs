/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::test::StringTestBase.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "stringtest.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/fileredirector.h"

#include "testutils/refdata.h"
#include "testutils/testexceptions.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

namespace
{
//! Stores the -stdout flag value to print out values instead of checking them.
bool g_bWriteToStdOut = false;

/*! \brief
 * Helper for checking a block of text, e.g., implementing the `-stdout`
 * option.
 *
 * \ingroup module_testutils
 */
void checkTextImpl(TestReferenceChecker *checker, const std::string &text,
                   const char *id)
{
    if (g_bWriteToStdOut)
    {
        printf("%s:\n", id);
        printf("%s[END]\n", text.c_str());
    }
    else
    {
        checker->checkStringBlock(text, id);
    }
}

}

// TODO: Only add this option to those test binaries that actually need it
// (depending on the linker, it may or may not appear right now),
// or replace by a generic mechanism in TestReferenceData.
//! \cond
GMX_TEST_OPTIONS(StringTestOptions, options)
{
    options->addOption(
            BooleanOption("stdout")
                .store(&g_bWriteToStdOut)
                .description("Print the test string to stdout instead of checking against reference data"));
}
//! \endcond

/********************************************************************
 * TestFileOutputRedirector
 */

/*! \internal
 * \brief
 * Implementation of FileOutputRedirectorInterface for tests.
 *
 * This class redirects all output files to temporary files managed by a
 * TestFileManager, and supports checking the contents of these files using the
 * reference data framework.
 *
 * \ingroup module_testutils
 */
class TestFileOutputRedirector : public FileOutputRedirectorInterface
{
    public:
        //! Initializes the redirector with the given file manager.
        explicit TestFileOutputRedirector(TestFileManager *fileManager)
            : fileManager_(*fileManager)
        {
        }

        virtual File &standardOutput()
        {
            if (!stdoutFile_)
            {
                const std::string path = fileManager_.getTemporaryFilePath("stdout.txt");
                stdoutFile_.reset(new File(path, "w"));
                fileList_.push_back(FileListEntry("<stdout>", path));
            }
            return *stdoutFile_;
        }
        virtual FileInitializer openFileForWriting(const char *filename)
        {
            std::string       suffix = filename;
            std::replace(suffix.begin(), suffix.end(), '/', '_');
            const std::string path = fileManager_.getTemporaryFilePath(suffix);
            fileList_.push_back(FileListEntry(filename, path));
            return FileInitializer(fileList_.back().second.c_str(), "w");
        }

        /*! \brief
         * Checks the contents of all redirected files.
         */
        void checkRedirectedFiles(TestReferenceChecker *checker)
        {
            if (stdoutFile_)
            {
                stdoutFile_->close();
                stdoutFile_.reset();
            }
            std::vector<FileListEntry>::const_iterator i;
            for (i = fileList_.begin(); i != fileList_.end(); ++i)
            {
                const std::string text = File::readToString(i->second);
                checkTextImpl(checker, text, i->first.c_str());
            }
        }

    private:
        typedef std::pair<std::string, std::string> FileListEntry;

        TestFileManager            &fileManager_;
        boost::scoped_ptr<File>     stdoutFile_;
        std::vector<FileListEntry>  fileList_;
};

/********************************************************************
 * StringTestBase::Impl
 */

class StringTestBase::Impl
{
    public:
        TestReferenceData                           data_;
        boost::scoped_ptr<TestReferenceChecker>     checker_;
        boost::scoped_ptr<TestFileOutputRedirector> redirector_;
};

/********************************************************************
 * StringTestBase
 */

StringTestBase::StringTestBase()
    : impl_(new Impl)
{
}

StringTestBase::~StringTestBase()
{
}

FileOutputRedirectorInterface &
StringTestBase::initOutputRedirector(TestFileManager *fileManager)
{
    if (impl_->redirector_)
    {
        GMX_THROW(TestException("initOutputRedirector() called more than once"));
    }
    impl_->redirector_.reset(new TestFileOutputRedirector(fileManager));
    return *impl_->redirector_;
}

TestReferenceChecker &
StringTestBase::checker()
{
    if (!impl_->checker_)
    {
        impl_->checker_.reset(new TestReferenceChecker(impl_->data_.rootChecker()));
    }
    return *impl_->checker_;
}

void
StringTestBase::checkText(const std::string &text, const char *id)
{
    checkTextImpl(&checker(), text, id);
}

void
StringTestBase::checkFileContents(const std::string &filename, const char *id)
{
    const std::string text = File::readToString(filename);
    checkText(text, id);
}

void
StringTestBase::checkRedirectedOutputFiles()
{
    if (!impl_->redirector_)
    {
        GMX_THROW(TestException("initOutputRedirector() not called"));
    }
    impl_->redirector_->checkRedirectedFiles(&checker());
}

} // namespace test
} // namespace gmx
