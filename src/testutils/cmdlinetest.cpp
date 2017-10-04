/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements classes from cmdlinetest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "cmdlinetest.h"

#include <cstdlib>
#include <cstring>

#include <memory>
#include <new>
#include <sstream>
#include <vector>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * CommandLine::Impl
 */

class CommandLine::Impl
{
    public:
        Impl(const char *const cmdline[], size_t count);
        ~Impl();

        std::vector<char *>     args_;
        std::vector<char *>     argv_;
        int                     argc_;
};

CommandLine::Impl::Impl(const char *const cmdline[], size_t count)
{
    args_.reserve(count);
    argv_.reserve(count + 1);
    argc_ = static_cast<int>(count);
    for (size_t i = 0; i < count; ++i)
    {
        char *arg = strdup(cmdline[i]);
        if (arg == nullptr)
        {
            throw std::bad_alloc();
        }
        args_.push_back(arg);
        argv_.push_back(arg);
    }
    argv_.push_back(nullptr);
}

CommandLine::Impl::~Impl()
{
    for (size_t i = 0; i < args_.size(); ++i)
    {
        std::free(args_[i]);
    }
}

/********************************************************************
 * CommandLine
 */

CommandLine::CommandLine()
    : impl_(new Impl(nullptr, 0))
{
}

CommandLine::CommandLine(const ArrayRef<const char *const> &cmdline)
    : impl_(new Impl(cmdline.data(), cmdline.size()))
{
}

CommandLine::CommandLine(const CommandLine &other)
    : impl_(new Impl(other.argv(), other.argc()))
{
}

CommandLine::~CommandLine()
{
}

void CommandLine::initFromArray(const ArrayRef<const char *const> &cmdline)
{
    impl_.reset(new Impl(cmdline.data(), cmdline.size()));
}

void CommandLine::append(const char *arg)
{
    GMX_RELEASE_ASSERT(impl_->argc_ == static_cast<int>(impl_->args_.size()),
                       "Command-line has been modified externally");
    size_t newSize = impl_->args_.size() + 1;
    impl_->args_.reserve(newSize);
    impl_->argv_.reserve(newSize + 1);
    char *newArg = strdup(arg);
    if (newArg == nullptr)
    {
        throw std::bad_alloc();
    }
    impl_->args_.push_back(newArg);
    impl_->argv_.pop_back(); // Remove the trailing NULL.
    impl_->argv_.push_back(newArg);
    impl_->argv_.push_back(nullptr);
    impl_->argc_ = static_cast<int>(newSize);
}

void CommandLine::addOption(const char *name)
{
    append(name);
}

void CommandLine::addOption(const char *name, const char *value)
{
    append(name);
    append(value);
}

void CommandLine::addOption(const char *name, const std::string &value)
{
    addOption(name, value.c_str());
}

void CommandLine::addOption(const char *name, int value)
{
    append(name);
    append(gmx::toString(value));
}

void CommandLine::addOption(const char *name, double value)
{
    append(name);
    append(gmx::toString(value));
}

void CommandLine::merge(const CommandLine &args)
{
    if (args.argc() > 0)
    {
        // Skip first argument if it is the module name.
        const int firstArg = (args.arg(0)[0] == '-' ? 0 : 1);
        for (int i = firstArg; i < args.argc(); ++i)
        {
            append(args.arg(i));
        }
    }
}

int &CommandLine::argc()
{
    return impl_->argc_;
}
char **CommandLine::argv()
{
    return &impl_->argv_[0];
}
int CommandLine::argc() const
{
    return impl_->argc_;
}
const char *const *CommandLine::argv() const
{
    return &impl_->argv_[0];
}
const char *CommandLine::arg(int i) const
{
    return impl_->argv_[i];
}

std::string CommandLine::toString() const
{
    return CommandLineProgramContext(argc(), argv()).commandLine();
}

bool CommandLine::contains(const char *name) const
{
    for (int i = 0; i < impl_->argc_; ++i)
    {
        if (std::strcmp(arg(i), name) == 0)
        {
            return true;
        }
    }
    return false;
}

/********************************************************************
 * CommandLineTestHelper::Impl
 */

class CommandLineTestHelper::Impl
{
    public:
        struct OutputFileInfo
        {
            OutputFileInfo(const char *option, const std::string &path,
                           FileMatcherPointer matcher)
                : option(option), path(path), matcher(move(matcher))
            {
            }

            std::string              option;
            std::string              path;
            FileMatcherPointer       matcher;
        };

        typedef std::vector<OutputFileInfo>        OutputFileList;

        explicit Impl(TestFileManager *fileManager)
            : fileManager_(*fileManager)
        {
        }

        TestFileManager &fileManager_;
        OutputFileList   outputFiles_;
};

/********************************************************************
 * CommandLineTestHelper
 */

// static
int CommandLineTestHelper::runModuleDirect(
        ICommandLineModule *module, CommandLine *commandLine)
{
    CommandLineModuleSettings settings;
    module->init(&settings);
    return module->run(commandLine->argc(), commandLine->argv());
}

// static
int CommandLineTestHelper::runModuleDirect(
        std::unique_ptr<ICommandLineOptionsModule> module, CommandLine *commandLine)
{
    // The name and description are not used in the tests, so they can be NULL.
    const std::unique_ptr<ICommandLineModule> wrapperModule(
            ICommandLineOptionsModule::createModule(nullptr, nullptr, std::move(module)));
    return runModuleDirect(wrapperModule.get(), commandLine);
}

// static
int CommandLineTestHelper::runModuleFactory(
        std::function<std::unique_ptr<ICommandLineOptionsModule>()>  factory,
        CommandLine                                                 *commandLine)
{
    return runModuleDirect(factory(), commandLine);
}

CommandLineTestHelper::CommandLineTestHelper(TestFileManager *fileManager)
    : impl_(new Impl(fileManager))
{
}

CommandLineTestHelper::~CommandLineTestHelper()
{
}

void CommandLineTestHelper::setInputFileContents(
        CommandLine *args, const char *option, const char *extension,
        const std::string &contents)
{
    GMX_ASSERT(extension[0] != '.', "Extension should not contain a dot");
    std::string fullFilename = impl_->fileManager_.getTemporaryFilePath(
                formatString("%d.%s", args->argc(), extension));
    TextWriter::writeFileFromString(fullFilename, contents);
    args->addOption(option, fullFilename);
}

void CommandLineTestHelper::setInputFileContents(
        CommandLine *args, const char *option, const char *extension,
        const ArrayRef<const char *const> &contents)
{
    GMX_ASSERT(extension[0] != '.', "Extension should not contain a dot");
    std::string fullFilename = impl_->fileManager_.getTemporaryFilePath(
                formatString("%d.%s", args->argc(), extension));
    TextWriter  file(fullFilename);
    ArrayRef<const char *const>::const_iterator i;
    for (i = contents.begin(); i != contents.end(); ++i)
    {
        file.writeLine(*i);
    }
    file.close();
    args->addOption(option, fullFilename);
}

void CommandLineTestHelper::setOutputFile(
        CommandLine *args, const char *option, const char *filename,
        const ITextBlockMatcherSettings &matcher)
{
    setOutputFile(args, option, filename, TextFileMatch(matcher));
}

void CommandLineTestHelper::setOutputFile(
        CommandLine *args, const char *option, const char *filename,
        const IFileMatcherSettings &matcher)
{
    std::string suffix(filename);
    if (startsWith(filename, "."))
    {
        suffix = formatString("%d.%s", args->argc(), filename);
    }
    std::string fullFilename = impl_->fileManager_.getTemporaryFilePath(suffix);
    args->addOption(option, fullFilename);
    impl_->outputFiles_.emplace_back(option, fullFilename, matcher.createFileMatcher());
}

void CommandLineTestHelper::checkOutputFiles(TestReferenceChecker checker) const
{
    if (!impl_->outputFiles_.empty())
    {
        TestReferenceChecker                 outputChecker(
                checker.checkCompound("OutputFiles", "Files"));
        Impl::OutputFileList::const_iterator outfile;
        for (const auto &outfile : impl_->outputFiles_)
        {
            TestReferenceChecker fileChecker(
                    outputChecker.checkCompound("File", outfile.option.c_str()));
            outfile.matcher->checkFile(outfile.path, &fileChecker);
        }
    }
}

/********************************************************************
 * CommandLineTestBase::Impl
 */

class CommandLineTestBase::Impl
{
    public:
        Impl() : helper_(&tempFiles_)
        {
            cmdline_.append("module");
        }

        TestReferenceData     data_;
        TestFileManager       tempFiles_;
        CommandLineTestHelper helper_;
        CommandLine           cmdline_;
};

/********************************************************************
 * CommandLineTestBase
 */

CommandLineTestBase::CommandLineTestBase()
    : impl_(new Impl)
{
}

CommandLineTestBase::~CommandLineTestBase()
{
}

void CommandLineTestBase::setInputFile(
        const char *option, const char *filename)
{
    impl_->cmdline_.addOption(option, TestFileManager::getInputFilePath(filename));
}

void CommandLineTestBase::setInputFileContents(
        const char *option, const char *extension, const std::string &contents)
{
    impl_->helper_.setInputFileContents(&impl_->cmdline_, option, extension,
                                        contents);
}

void CommandLineTestBase::setInputFileContents(
        const char *option, const char *extension,
        const ArrayRef<const char *const> &contents)
{
    impl_->helper_.setInputFileContents(&impl_->cmdline_, option, extension,
                                        contents);
}

void CommandLineTestBase::setOutputFile(
        const char *option, const char *filename,
        const ITextBlockMatcherSettings &matcher)
{
    impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
}

void CommandLineTestBase::setOutputFile(
        const char *option, const char *filename,
        const IFileMatcherSettings &matcher)
{
    impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
}

CommandLine &CommandLineTestBase::commandLine()
{
    return impl_->cmdline_;
}

TestFileManager &CommandLineTestBase::fileManager()
{
    return impl_->tempFiles_;
}

TestReferenceChecker CommandLineTestBase::rootChecker()
{
    return impl_->data_.rootChecker();
}

void CommandLineTestBase::testWriteHelp(ICommandLineModule *module)
{
    StringOutputStream     stream;
    TextWriter             writer(&stream);
    CommandLineHelpContext context(&writer, eHelpOutputFormat_Console, nullptr, "test");
    context.setModuleDisplayName(formatString("%s %s", "test", module->name()));
    module->writeHelp(context);
    TestReferenceChecker   checker(rootChecker());
    checker.checkTextBlock(stream.toString(), "HelpOutput");
}

void CommandLineTestBase::checkOutputFiles()
{
    impl_->helper_.checkOutputFiles(rootChecker());
}

} // namespace test
} // namespace gmx
