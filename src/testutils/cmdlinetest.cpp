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
/*! \internal \file
 * \brief
 * Implements classes from cmdlinetest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/cmdlinetest.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/futil.h"
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
class FloatingPointTolerance;
class ITextBlockMatcherSettings;

/********************************************************************
 * CommandLine::Impl
 */

class CommandLine::Impl
{
public:
    Impl(const ArrayRef<const char* const>& cmdline);
    Impl(const ArrayRef<const std::string>& cmdline);
    ~Impl();

    std::vector<char*> args_;
    std::vector<char*> argv_;
    int                argc_;
};

CommandLine::Impl::Impl(const ArrayRef<const char* const>& cmdline)
{
    args_.reserve(cmdline.size());
    argv_.reserve(cmdline.size() + 1);
    argc_ = ssize(cmdline);
    for (const auto& arg : cmdline)
    {
        char* argCopy = strdup(arg);
        if (argCopy == nullptr)
        {
            throw std::bad_alloc();
        }
        args_.push_back(argCopy);
        argv_.push_back(argCopy);
    }
    argv_.push_back(nullptr);
}

namespace
{

//! Helper function so we can delegate from the std::string constructor to the const char * one.
std::vector<const char*> convertFromStringArrayRef(const ArrayRef<const std::string>& cmdline)
{
    std::vector<const char*> v(cmdline.size());
    std::transform(
            cmdline.begin(), cmdline.end(), v.begin(), [](const std::string& s) { return s.c_str(); });
    return v;
}

} // namespace

// This makes a new temporary vector of views of the const char * in
// the view passed in. Those are then deep copied in the constructor
// delegated to.
CommandLine::Impl::Impl(const ArrayRef<const std::string>& cmdline) :
    Impl(convertFromStringArrayRef(cmdline))
{
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

CommandLine::CommandLine() : impl_(new Impl(ArrayRef<const char*>{})) {}

CommandLine::CommandLine(const ArrayRef<const char* const>& cmdline) : impl_(new Impl(cmdline)) {}

CommandLine::CommandLine(const ArrayRef<const std::string>& cmdline) : impl_(new Impl(cmdline)) {}

CommandLine::CommandLine(const CommandLine& other) :
    impl_(new Impl(arrayRefFromArray(other.argv(), other.argc())))
{
}

CommandLine::~CommandLine() {}

void CommandLine::initFromArray(const ArrayRef<const char* const>& cmdline)
{
    impl_ = std::make_unique<Impl>(cmdline);
}

void CommandLine::append(const char* arg)
{
    GMX_RELEASE_ASSERT(impl_->argc_ == gmx::ssize(impl_->args_),
                       "Command-line has been modified externally");
    size_t newSize = impl_->args_.size() + 1;
    impl_->args_.reserve(newSize);
    impl_->argv_.reserve(newSize + 1);
    char* newArg = strdup(arg);
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

void CommandLine::addOption(const char* name)
{
    append(name);
}

void CommandLine::addOption(const char* name, const char* value)
{
    append(name);
    append(value);
}

void CommandLine::addOption(const char* name, const std::string& value)
{
    addOption(name, value.c_str());
}

void CommandLine::addOption(const char* name, int value)
{
    append(name);
    append(gmx::toString(value));
}

void CommandLine::addOption(const char* name, double value)
{
    append(name);
    append(gmx::toString(value));
}

void CommandLine::merge(const CommandLine& args)
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

int& CommandLine::argc()
{
    return impl_->argc_;
}
char** CommandLine::argv()
{
    return &impl_->argv_[0];
}
int CommandLine::argc() const
{
    return impl_->argc_;
}
const char* const* CommandLine::argv() const
{
    return &impl_->argv_[0];
}
const char* CommandLine::arg(int i) const
{
    return impl_->argv_[i];
}

std::string CommandLine::toString() const
{
    return CommandLineProgramContext(argc(), argv()).commandLine();
}

namespace
{

//! Function object that returns whether the view matches the C-style string
struct ViewMatchesCString
{
    const std::string_view view_;
    explicit ViewMatchesCString(const std::string_view view) : view_(view) {}

    //! The operation
    bool operator()(const char* cString)
    {
        if (cString == nullptr)
        {
            return view_.empty();
        }
        return view_.compare(cString) == 0;
    }
};

} // namespace

bool CommandLine::contains(const std::string_view name) const
{
    return std::find_if(impl_->argv_.begin(), impl_->argv_.end(), ViewMatchesCString(name))
           != impl_->argv_.end();
}

std::optional<std::string_view> CommandLine::argumentOf(const std::string_view name) const
{
    auto foundIt = std::find_if(impl_->argv_.begin(), impl_->argv_.end(), ViewMatchesCString(name));
    if (foundIt == impl_->argv_.end())
    {
        // The named option was not found
        return std::nullopt;
    }
    return std::make_optional<std::string_view>(*++foundIt);
}

/********************************************************************
 * CommandLineTestHelper::Impl
 */

class CommandLineTestHelper::Impl
{
public:
    struct OutputFileInfo
    {
        OutputFileInfo(const char* option, const std::string& path, FileMatcherPointer matcher) :
            option_(option), path_(path), matcher_(std::move(matcher))
        {
        }

        std::string        option_;
        std::string        path_;
        FileMatcherPointer matcher_;
    };

    typedef std::vector<OutputFileInfo> OutputFileList;

    explicit Impl(TestFileManager* fileManager) : fileManager_(*fileManager) {}

    TestFileManager& fileManager_;
    OutputFileList   outputFiles_;
};

/********************************************************************
 * CommandLineTestHelper
 */

// static
int CommandLineTestHelper::runModuleDirect(ICommandLineModule* module, CommandLine* commandLine)
{
    CommandLineModuleSettings settings;
    module->init(&settings);
    return module->run(commandLine->argc(), commandLine->argv());
}

// static
int CommandLineTestHelper::runModuleDirect(std::unique_ptr<ICommandLineOptionsModule> module,
                                           CommandLine*                               commandLine)
{
    // The name and description are not used in the tests, so they can be NULL.
    const std::unique_ptr<ICommandLineModule> wrapperModule(
            ICommandLineOptionsModule::createModule(nullptr, nullptr, std::move(module)));
    return runModuleDirect(wrapperModule.get(), commandLine);
}

// static
int CommandLineTestHelper::runModuleFactory(
        const std::function<std::unique_ptr<ICommandLineOptionsModule>()>& factory,
        CommandLine*                                                       commandLine)
{
    return runModuleDirect(factory(), commandLine);
}

CommandLineTestHelper::CommandLineTestHelper(TestFileManager* fileManager) :
    impl_(new Impl(fileManager))
{
}

CommandLineTestHelper::~CommandLineTestHelper() {}

void CommandLineTestHelper::setInputFileContents(CommandLine*       args,
                                                 const char*        option,
                                                 const char*        extension,
                                                 const std::string& contents)
{
    GMX_ASSERT(extension[0] != '.', "Extension should not contain a dot");
    std::string fullFilename =
            impl_->fileManager_.getTemporaryFilePath(formatString("%d.%s", args->argc(), extension)).string();
    TextWriter::writeFileFromString(fullFilename, contents);
    args->addOption(option, fullFilename);
}

void CommandLineTestHelper::setInputFileContents(CommandLine*                       args,
                                                 const char*                        option,
                                                 const char*                        extension,
                                                 const ArrayRef<const char* const>& contents)
{
    GMX_ASSERT(extension[0] != '.', "Extension should not contain a dot");
    std::string fullFilename =
            impl_->fileManager_.getTemporaryFilePath(formatString("%d.%s", args->argc(), extension)).string();
    TextWriter                                  file(fullFilename);
    ArrayRef<const char* const>::const_iterator i;
    for (i = contents.begin(); i != contents.end(); ++i)
    {
        file.writeLine(*i);
    }
    file.close();
    args->addOption(option, fullFilename);
}

std::string CommandLineTestHelper::setOutputFile(CommandLine*                     args,
                                                 const char*                      option,
                                                 const char*                      filename,
                                                 const ITextBlockMatcherSettings& matcher)
{
    return setOutputFile(args, option, filename, TextFileMatch(matcher));
}

std::string CommandLineTestHelper::setOutputFile(CommandLine*                args,
                                                 const char*                 option,
                                                 const char*                 filename,
                                                 const IFileMatcherSettings& matcher)
{
    std::string suffix(filename);
    if (startsWith(filename, "."))
    {
        suffix = formatString("%d.%s", args->argc(), filename);
    }
    std::string fullFilename = impl_->fileManager_.getTemporaryFilePath(suffix).string();
    args->addOption(option, fullFilename);
    impl_->outputFiles_.emplace_back(option, fullFilename, matcher.createFileMatcher());
    return fullFilename;
}

std::string CommandLineTestHelper::setOutputFileWithGeneratedName(const char* filename,
                                                                  const ITextBlockMatcherSettings& matcher)
{
    return setOutputFileWithGeneratedName(std::string(filename), TextFileMatch(matcher));
}

std::string CommandLineTestHelper::setOutputFileWithGeneratedName(std::string&& filename,
                                                                  const ITextBlockMatcherSettings& matcher)
{
    return setOutputFileWithGeneratedName(std::move(filename), TextFileMatch(matcher));
}

std::string CommandLineTestHelper::setOutputFileWithGeneratedName(const char* filename,
                                                                  const IFileMatcherSettings& matcher)
{
    return setOutputFileWithGeneratedName(std::string(filename), matcher);
}

std::string CommandLineTestHelper::setOutputFileWithGeneratedName(std::string&& filename,
                                                                  const IFileMatcherSettings& matcher)
{
    impl_->outputFiles_.emplace_back(filename.c_str(), filename, matcher.createFileMatcher());
    std::string filenameToReturn(filename);
    impl_->fileManager_.manageGeneratedOutputFile(std::move(filename));
    return filenameToReturn;
}

void CommandLineTestHelper::checkOutputFiles(TestReferenceChecker checker) const
{
    if (!impl_->outputFiles_.empty())
    {
        TestReferenceChecker outputChecker(checker.checkCompound("OutputFiles", "Files"));
        for (const auto& outfile : impl_->outputFiles_)
        {
            TestReferenceChecker fileChecker(outputChecker.checkCompound("File", outfile.option_.c_str()));
            outfile.matcher_->checkFile(outfile.path_, &fileChecker);
        }
    }
}

/********************************************************************
 * CommandLineTestBase::Impl
 */

class CommandLineTestBase::Impl
{
public:
    Impl() : helper_(&tempFiles_) { cmdline_.append("module"); }

    TestReferenceData     data_;
    TestFileManager       tempFiles_;
    CommandLineTestHelper helper_;
    CommandLine           cmdline_;
};

/********************************************************************
 * CommandLineTestBase
 */

CommandLineTestBase::CommandLineTestBase() : impl_(new Impl) {}

CommandLineTestBase::~CommandLineTestBase() {}

void CommandLineTestBase::setInputFile(const char* option, const char* filename)
{
    impl_->cmdline_.addOption(option, TestFileManager::getInputFilePath(filename).string());
}

void CommandLineTestBase::setInputFile(const char* option, const std::string& filename)
{
    setInputFile(option, filename.c_str());
}

void CommandLineTestBase::setModifiableInputFile(const char* option, const std::string& filename)
{
    setModifiableInputFile(option, filename.c_str());
}

void CommandLineTestBase::setModifiableInputFile(const char* option, const char* filename)
{
    std::string originalFileName = gmx::test::TestFileManager::getInputFilePath(filename).string();
    std::string modifiableFileName = fileManager().getTemporaryFilePath(filename).string();
    gmx_file_copy(originalFileName.c_str(), modifiableFileName.c_str(), true);
    impl_->cmdline_.addOption(option, modifiableFileName);
}

void CommandLineTestBase::setInputFileContents(const char*        option,
                                               const char*        extension,
                                               const std::string& contents)
{
    impl_->helper_.setInputFileContents(&impl_->cmdline_, option, extension, contents);
}

void CommandLineTestBase::setInputFileContents(const char*                        option,
                                               const char*                        extension,
                                               const ArrayRef<const char* const>& contents)
{
    impl_->helper_.setInputFileContents(&impl_->cmdline_, option, extension, contents);
}

std::string CommandLineTestBase::setOutputFile(const char*                      option,
                                               const char*                      filename,
                                               const ITextBlockMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
}

std::string CommandLineTestBase::setOutputFile(const char*                 option,
                                               const char*                 filename,
                                               const IFileMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
}

std::string CommandLineTestBase::setOutputFileWithGeneratedName(const char* filename,
                                                                const ITextBlockMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFileWithGeneratedName(std::string(filename), matcher);
}

std::string CommandLineTestBase::setOutputFileWithGeneratedName(std::string&& filename,
                                                                const ITextBlockMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFileWithGeneratedName(std::move(filename), matcher);
}

std::string CommandLineTestBase::setOutputFileWithGeneratedName(const char* filename,
                                                                const IFileMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFileWithGeneratedName(std::string(filename), matcher);
}

std::string CommandLineTestBase::setOutputFileWithGeneratedName(std::string&& filename,
                                                                const IFileMatcherSettings& matcher)
{
    return impl_->helper_.setOutputFileWithGeneratedName(std::move(filename), matcher);
}

std::string CommandLineTestBase::setInputAndOutputFile(const char*                      option,
                                                       const char*                      filename,
                                                       const ITextBlockMatcherSettings& matcher)
{
    std::string originalFileName = gmx::test::TestFileManager::getInputFilePath(filename).string();
    std::string modifiableFileName = fileManager().getTemporaryFilePath(filename).string();
    gmx_file_copy(originalFileName.c_str(), modifiableFileName.c_str(), true);
    impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
    return modifiableFileName;
}

std::string CommandLineTestBase::setInputAndOutputFile(const char*                 option,
                                                       const char*                 filename,
                                                       const IFileMatcherSettings& matcher)
{
    std::string originalFileName = gmx::test::TestFileManager::getInputFilePath(filename).string();
    std::string modifiableFileName = fileManager().getTemporaryFilePath(filename).string();
    gmx_file_copy(originalFileName.c_str(), modifiableFileName.c_str(), true);
    impl_->helper_.setOutputFile(&impl_->cmdline_, option, filename, matcher);
    return modifiableFileName;
}

CommandLine& CommandLineTestBase::commandLine()
{
    return impl_->cmdline_;
}

TestFileManager& CommandLineTestBase::fileManager()
{
    return impl_->tempFiles_;
}

TestReferenceChecker CommandLineTestBase::rootChecker()
{
    return impl_->data_.rootChecker();
}

void CommandLineTestBase::setDefaultTolerance(const FloatingPointTolerance& tolerance)
{
    impl_->data_.rootChecker().setDefaultTolerance(tolerance);
}

void CommandLineTestBase::testWriteHelp(ICommandLineModule* module)
{
    StringOutputStream     stream;
    TextWriter             writer(&stream);
    CommandLineHelpContext context(&writer, eHelpOutputFormat_Console, nullptr, "test");
    context.setModuleDisplayName(formatString("%s %s", "test", module->name()));
    module->writeHelp(context);
    TestReferenceChecker checker(rootChecker());
    checker.checkTextBlock(stream.toString(), "HelpOutput");
}

void CommandLineTestBase::checkOutputFiles()
{
    impl_->helper_.checkOutputFiles(rootChecker());
}

} // namespace test
} // namespace gmx
