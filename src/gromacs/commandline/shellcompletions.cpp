/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * Implements gmx::ShellCompletionWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "shellcompletions.h"

#include <cstdio>

#include <algorithm>
#include <memory>
#include <string>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/options/abstractoption.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{
class OptionSectionInfo;

namespace
{

class OptionsListWriter : public OptionsVisitor
{
public:
    const std::string& optionList() const { return optionList_; }

    void visitSection(const OptionSectionInfo& section) override
    {
        OptionsIterator iterator(section);
        iterator.acceptSections(this);
        iterator.acceptOptions(this);
    }
    void visitOption(const OptionInfo& option) override
    {
        if (option.isHidden())
        {
            return;
        }
        if (!optionList_.empty())
        {
            optionList_.append("\\n");
        }
        optionList_.append("-");
        const BooleanOptionInfo* booleanOption = option.toType<BooleanOptionInfo>();
        if (booleanOption != nullptr && booleanOption->defaultValue())
        {
            optionList_.append("no");
        }
        optionList_.append(option.name());
    }

private:
    std::string optionList_;
};

class OptionCompletionWriter : public OptionsVisitor
{
public:
    explicit OptionCompletionWriter(TextWriter* out) : out_(*out) {}

    void visitSection(const OptionSectionInfo& section) override
    {
        OptionsIterator iterator(section);
        iterator.acceptSections(this);
        iterator.acceptOptions(this);
    }
    void visitOption(const OptionInfo& option) override;

private:
    void writeOptionCompletion(const OptionInfo& option, const std::string& completion);

    TextWriter& out_;
};

void OptionCompletionWriter::visitOption(const OptionInfo& option)
{
    if (option.isHidden())
    {
        return;
    }
    const FileNameOptionInfo* fileOption = option.toType<FileNameOptionInfo>();
    if (fileOption != nullptr)
    {
        if (fileOption->isDirectoryOption())
        {
            writeOptionCompletion(option, "compgen -S ' ' -d $c");
            return;
        }
        const FileNameOptionInfo::ExtensionList& extensionList = fileOption->extensions();
        if (extensionList.empty())
        {
            return;
        }
        std::string completion("compgen -S ' ' -X '!*");
        std::string extensions(joinStrings(extensionList, "|"));
        if (extensionList.size() > 1)
        {
            extensions = "@(" + extensions + ")";
        }
        completion.append(extensions);
        // TODO: Don't duplicate this from filenm.c/futil.c.
        completion.append("?(.gz|.Z)' -f -- $c ; compgen -S '/' -d $c");
        writeOptionCompletion(option, completion);
        return;
    }
    const StringOptionInfo* stringOption = option.toType<StringOptionInfo>();
    if (stringOption != nullptr && stringOption->isEnumerated())
    {
        std::string completion("compgen -S ' ' -W $'");
        completion.append(joinStrings(stringOption->allowedValues(), "\\n"));
        completion.append("' -- $c");
        writeOptionCompletion(option, completion);
        return;
    }
}

void OptionCompletionWriter::writeOptionCompletion(const OptionInfo& option, const std::string& completion)
{
    std::string result(formatString("-%s) ", option.name().c_str()));
    if (option.maxValueCount() >= 0)
    {
        result.append(formatString("(( $n <= %d )) && ", option.maxValueCount()));
    }
    result.append("COMPREPLY=( $(");
    result.append(completion);
    result.append("));;");
    out_.writeLine(result);
}

} // namespace

class ShellCompletionWriter::Impl
{
public:
    Impl(const std::string& binaryName, ShellCompletionFormat /*format*/) : binaryName_(binaryName)
    {
    }

    std::string completionFunctionName(const char* moduleName) const
    {
        std::string result = formatString("_%s_%s_compl", binaryName_.c_str(), moduleName);
        std::replace(result.begin(), result.end(), '-', '_');
        return result;
    }

    std::string binaryName_;
    // Never releases ownership.
    std::unique_ptr<TextWriter> file_;
};

ShellCompletionWriter::ShellCompletionWriter(const std::string& binaryName, ShellCompletionFormat format) :
    impl_(new Impl(binaryName, format))
{
}

ShellCompletionWriter::~ShellCompletionWriter() {}

TextWriter& ShellCompletionWriter::outputWriter()
{
    return *impl_->file_;
}

void ShellCompletionWriter::startCompletions()
{
    impl_->file_ = std::make_unique<TextWriter>(impl_->binaryName_ + "-completion.bash");
}

void ShellCompletionWriter::writeModuleCompletions(const char* moduleName, const Options& options)
{
    TextWriter& out = *impl_->file_;
    out.writeLine(formatString("%s() {", impl_->completionFunctionName(moduleName).c_str()));
    out.writeLine("local IFS=$'\\n'");
    out.writeLine("local c=${COMP_WORDS[COMP_CWORD]}");
    out.writeLine("local n");
    out.writeLine(
            "for ((n=1;n<COMP_CWORD;++n)) ; do [[ \"${COMP_WORDS[COMP_CWORD-n]}\" == -* ]] && "
            "break ; done");
    out.writeLine("local p=${COMP_WORDS[COMP_CWORD-n]}");
    out.writeLine("COMPREPLY=()");

    OptionsListWriter listWriter;
    listWriter.visitSection(options.rootSection());
    out.writeLine(
            formatString("if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen "
                         "-S ' '  -W $'%s' -- $c)); return 0; fi",
                         listWriter.optionList().c_str()));

    out.writeLine("case \"$p\" in");
    OptionCompletionWriter optionWriter(&out);
    optionWriter.visitSection(options.rootSection());
    out.writeLine("esac }");
}

void ShellCompletionWriter::writeWrapperCompletions(const ModuleNameList& modules, const Options& options)
{
    impl_->file_->writeLine("_" + impl_->binaryName_ + "_compl() {");
    impl_->file_->writeLine("local i c m");
    impl_->file_->writeLine("local IFS=$'\\n'\n");
    impl_->file_->writeLine("COMPREPLY=()");
    impl_->file_->writeLine("unset COMP_WORDS[0]");
    impl_->file_->writeLine("for ((i=1;i<COMP_CWORD;++i)) ; do");
    impl_->file_->writeLine("[[ \"${COMP_WORDS[i]}\" != -* ]] && break");
    impl_->file_->writeLine("unset COMP_WORDS[i]");
    impl_->file_->writeLine("done");
    impl_->file_->writeLine("if (( i == COMP_CWORD )); then");
    impl_->file_->writeLine("c=${COMP_WORDS[COMP_CWORD]}");
    OptionsListWriter lister;
    lister.visitSection(options.rootSection());
    std::string completions(lister.optionList());
    for (const auto& moduleName : modules)
    {
        completions.append("\\n");
        completions.append(moduleName);
    }
    impl_->file_->writeLine("COMPREPLY=( $(compgen -S ' ' -W $'" + completions + "' -- $c) )");
    impl_->file_->writeLine("return 0");
    impl_->file_->writeLine("fi");
    impl_->file_->writeLine("m=${COMP_WORDS[i]}");
    impl_->file_->writeLine("COMP_WORDS=( \"${COMP_WORDS[@]}\" )");
    impl_->file_->writeLine("COMP_CWORD=$((COMP_CWORD-i))");
    impl_->file_->writeLine("case \"$m\" in");
    for (const auto& moduleName : modules)
    {
        const char* const name = moduleName.c_str();
        impl_->file_->writeLine(
                formatString("%s) %s ;;", name, impl_->completionFunctionName(name).c_str()));
    }
    impl_->file_->writeLine("esac }");
}

void ShellCompletionWriter::finishCompletions()
{
    impl_->file_->close();
    impl_->file_.reset();
}

} // namespace gmx
