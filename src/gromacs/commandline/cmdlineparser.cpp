/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2017,2018, by the GROMACS development team, led by
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
 * Implements gmx::CommandLineParser.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlineparser.h"

#include <cstdlib>

#include <string>
#include <vector>

#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/********************************************************************
 * CommandLineParser::Impl
 */

/*! \internal \brief
 * Private implementation class for CommandLineParser.
 *
 * \ingroup module_commandline
 */
class CommandLineParser::Impl
{
    public:
        //! Sets the options object to parse to.
        explicit Impl(Options *options);

        /*! \brief
         * Determines whether a cmdline parameter starts an option and the name
         * of that option.
         *
         * \param[in] arg  Individual argument from \c argv.
         * \returns The beginning of the option name in \p arg, or NULL if
         *     \p arg does not look like an option.
         */
        const char *toOptionName(const char *arg) const;

        //! Helper object for assigning the options.
        OptionsAssigner         assigner_;
        //! Whether to allow and skip unknown options.
        bool                    bSkipUnknown_;
        /*! \brief Whether to allow positional arguments
         *
         * These are not options (no leading hyphen), and come before
         * all options. */
        bool                    bAllowPositionalArguments_;
};

CommandLineParser::Impl::Impl(Options *options)
    : assigner_(options), bSkipUnknown_(false), bAllowPositionalArguments_(false)
{
    assigner_.setAcceptBooleanNoPrefix(true);
}

const char *CommandLineParser::Impl::toOptionName(const char *arg) const
{
    // Lone '-' or '--' is not an option.
    if (arg[0] != '-' || arg[1] == '\0' || (arg[1] == '-' && arg[2] == '\0'))
    {
        return nullptr;
    }
    // Something starting with '--' is always an option.
    if (arg[1] == '-')
    {
        return arg + 2;
    }
    // Don't return numbers as option names.
    char *endptr;
    // We are only interested in endptr, not in the actual value.
    GMX_IGNORE_RETURN_VALUE(std::strtod(arg, &endptr));
    if (*endptr == '\0')
    {
        return nullptr;
    }
    return arg + 1;
}

/********************************************************************
 * CommandLineParser
 */

CommandLineParser::CommandLineParser(Options *options)
    : impl_(new Impl(options))
{
}

CommandLineParser::~CommandLineParser()
{
}

CommandLineParser &CommandLineParser::skipUnknown(bool bEnabled)
{
    impl_->bSkipUnknown_ = bEnabled;
    return *this;
}

CommandLineParser &CommandLineParser::allowPositionalArguments(bool bEnabled)
{
    impl_->bAllowPositionalArguments_ = bEnabled;
    return *this;
}

void CommandLineParser::parse(int *argc, char *argv[])
{
    ExceptionInitializer errors("Invalid command-line options");
    std::string          currentContext;
    bool                 bInOption = false;

    // Note that this function gets called multiple times in typical
    // cases of calling gmx. Command lines like "gmx -hidden mdrun -h"
    // work because the first call has argv[0] == "gmx" and skips
    // unknown things, and the second has argv[0] == "mdrun".
    int i = 1, newi = 1;

    // First, process any permitted leading positional arguments.
    for (; i < *argc; ++i)
    {
        const char *const arg        = argv[i];
        if (impl_->toOptionName(arg) != nullptr)
        {
            // If we find an option, no more positional arguments
            // can be handled.
            break;
        }

        if (!impl_->bAllowPositionalArguments_)
        {
            GMX_THROW(InvalidInputError
                          ("Positional argument '" + std::string(arg) + "' cannot be accepted. "
                          "Perhaps you forgot to put a hyphen before an option name."));
        }
        // argv[i] is not an option, so preserve it in the argument list
        // by incrementing newi. There's no need to copy argv contents
        // because they cannot have changed yet.
        ++newi;
    }

    // Now handle the option arguments.
    impl_->assigner_.start();
    for (; i < *argc; ++i)
    {
        const char *const arg        = argv[i];
        const char *const optionName = impl_->toOptionName(arg);
        if (optionName != nullptr)
        {
            if (bInOption)
            {
                try
                {
                    impl_->assigner_.finishOption();
                }
                catch (UserInputError &ex)
                {
                    ex.prependContext(currentContext);
                    errors.addCurrentExceptionAsNested();
                }
            }
            currentContext = "In command-line option " + std::string(arg);
            try
            {
                bInOption = impl_->assigner_.tryStartOption(optionName);
                if (!bInOption)
                {
                    currentContext.clear();
                    if (!impl_->bSkipUnknown_)
                    {
                        std::string message =
                            "Unknown command-line option " + std::string(arg);
                        GMX_THROW(InvalidInputError(message));
                    }
                }
            }
            catch (UserInputError &ex)
            {
                // If tryStartOption() throws, make sure that the rest gets
                // ignored.
                // TODO: Consider whether we should remove the option from the
                // command line nonetheless, as it is recognized, but just
                // invalid.
                bInOption = false;
                ex.prependContext(currentContext);
                errors.addCurrentExceptionAsNested();
                currentContext.clear();
            }
        }
        else if (bInOption)
        {
            try
            {
                impl_->assigner_.appendValue(arg);
            }
            // TODO: Consider if some types of exceptions would be better left
            // unhandled.
            catch (GromacsException &ex)
            {
                ex.prependContext(currentContext);
                errors.addCurrentExceptionAsNested();
            }
        }
        // Retain unrecognized options if applicable.
        if (!bInOption && impl_->bSkipUnknown_)
        {
            argv[newi] = argv[i];
            ++newi;
        }
    }
    // Update the args if argv was modified.
    if (impl_->bSkipUnknown_)
    {
        *argc      = newi;
        argv[newi] = nullptr;
    }

    // Finish the last option.
    if (bInOption)
    {
        try
        {
            impl_->assigner_.finishOption();
        }
        catch (UserInputError &ex)
        {
            ex.prependContext(currentContext);
            errors.addCurrentExceptionAsNested();
        }
    }
    impl_->assigner_.finish();
    if (errors.hasNestedExceptions())
    {
        // TODO: This exception type may not always be appropriate.
        GMX_THROW(InvalidInputError(errors));
    }
}

} // namespace gmx
