/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::CommandLineParser.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
#include "cmdlineparser.h"

#include <cctype>

#include <string>
#include <vector>

#include "gromacs/options/optionsassigner.h"
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

        //! Helper object for assigning the options.
        OptionsAssigner         assigner_;
};

CommandLineParser::Impl::Impl(Options *options)
    : assigner_(options)
{
    assigner_.setAcceptBooleanNoPrefix(true);
    assigner_.setNoStrictSectioning(true);
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

void CommandLineParser::parse(int *argc, char *argv[])
{
    std::vector<std::string> commandLine;
    for (int i = 0; i < *argc; ++i)
    {
        commandLine.push_back(argv[i]);
    }
    parse(&commandLine);
}

void CommandLineParser::parse(std::vector<std::string> *commandLine)
{
    ExceptionInitializer errors("Invalid command-line options");
    std::string          currentContext;
    // Start in the discard phase to skip options that can't be understood.
    bool                 bDiscard = true;

    impl_->assigner_.start();
    std::vector<std::string>::const_iterator arg;
    for (arg = commandLine->begin() + 1; arg != commandLine->end(); ++arg)
    {
        // Lone '-' and numbers are passed as values.
        if ((*arg)[0] == '-' && std::isalpha((*arg)[1]))
        {
            if (!bDiscard)
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
                currentContext.clear();
            }
            currentContext = "In command-line option " + *arg;
            bDiscard       = false;
            try
            {
                const char *name = arg->c_str() + 1;
                impl_->assigner_.startOption(name);
            }
            catch (UserInputError &ex)
            {
                bDiscard = true;
                ex.prependContext(currentContext);
                errors.addCurrentExceptionAsNested();
                currentContext.clear();
            }
        }
        else if (!bDiscard)
        {
            try
            {
                impl_->assigner_.appendValue(*arg);
            }
            catch (UserInputError &ex)
            {
                ex.prependContext(currentContext);
                errors.addCurrentExceptionAsNested();
            }
        }
    }
    if (!bDiscard)
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
