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
 * Declares private implementation class for gmx::CommandLineHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_CMDLINEHELPWRITER_IMPL_H
#define GMX_OPTIONS_CMDLINEHELPWRITER_IMPL_H

#include "cmdlinehelpwriter.h"

#include "basicoptioninfo.h"
#include "optionsvisitor.h"
#include "../utility/flags.h"

namespace gmx
{

class Options;

/*! \internal \brief
 * Private implementation class for CommandLineHelpWriter.
 *
 * \ingroup module_options
 */
class CommandLineHelpWriter::Impl
{
    public:
        /*! \brief
         * Possible flags for controlling output.
         */
        enum Flag
        {
            efShowDescriptions  = 1<<0,
            efShowHidden        = 1<<1
        };

        //! Sets the Options object to use for generating help.
        explicit Impl(const Options &options);

        //! Options object to use for generating help.
        const Options          &_options;
        //! Flags for controlling what to write.
        FlagsTemplate<Flag>     _flags;
};

/*! \internal \brief
 * Helper object for writing section descriptions to help.
 *
 * \ingroup module_options
 */
class AsciiDescriptionWriter : public OptionsVisitor
{
    public:
        //! Creates a helper object for writing section descriptions.
        explicit AsciiDescriptionWriter(FILE *fp) : _fp(fp) {}

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo & /*option*/) { }

    private:
        FILE                   *_fp;
};

/*! \internal \brief
 * Helper object for writing help for file parameters.
 *
 * \ingroup module_options
 */
class AsciiFileParameterWriter : public OptionsTypeVisitor<FileNameOptionInfo>
{
    public:
        //! Creates a helper object for writing file parameters.
        explicit AsciiFileParameterWriter(FILE *fp)
            : _fp(fp), _bFirst(true)
        {
        }

        //! Returns true if anything was written out.
        bool didOutput() const { return !_bFirst; }

        virtual void visitSubSection(const Options &section);
        virtual void visitOptionType(const FileNameOptionInfo &option);

    private:
        FILE                   *_fp;
        bool                    _bFirst;
};

/*! \internal \brief
 * Helper object for writing help for non-file parameters.
 *
 * \ingroup module_options
 */
class AsciiParameterWriter : public OptionsVisitor
{
    public:
        //! Creates a helper object for writing non-file parameters.
        explicit AsciiParameterWriter(FILE *fp)
            : _fp(fp), _bFirst(true), _bShowHidden(false)
        {
        }

        //! Sets the writer to show hidden options.
        void setShowHidden(bool bSet)
        {
            _bShowHidden = bSet;
        }
        //! Returns true if anything was written out.
        bool didOutput() const { return !_bFirst; }

        virtual void visitSubSection(const Options &section);
        virtual void visitOption(const OptionInfo &option);

    private:
        FILE                   *_fp;
        bool                    _bFirst;
        bool                    _bShowHidden;
};

} // namespace gmx

#endif
