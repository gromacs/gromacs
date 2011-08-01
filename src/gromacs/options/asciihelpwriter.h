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
/*! \file
 * \brief
 * Declares gmx::AsciiHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ASCIIHELPWRITER_H
#define GMX_OPTIONS_ASCIIHELPWRITER_H

#include <cstdio>

namespace gmx
{

class Options;

/*! \brief
 * Writes help information for Options in ascii format.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class AsciiHelpWriter
{
    public:
        /*! \brief
         * Creates an object that writer ascii-formatted help for Options.
         *
         * \param[in] options  Options for which help should be printed.
         */
        explicit AsciiHelpWriter(const Options &options);
        ~AsciiHelpWriter();

        /*! \brief
         * Sets whether hidden options are shown in the help.
         */
        AsciiHelpWriter &setShowHidden(bool bShow);
        /*! \brief
         * Sets whether long descriptions for sections are shown in the help.
         */
        AsciiHelpWriter &setShowDescriptions(bool bShow);

        /*! \brief
         * Writes the help.
         *
         * \param[in] fp  File to write the help to.
         */
        void writeHelp(FILE *fp);

    private:
        class Impl;

        Impl                   *_impl;

        // Disallow copy and assign.
        AsciiHelpWriter(const AsciiHelpWriter &);
        void operator =(const AsciiHelpWriter &);
};

} // namespace gmx

#endif
