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
 * Declares option info objects for basic option types.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONINFO_H
#define GMX_OPTIONS_BASICOPTIONINFO_H

#include "optioninfo.h"

namespace gmx
{

class BooleanOptionStorage;
class IntegerOptionStorage;
class DoubleOptionStorage;
class StringOptionStorage;
class FileNameOptionStorage;

/*! \addtogroup module_options
 * \{
 */

/*! \brief
 * Wrapper class for accessing boolean option information.
 *
 * \inpublicapi
 */
class BooleanOptionInfo : public OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit BooleanOptionInfo(BooleanOptionStorage *option);
};

/*! \brief
 * Wrapper class for accessing integer option information.
 *
 * \inpublicapi
 */
class IntegerOptionInfo : public OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit IntegerOptionInfo(IntegerOptionStorage *option);
};

/*! \brief
 * Wrapper class for accessing floating-point option information.
 *
 * \inpublicapi
 */
class DoubleOptionInfo : public OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit DoubleOptionInfo(DoubleOptionStorage *option);

        //! Whether the option specifies a time value.
        bool isTime() const;

        /*! \brief
         * Sets a scale factor for user-provided values.
         *
         * Any user-provided value is scaled by the provided factor.
         * Programmatically set default values are not scaled.
         * If called multiple times, later calls override the previously set
         * value.  In other words, the scaling is not cumulative.
         */
        void setScaleFactor(double factor);

    private:
        DoubleOptionStorage &option();
        const DoubleOptionStorage &option() const;
};

/*! \brief
 * Wrapper class for accessing string option information.
 *
 * \inpublicapi
 */
class StringOptionInfo : public OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit StringOptionInfo(StringOptionStorage *option);
};

/*! \brief
 * Wrapper class for accessing file name option information.
 *
 * \inpublicapi
 */
class FileNameOptionInfo : public OptionInfo
{
    public:
        //! Creates an option info object for the given option.
        explicit FileNameOptionInfo(FileNameOptionStorage *option);

        //! Whether the option specifies an input file.
        bool isInputFile() const;
        //! Whether the option specifies an output file.
        bool isOutputFile() const;
        //! Whether the option specifies a file used for both input and output.
        bool isInputOutputFile() const;
        /*! \brief
         * Whether the option specifies a library file.
         *
         * \see FileNameOption::libraryFile()
         */
        bool isLibraryFile() const;

    private:
        const FileNameOptionStorage &option() const;
};

/*!\}*/

} // namespace gmx

#endif
