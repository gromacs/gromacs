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
 * Declares gmx::OptionInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONINFO_H
#define GMX_OPTIONS_OPTIONINFO_H

#include <cstddef>

#include <string>

#include "../utility/common.h"

namespace gmx
{

class AbstractOptionStorage;

/*! \brief
 * Wrapper class for accessing option information.
 *
 * This class isolates the details of the internal option implementation
 * from option visitors.  Non-const methods in this class or in derived classes
 * also allow modifying the underlying option after its initial creation with
 * Options::addOption().
 *
 * \see OptionsVisitor
 * \see OptionsModifyingVisitor
 *
 * \inpublicapi
 * \ingroup module_options
 */
class OptionInfo
{
    public:
        virtual ~OptionInfo();

        /*! \brief
         * Test whether the option is of a particular type.
         *
         * \tparam InfoType  Option type to test for. Should be a class derived
         *      from OptionInfo.
         */
        template <class InfoType>
        bool isType() const
        {
            return toType<InfoType>() != NULL;
        }
        /*! \brief
         * Convert the info object to a particular type if the type is correct.
         *
         * \tparam InfoType  Option type to convert to. Should be a class
         *      derived from OptionInfo.
         * \retval this converted to a pointer to \p InfoType, or NULL if the
         *      conversion is not possible.
         */
        template <class InfoType>
        InfoType *toType()
        {
            return dynamic_cast<InfoType *>(this);
        }
        //! \copydoc toType()
        template <class InfoType>
        const InfoType *toType() const
        {
            return dynamic_cast<const InfoType *>(this);
        }

        //! Returns true if the option has been set.
        bool isSet() const;
        //! Returns true if the option is a hidden option.
        bool isHidden() const;
        //! Returns the name of the option.
        const std::string &name() const;
        //! Returns the description of the option.
        const std::string &description() const;
        //! Returns the type of the option as a string.
        const char *type() const;
        //! Returns the number of values given for the option.
        int valueCount() const;
        //! Returns the i'th value of the option as a string.
        std::string formatValue(int i) const;
        //! Returns all the values of the option as a single string.
        std::string formatValues() const;

    protected:
        /*! \cond libapi */
        /*! \brief
         * Wraps a given option object.
         *
         * Does not throw.
         */
        explicit OptionInfo(AbstractOptionStorage *option);

        //! Returns the wrapped option storage object.
        AbstractOptionStorage &option() { return _option; }
        //! Returns the wrapped option storage object.
        const AbstractOptionStorage &option() const { return _option; }
        //! \endcond

    private:
        //! The wrapped option.
        AbstractOptionStorage  &_option;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionInfo);
};

} // namespace gmx

#endif
