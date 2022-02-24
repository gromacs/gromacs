/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Provides gmx::OptionValueConverterSimple.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_VALUECONVERTER_H
#define GMX_OPTIONS_VALUECONVERTER_H

#include <functional>
#include <map>
#include <typeindex>

#include "gromacs/utility/any.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \libinternal \brief
 * Helper for converting from Any to a given type.
 *
 * \tparam OutType  Type this converter converts to.
 *
 * Default-constructed converter only supports identity mapping from the a
 * Any holding `OutType`.  To add support for additional input types,
 * provide conversion functions with addConverter().  To use a non-identity
 * mapping for an `OutType` -> `OutType` conversion, provide an alternative
 * conversion from `OutType` with addConverter().
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template<typename OutType>
class OptionValueConverterSimple
{
public:
    /*! \brief
     * Converts a Any value to the output type.
     *
     * \returns  Converted value.
     * \throws InvalidInputError If the input Any has a type that is
     *     not recognized by any conversion.
     */
    OutType convert(const Any& value) const
    {
        std::type_index type(value.type());
        auto            iter = converters_.find(type);
        if (iter == converters_.end())
        {
            if (value.isType<OutType>())
            {
                return value.cast<OutType>();
            }
            GMX_THROW(InvalidInputError("Invalid type of value"));
        }
        return iter->second(value);
    }

    /*! \brief
     * Adds a supported conversion.
     *
     * \tparam InType  Type to convert from.
     * \param  func    Function to convert from `InType` to `OutType`.
     */
    template<typename InType>
    void addConverter(std::function<OutType(const InType&)> func)
    {
        converters_[std::type_index(typeid(InType))] = [func](const Any& value) {
            return func(value.cast<InType>());
        };
    }
    /*! \brief
     * Adds a supported conversion from a type that can be directly cast.
     *
     * \tparam InType  Type to convert from with a simple cast.
     */
    template<typename InType>
    void addCastConversion()
    {
        converters_[std::type_index(typeid(InType))] = [](const Any& value) {
            return static_cast<OutType>(value.cast<InType>());
        };
    }

private:
    typedef std::function<OutType(const Any& value)> ConversionFunction;

    std::map<std::type_index, ConversionFunction> converters_;
};

} // namespace gmx

#endif
