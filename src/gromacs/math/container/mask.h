/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#ifndef GMX_MATH_MASK_H
#define GMX_MATH_MASK_H

#include <cmath>

#include <vector>

#endif /* end of include guard: GMX_MATH_MASK_H */

namespace gmx
{

/*! \brief Generate and apply masking to containers size.
 *
 * Container values are kept where mask is true.
 */
class Mask
{
    public:

        Mask() = default;

        template<typename T>
        Mask(const T &container, typename T::value_type threshold)
        {
            for (const auto &value : container)
            {
                mask_.push_back(value >= threshold);
            }
        }

        bool operator[](size_t i) const
        {
            return ((i >= mask_.size()) || mask_[i]);
        }

        template<typename T>
        void maskToValue(T* container, typename T::value_type value) const
        {
            auto maskIt = mask_.begin();
            for (auto containerIt = container->begin(); containerIt != container->end() && maskIt != mask_.end(); ++containerIt, ++maskIt)
            {
                if (!(*maskIt))
                {
                    *containerIt =  value;
                }
            }
        }

    private:
        std::vector<bool> mask_;
};

}
