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
//
// Created by Eric Irrgang on 11/14/17.
//

#ifndef GMXAPI_STATUS_H
#define GMXAPI_STATUS_H
/*! \file
 * \brief Declare the base class for Status returned by API calls.
 *
 * \ingroup gmxapi
 */

#include <memory>

namespace gmxapi
{
/*! \brief Container for results of API operations.
 *
 * \internal
 * I'm leaning towards saying this should not be derivable, but that it
 * may contain one or more additional objects, perhaps including exceptions
 * or chains of status / exceptions. Maybe it is a stack. Maybe all
 * API objects should have a Status member that can accumulate Status
 * objects of child objects/operations.
 */
class Status final
{
    public:
        /*!
         * \brief Default constructor.
         */
        Status();
        /*!
         * \brief Copy constructor
         * \param status
         */
        Status(const Status &status);
        /*!
         * \brief Move constructor.
         * \param status
         */
        Status(Status && status) noexcept;
        /*!
         * \brief Copy assignment operator.
         * \param status
         * \return reference to lhs.
         */
        Status &operator=(const Status &status);
        /*!
         * \brief Move assignment operator.
         * \param status
         * \return reference to lhs
         */
        Status &operator=(Status &&status) noexcept;
        /*!
         * \brief Converting assignment operator.
         * \param success
         * \return reference to lhs
         */
        Status &operator=(bool success);

        /*!
         * \brief Converting constructor.
         * \param success
         */
        explicit Status(bool success);

        /*!
         * \brief non-virtual destructor
         *
         * Do not inherit from this class.
         */
        ~Status();
        /*
         * \brief Check success status.
         *
         * \return true if the operation described was successful.
         */
        bool success() const;
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

}      // end namespace gmxapi

#endif //GMXAPI_STATUS_H
