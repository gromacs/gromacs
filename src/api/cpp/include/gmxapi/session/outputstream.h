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
// Created by Eric Irrgang on 5/21/18.
//

#ifndef GMXAPI_SESSION_OUTPUTSTREAM_H
#define GMXAPI_SESSION_OUTPUTSTREAM_H

#include <memory>
#include <string>

#include "gmxapi/status.h"

namespace gmxapi
{
namespace session
{

/*!
 * \brief Extensible adapter to output data streams for gmxapi client code running in Sessions.
 *
 * Code that is launched as part of a gmxapi Session receives a bundle of resources for the Session. Input and output
 * "ports" registered by the client code appear as collections of publish and visit functions with two arguments:
 * one naming the port, and one providing the data to publish or the memory into which to receive data. The type of the
 * data argument must match the registered type of the named port.
 *
 * Objects of this type are not created by API client code, but are received during session launch.
 *
 * /internal
 * I don't really want to do run-time type checking, so I think we need to just decide that the data types of graph
 * edges are specified in the API. We can provide one or two more free-form data types to allow flexibilty, and users
 * can request officially-supported data types if it makes sense. I think that a robust matrix class like Eigen provides
 * should be sufficient for binary data and we can provide a string type for arbitrary serializeable data. Key-value
 * maps and simple scalars also make sense.
 */
class OutputStream final
{
    public:
        OutputStream(const OutputStream&) = delete;
        OutputStream(OutputStream &&)     = default;
        OutputStream                   &operator=(const OutputStream &) = delete;
        OutputStream                   &operator=(OutputStream &&)      = default;
        ~OutputStream();

        // Opaque implementation class.
        class Impl;

        /*!
         * \brief Set data for a registered output stream.
         *
         * \param outputName Registered name of the output port
         * \param data data to set with the registered output handler.
         *
         * We should not use a template here to handle the different data types because the template might be expressed
         * with different munged symbol names by different compilers. But we want this interface to be extensible, so
         * we need to consider how to allow flexible types. We could wrap all data in a gmxapi::Data wrapper or something,
         * but that makes calling the set() method more cumbersome in the client code.
         *
         * What we could do, though, is to provide a template function as a helper that is compiled in the client code
         * and just makes it easier to make the correct call here. Then a gmxapi::Data wrapper wouldn't be cumbersome.
         */
        gmxapi::Status set(const std::string &outputName,
                           bool               data);
        gmxapi::Status set(const std::string &outputName,
                           double             data);

        /*!
         * \brief Wrapper to accept const C-string arguments
         *
         * \tparam T output stream data type.
         * \param outputName registered output port name.
         * \param data data to set with the registered output handler.
         */
        template<class T>
        void set(const char* outputName, T data)
        {
            this->set(std::string(outputName), data);
        }

    private:
        // opaque pointer to implementation.
        std::unique_ptr<Impl> impl_;

        // Private constructor. Objects of this type are provided by the framework and are a detail of the Context
        // implementation.
        explicit OutputStream(std::unique_ptr<Impl> &&implementation);
};

} // end namespace gmxapi::session
} //end namespace gmxapi


#endif //GMXAPI_SESSION_OUTPUTSTREAM_H
