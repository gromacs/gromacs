/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#ifndef GMXPY_MDPARAMS_H
#define GMXPY_MDPARAMS_H

/*! \file
 * \brief Compatibility header for functionality differences in gmxapi releases.
 *
 * Also handle the transitioning installed headers from GROMACS 2019 moving forward.
 *
 * \todo Configure for gmxapi 0.0.7, 0.0.8, GROMACS 2019, GROMACS master...
 *
 * \defgroup gmxapi_compat
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#include <map>
#include <memory>
#include <string>

#include "compat_exceptions.h"

struct t_inputrec;

/*!
 * \brief Compatibility code for features that may not be in gmxapi yet.
 */
namespace gmxapicompat
{


/*!
 * \brief Label the types recognized by gmxapi.
 *
 * Provide an enumeration to aid in translating data between languages, APIs,
 * and storage formats.
 *
 * \todo The spec should explicitly map these to types in APIs already used.
 * e.g. MPI, Python, numpy, GROMACS, JSON, etc.
 * \todo Actually check the size of the types.
 *
 * \see https://redmine.gromacs.org/issues/2993 for discussion.
 */
enum class GmxapiType
{
    NULLTYPE, //! Reserved
    MAP,      //! Mapping of key name (string) to a value of some MdParamType
    BOOL,     //! Boolean logical type
    INT64,    //! 64-bit integer type
    FLOAT64,  //! 64-bit float type
    STRING,   //! string with metadata
    NDARRAY,  //! multi-dimensional array with metadata
};


/*!
 * \brief Static map of GROMACS MDP user input to normalized "type".
 *
 * Note that only fields present in the TPR file are named. Additional names
 * may be accepted as mdp file entries, but we cannot discern which parameter
 * name was used from inspection of the TPR file and this is an interim solution
 * that does not need to support a complete MDP file converter.
 */
const std::map<std::string, GmxapiType> simulationParameterTypeMap();

const std::map<std::string, bool t_inputrec::*> boolParams();
const std::map<std::string, int t_inputrec::*> int32Params();
const std::map<std::string, float t_inputrec::*> float32Params();
const std::map<std::string, double t_inputrec::*> float64Params();
const std::map<std::string, int64_t t_inputrec::*> int64Params();

/*!
 * \brief Static mapping of parameter names to gmxapi types for GROMACS.
 *
 * \param name MDP entry name.
 * \return enumeration value for known parameters.
 *
 * \throws gmxapi_compat::ValueError for parameters with no mapping.
 */
GmxapiType mdParamToType(const std::string &name);

// Forward declaration for private implementation class for GmxMdParams
class GmxMdParamsImpl;

/*!
 * \brief Handle / manager for GROMACS molecular computation input parameters.
 *
 * Interface should be consistent with MDP file entries, but data maps to TPR
 * file interface. For type safety and simplicity, we don't have generic operator
 * accessors. Instead, we have templated accessors that throw exceptions when
 * there is trouble.
 *
 * When MDP input is entirely stored in a key-value tree, this class can be a
 * simple adapter or wrapper. Until then, we need a manually maintained mapping
 * of MDP entries to TPR data.
 *
 * Alternatively, we could update the infrastructure used by list_tpx to provide
 * more generic output, but our efforts may be better spent in updating the
 * infrastructure for the key-value tree input system.
 */
class GmxMdParams
{
    public:
        GmxMdParams();
        ~GmxMdParams();
        GmxMdParams(const GmxMdParams &)            = delete;
        GmxMdParams &operator=(const GmxMdParams &) = delete;
        GmxMdParams(GmxMdParams &&) noexcept;
        GmxMdParams &operator=(GmxMdParams &&) noexcept;

        std::unique_ptr<GmxMdParamsImpl> params_;
};

/*!
 * \brief A set of overloaded functions to fetch parameters of the indicated type, if possible.
 *
 * \param params Handle to a parameters structure from which to extract.
 * \param name Parameter name
 * \param (tag) type for dispatch
 *
 * Could be used for dispatch and/or some sort of templating in the future, but
 * invoked directly for now.
 */
int extractParam(const gmxapicompat::GmxMdParams &params, const std::string &name, int);
int64_t extractParam(const gmxapicompat::GmxMdParams& params, const std::string& name, int64_t);
float extractParam(const gmxapicompat::GmxMdParams &params, const std::string &name, float);
double extractParam(const gmxapicompat::GmxMdParams &params, const std::string &name, double);

void setParam(gmxapicompat::GmxMdParams* params, const std::string &name, double value);
void setParam(gmxapicompat::GmxMdParams* params, const std::string &name, int64_t value);
// TODO: unsetParam


// Anonymous namespace to confine helper function definitions to file scope.
namespace
{

bool isFloat(GmxapiType dataType)
{
    return (dataType == GmxapiType::FLOAT64);
}

bool isInt(GmxapiType dataType)
{
    return (dataType == GmxapiType::INT64);
}

}      // end anonymous namespace

}      // end namespace gmxapicompat

#endif //GMXPY_MDPARAMS_H
