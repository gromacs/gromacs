/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

/*! \file
 * \brief Compatibility header for functionality differences in gmxapi releases.
 *
 * Also handle the transitioning installed headers from GROMACS 2019 moving forward.
 *
 * \todo Configure for gmxapi 0.0.7, 0.0.8, GROMACS 2019, GROMACS main...
 *
 * \defgroup gmxapi_compat
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#ifndef GMXAPICOMPAT_H
#define GMXAPICOMPAT_H

#include <map>
#include <string>

#include "gmxapi/exceptions.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/gromacsfwd.h"

/*!
 * \brief Compatibility code for features that may not be in the gmxapi specification yet.
 *
 * \ingroup gmxapi_compat
 */
namespace gmxapicompat
{

/*!
 * \brief The key name provided for a key-value mapping is invalid.
 */
class KeyError : public gmxapi::BasicException<KeyError>
{
    using gmxapi::BasicException<KeyError>::BasicException;
};

/*!
 * \brief The value provided for a key-value mapping is invalid.
 */
class ValueError : public gmxapi::BasicException<ValueError>
{
    using gmxapi::BasicException<ValueError>::BasicException;
};

/*!
 * \brief Value provided for a key-value mapping is of an incompatible type.
 */
class TypeError : public gmxapi::BasicException<TypeError>
{
    using gmxapi::BasicException<TypeError>::BasicException;
};

/*!
 * \brief Static map of GROMACS MDP user input to normalized "type".
 *
 * Note that only fields present in the TPR file are named. Additional names
 * may be accepted as mdp file entries, but we cannot discern which parameter
 * name was used from inspection of the TPR file and this is an interim solution
 * that does not need to support a complete MDP file converter.
 */
std::map<std::string, gmxapi::GmxapiType> simulationParameterTypeMap();

std::map<std::string, bool t_inputrec::*>    boolParams();
std::map<std::string, int t_inputrec::*>     int32Params();
std::map<std::string, float t_inputrec::*>   float32Params();
std::map<std::string, double t_inputrec::*>  float64Params();
std::map<std::string, int64_t t_inputrec::*> int64Params();

/*!
 * \brief Static mapping of parameter names to gmxapi types for GROMACS.
 *
 * \param name MDP entry name.
 * \return enumeration value for known parameters.
 *
 * \throws gmxapi_compat::ValueError for parameters with no mapping.
 */
gmxapi::GmxapiType mdParamToType(const std::string& name);

/*!
 * \brief Facade for objects that can provide atomic data for a configuration.
 */
class StructureSource;

/*!
 * \brief Facade for objects that can provide molecular topology information for a structure.
 */
class TopologySource;

/*!
 * \brief Proxy to simulation state data.
 */
class SimulationState;

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
class GmxMdParams;

/*!
 * \brief Handle for a TPR data resource.
 *
 * Can provide StructureSource, TopologySource, GmxMdParams, and SimulationState.
 *
 * This is the type of object we allow Python clients to hold references to, though
 * we don't expose any methods to Python. Python clients should acquire access
 * to TPR file contents with read_tpr().
 *
 * \todo gmxapi C++ API should provide mechanisms for subscribing to simulation
 *       input data from various sources.
 */
class TprReadHandle;

/*!
 * \brief Open a TPR file and retrieve a handle.
 *
 * \param filename Path of file to read.
 * \return handle that may share ownership of TPR file resource.
 */
std::unique_ptr<TprReadHandle> readTprFile(const std::string& filename);

/*!
 * \brief Write a new TPR file to the filesystem with the provided contents.
 *
 * \param filename output file path
 * \param params simulation parameters
 * \param structure system structure (atomic configuration)
 * \param state simulation state
 * \param topology molecular topology
 *
 * \throws ValueError for invalid or irreconcilable input.
 */
void writeTprFile(const std::string&     filename,
                  const GmxMdParams&     params,
                  const StructureSource& structure,
                  const SimulationState& state,
                  const TopologySource&  topology);

/*!
 * \brief Get a topology source from the TPR contents collection.
 * \param handle
 * \return
 *
 * \todo replace with a helper template on T::topologySource() member function existence.
 */

std::unique_ptr<TopologySource> getTopologySource(const TprReadHandle& handle);

/*!
 * \brief Get a source of simulation state from the TPR contents collection.
 * \param handle
 * \return
 *
 * \todo template on T::simulationState() member function existence.
 */
std::unique_ptr<SimulationState> getSimulationState(const TprReadHandle& handle);

/*!
 * \brief Get a source of atomic structure from the TPR contents collection.
 * \param handle
 * \return
 */
std::unique_ptr<StructureSource> getStructureSource(const TprReadHandle& handle);

/*!
 * \brief Get an initialized parameters structure.
 * \param handle
 * \return
 */
std::unique_ptr<GmxMdParams> getMdParams(const TprReadHandle& handle);

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
int     extractParam(const GmxMdParams& params, const std::string& name, int /*unused*/);
int64_t extractParam(const GmxMdParams& params, const std::string& name, int64_t /*unused*/);
float   extractParam(const GmxMdParams& params, const std::string& name, float /*unused*/);
double  extractParam(const GmxMdParams& params, const std::string& name, double /*unused*/);

void setParam(GmxMdParams* params, const std::string& name, double value);
void setParam(GmxMdParams* params, const std::string& name, int64_t value);
// TODO: unsetParam

} // end namespace gmxapicompat

#endif // GMXAPICOMPAT_H
