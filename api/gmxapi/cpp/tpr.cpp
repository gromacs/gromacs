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
 * \brief Helper code for TPR file access.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#include "gmxapi/compat/tpr.h"

#include <cassert>
#include <cmath>
#include <cstdint>

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/real.h"

#include "gmxapi/compat/mdparams.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/gmxapicompat.h"

struct gmx_output_env_t;

using gmxapi::GmxapiType;

namespace gmxapicompat
{

class TprContents
{
public:
    explicit TprContents(const std::string& infile) :
        irInstance_{ std::make_unique<t_inputrec>() },
        mtop_{ std::make_unique<gmx_mtop_t>() },
        state_{ std::make_unique<t_state>() }
    {
        read_tpx_state(infile, irInstance_.get(), state_.get(), mtop_.get());
    }
    ~TprContents()                             = default;
    TprContents(TprContents&& source) noexcept = default;
    TprContents& operator=(TprContents&&) noexcept = default;

    /*!
     * \brief Get a reference to the input record in the TPR file.
     *
     * Note that this implementation allows different objects to share ownership
     * of the TprFile and does not provide access restrictions to prevent multiple
     * code blocks writing to the input record. This should be resolved with a
     * combination of managed access controlled handles and through better
     * management of the data structures in the TPR file. I.e. the t_inputrec is
     * not copyable, moveable, nor default constructable (at least, to produce a
     * valid record), and it does not necessarily make sense to map the library
     * data structure to the file data structure (except that we don't have another
     * way of constructing a complete and valid input record).
     *
     * \todo We can't play fast and loose with the irInstance for long...
     *
     * \return
     */
    t_inputrec& inputRecord() const
    {
        assert(irInstance_);
        return *irInstance_;
    }

    gmx_mtop_t& molecularTopology() const
    {
        assert(mtop_);
        return *mtop_;
    }

    t_state& state() const
    {
        assert(state_);
        return *state_;
    }

private:
    // These types are not moveable in GROMACS 2019, so we use unique_ptr as a
    // moveable wrapper to let TprContents be moveable.
    std::unique_ptr<t_inputrec> irInstance_;
    std::unique_ptr<gmx_mtop_t> mtop_;
    std::unique_ptr<t_state>    state_;
};

// Note: This mapping is incomplete. Hopefully we can replace it before more mapping is necessary.
// TODO: (#2993) Replace with GROMACS library header resources when available.
std::map<std::string, GmxapiType> simulationParameterTypeMap()
{
    return {
        { "integrator", GmxapiType::STRING },
        { "tinit", GmxapiType::FLOAT64 },
        { "dt", GmxapiType::FLOAT64 },
        { "nsteps", GmxapiType::INT64 },
        { "init-step", GmxapiType::INT64 },
        { "simulation-part", GmxapiType::INT64 },
        { "comm-mode", GmxapiType::STRING },
        { "nstcomm", GmxapiType::INT64 },
        { "comm-grps", GmxapiType::NDARRAY }, // Note: we do not have processing for this yet.
        { "bd-fric", GmxapiType::FLOAT64 },
        { "ld-seed", GmxapiType::INT64 },
        { "emtol", GmxapiType::FLOAT64 },
        { "emstep", GmxapiType::FLOAT64 },
        { "niter", GmxapiType::INT64 },
        { "fcstep", GmxapiType::FLOAT64 },
        { "nstcgsteep", GmxapiType::INT64 },
        { "nbfgscorr", GmxapiType::INT64 },
        { "rtpi", GmxapiType::FLOAT64 },
        { "nstxout", GmxapiType::INT64 },
        { "nstvout", GmxapiType::INT64 },
        { "nstfout", GmxapiType::INT64 },
        { "nstlog", GmxapiType::INT64 },
        { "nstcalcenergy", GmxapiType::INT64 },
        { "nstenergy", GmxapiType::INT64 },
        { "nstxout-compressed", GmxapiType::INT64 },
        { "compressed-x-precision", GmxapiType::FLOAT64 },
        { "cutoff-scheme", GmxapiType::STRING },
        { "nstlist", GmxapiType::INT64 },
        { "pbc", GmxapiType::STRING },
        { "periodic-molecules", GmxapiType::BOOL },
        //            TBD

    };
}

/*
 * TODO: Visitor for predetermined known types.
 *
 * Development sequence:
 * 1. map pointers
 * 2. map setters ()
 * 3. template the Visitor setter for compile-time extensibility of type and to prune incompatible types.
 * 4. switch to Variant type for handling (setter templated on caller input)
 * 5. switch to Variant type for input as well? (Variant in public API?)
 */

std::map<std::string, bool t_inputrec::*> boolParams()
{
    return {
        { "periodic-molecules", &t_inputrec::bPeriodicMols },
        //            ...
    };
}

std::map<std::string, int t_inputrec::*> int32Params()
{
    return {
        { "simulation-part", &t_inputrec::simulation_part },
        { "nstcomm", &t_inputrec::nstcomm },
        { "niter", &t_inputrec::niter },
        { "nstcgsteep", &t_inputrec::nstcgsteep },
        { "nbfgscorr", &t_inputrec::nbfgscorr },
        { "nstxout", &t_inputrec::nstxout },
        { "nstvout", &t_inputrec::nstvout },
        { "nstfout", &t_inputrec::nstfout },
        { "nstlog", &t_inputrec::nstlog },
        { "nstcalcenergy", &t_inputrec::nstcalcenergy },
        { "nstenergy", &t_inputrec::nstenergy },
        { "nstxout-compressed", &t_inputrec::nstxout_compressed },
        { "nstlist", &t_inputrec::nstlist },
        //            ...
    };
}

namespace
{

// Provide a helper to disambiguate `real` typed inputrec values.

// Primary template returns empty map.
template<typename RealT>
std::map<std::string, RealT t_inputrec::*> compatibleRealParams()
{
    return {};
}

// Specialize for the compile-time typedef of `real` to get the inputrec parameters
// compatible with the template parameter whose type is not known until compile time.
template<>
std::map<std::string, real t_inputrec::*> compatibleRealParams()
{
    return {
        { "bd-fric", &t_inputrec::bd_fric },
        { "emtol", &t_inputrec::em_tol },
        { "emstep", &t_inputrec::em_stepsize },
        { "fcstep", &t_inputrec::fc_stepsize },
        { "rtpi", &t_inputrec::rtpi },
        { "compressed-x-precision", &t_inputrec::x_compression_precision },
        //            ...

    };
}

} // namespace

std::map<std::string, float t_inputrec::*> float32Params()
{
    return compatibleRealParams<float>();
}

std::map<std::string, double t_inputrec::*> float64Params()
{
    static const std::map<std::string, double t_inputrec::*> explicitDoubles = {
        { "dt", &t_inputrec::delta_t }, { "tinit", &t_inputrec::init_t },
        //            ...

    };

    // Initialize map to be returned with any compile-time-only doubles.
    auto fullMap = compatibleRealParams<double>();

    // Get the explicitly `double` parameters.
    for (const auto& item : explicitDoubles)
    {
        fullMap.emplace(item);
    }

    return fullMap;
}

std::map<std::string, int64_t t_inputrec::*> int64Params()
{
    return {
        { "nsteps", &t_inputrec::nsteps },
        { "init-step", &t_inputrec::init_step },
        { "ld-seed", &t_inputrec::ld_seed },
        //            ...

    };
}

/*!
 * \brief Static mapping of parameter names to gmxapi types for GROMACS 2019.
 *
 * \param name MDP entry name.
 * \return enumeration value for known parameters.
 *
 * \throws gmxapi_compat::ValueError for parameters with no mapping.
 */
GmxapiType mdParamToType(const std::string& name)
{
    const auto staticMap = simulationParameterTypeMap();
    auto       entry     = staticMap.find(name);
    if (entry == staticMap.end())
    {
        throw ValueError("Named parameter has unknown type mapping.");
    }
    return entry->second;
}


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
class GmxMdParamsImpl final
{
public:
    /*!
     * \brief Create an initialized but empty parameters structure.
     *
     * Parameter keys are set at construction, but all values are empty. This
     * allows the caller to check for valid parameter names or their types,
     * while allowing the consuming code to know which parameters were explicitly
     * set by the caller.
     *
     * To load values from a TPR file, see getMdParams().
     */
    GmxMdParamsImpl();

    explicit GmxMdParamsImpl(std::shared_ptr<TprContents> tprContents);

    /*!
     * \brief Get the current list of keys.
     *
     * \return
     */
    std::vector<std::string> keys() const
    {
        std::vector<std::string> keyList;
        for (auto&& entry : int64Params_)
        {
            keyList.emplace_back(entry.first);
        }
        for (auto&& entry : intParams_)
        {
            keyList.emplace_back(entry.first);
        }
        for (auto&& entry : floatParams_)
        {
            keyList.emplace_back(entry.first);
        }
        for (auto&& entry : float64Params_)
        {
            keyList.emplace_back(entry.first);
        }
        return keyList;
    };

    template<typename T>
    T extract(const std::string& /* key */) const
    {
        auto value = T();
        // should be an APIError
        throw TypeError("unhandled type");
    }

    void set(const std::string& key, const int64_t& value)
    {
        if (int64Params_.find(key) != int64Params_.end())
        {
            int64Params_[key] = std::make_pair(value, true);

            if (source_)
            {
                auto memberPointer                    = int64Params().at(key);
                source_->inputRecord().*memberPointer = value;
            }
        }
        else if (intParams_.find(key) != intParams_.end())
        {
            // TODO: check whether value is too large?
            intParams_[key] = std::make_pair(static_cast<int>(value), true);

            if (source_)
            {
                auto memberPointer                    = int32Params().at(key);
                source_->inputRecord().*memberPointer = value;
            }
        }
        else
        {
            throw KeyError("Named parameter is incompatible with integer type value.");
        }
    };

    void set(const std::string& key, const double& value)
    {
        if (float64Params_.find(key) != float64Params_.end())
        {
            float64Params_[key] = std::make_pair(value, true);

            if (source_)
            {
                auto memberPointer                    = float64Params().at(key);
                source_->inputRecord().*memberPointer = value;
            }
        }
        else if (floatParams_.find(key) != floatParams_.end())
        {
            // TODO: check whether value is too large?
            floatParams_[key] = std::make_pair(static_cast<float>(value), true);

            if (source_)
            {
                auto memberPointer                    = float32Params().at(key);
                source_->inputRecord().*memberPointer = static_cast<float>(value);
            }
        }
        else
        {
            throw KeyError("Named parameter is incompatible with floating point type value.");
        }
    };

    TprReadHandle getSource() const
    {
        // Note: might return a null handle. Need to decide what that means and how to address it.
        return TprReadHandle(source_);
    }

private:
    // Hold the settable parameters and whether or not they have been set.
    // TODO: update to gmxapi named types?
    // TODO: update to gmx::compat::optional now that this file is in the GROMACS source.
    std::map<std::string, std::pair<int64_t, bool>> int64Params_;
    std::map<std::string, std::pair<int, bool>>     intParams_;
    std::map<std::string, std::pair<float, bool>>   floatParams_;
    std::map<std::string, std::pair<double, bool>>  float64Params_;

    /*! \brief Shared ownership of a pack of TPR data.
     *
     * This is a non-normative way to retain access to gmxapi resources.
     * \todo Subscribe to a Context-managed resource.
     */
    std::shared_ptr<TprContents> source_;
};

void setParam(gmxapicompat::GmxMdParams* params, const std::string& name, double value)
{
    assert(params != nullptr);
    assert(params->params_ != nullptr);
    params->params_->set(name, value);
}

void setParam(gmxapicompat::GmxMdParams* params, const std::string& name, int64_t value)
{
    assert(params != nullptr);
    assert(params->params_ != nullptr);
    params->params_->set(name, value);
}

template<typename ParamsContainerT, typename Mapping>
static void updateParamsContainer(ParamsContainerT* params, const TprContents& source, const Mapping& map)
{
    for (const auto& definition : map)
    {
        const auto& key           = definition.first;
        auto        memberPointer = definition.second;
        auto&       irInstance    = source.inputRecord();
        auto        fileValue     = irInstance.*memberPointer;
        (*params)[key]            = std::make_pair(fileValue, true);
    }
}

/*!
 * \brief A GmxMdParams implementation that depends on TPR files.
 *
 * \param tprContents
 */
GmxMdParamsImpl::GmxMdParamsImpl(std::shared_ptr<gmxapicompat::TprContents> tprContents) :
    source_{ std::move(tprContents) }
{
    if (source_)
    {
        updateParamsContainer(&int64Params_, *source_, int64Params());
        updateParamsContainer(&intParams_, *source_, int32Params());
        updateParamsContainer(&floatParams_, *source_, float32Params());
        updateParamsContainer(&float64Params_, *source_, float64Params());
    }
}

GmxMdParamsImpl::GmxMdParamsImpl() : GmxMdParamsImpl(nullptr) {}

template<>
int GmxMdParamsImpl::extract<int>(const std::string& key) const
{
    const auto& params = intParams_;
    const auto& entry  = params.find(key);
    if (entry == params.cend())
    {
        throw KeyError("Parameter of the requested name and type not defined.");
    }
    else if (!entry->second.second)
    {
        // TODO: handle invalid and unset parameters differently.
        throw KeyError("Parameter of the requested name not set.");
    }
    else
    {
        return entry->second.first;
    }
}

template<>
int64_t GmxMdParamsImpl::extract<int64_t>(const std::string& key) const
{
    const auto& params = int64Params_;
    const auto& entry  = params.find(key);
    if (entry == params.cend())
    {
        throw KeyError("Parameter of the requested name and type not defined.");
    }
    else if (!entry->second.second)
    {
        // TODO: handle invalid and unset parameters differently.
        throw KeyError("Parameter of the requested name not set.");
    }
    else
    {
        return entry->second.first;
    }
}
template<>
float GmxMdParamsImpl::extract<float>(const std::string& key) const
{
    const auto& params = floatParams_;
    const auto& entry  = params.find(key);
    if (entry == params.cend())
    {
        throw KeyError("Parameter of the requested name and type not defined.");
    }
    else if (!entry->second.second)
    {
        // TODO: handle invalid and unset parameters differently.
        throw KeyError("Parameter of the requested name not set.");
    }
    else
    {
        return entry->second.first;
    }
}
template<>
double GmxMdParamsImpl::extract<double>(const std::string& key) const
{
    const auto& params = float64Params_;
    const auto& entry  = params.find(key);
    if (entry == params.cend())
    {
        throw KeyError("Parameter of the requested name and type not defined.");
    }
    else if (!entry->second.second)
    {
        // TODO: handle invalid and unset parameters differently.
        throw KeyError("Parameter of the requested name not set.");
    }
    else
    {
        return entry->second.first;
    }
}


int extractParam(const GmxMdParams& params, const std::string& name, int /*unused*/)
{
    assert(params.params_);
    return params.params_->extract<int>(name);
}

int64_t extractParam(const GmxMdParams& params, const std::string& name, int64_t /*unused*/)
{
    assert(params.params_);
    int64_t value{};
    // Allow fetching both known integer types.
    try
    {
        value = params.params_->extract<int>(name);
    }
    catch (const KeyError& error)
    {
        // If not found as a regular int, check for int64.
        try
        {
            value = params.params_->extract<int64_t>(name);
        }
        catch (const KeyError& error64)
        {
            throw KeyError("Parameter of the requested name not set.");
        }
    }
    // Any other exceptions propagate out.
    return value;
}

float extractParam(const GmxMdParams& params, const std::string& name, float /*unused*/)
{
    assert(params.params_);
    return params.params_->extract<float>(name);
}

double extractParam(const GmxMdParams& params, const std::string& name, double /*unused*/)
{
    assert(params.params_);
    double value{};
    // Allow fetching both single and double precision.
    try
    {
        value = params.params_->extract<double>(name);
    }
    catch (const KeyError& errorDouble)
    {
        // If not found as a double precision value, check for single-precision.
        try
        {
            value = params.params_->extract<float>(name);
        }
        catch (const KeyError& errorFloat)
        {
            throw KeyError("Parameter of the requested name not set.");
        }
    }
    // Any other exceptions propagate out.
    return value;
}

std::vector<std::string> keys(const GmxMdParams& params)
{
    return params.params_->keys();
}


std::unique_ptr<TprReadHandle> readTprFile(const std::string& filename)
{
    auto tprfile = gmxapicompat::TprContents(filename);
    auto handle  = std::make_unique<gmxapicompat::TprReadHandle>(std::move(tprfile));
    return handle;
}

std::unique_ptr<GmxMdParams> getMdParams(const TprReadHandle& handle)
{
    auto tprfile = handle.get();
    // TODO: convert to exception / decide whether null handles are allowed.
    assert(tprfile);
    auto params_impl = std::make_unique<GmxMdParamsImpl>(tprfile);
    auto params      = std::make_unique<GmxMdParams>(std::move(params_impl));
    return params;
}

std::unique_ptr<TopologySource> getTopologySource(const TprReadHandle& handle)
{
    auto source      = std::make_unique<TopologySource>();
    source->tprFile_ = handle.get();
    return source;
}

std::unique_ptr<SimulationState> getSimulationState(const TprReadHandle& handle)
{
    auto source      = std::make_unique<SimulationState>();
    source->tprFile_ = handle.get();
    return source;
}

std::unique_ptr<StructureSource> getStructureSource(const TprReadHandle& handle)
{
    auto source      = std::make_unique<StructureSource>();
    source->tprFile_ = handle.get();
    return source;
}

TprReadHandle::TprReadHandle(std::shared_ptr<TprContents> tprFile) :
    tprContents_{ std::move(tprFile) }
{
}

TprReadHandle getSourceFileHandle(const GmxMdParams& params)
{
    return params.params_->getSource();
}

void writeTprFile(const std::string&     filename,
                  const GmxMdParams&     params,
                  const StructureSource& structure,
                  const SimulationState& state,
                  const TopologySource&  topology)
{
    assert(params.params_);
    // The only way we can check for consistent input right now is to make sure
    // it all comes from the same file.
    if (structure.tprFile_.get() != state.tprFile_.get()
        || state.tprFile_.get() != topology.tprFile_.get()
        || topology.tprFile_.get() != params.params_->getSource().get().get()
        || params.params_->getSource().get().get() != structure.tprFile_.get())
    {
        throw ValueError(
                "writeTprFile does not yet know how to reconcile data from different TPR file "
                "sources.");
    }

    const auto tprFileHandle = params.params_->getSource();
    const auto tprFile       = tprFileHandle.get();
    assert(tprFile);
    const auto& inputRecord   = tprFile->inputRecord();
    const auto& writeState    = tprFile->state();
    const auto& writeTopology = tprFile->molecularTopology();
    write_tpx_state(filename.c_str(), &inputRecord, &writeState, writeTopology);
}

TprReadHandle::TprReadHandle(TprContents&& tprFile) :
    TprReadHandle{ std::make_shared<TprContents>(std::move(tprFile)) }
{
}

std::shared_ptr<TprContents> TprReadHandle::get() const
{
    return tprContents_;
}

// defaulted here to delay definition until after member types are defined.
TprReadHandle::~TprReadHandle() = default;

GmxMdParams::~GmxMdParams() = default;

GmxMdParams::GmxMdParams() : params_{ std::make_unique<GmxMdParamsImpl>() } {}

GmxMdParams::GmxMdParams(GmxMdParams&&) noexcept = default;

GmxMdParams& GmxMdParams::operator=(GmxMdParams&&) noexcept = default;

GmxMdParams::GmxMdParams(std::unique_ptr<GmxMdParamsImpl>&& impl)
{
    // We use swap instead of move construction so that we don't have
    // to worry about the restrictions on Deleters.
    // Ref: https://en.cppreference.com/w/cpp/memory/unique_ptr/unique_ptr
    params_.swap(impl);
}

// maybe this should return a handle to the new file?
bool copy_tprfile(const gmxapicompat::TprReadHandle& input, const std::string& outFile)
{
    if (!input.get())
    {
        return false;
    }
    gmxapicompat::writeTprFile(outFile,
                               *gmxapicompat::getMdParams(input),
                               *gmxapicompat::getStructureSource(input),
                               *gmxapicompat::getSimulationState(input),
                               *gmxapicompat::getTopologySource(input));
    return true;
}

bool rewrite_tprfile(const std::string& inFile, const std::string& outFile, double endTime)
{
    t_inputrec irInstance;
    gmx_mtop_t mtop;
    t_state    state;
    read_tpx_state(inFile, &irInstance, &state, &mtop);

    /* set program name, command line, and default values for output options */
    gmx_output_env_t* oenv;
    gmx::TimeUnit     timeUnit = gmx::TimeUnit::Default;
    bool              bView{ false }; // argument that says we don't want to view graphs.
    output_env_init(&oenv, gmx::getProgramContext(), timeUnit, bView, XvgFormat::Xmgrace, 0);

    double run_t = irInstance.init_step * irInstance.delta_t + irInstance.init_t;

    irInstance.nsteps = std::lround((endTime - run_t) / irInstance.delta_t);

    write_tpx_state(outFile.c_str(), &irInstance, &state, mtop);

    return true;
}

} // end namespace gmxapicompat
