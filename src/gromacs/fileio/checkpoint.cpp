/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2008- The GROMACS Authors
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

/* The source code in this file should be thread-safe.
   Please keep it that way. */
#include "gmxpre.h"

#include "checkpoint.h"

#include <cerrno>
#include <cinttypes>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <type_traits>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/xdr_serializer.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/awh_correlation_history.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/edsamhistory.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/swaphistory.h"
#include "gromacs/modularsimulator/modularsimulator.h"
#include "gromacs/serialization/iserializer.h"
#include "gromacs/serialization/keyvaluetreeserializer.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/int64_to_int.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/txtdump.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vecdump.h"
#include "gromacs/utility/vectypes.h"

#include "buildinfo.h"


enum class ChannelHistory : int;

#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

using gmx::ISerializer;
using gmx::XdrSerializer;

namespace gmx
{

template<typename ValueType>
void readKvtCheckpointValue(compat::not_null<ValueType*> value,
                            std::string_view             name,
                            std::string_view             identifier,
                            const KeyValueTreeObject&    kvt)
{
    std::string key(identifier);
    key += "-";
    key += name;
    if (!kvt.keyExists(key))
    {
        std::string errorMessage = "Cannot read requested checkpoint value " + key + " .";
        GMX_THROW(InternalError(errorMessage));
    }
    *value = kvt[key].cast<ValueType>();
}

template void readKvtCheckpointValue(compat::not_null<std::int64_t*> value,
                                     std::string_view                name,
                                     std::string_view                identifier,
                                     const KeyValueTreeObject&       kvt);
template void readKvtCheckpointValue(compat::not_null<real*>   value,
                                     std::string_view          name,
                                     std::string_view          identifier,
                                     const KeyValueTreeObject& kvt);

template<typename ValueType>
void writeKvtCheckpointValue(const ValueType&          value,
                             std::string_view          name,
                             std::string_view          identifier,
                             KeyValueTreeObjectBuilder kvtBuilder)
{
    std::string key(identifier);
    key += "-";
    key += name;
    kvtBuilder.addValue<ValueType>(key, value);
}

template void writeKvtCheckpointValue(const std::int64_t&       value,
                                      std::string_view          name,
                                      std::string_view          identifier,
                                      KeyValueTreeObjectBuilder kvtBuilder);
template void writeKvtCheckpointValue(const real&               value,
                                      std::string_view          name,
                                      std::string_view          identifier,
                                      KeyValueTreeObjectBuilder kvtBuilder);


} // namespace gmx

/*! \brief Version number of the file format written to checkpoint
 * files by this version of the code.
 *
 * cpt_version should normally only be changed, via adding a new field
 * to cptv enumeration, when the header or footer format changes.
 *
 * The state data format itself is backward and forward compatible.
 * But old code can not read a new entry that is present in the file
 * (but can read a new format when new entries are not present).
 *
 * The cpt_version increases whenever the file format in the main
 * development branch changes, due to an extension of the CheckPointVersion
 * enum (see src/gromacs/fileio/checkpoint.h).
 * Backward compatibility for reading old run input files is maintained
 * by checking this version number against that of the file and then using
 * the correct code path. */
static constexpr CheckPointVersion cpt_version = CheckPointVersion::CurrentVersion;

const char* enumValueToString(StateEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StateEntry, const char*> stateEntryNames = {
        "FE-lambda",
        "box",
        "box-rel",
        "box-v",
        "pres_prev",
        "nosehoover-xi",
        "thermostat-integral",
        "x",
        "v",
        "sdx-unsupported",
        "CGp",
        "LD-rng-unsupported",
        "LD-rng-i-unsupported",
        "disre_initf",
        "disre_rm3tav",
        "orire_initf",
        "orire_Dtav",
        "svir_prev",
        "nosehoover-vxi",
        "v_eta",
        "vol0",
        "nhpres_xi",
        "nhpres_vxi",
        "fvir_prev",
        "fep_state",
        "MC-rng-unsupported",
        "MC-rng-i-unsupported",
        "barostat-integral"
    };
    return stateEntryNames[enumValue];
}

enum class StateKineticEntry : int
{
    EkinNumber,
    EkinHalfStep,
    DEkinDLambda,
    Mvcos,
    EkinFullStep,
    EkinHalfStepOld,
    EkinNoseHooverScaleFullStep,
    EkinNoseHooverScaleHalfStep,
    VelocityScale,
    EkinTotal,
    Count
};

static const char* enumValueToString(StateKineticEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StateKineticEntry, const char*> stateKineticEntryNames = {
        "Ekin_n",    "Ekinh",          "dEkindlambda",   "mv_cos",     "Ekinf",
        "Ekinh_old", "EkinScaleF_NHC", "EkinScaleH_NHC", "Vscale_NHC", "Ekin_Total"
    };
    return stateKineticEntryNames[enumValue];
}

enum class StateEnergyEntry : int
{
    N,
    Aver,
    Sum,
    NumSum,
    SumSim,
    NumSumSim,
    NumSteps,
    NumStepsSim,
    DeltaHNN,
    DeltaHList,
    DeltaHStartTime,
    DeltaHStartLambda,
    Count
};

enum class StatePullEntry : int
{
    NumCoordinates,
    NumGroups,
    NumValuesInXSum,
    NumValuesInFSum,
    Count
};

enum class StatePullCoordEntry : int
{
    ValueReferenceSum,
    ValueSum,
    DR01Sum,
    DR23Sum,
    DR45Sum,
    FScalarSum,
    DynaxSum,
    Count
};

static const char* enumValueToString(StatePullCoordEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StatePullCoordEntry, const char*> statePullCoordEntryNames = {
        "reference-sum", "sum", "dr01-sum", "dr23-sum", "dr45-sum", "fscal-sum", "dynax-sum"
    };
    return statePullCoordEntryNames[enumValue];
}

enum class StatePullGroupEntry : int
{
    XSum,
    Count
};

static const char* enumValueToString(StatePullGroupEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StatePullGroupEntry, const char*> statePullGroupEntryNames = {
        "coordinate-sum"
    };
    return statePullGroupEntryNames[enumValue];
}

static const char* enumValueToString(StateEnergyEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StateEnergyEntry, const char*> stateEnergyEntryNames = {
        "energy_n",
        "energy_aver",
        "energy_sum",
        "energy_nsum",
        "energy_sum_sim",
        "energy_nsum_sim",
        "energy_nsteps",
        "energy_nsteps_sim",
        "energy_delta_h_nn",
        "energy_delta_h_list",
        "energy_delta_h_start_time",
        "energy_delta_h_start_lambda"
    };
    return stateEnergyEntryNames[enumValue];
}

static const char* enumValueToString(StatePullEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StatePullEntry, const char*> statePullEntryNames = {
        "pullhistory_numcoordinates",
        "pullhistory_numgroups",
        "pullhistory_numvaluesinxsum",
        "pullhistory_numvaluesinfsum"
    };
    return statePullEntryNames[enumValue];
}

/* free energy history variables -- need to be preserved over checkpoint */
enum class StateFepEntry : int
{
    IsEquilibrated,
    NumAtLambdaStats,
    NumAtLambdaEquil,
    WangLandauHistogram,
    WangLandauDelta,
    SumWeights,
    SumDG,
    SumMinVar,
    SumVar,
    Accump,
    Accumm,
    Accump2,
    Accumm2,
    Tij,
    TijEmp,
    Count
};

//! free energy history names
static const char* enumValueToString(StateFepEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StateFepEntry, const char*> stateFepEntryNames = {
        "bEquilibrated",
        "Number at State Stats",
        "Number at State Equil",
        "Wang-Landau Histogram",
        "Wang-Landau Delta",
        "Weights",
        "Free Energies",
        "minvar",
        "variance",
        "accumulated_plus",
        "accumulated_minus",
        "accumulated_plus_2",
        "accumulated_minus_2",
        "Tij",
        "Tij_empirical"
    };
    return stateFepEntryNames[enumValue];
}

//! AWH biasing history variables
enum class StateAwhEntry : int
{
    InInitial,
    EquilibrateHistogram,
    HistogramSize,
    NumPoints,
    CoordPoint,
    UmbrellaGridPoint,
    UpdateList,
    LogScaledSampleWeight,
    NumUpdates,
    ForceCorrelationGrid,
    Count
};

static const char* enumValueToString(StateAwhEntry enumValue)
{
    static constexpr gmx::EnumerationArray<StateAwhEntry, const char*> stateAwhEntryNames = {
        "awh_in_initial", "awh_equilibrateHistogram", "awh_histsize",   "awh_npoints",
        "awh_coordpoint", "awh_umbrellaGridpoint",    "awh_updatelist", "awh_logScaledSampleWeight",
        "awh_numupdates", "awh_forceCorrelationGrid"
    };
    return stateAwhEntryNames[enumValue];
}

enum class StatePullCommunicationEntry : int
{
    PreviousStepCom,
    Count
};

//! Higher level vector element type, only used for formatting checkpoint dumps
enum class CptElementType
{
    integer,  //!< integer
    realnum,  //!< float or double, not linked to precision of type real
    real3,    //!< float[3] or double[3], not linked to precision of type real
    matrix3x3 //!< float[3][3] or double[3][3], not linked to precision of type real
};

static void cp_warning(FILE* fp)
{
    fprintf(fp, "\nWARNING: Checkpoint file is corrupted or truncated\n\n");
}

[[noreturn]] static void cp_error()
{
    gmx_fatal(FARGS, "Checkpoint file corrupted/truncated, or maybe you are out of disk space?");
}

static void do_cpt_string_err(XdrSerializer* serializer, const char* desc, gmx::ArrayRef<char> s, FILE* list)
{
    char* data = s.data();
    if (xdr_string(serializer->xdr(), &data, s.size()) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, data);
    }
}

static int do_cpt_int(XdrSerializer* serializer, const char* desc, int* i, FILE* list)
{
    if (xdr_int(serializer->xdr(), i) == 0)
    {
        return -1;
    }
    if (list)
    {
        fprintf(list, "%s = %d\n", desc, *i);
    }
    return 0;
}

static int do_cpt_u_chars(XdrSerializer* serializer, const char* desc, int n, unsigned char* i, FILE* list)
{
    if (list)
    {
        fprintf(list, "%s = ", desc);
    }
    bool_t res = 1;
    for (int j = 0; j < n && res; j++)
    {
        res &= xdr_u_char(serializer->xdr(), &i[j]);
        if (list)
        {
            fprintf(list, "%02x", i[j]);
        }
    }
    if (list)
    {
        fprintf(list, "\n");
    }
    if (res == 0)
    {
        return -1;
    }

    return 0;
}

template<typename EnumType>
static int do_cpt_enum_as_int(XdrSerializer* serializer, const char* desc, EnumType* enumValue, FILE* list)
{
    static_assert(std::is_same_v<std::underlying_type_t<EnumType>, int>,
                  "Only enums with underlying type int are supported.");
    auto castedValue = static_cast<int>(*enumValue);
    if (xdr_int(serializer->xdr(), &castedValue) == 0)
    {
        return -1;
    }
    *enumValue = static_cast<EnumType>(castedValue);
    if (list)
    {
        fprintf(list, "%s = %d\n", desc, castedValue);
    }
    return 0;
}

template<typename EnumType>
static int do_cpt_n_enum_as_int(XdrSerializer* serializer, const char* desc, int n, EnumType* enumValue, FILE* list)
{
    bool_t res = 1;
    for (int j = 0; j < n && res; j++)
    {
        res &= do_cpt_enum_as_int<EnumType>(serializer, desc, &enumValue[j], list);
    }
    if (res == 0)
    {
        return -1;
    }

    return 0;
}

static void do_cpt_int_err(XdrSerializer* serializer, const char* desc, int* i, FILE* list)
{
    if (do_cpt_int(serializer, desc, i, list) < 0)
    {
        cp_error();
    }
}

static void do_cpt_bool_err(XdrSerializer* serializer, const char* desc, bool* b, FILE* list)
{
    int i = static_cast<int>(*b);

    if (do_cpt_int(serializer, desc, &i, list) < 0)
    {
        cp_error();
    }

    *b = (i != 0);
}

static void do_cpt_step_err(XdrSerializer* serializer, const char* desc, int64_t* i, FILE* list)
{
    char buf[STEPSTRSIZE];

    if (xdr_int64(serializer->xdr(), i) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, gmx_step_str(*i, buf));
    }
}

static void do_cpt_double_err(XdrSerializer* serializer, const char* desc, double* f, FILE* list)
{
    if (xdr_double(serializer->xdr(), f) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %f\n", desc, *f);
    }
}

static void do_cpt_real_err(XdrSerializer* serializer, real* f)
{
#if GMX_DOUBLE
    bool_t res = xdr_double(serializer->xdr(), f);
#else
    bool_t res = xdr_float(serializer->xdr(), f);
#endif
    if (res == 0)
    {
        cp_error();
    }
}

static void do_cpt_n_rvecs_err(XdrSerializer* serializer, const char* desc, int n, rvec f[], FILE* list)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            do_cpt_real_err(serializer, &f[i][j]);
        }
    }

    if (list)
    {
        pr_rvecs(list, 0, desc, f, n);
    }
}

template<typename T>
struct xdr_type
{
};

template<>
struct xdr_type<int>
{
    static const XdrDataType value = XdrDataType::Int;
};

template<>
struct xdr_type<float>
{
    static const XdrDataType value = XdrDataType::Float;
};

template<>
struct xdr_type<double>
{
    static const XdrDataType value = XdrDataType::Double;
};

//! \brief Returns size in byte of an XdrDataType
static inline unsigned int sizeOfXdrType(XdrDataType xdrType)
{
    switch (xdrType)
    {
        case XdrDataType::Int: return sizeof(int);
        case XdrDataType::Float: return sizeof(float);
        case XdrDataType::Double: return sizeof(double);
        default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
    }

    return 0;
}

//! \brief Returns the XDR process function for i/o of an XDR type
static inline xdrproc_t xdrProc(XdrDataType xdrType)
{
    switch (xdrType)
    {
        case XdrDataType::Int: return reinterpret_cast<xdrproc_t>(xdr_int);
        case XdrDataType::Float: return reinterpret_cast<xdrproc_t>(xdr_float);
        case XdrDataType::Double: return reinterpret_cast<xdrproc_t>(xdr_double);
        default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
    }

    return nullptr;
}

/*! \brief Lists or only reads an xdr vector from checkpoint file
 *
 * When list!=NULL reads and lists the \p nf vector elements of type \p xdrType.
 * The header for the print is set by \p part and \p ecpt.
 * The formatting of the printing is set by \p cptElementType.
 * When list==NULL only reads the elements.
 */
template<typename Enum>
static bool_t listXdrVector(XdrSerializer* serializer, Enum ecpt, int nf, XdrDataType xdrType, FILE* list, CptElementType cptElementType)
{
    bool_t res = 0;

    const unsigned int elemSize = sizeOfXdrType(xdrType);
    std::vector<char>  data(nf * elemSize);
    res = xdr_vector(serializer->xdr(), data.data(), nf, elemSize, xdrProc(xdrType));

    if (list != nullptr)
    {
        switch (xdrType)
        {
            case XdrDataType::Int:
                pr_ivec(list, 0, enumValueToString(ecpt), reinterpret_cast<const int*>(data.data()), nf, TRUE);
                break;
            case XdrDataType::Float:
#if !GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    pr_rvecs(list, 0, enumValueToString(ecpt), reinterpret_cast<const rvec*>(data.data()), nf / 3);
                }
                else
#endif
                {
                    /* Note: With double precision code dumping a single precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    pr_fvec(list, 0, enumValueToString(ecpt), reinterpret_cast<const float*>(data.data()), nf, TRUE);
                }
                break;
            case XdrDataType::Double:
#if GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    pr_rvecs(list, 0, enumValueToString(ecpt), reinterpret_cast<const rvec*>(data.data()), nf / 3);
                }
                else
#endif
                {
                    /* Note: With single precision code dumping a double precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    pr_dvec(list, 0, enumValueToString(ecpt), reinterpret_cast<const double*>(data.data()), nf, TRUE);
                }
                break;
            default: GMX_RELEASE_ASSERT(false, "Data type not implemented for listing");
        }
    }

    return res;
}

//! \brief Convert a double array, typed char*, to float
gmx_unused static void convertArrayRealPrecision(const char* c, float* v, int n)
{
    const double* d = reinterpret_cast<const double*>(c);
    for (int i = 0; i < n; i++)
    {
        v[i] = static_cast<float>(d[i]);
    }
}

//! \brief Convert a float array, typed char*, to double
static void convertArrayRealPrecision(const char* c, double* v, int n)
{
    const float* f = reinterpret_cast<const float*>(c);
    for (int i = 0; i < n; i++)
    {
        v[i] = static_cast<double>(f[i]);
    }
}

//! \brief Generate an error for trying to convert to integer
static void convertArrayRealPrecision(const char gmx_unused* c, int gmx_unused* v, int gmx_unused n)
{
    GMX_RELEASE_ASSERT(false,
                       "We only expect type mismatches between float and double, not integer");
}

/*! \brief Low-level routine for reading/writing a vector of reals from/to file.
 *
 * This is the only routine that does the actually i/o of real vector,
 * all other routines are intermediate level routines for specific real
 * data types, calling this routine.
 * Currently this routine is (too) complex, since it handles both real *
 * and std::vector<real>. Using real * is deprecated and this routine
 * will simplify a lot when only std::vector needs to be supported.
 *
 * The routine is generic to vectors with different allocators,
 * because as part of reading a checkpoint there are vectors whose
 * size is not known until reading has progressed far enough, so a
 * resize method must be called.
 *
 * When not listing, we use either v or vector, depending on which is !=NULL.
 * If nval >= 0, nval is used; on read this should match the passed value.
 * If nval n<0, vector->size() is used.
 */
template<typename T, typename AllocatorType, typename Enum>
static int doVectorLow(XdrSerializer*                 serializer,
                       Enum                           ecpt,
                       const int64_t                  nval,
                       T**                            v,
                       std::vector<T, AllocatorType>* vector,
                       FILE*                          list,
                       CptElementType                 cptElementType)
{
    GMX_RELEASE_ASSERT(list != nullptr || (v != nullptr && vector == nullptr)
                               || (v == nullptr && vector != nullptr),
                       "Without list, we should have exactly one of v and vector != NULL");

    bool_t res = 0;

    unsigned int numElemInTheFile;
    if (list == nullptr)
    {
        if (nval >= 0)
        {
            // Remove this check when we store int64_t in the file
            GMX_RELEASE_ASSERT(nval <= std::numeric_limits<unsigned int>::max(),
                               "Vector size in checkpoint beyond max uint");

            numElemInTheFile = nval;
        }
        else
        {
            GMX_RELEASE_ASSERT(v == nullptr, "With nval<0 we should have v=nullptr");
            // Remove this check when we store int64_t in the file
            GMX_RELEASE_ASSERT(
                    vector->size() <= static_cast<std::size_t>(std::numeric_limits<unsigned int>::max()),
                    "Vector size in checkpoint beyond max uint");

            numElemInTheFile = vector->size();
        }
    }
    /* Read/write the vector element count */
    /* We store an unsigned int as a signed int to avoid changing the file format */
    int* numElemInTheFileIntPtr = reinterpret_cast<int*>(&numElemInTheFile);
    res                         = xdr_int(serializer->xdr(), numElemInTheFileIntPtr);
    if (res == 0)
    {
        return -1;
    }
    /* Read/write the element data type */
    constexpr XdrDataType xdrTypeInTheCode      = xdr_type<T>::value;
    XdrDataType           xdrTypeInTheFile      = xdrTypeInTheCode;
    int                   xdrTypeInTheFileAsInt = static_cast<int>(xdrTypeInTheFile);
    res              = xdr_int(serializer->xdr(), &xdrTypeInTheFileAsInt);
    xdrTypeInTheFile = static_cast<XdrDataType>(xdrTypeInTheFileAsInt);
    if (res == 0)
    {
        return -1;
    }

    if (list == nullptr)
    {
        if (nval >= 0)
        {
            if (numElemInTheFile != nval)
            {
                gmx_fatal(FARGS,
                          "Count mismatch for state entry %s, code count is %" PRId64
                          ", file count is %u\n",
                          enumValueToString(ecpt),
                          nval,
                          numElemInTheFile);
            }
        }

        bool typesMatch = (xdrTypeInTheFile == xdrTypeInTheCode);
        if (!typesMatch)
        {
            char buf[STRLEN];
            sprintf(buf,
                    "mismatch for state entry %s, code precision is %s, file precision is %s",
                    enumValueToString(ecpt),
                    enumValueToString(xdrTypeInTheCode),
                    enumValueToString(xdrTypeInTheFile));

            /* Matching int and real should never occur, but check anyhow */
            if (xdrTypeInTheFile == XdrDataType::Int || xdrTypeInTheCode == XdrDataType::Int)
            {
                gmx_fatal(FARGS,
                          "Type %s: incompatible checkpoint formats or corrupted checkpoint file.",
                          buf);
            }
        }

        T* vp;
        if (v != nullptr)
        {
            if (*v == nullptr)
            {
                snew(*v, numElemInTheFile);
            }
            vp = *v;
        }
        else
        {
            GMX_RELEASE_ASSERT(vector != nullptr, "Without list or v, vector should be supplied");
            /* This conditional ensures that we don't resize on write.
             * In particular in the state where this code was written
             * vector has a size of numElemInThefile and we
             * don't want to lose that padding here.
             */
            if (vector->size() < static_cast<unsigned int>(numElemInTheFile))
            {
                vector->resize(numElemInTheFile);
            }
            vp = vector->data();
        }

        char* vChar;
        if (typesMatch)
        {
            vChar = reinterpret_cast<char*>(vp);
        }
        else
        {
            snew(vChar, numElemInTheFile * sizeOfXdrType(xdrTypeInTheFile));
        }
        res = xdr_vector(serializer->xdr(),
                         vChar,
                         numElemInTheFile,
                         sizeOfXdrType(xdrTypeInTheFile),
                         xdrProc(xdrTypeInTheFile));
        if (res == 0)
        {
            return -1;
        }

        if (!typesMatch)
        {
            /* In the old code float-double conversion came for free.
             * In the new code we still support it, mainly because
             * the tip4p_continue regression test makes use of this.
             * It's an open question if we do or don't want to allow this.
             */
            convertArrayRealPrecision(vChar, vp, numElemInTheFile);
            sfree(vChar);
        }
    }
    else
    {
        res = listXdrVector(serializer, ecpt, numElemInTheFile, xdrTypeInTheFile, list, cptElementType);
    }

    return 0;
}

//! \brief Read/Write a std::vector, on read checks the number of elements matches \p numElements, if specified.
template<typename T, typename Enum>
static int doVector(XdrSerializer* serializer, Enum ecpt, std::vector<T>* vector, FILE* list, int numElements = -1)
{
    return doVectorLow<T>(serializer, ecpt, numElements, nullptr, vector, list, CptElementType::realnum);
}

//! \brief Read/Write an ArrayRef<real>.
template<typename Enum>
static int doRealArrayRef(XdrSerializer* serializer, Enum ecpt, gmx::ArrayRef<real> vector, FILE* list)
{
    real* v_real = vector.data();
    return doVectorLow<real, std::allocator<real>>(
            serializer, ecpt, vector.size(), &v_real, nullptr, list, CptElementType::realnum);
}

//! Convert from view of RVec to view of real.
static gmx::ArrayRef<real> realArrayRefFromRVecArrayRef(gmx::ArrayRef<gmx::RVec> ofRVecs)
{
    return gmx::arrayRefFromArray<real>(reinterpret_cast<real*>(ofRVecs.data()), ofRVecs.size() * DIM);
}

//! \brief Read/Write a PaddedVector whose value_type is RVec.
template<typename PaddedVectorOfRVecType, typename Enum>
static int doRvecVector(XdrSerializer* serializer, Enum ecpt, PaddedVectorOfRVecType* v, int numAtoms, FILE* list)
{
    const int numReals = numAtoms * DIM;

    if (list == nullptr)
    {
        GMX_RELEASE_ASSERT(v->size() == numAtoms, "v should have sufficient size for numAtoms");

        return doRealArrayRef(serializer, ecpt, realArrayRefFromRVecArrayRef(makeArrayRef(*v)), list);
    }
    else
    {
        // Use the rebind facility to change the value_type of the
        // allocator from RVec to real.
        using realAllocator =
                typename std::allocator_traits<typename PaddedVectorOfRVecType::allocator_type>::template rebind_alloc<real>;
        return doVectorLow<real, realAllocator>(
                serializer, ecpt, numReals, nullptr, nullptr, list, CptElementType::realnum);
    }
}

/* This function stores n along with the reals for reading,
 * but on reading it assumes that n matches the value in the checkpoint file,
 * a fatal error is generated when this is not the case.
 */
template<typename Enum>
static int do_cpte_reals(XdrSerializer* serializer, Enum ecpt, int n, real** v, FILE* list)
{
    return doVectorLow<real, std::allocator<real>>(
            serializer, ecpt, n, v, nullptr, list, CptElementType::realnum);
}

template<typename Enum>
static int do_cpte_real(XdrSerializer* serializer, Enum ecpt, real* r, FILE* list)
{
    return doVectorLow<real, std::allocator<real>>(
            serializer, ecpt, 1, &r, nullptr, list, CptElementType::realnum);
}

template<typename Enum>
static int do_cpte_ints(XdrSerializer* serializer, Enum ecpt, int n, int** v, FILE* list)
{
    return doVectorLow<int, std::allocator<int>>(
            serializer, ecpt, n, v, nullptr, list, CptElementType::integer);
}

template<typename Enum>
static int do_cpte_int(XdrSerializer* serializer, Enum ecpt, int* i, FILE* list)
{
    return do_cpte_ints(serializer, ecpt, 1, &i, list);
}

template<typename Enum>
static int do_cpte_bool(XdrSerializer* serializer, Enum ecpt, bool* b, FILE* list)
{
    int i   = static_cast<int>(*b);
    int ret = do_cpte_int(serializer, ecpt, &i, list);
    *b      = (i != 0);
    return ret;
}

template<typename Enum>
static int do_cpte_doubles(XdrSerializer* serializer, Enum ecpt, int n, double** v, FILE* list)
{
    return doVectorLow<double, std::allocator<double>>(
            serializer, ecpt, n, v, nullptr, list, CptElementType::realnum);
}

template<typename Enum>
static int do_cpte_double(XdrSerializer* serializer, Enum ecpt, double* r, FILE* list)
{
    return do_cpte_doubles(serializer, ecpt, 1, &r, list);
}

template<typename Enum>
static int do_cpte_matrix(XdrSerializer* serializer, Enum ecpt, matrix v, FILE* list)
{
    real* vr;
    int   ret;

    vr  = &(v[0][0]);
    ret = doVectorLow<real, std::allocator<real>>(
            serializer, ecpt, DIM * DIM, &vr, nullptr, nullptr, CptElementType::matrix3x3);

    if (list && ret == 0)
    {
        pr_rvecs(list, 0, enumValueToString(ecpt), v, DIM);
    }

    return ret;
}

template<typename Enum>
static int do_cpte_nmatrix(XdrSerializer* serializer, Enum ecpt, int n, real** v, FILE* list)
{
    int  i;
    int  ret, reti;
    char name[CPTSTRLEN];

    ret = 0;
    if (v == nullptr)
    {
        snew(v, n);
    }
    for (i = 0; i < n; i++)
    {
        reti = doVectorLow<real, std::allocator<real>>(
                serializer, ecpt, n, &(v[i]), nullptr, nullptr, CptElementType::matrix3x3);
        if (list && reti == 0)
        {
            sprintf(name, "%s[%d]", enumValueToString(ecpt), i);
            pr_reals(list, 0, name, v[i], n);
        }
        if (reti != 0)
        {
            ret = reti;
        }
    }
    return ret;
}

template<typename Enum>
static int do_cpte_matrices(XdrSerializer* serializer, Enum ecpt, int n, matrix** v, FILE* list)
{
    bool_t  res = 0;
    matrix *vp, *va = nullptr;
    real*   vr;
    int     nf, i, j, k;
    int     ret;

    nf  = n;
    res = xdr_int(serializer->xdr(), &nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == nullptr && nf != n)
    {
        gmx_fatal(FARGS,
                  "Count mismatch for state entry %s, code count is %d, file count is %d\n",
                  enumValueToString(ecpt),
                  n,
                  nf);
    }
    if (list)
    {
        snew(va, nf);
        vp = va;
    }
    else
    {
        if (*v == nullptr)
        {
            snew(*v, nf);
        }
        vp = *v;
    }
    snew(vr, nf * DIM * DIM);
    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                vr[(i * DIM + j) * DIM + k] = vp[i][j][k];
            }
        }
    }
    ret = doVectorLow<real, std::allocator<real>>(
            serializer, ecpt, nf * DIM * DIM, &vr, nullptr, nullptr, CptElementType::matrix3x3);
    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                vp[i][j][k] = vr[(i * DIM + j) * DIM + k];
            }
        }
    }
    sfree(vr);

    if (list && ret == 0)
    {
        for (i = 0; i < nf; i++)
        {
            pr_rvecs(list, 0, enumValueToString(ecpt), vp[i], DIM);
        }
    }
    if (va)
    {
        sfree(va);
    }

    return ret;
}

static void do_cpt_header(XdrSerializer* serializer, gmx_bool bRead, FILE* list, CheckpointHeaderContents* contents)
{
    bool_t res = 0;
    int    magic;

    if (bRead)
    {
        magic = -1;
    }
    else
    {
        magic = CPT_MAGIC1;
    }
    res = xdr_int(serializer->xdr(), &magic);
    if (res == 0)
    {
        gmx_fatal(FARGS,
                  "The checkpoint file is empty/corrupted, or maybe you are out of disk space?");
    }
    if (magic != CPT_MAGIC1)
    {
        gmx_fatal(FARGS,
                  "Start of file magic number mismatch, checkpoint file has %d, should be %d\n"
                  "The checkpoint file is corrupted or not a checkpoint file",
                  magic,
                  CPT_MAGIC1);
    }
    char fhost[255];
    if (!bRead)
    {
        gmx_gethostname(fhost, 255);
    }
    do_cpt_string_err(serializer, "GROMACS version", contents->version, list);
    // The following fields are no longer ever written with meaningful
    // content, but because they precede the file version, there is no
    // good way for new code to read the old and new formats, nor a
    // good way for old code to avoid giving an error while reading a
    // new format. So we read and write a field that no longer has a
    // purpose.
    do_cpt_string_err(serializer, "GROMACS build time UNUSED", contents->btime_UNUSED, list);
    do_cpt_string_err(serializer, "GROMACS build user UNUSED", contents->buser_UNUSED, list);
    do_cpt_string_err(serializer, "GROMACS build host UNUSED", contents->bhost_UNUSED, list);
    do_cpt_string_err(serializer, "generating program", contents->fprog, list);
    do_cpt_string_err(serializer, "generation time", contents->ftime, list);
    contents->file_version = cpt_version;
    do_cpt_enum_as_int<CheckPointVersion>(
            serializer, "checkpoint file version", &contents->file_version, list);
    if (contents->file_version > cpt_version)
    {
        gmx_fatal(FARGS,
                  "Attempting to read a checkpoint file of version %d with code of version %d\n",
                  static_cast<int>(contents->file_version),
                  static_cast<int>(cpt_version));
    }
    if (contents->file_version >= CheckPointVersion::DoublePrecisionBuild)
    {
        do_cpt_int_err(serializer, "GROMACS double precision", &contents->double_prec, list);
    }
    else
    {
        contents->double_prec = -1;
    }
    if (contents->file_version >= CheckPointVersion::HostInformation)
    {
        do_cpt_string_err(serializer, "generating host", fhost, list);
    }
    do_cpt_int_err(serializer, "#atoms", &contents->natoms, list);
    do_cpt_int_err(serializer, "#T-coupling groups", &contents->ngtc, list);
    if (contents->file_version >= CheckPointVersion::NoseHooverThermostat)
    {
        do_cpt_int_err(serializer, "#Nose-Hoover T-chains", &contents->nhchainlength, list);
    }
    else
    {
        contents->nhchainlength = 1;
    }
    if (contents->file_version >= CheckPointVersion::NoseHooverBarostat)
    {
        do_cpt_int_err(serializer, "#Nose-Hoover T-chains for barostat ", &contents->nnhpres, list);
    }
    else
    {
        contents->nnhpres = 0;
    }
    if (contents->file_version >= CheckPointVersion::LambdaStateAndHistory)
    {
        do_cpt_int_err(serializer, "# of total lambda states ", &contents->nlambda, list);
    }
    else
    {
        contents->nlambda = 0;
    }
    {
        int integrator = static_cast<int>(contents->eIntegrator);
        do_cpt_int_err(serializer, "integrator", &integrator, list);
        if (bRead)
        {
            contents->eIntegrator = static_cast<IntegrationAlgorithm>(integrator);
        }
    }
    if (contents->file_version >= CheckPointVersion::SafeSimulationPart)
    {
        do_cpt_int_err(serializer, "simulation part #", &contents->simulation_part, list);
    }
    else
    {
        contents->simulation_part = 1;
    }
    if (contents->file_version >= CheckPointVersion::SafeSteps)
    {
        do_cpt_step_err(serializer, "step", &contents->step, list);
    }
    else
    {
        int idum = 0;
        do_cpt_int_err(serializer, "step", &idum, list);
        contents->step = static_cast<int64_t>(idum);
    }
    do_cpt_double_err(serializer, "t", &contents->t, list);
    do_cpt_int_err(serializer, "#PP-ranks", &contents->nnodes, list);
    do_cpt_int_err(serializer, "dd_nc[x]", &contents->dd_nc[XX], list);
    do_cpt_int_err(serializer, "dd_nc[y]", &contents->dd_nc[YY], list);
    do_cpt_int_err(serializer, "dd_nc[z]", &contents->dd_nc[ZZ], list);
    do_cpt_int_err(serializer, "#PME-only ranks", &contents->npme, list);
    do_cpt_int_err(serializer, "state flags", &contents->flags_state, list);
    if (contents->file_version >= CheckPointVersion::EkinDataAndFlags)
    {
        do_cpt_int_err(serializer, "ekin data flags", &contents->flags_eks, list);
        do_cpt_int_err(serializer, "energy history flags", &contents->flags_enh, list);
    }
    else
    {
        contents->flags_eks = 0;
        contents->flags_enh = (contents->flags_state >> (static_cast<int>(StateEntry::OrireDtav) + 1));
        contents->flags_state = (contents->flags_state
                                 & ~((1 << (static_cast<int>(StateEntry::OrireDtav) + 1))
                                     | (1 << (static_cast<int>(StateEntry::OrireDtav) + 2))
                                     | (1 << (static_cast<int>(StateEntry::OrireDtav) + 3))));
    }
    if (contents->file_version >= CheckPointVersion::LambdaStateAndHistory)
    {
        do_cpt_int_err(serializer, "df history flags", &contents->flags_dfh, list);
    }
    else
    {
        contents->flags_dfh = 0;
    }

    if (contents->file_version >= CheckPointVersion::EssentialDynamics)
    {
        do_cpt_int_err(serializer, "ED data sets", &contents->nED, list);
    }
    else
    {
        contents->nED = 0;
    }

    if (contents->file_version >= CheckPointVersion::SwapState)
    {
        int swapState = static_cast<int>(contents->eSwapCoords);
        do_cpt_int_err(serializer, "swap", &swapState, list);
        if (bRead)
        {
            contents->eSwapCoords = static_cast<SwapType>(swapState);
        }
    }
    else
    {
        contents->eSwapCoords = SwapType::No;
    }

    if (contents->file_version >= CheckPointVersion::AwhHistoryFlags)
    {
        do_cpt_int_err(serializer, "AWH history flags", &contents->flags_awhh, list);
    }
    else
    {
        contents->flags_awhh = 0;
    }

    if (contents->file_version >= CheckPointVersion::RemoveBuildMachineInformation)
    {
        do_cpt_int_err(serializer, "pull history flags", &contents->flagsPullHistory, list);
    }
    else
    {
        contents->flagsPullHistory = 0;
    }

    if (contents->file_version >= CheckPointVersion::ModularSimulator)
    {
        do_cpt_bool_err(
                serializer, "Is modular simulator checkpoint", &contents->isModularSimulatorCheckpoint, list);
    }
    else
    {
        contents->isModularSimulatorCheckpoint = false;
    }
}

static int do_cpt_footer(XdrSerializer* serializer, CheckPointVersion file_version)
{
    bool_t res = 0;
    int    magic;

    if (file_version >= CheckPointVersion::AddMagicNumber)
    {
        magic = CPT_MAGIC2;
        res   = xdr_int(serializer->xdr(), &magic);
        if (res == 0)
        {
            cp_error();
        }
        if (magic != CPT_MAGIC2)
        {
            return -1;
        }
    }

    return 0;
}

static int do_cpt_state(XdrSerializer* serializer, t_state* state, FILE* list)
{
    GMX_RELEASE_ASSERT(static_cast<unsigned int>(state->numAtoms())
                               <= std::numeric_limits<unsigned int>::max() / 3,
                       "Can not write more than max_int/3 atoms to checkpoint");

    int       ret    = 0;
    const int sflags = state->flags();
    using StateFlags = gmx::EnumerationArray<StateEntry, bool>;
    for (auto i = StateFlags::keys().begin(); (i != StateFlags::keys().end() && ret == 0); i++)
    {
        if (sflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {
                case StateEntry::Lambda:
                    ret = doRealArrayRef(serializer, *i, state->lambda, list);
                    break;
                case StateEntry::FepState:
                    ret = do_cpte_int(serializer, *i, &state->fep_state, list);
                    break;
                case StateEntry::Box: ret = do_cpte_matrix(serializer, *i, state->box, list); break;
                case StateEntry::BoxRel:
                    ret = do_cpte_matrix(serializer, *i, state->box_rel, list);
                    break;
                case StateEntry::BoxV:
                    ret = do_cpte_matrix(serializer, *i, state->boxv, list);
                    break;
                case StateEntry::PressurePrevious:
                    ret = do_cpte_matrix(serializer, *i, state->pres_prev, list);
                    break;
                case StateEntry::SVirPrev:
                    ret = do_cpte_matrix(serializer, *i, state->svir_prev, list);
                    break;
                case StateEntry::FVirPrev:
                    ret = do_cpte_matrix(serializer, *i, state->fvir_prev, list);
                    break;
                case StateEntry::Nhxi:
                    ret = doVector<double>(serializer, *i, &state->nosehoover_xi, list);
                    break;
                case StateEntry::Nhvxi:
                    ret = doVector<double>(serializer, *i, &state->nosehoover_vxi, list);
                    break;
                case StateEntry::Nhpresxi:
                    ret = doVector<double>(serializer, *i, &state->nhpres_xi, list);
                    break;
                case StateEntry::Nhpresvxi:
                    ret = doVector<double>(serializer, *i, &state->nhpres_vxi, list);
                    break;
                case StateEntry::ThermInt:
                    ret = doVector<double>(serializer, *i, &state->therm_integral, list);
                    break;
                case StateEntry::BarosInt:
                    ret = do_cpte_double(serializer, *i, &state->baros_integral, list);
                    break;
                case StateEntry::Veta:
                    ret = do_cpte_real(serializer, *i, &state->veta, list);
                    break;
                case StateEntry::Vol0:
                    ret = do_cpte_real(serializer, *i, &state->vol0, list);
                    break;
                case StateEntry::X:
                    ret = doRvecVector(serializer, *i, &state->x, state->numAtoms(), list);
                    break;
                case StateEntry::V:
                    ret = doRvecVector(serializer, *i, &state->v, state->numAtoms(), list);
                    break;
                /* The RNG entries are no longer written,
                 * the next 4 lines are only for reading old files.
                 * It's OK that three case statements fall through.
                 */
                case StateEntry::LDRngNotSupported:
                case StateEntry::LDRngINotSupported:
                case StateEntry::MCRngNotSupported:
                case StateEntry::MCRngINotSupported:
                    ret = do_cpte_ints(serializer, *i, 0, nullptr, list);
                    break;
                case StateEntry::DisreInitF:
                    ret = do_cpte_real(serializer, *i, &state->hist.disre_initf, list);
                    break;
                case StateEntry::DisreRm3Tav:
                    ret = doVector<real>(serializer, *i, &state->hist.disre_rm3tav, list);
                    break;
                case StateEntry::OrireInitF:
                    ret = do_cpte_real(serializer, *i, &state->hist.orire_initf, list);
                    break;
                case StateEntry::OrireDtav:
                    ret = doVector<real>(serializer, *i, &state->hist.orire_Dtav, list);
                    break;
                case StateEntry::PullComPrevStep:
                    ret = doVector<double>(serializer, *i, &state->pull_com_prev_step, list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown state entry %d\n"
                              "You are reading a checkpoint file written by different code, which "
                              "is not supported",
                              enumValueToBitMask(*i));
            }
        }
    }
    return ret;
}

static int do_cpt_ekinstate(XdrSerializer* serializer, int fflags, ekinstate_t* ekins, FILE* list)
{
    int ret = 0;

    using StateFlags = gmx::EnumerationArray<StateKineticEntry, bool>;
    for (auto i = StateFlags::keys().begin(); (i != StateFlags::keys().end() && ret == 0); i++)
    {
        if (fflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {

                case StateKineticEntry::EkinNumber:
                    ret = do_cpte_int(serializer, *i, &ekins->ekin_n, list);
                    break;
                case StateKineticEntry::EkinHalfStep:
                    ret = do_cpte_matrices(serializer, *i, ekins->ekin_n, &ekins->ekinh, list);
                    break;
                case StateKineticEntry::EkinFullStep:
                    ret = do_cpte_matrices(serializer, *i, ekins->ekin_n, &ekins->ekinf, list);
                    break;
                case StateKineticEntry::EkinHalfStepOld:
                    ret = do_cpte_matrices(serializer, *i, ekins->ekin_n, &ekins->ekinh_old, list);
                    break;
                case StateKineticEntry::EkinTotal:
                    ret = do_cpte_matrix(serializer, *i, ekins->ekin_total, list);
                    break;
                case StateKineticEntry::EkinNoseHooverScaleFullStep:
                    ret = doVector<double>(serializer, *i, &ekins->ekinscalef_nhc, list);
                    break;
                case StateKineticEntry::VelocityScale:
                    ret = doVector<double>(serializer, *i, &ekins->vscale_nhc, list);
                    break;
                case StateKineticEntry::EkinNoseHooverScaleHalfStep:
                    ret = doVector<double>(serializer, *i, &ekins->ekinscaleh_nhc, list);
                    break;
                case StateKineticEntry::DEkinDLambda:
                    ret = do_cpte_real(serializer, *i, &ekins->dekindl, list);
                    break;
                case StateKineticEntry::Mvcos:
                    ret = do_cpte_real(serializer, *i, &ekins->mvcos, list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown ekin data state entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              enumValueToBitMask(*i));
            }
        }
    }

    return ret;
}


static int do_cpt_swapstate(XdrSerializer* serializer,
                            gmx_bool       bRead,
                            SwapType       eSwapCoords,
                            swaphistory_t* swapstate,
                            FILE*          list)
{
    int swap_cpt_version = 2;

    if (eSwapCoords == SwapType::No)
    {
        return 0;
    }

    swapstate->bFromCpt    = bRead;
    swapstate->eSwapCoords = eSwapCoords;

    do_cpt_int_err(serializer, "swap checkpoint version", &swap_cpt_version, list);
    if (bRead && swap_cpt_version < 2)
    {
        gmx_fatal(FARGS,
                  "Cannot read checkpoint files that were written with old versions"
                  "of the ion/water position swapping protocol.\n");
    }

    do_cpt_int_err(serializer, "swap coupling steps", &swapstate->nAverage, list);

    /* When reading, init_swapcoords has not been called yet,
     * so we have to allocate memory first. */
    int numIonTypes = swapstate->ionType.size();
    do_cpt_int_err(serializer, "number of ion types", &numIonTypes, list);
    if (bRead)
    {
        swapstate->ionType.resize(numIonTypes);
    }

    for (auto ic : gmx::EnumerationWrapper<Compartment>{})
    {
        for (swapstateIons_t& gs : swapstate->ionType)
        {
            if (bRead)
            {
                do_cpt_int_err(serializer, "swap requested atoms", &gs.nMolReq[ic], list);
            }
            else
            {
                do_cpt_int_err(serializer, "swap requested atoms p", gs.nMolReq_p[ic], list);
            }

            if (bRead)
            {
                do_cpt_int_err(serializer, "swap influx net", &gs.inflow_net[ic], list);
            }
            else
            {
                do_cpt_int_err(serializer, "swap influx net p", gs.inflow_net_p[ic], list);
            }

            if (bRead && (nullptr == gs.nMolPast[ic]))
            {
                snew(gs.nMolPast[ic], swapstate->nAverage);
            }

            for (int j = 0; j < swapstate->nAverage; j++)
            {
                if (bRead)
                {
                    do_cpt_int_err(serializer, "swap past atom counts", &gs.nMolPast[ic][j], list);
                }
                else
                {
                    do_cpt_int_err(serializer, "swap past atom counts p", &gs.nMolPast_p[ic][j], list);
                }
            }
        }
    }

    /* Ion flux per channel */
    for (auto ic : gmx::EnumerationWrapper<Channel>{})
    {
        for (swapstateIons_t& gs : swapstate->ionType)
        {
            if (bRead)
            {
                do_cpt_int_err(serializer, "channel flux A->B", &gs.fluxfromAtoB[ic], list);
            }
            else
            {
                do_cpt_int_err(serializer, "channel flux A->B p", gs.fluxfromAtoB_p[ic], list);
            }
        }
    }

    /* Ion flux leakage */
    if (bRead)
    {
        do_cpt_int_err(serializer, "flux leakage", &swapstate->fluxleak, list);
    }
    else
    {
        do_cpt_int_err(serializer, "flux leakage", swapstate->fluxleak_p, list);
    }

    /* Ion history */
    for (swapstateIons_t& gs : swapstate->ionType)
    {
        int numMolecules;
        if (!bRead)
        {
            numMolecules = gs.channel_label->size();
        }
        do_cpt_int_err(serializer, "number of ions", &numMolecules, list);

        if (bRead)
        {
            gs.channel_label = std::make_shared<std::vector<ChannelHistory>>(numMolecules);
            gs.comp_from     = std::make_shared<std::vector<Domain>>(numMolecules);
        }

        do_cpt_n_enum_as_int<ChannelHistory>(
                serializer, "channel history", numMolecules, gs.channel_label->data(), list);
        do_cpt_n_enum_as_int<Domain>(serializer, "domain history", numMolecules, gs.comp_from->data(), list);
    }

    /* Save the last known whole positions to checkpoint
     * file to be able to also make multimeric channels whole in PBC */
    int numAtoms;
    if (!bRead)
    {
        numAtoms = swapstate->xc_old_whole_p[Channel::Zero]->size();
    }
    do_cpt_int_err(serializer, "Ch0 atoms", &numAtoms, list);
    if (bRead)
    {
        swapstate->xc_old_whole[Channel::Zero].resize(numAtoms);
    }
    do_cpt_int_err(serializer, "Ch1 atoms", &numAtoms, list);
    if (bRead)
    {
        swapstate->xc_old_whole[Channel::One].resize(numAtoms);
        do_cpt_n_rvecs_err(serializer,
                           "Ch0 whole x",
                           swapstate->xc_old_whole[Channel::Zero].size(),
                           as_rvec_array(swapstate->xc_old_whole[Channel::Zero].data()),
                           list);
        do_cpt_n_rvecs_err(serializer,
                           "Ch1 whole x",
                           swapstate->xc_old_whole[Channel::One].size(),
                           as_rvec_array(swapstate->xc_old_whole[Channel::One].data()),
                           list);
    }
    else
    {
        do_cpt_n_rvecs_err(serializer,
                           "Ch0 whole x",
                           swapstate->xc_old_whole_p[Channel::Zero]->size(),
                           as_rvec_array(swapstate->xc_old_whole_p[Channel::Zero]->data()),
                           list);
        do_cpt_n_rvecs_err(serializer,
                           "Ch1 whole x",
                           swapstate->xc_old_whole_p[Channel::One]->size(),
                           as_rvec_array(swapstate->xc_old_whole_p[Channel::One]->data()),
                           list);
    }

    return 0;
}


static int do_cpt_enerhist(XdrSerializer* serializer, gmx_bool bRead, int fflags, energyhistory_t* enerhist, FILE* list)
{
    int ret = 0;

    if (fflags == 0)
    {
        return ret;
    }

    GMX_RELEASE_ASSERT(enerhist != nullptr, "With energy history, we need a valid enerhist pointer");

    /* This is stored/read for backward compatibility */
    int energyHistoryNumEnergies = 0;
    if (bRead)
    {
        enerhist->nsteps     = 0;
        enerhist->nsum       = 0;
        enerhist->nsteps_sim = 0;
        enerhist->nsum_sim   = 0;
    }
    else if (enerhist != nullptr)
    {
        energyHistoryNumEnergies = enerhist->ener_sum_sim.size();
    }

    delta_h_history_t* deltaH = enerhist->deltaHForeignLambdas.get();
    using StateFlags          = gmx::EnumerationArray<StateEnergyEntry, bool>;
    for (auto i = StateFlags::keys().begin(); (i != StateFlags::keys().end() && ret == 0); i++)
    {
        if (fflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {
                case StateEnergyEntry::N:
                    ret = do_cpte_int(serializer, *i, &energyHistoryNumEnergies, list);
                    break;
                case StateEnergyEntry::Aver:
                    ret = doVector<double>(serializer, *i, &enerhist->ener_ave, list);
                    break;
                case StateEnergyEntry::Sum:
                    ret = doVector<double>(serializer, *i, &enerhist->ener_sum, list);
                    break;
                case StateEnergyEntry::NumSum:
                    do_cpt_step_err(serializer, enumValueToString(*i), &enerhist->nsum, list);
                    break;
                case StateEnergyEntry::SumSim:
                    ret = doVector<double>(serializer, *i, &enerhist->ener_sum_sim, list);
                    break;
                case StateEnergyEntry::NumSumSim:
                    do_cpt_step_err(serializer, enumValueToString(*i), &enerhist->nsum_sim, list);
                    break;
                case StateEnergyEntry::NumSteps:
                    do_cpt_step_err(serializer, enumValueToString(*i), &enerhist->nsteps, list);
                    break;
                case StateEnergyEntry::NumStepsSim:
                    do_cpt_step_err(serializer, enumValueToString(*i), &enerhist->nsteps_sim, list);
                    break;
                case StateEnergyEntry::DeltaHNN:
                {
                    int numDeltaH = 0;
                    if (!bRead && deltaH != nullptr)
                    {
                        numDeltaH = deltaH->dh.size();
                    }
                    do_cpt_int_err(serializer, enumValueToString(*i), &numDeltaH, list);
                    if (bRead)
                    {
                        if (deltaH == nullptr)
                        {
                            enerhist->deltaHForeignLambdas = std::make_unique<delta_h_history_t>();
                            deltaH                         = enerhist->deltaHForeignLambdas.get();
                        }
                        deltaH->dh.resize(numDeltaH);
                        deltaH->start_lambda_set = FALSE;
                    }
                    break;
                }
                case StateEnergyEntry::DeltaHList:
                    for (auto dh : deltaH->dh)
                    {
                        ret = doVector<real>(serializer, *i, &dh, list);
                    }
                    break;
                case StateEnergyEntry::DeltaHStartTime:
                    ret = do_cpte_double(serializer, *i, &(deltaH->start_time), list);
                    break;
                case StateEnergyEntry::DeltaHStartLambda:
                    ret = do_cpte_double(serializer, *i, &(deltaH->start_lambda), list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown energy history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              enumValueToBitMask(*i));
            }
        }
    }

    if ((fflags & enumValueToBitMask(StateEnergyEntry::Sum))
        && !(fflags & enumValueToBitMask(StateEnergyEntry::SumSim)))
    {
        /* Assume we have an old file format and copy sum to sum_sim */
        enerhist->ener_sum_sim = enerhist->ener_sum;
    }

    if ((fflags & enumValueToBitMask(StateEnergyEntry::NumSum))
        && !(fflags & enumValueToBitMask(StateEnergyEntry::NumSteps)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps = enerhist->nsum;
    }
    if ((fflags & enumValueToBitMask(StateEnergyEntry::NumSumSim))
        && !(fflags & enumValueToBitMask(StateEnergyEntry::NumStepsSim)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps_sim = enerhist->nsum_sim;
    }

    return ret;
}

static int doCptPullCoordHist(XdrSerializer* serializer, PullCoordinateHistory* pullCoordHist, FILE* list)
{
    if (do_cpte_double(serializer, StatePullCoordEntry::ValueReferenceSum, &(pullCoordHist->valueRef), list)
        != 0)
    {
        return -1;
    }
    if (do_cpte_double(serializer, StatePullCoordEntry::ValueSum, &(pullCoordHist->value), list) != 0)
    {
        return -1;
    }
    for (int j = 0; j < DIM; j++)
    {
        if (do_cpte_double(serializer, StatePullCoordEntry::DR01Sum, &(pullCoordHist->dr01[j]), list) != 0)
        {
            return -1;
        }
    }
    for (int j = 0; j < DIM; j++)
    {
        if (do_cpte_double(serializer, StatePullCoordEntry::DR23Sum, &(pullCoordHist->dr23[j]), list) != 0)
        {
            return -1;
        }
    }
    for (int j = 0; j < DIM; j++)
    {
        if (do_cpte_double(serializer, StatePullCoordEntry::DR45Sum, &(pullCoordHist->dr45[j]), list) != 0)
        {
            return -1;
        }
    }
    if (do_cpte_double(serializer, StatePullCoordEntry::FScalarSum, &(pullCoordHist->scalarForce), list) != 0)
    {
        return -1;
    }
    for (int j = 0; j < DIM; j++)
    {
        if (do_cpte_double(serializer, StatePullCoordEntry::DynaxSum, &(pullCoordHist->dynaX[j]), list) != 0)
        {
            return -1;
        }
    }

    return 0;
}

static int doCptPullGroupHist(XdrSerializer* serializer, PullGroupHistory* pullGroupHist, FILE* list)
{
    for (int j = 0; j < DIM; j++)
    {
        if (do_cpte_double(serializer, StatePullGroupEntry::XSum, &(pullGroupHist->x[j]), list) != 0)
        {
            return -1;
        }
    }

    return 0;
}


static int doCptPullHist(XdrSerializer* serializer, gmx_bool bRead, int fflags, PullHistory* pullHist, FILE* list)
{
    int ret                       = 0;
    int pullHistoryNumCoordinates = 0;
    int pullHistoryNumGroups      = 0;

    /* Retain the number of terms in the sum and the number of coordinates (used for writing
     * average pull forces and coordinates) in the pull history, in temporary variables,
     * in case they cannot be read from the checkpoint, in order to have backward compatibility */
    if (bRead)
    {
        pullHist->numValuesInXSum = 0;
        pullHist->numValuesInFSum = 0;
    }
    else if (pullHist != nullptr)
    {
        pullHistoryNumCoordinates = pullHist->pullCoordinateSums.size();
        pullHistoryNumGroups      = pullHist->pullGroupSums.size();
    }
    else
    {
        GMX_RELEASE_ASSERT(fflags == 0, "Without pull history, all flags should be off");
    }

    using StateFlags = gmx::EnumerationArray<StatePullEntry, bool>;
    for (auto i = StateFlags::keys().begin(); i != StateFlags::keys().end(); i++)
    {
        if (fflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {
                case StatePullEntry::NumCoordinates:
                    do_cpt_int_err(serializer, enumValueToString(*i), &pullHistoryNumCoordinates, list);
                    break;
                case StatePullEntry::NumGroups:
                    do_cpt_int_err(serializer, enumValueToString(*i), &pullHistoryNumGroups, list);
                    break;
                case StatePullEntry::NumValuesInXSum:
                    do_cpt_int_err(serializer, enumValueToString(*i), &pullHist->numValuesInXSum, list);
                    break;
                case StatePullEntry::NumValuesInFSum:
                    do_cpt_int_err(serializer, enumValueToString(*i), &pullHist->numValuesInFSum, list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown pull history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              enumValueToBitMask(*i));
            }
        }
    }
    if (bRead)
    {
        pullHist->pullCoordinateSums.resize(pullHistoryNumCoordinates);
        pullHist->pullGroupSums.resize(pullHistoryNumGroups);
    }
    if (pullHist->numValuesInXSum > 0 || pullHist->numValuesInFSum > 0)
    {
        for (size_t i = 0; i < pullHist->pullCoordinateSums.size() && ret == 0; i++)
        {
            ret = doCptPullCoordHist(serializer, &(pullHist->pullCoordinateSums[i]), list);
        }
        for (size_t i = 0; i < pullHist->pullGroupSums.size() && ret == 0; i++)
        {
            ret = doCptPullGroupHist(serializer, &(pullHist->pullGroupSums[i]), list);
        }
    }

    return ret;
}

static int do_cpt_df_hist(XdrSerializer* serializer, int fflags, int nlambda, df_history_t** dfhistPtr, FILE* list)
{
    int ret = 0;

    if (fflags == 0)
    {
        return 0;
    }

    std::unique_ptr<df_history_t> localDFHistory = nullptr;
    if (*dfhistPtr == nullptr)
    {
        localDFHistory        = std::make_unique<df_history_t>();
        *dfhistPtr            = localDFHistory.get();
        (*dfhistPtr)->nlambda = nlambda;
        init_df_history(*dfhistPtr, nlambda);
    }
    df_history_t* dfhist = *dfhistPtr;

    using StateFlags = gmx::EnumerationArray<StateFepEntry, bool>;
    for (auto i = StateFlags::keys().begin(); (i != StateFlags::keys().end() && ret == 0); i++)
    {
        if (fflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {
                case StateFepEntry::IsEquilibrated:
                    ret = do_cpte_bool(serializer, *i, &dfhist->bEquil, list);
                    break;
                case StateFepEntry::NumAtLambdaStats:
                    ret = do_cpte_ints(
                            serializer, *i, nlambda, &dfhist->numSamplesAtLambdaForStatistics, list);
                    break;
                case StateFepEntry::NumAtLambdaEquil:
                    ret = do_cpte_ints(
                            serializer, *i, nlambda, &dfhist->numSamplesAtLambdaForEquilibration, list);
                    break;
                case StateFepEntry::WangLandauHistogram:
                    ret = do_cpte_reals(serializer, *i, nlambda, &dfhist->wl_histo, list);
                    break;
                case StateFepEntry::WangLandauDelta:
                    ret = do_cpte_real(serializer, *i, &dfhist->wl_delta, list);
                    break;
                case StateFepEntry::SumWeights:
                    ret = do_cpte_reals(serializer, *i, nlambda, &dfhist->sum_weights, list);
                    break;
                case StateFepEntry::SumDG:
                    ret = do_cpte_reals(serializer, *i, nlambda, &dfhist->sum_dg, list);
                    break;
                case StateFepEntry::SumMinVar:
                    ret = do_cpte_reals(serializer, *i, nlambda, &dfhist->sum_minvar, list);
                    break;
                case StateFepEntry::SumVar:
                    ret = do_cpte_reals(serializer, *i, nlambda, &dfhist->sum_variance, list);
                    break;
                case StateFepEntry::Accump:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->accum_p, list);
                    break;
                case StateFepEntry::Accumm:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->accum_m, list);
                    break;
                case StateFepEntry::Accump2:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->accum_p2, list);
                    break;
                case StateFepEntry::Accumm2:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->accum_m2, list);
                    break;
                case StateFepEntry::Tij:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->Tij, list);
                    break;
                case StateFepEntry::TijEmp:
                    ret = do_cpte_nmatrix(serializer, *i, nlambda, dfhist->Tij_empirical, list);
                    break;

                default:
                    gmx_fatal(FARGS,
                              "Unknown df history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              enumValueToBitMask(*i));
            }
        }
    }

    return ret;
}


/* This function stores the last whole configuration of the reference and
 * average structure in the .cpt file
 */
static int do_cpt_EDstate(XdrSerializer* serializer, gmx_bool bRead, int nED, edsamhistory_t* EDstate, FILE* list)
{
    if (nED == 0)
    {
        return 0;
    }

    EDstate->bFromCpt = bRead;
    EDstate->nED      = nED;

    /* When reading, init_edsam has not been called yet,
     * so we have to allocate memory first. */
    if (bRead)
    {
        snew(EDstate->nref, EDstate->nED);
        snew(EDstate->old_sref, EDstate->nED);
        snew(EDstate->nav, EDstate->nED);
        snew(EDstate->old_sav, EDstate->nED);
    }

    /* Read/write the last whole conformation of SREF and SAV for each ED dataset (usually only one) */
    for (int i = 0; i < EDstate->nED; i++)
    {
        char buf[STRLEN];

        /* Reference structure SREF */
        sprintf(buf, "ED%d # of atoms in reference structure", i + 1);
        do_cpt_int_err(serializer, buf, &EDstate->nref[i], list);
        sprintf(buf, "ED%d x_ref", i + 1);
        if (bRead)
        {
            snew(EDstate->old_sref[i], EDstate->nref[i]);
            do_cpt_n_rvecs_err(serializer, buf, EDstate->nref[i], EDstate->old_sref[i], list);
        }
        else
        {
            do_cpt_n_rvecs_err(serializer, buf, EDstate->nref[i], EDstate->old_sref_p[i], list);
        }

        /* Average structure SAV */
        sprintf(buf, "ED%d # of atoms in average structure", i + 1);
        do_cpt_int_err(serializer, buf, &EDstate->nav[i], list);
        sprintf(buf, "ED%d x_av", i + 1);
        if (bRead)
        {
            snew(EDstate->old_sav[i], EDstate->nav[i]);
            do_cpt_n_rvecs_err(serializer, buf, EDstate->nav[i], EDstate->old_sav[i], list);
        }
        else
        {
            do_cpt_n_rvecs_err(serializer, buf, EDstate->nav[i], EDstate->old_sav_p[i], list);
        }
    }

    return 0;
}

static int do_cpt_correlation_grid(XdrSerializer*               serializer,
                                   gmx_bool                     bRead,
                                   gmx_unused int               fflags,
                                   gmx::CorrelationGridHistory* corrGrid,
                                   FILE*                        list,
                                   StateAwhEntry                eawhh)
{
    int ret = 0;

    do_cpt_int_err(serializer, enumValueToString(eawhh), &(corrGrid->numCorrelationTensors), list);
    do_cpt_int_err(serializer, enumValueToString(eawhh), &(corrGrid->tensorSize), list);
    do_cpt_int_err(serializer, enumValueToString(eawhh), &(corrGrid->blockDataListSize), list);

    if (bRead)
    {
        initCorrelationGridHistory(
                corrGrid, corrGrid->numCorrelationTensors, corrGrid->tensorSize, corrGrid->blockDataListSize);
    }

    for (gmx::CorrelationBlockDataHistory& blockData : corrGrid->blockDataBuffer)
    {
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.blockSumWeight), list);
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.blockSumSquareWeight), list);
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.blockSumWeightX), list);
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.blockSumWeightY), list);
        do_cpt_double_err(
                serializer, enumValueToString(eawhh), &(blockData.sumOverBlocksSquareBlockWeight), list);
        do_cpt_double_err(
                serializer, enumValueToString(eawhh), &(blockData.sumOverBlocksBlockSquareWeight), list);
        do_cpt_double_err(
                serializer, enumValueToString(eawhh), &(blockData.sumOverBlocksBlockWeightBlockWeightX), list);
        do_cpt_double_err(
                serializer, enumValueToString(eawhh), &(blockData.sumOverBlocksBlockWeightBlockWeightY), list);
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.blockLength), list);
        do_cpt_int_err(serializer, enumValueToString(eawhh), &(blockData.previousBlockIndex), list);
        do_cpt_double_err(serializer, enumValueToString(eawhh), &(blockData.correlationIntegral), list);
    }

    return ret;
}

static int do_cpt_awh_bias(XdrSerializer*          serializer,
                           gmx_bool                bRead,
                           int                     fflags,
                           gmx::AwhBiasHistory*    biasHistory,
                           FILE*                   list,
                           const CheckPointVersion fileVersion)
{
    int ret = 0;

    gmx::AwhBiasStateHistory* state = &biasHistory->state;
    using StateFlags                = gmx::EnumerationArray<StateAwhEntry, bool>;
    for (auto i = StateFlags::keys().begin(); (i != StateFlags::keys().end() && ret == 0); i++)
    {
        if (fflags & enumValueToBitMask(*i))
        {
            switch (*i)
            {
                case StateAwhEntry::InInitial:
                    do_cpt_bool_err(serializer, enumValueToString(*i), &state->in_initial, list);
                    break;
                case StateAwhEntry::EquilibrateHistogram:
                    do_cpt_bool_err(serializer, enumValueToString(*i), &state->equilibrateHistogram, list);
                    break;
                case StateAwhEntry::HistogramSize:
                    do_cpt_double_err(serializer, enumValueToString(*i), &state->histSize, list);
                    break;
                case StateAwhEntry::NumPoints:
                {
                    int numPoints;
                    if (!bRead)
                    {
                        numPoints = biasHistory->pointState.size();
                    }
                    do_cpt_int_err(serializer, enumValueToString(*i), &numPoints, list);
                    if (bRead)
                    {
                        biasHistory->pointState.resize(numPoints);
                    }
                }
                break;
                case StateAwhEntry::CoordPoint:
                    for (auto& psh : biasHistory->pointState)
                    {
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.target, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.free_energy, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.bias, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.weightsum_iteration, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.weightsum_covering, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.weightsum_tot, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.weightsum_ref, list);
                        do_cpt_step_err(serializer, enumValueToString(*i), &psh.last_update_index, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.log_pmfsum, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.visits_iteration, list);
                        do_cpt_double_err(serializer, enumValueToString(*i), &psh.visits_tot, list);
                        if (fileVersion >= CheckPointVersion::AwhLocalWeightSum)
                        {
                            do_cpt_double_err(serializer, enumValueToString(*i), &psh.localWeightSum, list);
                        }
                        else
                        {
                            psh.localWeightSum = 0;
                        }
                    }
                    break;
                case StateAwhEntry::UmbrellaGridPoint:
                    do_cpt_int_err(serializer, enumValueToString(*i), &(state->umbrellaGridpoint), list);
                    break;
                case StateAwhEntry::UpdateList:
                    do_cpt_int_err(
                            serializer, enumValueToString(*i), &(state->origin_index_updatelist), list);
                    do_cpt_int_err(serializer, enumValueToString(*i), &(state->end_index_updatelist), list);
                    break;
                case StateAwhEntry::LogScaledSampleWeight:
                    do_cpt_double_err(
                            serializer, enumValueToString(*i), &(state->logScaledSampleWeight), list);
                    do_cpt_double_err(
                            serializer, enumValueToString(*i), &(state->maxLogScaledSampleWeight), list);
                    break;
                case StateAwhEntry::NumUpdates:
                    do_cpt_step_err(serializer, enumValueToString(*i), &(state->numUpdates), list);
                    break;
                case StateAwhEntry::ForceCorrelationGrid:
                    ret = do_cpt_correlation_grid(
                            serializer, bRead, fflags, &biasHistory->forceCorrelationGrid, list, *i);
                    break;
                default: gmx_fatal(FARGS, "Unknown awh history entry %d\n", enumValueToBitMask(*i));
            }
        }
    }

    return ret;
}

static int do_cpt_awh(XdrSerializer*          serializer,
                      gmx_bool                bRead,
                      int                     fflags,
                      gmx::AwhHistory*        awhHistory,
                      FILE*                   list,
                      const CheckPointVersion fileVersion)
{
    int ret = 0;

    if (fflags != 0)
    {
        std::shared_ptr<gmx::AwhHistory> awhHistoryLocal;

        if (awhHistory == nullptr)
        {
            GMX_RELEASE_ASSERT(bRead,
                               "do_cpt_awh should not be called for writing without an AwhHistory");

            awhHistoryLocal = std::make_shared<gmx::AwhHistory>();
            awhHistory      = awhHistoryLocal.get();
        }

        /* To be able to read and write the AWH data several parameters determining the layout of the AWH structs need to be known
           (nbias, npoints, etc.). The best thing (?) would be to have these passed to this function. When writing to a checkpoint
           these parameters are available in awh_history (after calling init_awh_history). When reading from a checkpoint though, there
           is no initialized awh_history (it is initialized and set in this function). The AWH parameters have not always been read
           at the time when this function is called for reading so I don't know how to pass them as input. Here, this is solved by
           when writing a checkpoint, also storing parameters needed for future reading of the checkpoint.

           Another issue is that some variables that AWH checkpoints don't have a registered enum and string (e.g. nbias below).
           One difficulty is the multilevel structure of the data which would need to be represented somehow. */

        int numBias;
        if (!bRead)
        {
            numBias = awhHistory->bias.size();
        }
        do_cpt_int_err(serializer, "awh_nbias", &numBias, list);

        if (bRead)
        {
            awhHistory->bias.resize(numBias);
        }
        for (auto& bias : awhHistory->bias)
        {
            ret = do_cpt_awh_bias(serializer, bRead, fflags, &bias, list, fileVersion);
            if (ret)
            {
                return ret;
            }
        }
        do_cpt_double_err(serializer, "awh_potential_offset", &awhHistory->potentialOffset, list);
    }
    return ret;
}

static void do_cpt_mdmodules(CheckPointVersion              fileVersion,
                             ISerializer*                   serializer,
                             const gmx::MDModulesNotifiers& mdModulesNotifiers,
                             FILE*                          outputFile)
{
    if (fileVersion >= CheckPointVersion::MDModules)
    {
        gmx::KeyValueTreeObject mdModuleCheckpointParameterTree = gmx::deserializeKeyValueTree(serializer);
        if (outputFile)
        {
            gmx::TextWriter textWriter(outputFile);
            gmx::dumpKeyValueTree(&textWriter, mdModuleCheckpointParameterTree);
        }
        gmx::MDModulesCheckpointReadingDataOnMain mdModuleCheckpointReadingDataOnMain = {
            mdModuleCheckpointParameterTree
        };
        mdModulesNotifiers.checkpointingNotifier_.notify(mdModuleCheckpointReadingDataOnMain);
    }
}

static int do_cpt_files(XdrSerializer*                    serializer,
                        gmx_bool                          bRead,
                        std::vector<gmx_file_position_t>* outputfiles,
                        FILE*                             list,
                        CheckPointVersion                 file_version)
{
    gmx_off_t                   offset;
    gmx_off_t                   mask = 0xFFFFFFFFL;
    int                         offset_high, offset_low;
    std::array<char, CPTSTRLEN> buf;
    GMX_RELEASE_ASSERT(outputfiles, "Must have valid outputfiles");

    // Ensure that reading pre-allocates outputfiles, while writing
    // writes what is already there.
    int nfiles = outputfiles->size();
    if (do_cpt_int(serializer, "number of output files", &nfiles, list) != 0)
    {
        return -1;
    }
    if (bRead)
    {
        outputfiles->resize(nfiles);
    }

    for (auto& outputfile : *outputfiles)
    {
        /* 64-bit XDR numbers are not portable, so it is stored as separate high/low fractions */
        if (bRead)
        {
            do_cpt_string_err(serializer, "output filename", buf, list);
            std::copy(std::begin(buf), std::end(buf), std::begin(outputfile.filename));

            if (do_cpt_int(serializer, "file_offset_high", &offset_high, list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(serializer, "file_offset_low", &offset_low, list) != 0)
            {
                return -1;
            }
            outputfile.offset = (static_cast<gmx_off_t>(offset_high) << 32)
                                | (static_cast<gmx_off_t>(offset_low) & mask);
        }
        else
        {
            do_cpt_string_err(serializer, "output filename", outputfile.filename, list);
            /* writing */
            offset = outputfile.offset;
            if (offset == -1)
            {
                offset_low  = -1;
                offset_high = -1;
            }
            else
            {
                offset_low  = static_cast<int>(offset & mask);
                offset_high = static_cast<int>((offset >> 32) & mask);
            }
            if (do_cpt_int(serializer, "file_offset_high", &offset_high, list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(serializer, "file_offset_low", &offset_low, list) != 0)
            {
                return -1;
            }
        }
        if (file_version >= CheckPointVersion::FileChecksumAndSize)
        {
            if (do_cpt_int(serializer, "file_checksum_size", &outputfile.checksumSize, list) != 0)
            {
                return -1;
            }
            if (do_cpt_u_chars(
                        serializer, "file_checksum", outputfile.checksum.size(), outputfile.checksum.data(), list)
                != 0)
            {
                return -1;
            }
        }
        else
        {
            outputfile.checksumSize = -1;
        }
    }
    return 0;
}

void write_checkpoint_data(const std::filesystem::path&      filename,
                           CheckpointHeaderContents          headerContents,
                           gmx_bool                          bExpanded,
                           LambdaWeightCalculation           elamstats,
                           t_state*                          state,
                           ObservablesHistory*               observablesHistory,
                           const gmx::MDModulesNotifiers&    mdModulesNotifiers,
                           std::vector<gmx_file_position_t>* outputfiles,
                           gmx::WriteCheckpointDataHolder*   modularSimulatorCheckpointData)
{
    XdrSerializer serializer(filename, "w");
    headerContents.flags_eks = 0;
    if (state->ekinstate.bUpToDate)
    {
        // Likely only EkinNumber, EkinHalfStep, EkinFullStep and DEkinDLambda
        // are necessary and the rest can go
        headerContents.flags_eks = (enumValueToBitMask(StateKineticEntry::EkinNumber)
                                    | enumValueToBitMask(StateKineticEntry::EkinHalfStep)
                                    | enumValueToBitMask(StateKineticEntry::EkinFullStep)
                                    | enumValueToBitMask(StateKineticEntry::EkinNoseHooverScaleFullStep)
                                    | enumValueToBitMask(StateKineticEntry::EkinNoseHooverScaleHalfStep)
                                    | enumValueToBitMask(StateKineticEntry::VelocityScale)
                                    | enumValueToBitMask(StateKineticEntry::DEkinDLambda)
                                    | enumValueToBitMask(StateKineticEntry::Mvcos));
    }
    headerContents.isModularSimulatorCheckpoint = !modularSimulatorCheckpointData->empty();

    energyhistory_t* enerhist = observablesHistory->energyHistory.get();
    headerContents.flags_enh  = 0;
    if (enerhist != nullptr && (enerhist->nsum > 0 || enerhist->nsum_sim > 0))
    {
        headerContents.flags_enh |= enumValueToBitMask(StateEnergyEntry::N)
                                    | enumValueToBitMask(StateEnergyEntry::NumSteps)
                                    | enumValueToBitMask(StateEnergyEntry::NumStepsSim);
        if (enerhist->nsum > 0)
        {
            headerContents.flags_enh |= (enumValueToBitMask(StateEnergyEntry::Aver)
                                         | enumValueToBitMask(StateEnergyEntry::Sum)
                                         | enumValueToBitMask(StateEnergyEntry::NumSum));
        }
        if (enerhist->nsum_sim > 0)
        {
            headerContents.flags_enh |= (enumValueToBitMask(StateEnergyEntry::SumSim)
                                         | enumValueToBitMask(StateEnergyEntry::NumSumSim));
        }
        if (enerhist->deltaHForeignLambdas != nullptr)
        {
            headerContents.flags_enh |= (enumValueToBitMask(StateEnergyEntry::DeltaHNN)
                                         | enumValueToBitMask(StateEnergyEntry::DeltaHList)
                                         | enumValueToBitMask(StateEnergyEntry::DeltaHStartTime)
                                         | enumValueToBitMask(StateEnergyEntry::DeltaHStartLambda));
        }
    }

    PullHistory* pullHist           = observablesHistory->pullHistory.get();
    headerContents.flagsPullHistory = 0;
    if (pullHist != nullptr && (pullHist->numValuesInXSum > 0 || pullHist->numValuesInFSum > 0))
    {
        headerContents.flagsPullHistory |= enumValueToBitMask(StatePullEntry::NumCoordinates);
        headerContents.flagsPullHistory |= (enumValueToBitMask(StatePullEntry::NumGroups)
                                            | enumValueToBitMask(StatePullEntry::NumValuesInXSum)
                                            | enumValueToBitMask(StatePullEntry::NumValuesInFSum));
    }

    headerContents.flags_dfh = 0;
    // Modular simulator uses the modularSimulatorCheckpointData object to store the expanded ensemble history
    if (bExpanded && !headerContents.isModularSimulatorCheckpoint)
    {
        headerContents.flags_dfh =
                (enumValueToBitMask(StateFepEntry::IsEquilibrated)
                 | enumValueToBitMask(StateFepEntry::NumAtLambdaStats)
                 | enumValueToBitMask(StateFepEntry::NumAtLambdaEquil)
                 | enumValueToBitMask(StateFepEntry::SumWeights) | enumValueToBitMask(StateFepEntry::SumDG)
                 | enumValueToBitMask(StateFepEntry::Tij) | enumValueToBitMask(StateFepEntry::TijEmp));
        if (EWL(elamstats))
        {
            headerContents.flags_dfh |= (enumValueToBitMask(StateFepEntry::WangLandauDelta)
                                         | enumValueToBitMask(StateFepEntry::WangLandauHistogram));
        }
        if ((elamstats == LambdaWeightCalculation::Minvar) || (elamstats == LambdaWeightCalculation::Barker)
            || (elamstats == LambdaWeightCalculation::Metropolis))
        {
            headerContents.flags_dfh |= (enumValueToBitMask(StateFepEntry::Accump)
                                         | enumValueToBitMask(StateFepEntry::Accumm)
                                         | enumValueToBitMask(StateFepEntry::Accump2)
                                         | enumValueToBitMask(StateFepEntry::Accumm2)
                                         | enumValueToBitMask(StateFepEntry::SumMinVar)
                                         | enumValueToBitMask(StateFepEntry::SumVar));
        }
    }

    headerContents.flags_awhh = 0;
    if (state->awhHistory != nullptr && !state->awhHistory->bias.empty())
    {
        headerContents.flags_awhh |= (enumValueToBitMask(StateAwhEntry::InInitial)
                                      | enumValueToBitMask(StateAwhEntry::EquilibrateHistogram)
                                      | enumValueToBitMask(StateAwhEntry::HistogramSize)
                                      | enumValueToBitMask(StateAwhEntry::NumPoints)
                                      | enumValueToBitMask(StateAwhEntry::CoordPoint)
                                      | enumValueToBitMask(StateAwhEntry::UmbrellaGridPoint)
                                      | enumValueToBitMask(StateAwhEntry::UpdateList)
                                      | enumValueToBitMask(StateAwhEntry::LogScaledSampleWeight)
                                      | enumValueToBitMask(StateAwhEntry::NumUpdates)
                                      | enumValueToBitMask(StateAwhEntry::ForceCorrelationGrid));
    }

    do_cpt_header(&serializer, FALSE, nullptr, &headerContents);

    if ((do_cpt_state(&serializer, state, nullptr) < 0)
        || (do_cpt_ekinstate(&serializer, headerContents.flags_eks, &state->ekinstate, nullptr) < 0)
        || (do_cpt_enerhist(&serializer, FALSE, headerContents.flags_enh, enerhist, nullptr) < 0)
        || (doCptPullHist(&serializer, FALSE, headerContents.flagsPullHistory, pullHist, nullptr) < 0)
        || (do_cpt_df_hist(&serializer, headerContents.flags_dfh, headerContents.nlambda, &state->dfhist, nullptr)
            < 0)
        || (do_cpt_EDstate(&serializer, FALSE, headerContents.nED, observablesHistory->edsamHistory.get(), nullptr)
            < 0)
        || (do_cpt_awh(&serializer, FALSE, headerContents.flags_awhh, state->awhHistory.get(), nullptr, CheckPointVersion::CurrentVersion)
            < 0)
        || (do_cpt_swapstate(
                    &serializer, FALSE, headerContents.eSwapCoords, observablesHistory->swapHistory.get(), nullptr)
            < 0)
        || (do_cpt_files(&serializer, FALSE, outputfiles, nullptr, headerContents.file_version) < 0))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    // Checkpointing MDModules
    {
        gmx::KeyValueTreeBuilder          builder;
        gmx::MDModulesWriteCheckpointData mdModulesWriteCheckpoint = { builder.rootObject() };
        mdModulesNotifiers.checkpointingNotifier_.notify(mdModulesWriteCheckpoint);
        auto tree = builder.build();
        gmx::serializeKeyValueTree(tree, &serializer);
    }

    // Checkpointing modular simulator
    {
        modularSimulatorCheckpointData->serialize(&serializer);
    }

    do_cpt_footer(&serializer, headerContents.file_version);
#if GMX_FAHCORE
    /* Always FAH checkpoint immediately after a Gromacs checkpoint.
     *
     * Note that it is critical that we save a FAH checkpoint directly
     * after writing a Gromacs checkpoint.  If the program dies, either
     * by the machine powering off suddenly or the process being,
     * killed, FAH can recover files that have only appended data by
     * truncating them to the last recorded length.  The Gromacs
     * checkpoint does not just append data, it is fully rewritten each
     * time so a crash between moving the new Gromacs checkpoint file in
     * to place and writing a FAH checkpoint is not recoverable.  Thus
     * the time between these operations must be kept as short a
     * possible.
     */
    fcCheckpoint();
#endif
}

static void check_int(FILE* fplog, const char* type, int p, int f, gmx_bool* mm)
{
    bool foundMismatch = (p != f);
    if (!foundMismatch)
    {
        return;
    }
    *mm = TRUE;
    if (fplog)
    {
        fprintf(fplog, "  %s mismatch,\n", type);
        fprintf(fplog, "    current program: %d\n", p);
        fprintf(fplog, "    checkpoint file: %d\n", f);
        fprintf(fplog, "\n");
    }
}

static void check_string(FILE* fplog, const char* type, const char* p, const char* f, gmx_bool* mm)
{
    bool foundMismatch = (std::strcmp(p, f) != 0);
    if (!foundMismatch)
    {
        return;
    }
    *mm = TRUE;
    if (fplog)
    {
        fprintf(fplog, "  %s mismatch,\n", type);
        fprintf(fplog, "    current program: %s\n", p);
        fprintf(fplog, "    checkpoint file: %s\n", f);
        fprintf(fplog, "\n");
    }
}

static void check_match(FILE*                           fplog,
                        const t_commrec*                cr,
                        const ivec                      dd_nc,
                        const CheckpointHeaderContents& headerContents,
                        gmx_bool                        reproducibilityRequested)
{
    /* Note that this check_string on the version will also print a message
     * when only the minor version differs. But we only print a warning
     * message further down with reproducibilityRequested=TRUE.
     */
    gmx_bool versionDiffers = FALSE;
    check_string(fplog, "Version", gmx_version(), headerContents.version, &versionDiffers);

    gmx_bool precisionDiffers = FALSE;
    check_int(fplog, "Double prec.", GMX_DOUBLE, headerContents.double_prec, &precisionDiffers);
    if (precisionDiffers)
    {
        const char msg_precision_difference[] =
                "You are continuing a simulation with a different precision. Not matching\n"
                "mixed/double precision will lead to precision or performance loss.\n";
        if (fplog)
        {
            fprintf(fplog, "%s\n", msg_precision_difference);
        }
    }

    gmx_bool mm = (versionDiffers || precisionDiffers);

    if (reproducibilityRequested)
    {
        check_string(fplog,
                     "Program name",
                     gmx::getProgramContext().fullBinaryPath().string().c_str(),
                     headerContents.fprog,
                     &mm);

        check_int(fplog, "#ranks", cr->commMySim.size(), headerContents.nnodes, &mm);
    }

    if (cr->mpiDefaultCommunicator.size() > 1 && reproducibilityRequested)
    {
        // TODO: These checks are incorrect (see redmine #3309)
        const int numPmeOnlyRanks = (cr->dd ? cr->dd->numPmeOnlyRanks : 0);
        check_int(fplog, "#PME-ranks", numPmeOnlyRanks, headerContents.npme, &mm);

        int npp = cr->mpiDefaultCommunicator.size();
        if (numPmeOnlyRanks >= 0)
        {
            npp -= numPmeOnlyRanks;
        }
        int npp_f = headerContents.nnodes;
        if (headerContents.npme >= 0)
        {
            npp_f -= headerContents.npme;
        }
        if (npp == npp_f)
        {
            check_int(fplog, "#DD-cells[x]", dd_nc[XX], headerContents.dd_nc[XX], &mm);
            check_int(fplog, "#DD-cells[y]", dd_nc[YY], headerContents.dd_nc[YY], &mm);
            check_int(fplog, "#DD-cells[z]", dd_nc[ZZ], headerContents.dd_nc[ZZ], &mm);
        }
    }

    if (mm)
    {
        /* Gromacs should be able to continue from checkpoints between
         * different patch level versions, but we do not guarantee
         * compatibility between different major/minor versions - check this.
         */
        int gmx_major;
        int cpt_major;
        sscanf(gmx_version(), "%5d", &gmx_major);
        int      ret                 = sscanf(headerContents.version, "%5d", &cpt_major);
        gmx_bool majorVersionDiffers = (ret < 1 || gmx_major != cpt_major);

        const char msg_major_version_difference[] =
                "The current GROMACS major version is not identical to the one that\n"
                "generated the checkpoint file. In principle GROMACS does not support\n"
                "continuation from checkpoints between different versions, so we advise\n"
                "against this. If you still want to try your luck we recommend that you use\n"
                "the -noappend flag to keep your output files from the two versions separate.\n"
                "This might also work around errors where the output fields in the energy\n"
                "file have changed between the different versions.\n";

        const char msg_mismatch_notice[] =
                "GROMACS patchlevel, binary or parallel settings differ from previous run.\n"
                "Continuation is exact, but not guaranteed to be binary identical.\n";

        if (majorVersionDiffers)
        {
            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_major_version_difference);
            }
        }
        else if (reproducibilityRequested)
        {
            /* Major & minor versions match at least, but something is different. */
            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_mismatch_notice);
            }
        }
    }
}

static void read_checkpoint(const std::filesystem::path&   fn,
                            t_fileio*                      logfio,
                            const t_commrec*               cr,
                            const ivec                     dd_nc,
                            IntegrationAlgorithm           eIntegrator,
                            int*                           init_fep_state,
                            CheckpointHeaderContents*      headerContents,
                            t_state*                       state,
                            ObservablesHistory*            observablesHistory,
                            gmx_bool                       reproducibilityRequested,
                            const gmx::MDModulesNotifiers& mdModulesNotifiers,
                            gmx::ReadCheckpointDataHolder* modularSimulatorCheckpointData,
                            bool                           useModularSimulator)
{
    char buf[STEPSTRSIZE];
    int  ret;

    XdrSerializer serializer(fn, "r");
    do_cpt_header(&serializer, TRUE, nullptr, headerContents);

    // If we are appending, then we don't want write to the open log
    // file because we still need to compute a checksum for it. In
    // that case, the filehandle will be nullptr. Otherwise, we report
    // to the new log file about the checkpoint file that we are
    // reading from.
    FILE* fplog = gmx_fio_getfp(logfio);
    if (fplog)
    {
        fprintf(fplog, "\n");
        fprintf(fplog, "Reading checkpoint file %s\n", fn.string().c_str());
        fprintf(fplog, "  file generated by:     %s\n", headerContents->fprog);
        fprintf(fplog, "  file generated at:     %s\n", headerContents->ftime);
        fprintf(fplog, "  GROMACS double prec.:  %d\n", headerContents->double_prec);
        fprintf(fplog, "  simulation part #:     %d\n", headerContents->simulation_part);
        fprintf(fplog, "  step:                  %s\n", gmx_step_str(headerContents->step, buf));
        fprintf(fplog, "  time:                  %f\n", headerContents->t);
        fprintf(fplog, "\n");
    }

    if (headerContents->natoms != state->numAtoms())
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d atoms, while the current system consists "
                  "of %d atoms",
                  headerContents->natoms,
                  state->numAtoms());
    }
    if (headerContents->ngtc != state->ngtc)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d T-coupling groups, while the current "
                  "system consists of %d T-coupling groups",
                  headerContents->ngtc,
                  state->ngtc);
    }
    if (headerContents->nnhpres != state->nnhpres)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d NH-pressure-coupling variables, while the "
                  "current system consists of %d NH-pressure-coupling variables",
                  headerContents->nnhpres,
                  state->nnhpres);
    }

    int nlambdaHistory = (state->dfhist ? state->dfhist->nlambda : 0);
    if (headerContents->nlambda != nlambdaHistory)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system with %d lambda states, while the current system "
                  "consists of %d lambda states",
                  headerContents->nlambda,
                  nlambdaHistory);
    }

    init_gtc_state(state,
                   state->ngtc,
                   state->nnhpres,
                   headerContents->nhchainlength); /* need to keep this here to keep the tpr format working */
    /* write over whatever was read; we use the number of Nose-Hoover chains from the checkpoint */

    if (headerContents->eIntegrator != eIntegrator)
    {
        gmx_fatal(FARGS,
                  "Cannot change integrator during a checkpoint restart. Perhaps you should make a "
                  "new .tpr with grompp -f new.mdp -t %s",
                  fn.string().c_str());
    }

    // For modular simulator, no state object is populated, so we cannot do this check here!
    if (headerContents->flags_state != state->flags() && !useModularSimulator)
    {
        gmx_fatal(FARGS,
                  "Cannot change a simulation algorithm during a checkpoint restart. Perhaps you "
                  "should make a new .tpr with grompp -f new.mdp -t %s",
                  fn.string().c_str());
    }

    GMX_RELEASE_ASSERT(!(headerContents->isModularSimulatorCheckpoint && !useModularSimulator),
                       "Checkpoint file was written by modular simulator, but the current "
                       "simulation uses the legacy simulator.\n\n"
                       "Try the following steps:\n"
                       "1. Make sure the GMX_DISABLE_MODULAR_SIMULATOR environment variable is not "
                       "set to return to the default behavior. Retry running the simulation.\n"
                       "2. If the problem persists, set the environment variable "
                       "GMX_USE_MODULAR_SIMULATOR=ON to overwrite the default behavior and use "
                       "modular simulator for all implemented use cases.");
    GMX_RELEASE_ASSERT(!(!headerContents->isModularSimulatorCheckpoint && useModularSimulator),
                       "Checkpoint file was written by legacy simulator, but the current "
                       "simulation uses the modular simulator.\n\n"
                       "Try the following steps:\n"
                       "1. Make sure the GMX_USE_MODULAR_SIMULATOR environment variable is not set "
                       "to return to the default behavior. Retry running the simulation.\n"
                       "2. If the problem persists, set the environment variable "
                       "GMX_DISABLE_MODULAR_SIMULATOR=ON to overwrite the default behavior and use "
                       "legacy simulator for all implemented use cases.");

    if (MAIN(cr))
    {
        check_match(fplog, cr, dd_nc, *headerContents, reproducibilityRequested);
    }

    ret             = do_cpt_state(&serializer, state, nullptr);
    *init_fep_state = state->fep_state; /* there should be a better way to do this than setting it
                                           here. Investigate for 5.0. */
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(&serializer, headerContents->flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }
    state->ekinstate.hasReadEkinState =
            (((headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::EkinHalfStep)) != 0)
             || ((headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::EkinFullStep)) != 0)
             || ((headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::EkinHalfStepOld)) != 0)
             || (((headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::EkinNoseHooverScaleFullStep))
                  | (headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::EkinNoseHooverScaleHalfStep))
                  | (headerContents->flags_eks & enumValueToBitMask(StateKineticEntry::VelocityScale)))
                 != 0));

    if (headerContents->flags_enh && observablesHistory->energyHistory == nullptr)
    {
        observablesHistory->energyHistory = std::make_unique<energyhistory_t>();
    }
    ret = do_cpt_enerhist(
            &serializer, TRUE, headerContents->flags_enh, observablesHistory->energyHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->flagsPullHistory)
    {
        if (observablesHistory->pullHistory == nullptr)
        {
            observablesHistory->pullHistory = std::make_unique<PullHistory>();
        }
        ret = doCptPullHist(&serializer,
                            TRUE,
                            headerContents->flagsPullHistory,
                            observablesHistory->pullHistory.get(),
                            nullptr);
        if (ret)
        {
            cp_error();
        }
    }

    if (headerContents->file_version < CheckPointVersion::Version45)
    {
        gmx_fatal(FARGS,
                  "Continuing from checkpoint files written before GROMACS 4.5 is not supported");
    }

    ret = do_cpt_df_hist(
            &serializer, headerContents->flags_dfh, headerContents->nlambda, &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->nED > 0 && observablesHistory->edsamHistory == nullptr)
    {
        observablesHistory->edsamHistory = std::make_unique<edsamhistory_t>(edsamhistory_t{});
    }
    ret = do_cpt_EDstate(
            &serializer, TRUE, headerContents->nED, observablesHistory->edsamHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->flags_awhh != 0 && state->awhHistory == nullptr)
    {
        state->awhHistory = std::make_shared<gmx::AwhHistory>();
    }
    ret = do_cpt_awh(&serializer,
                     TRUE,
                     headerContents->flags_awhh,
                     state->awhHistory.get(),
                     nullptr,
                     headerContents->file_version);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->eSwapCoords != SwapType::No && observablesHistory->swapHistory == nullptr)
    {
        observablesHistory->swapHistory = std::make_unique<swaphistory_t>(swaphistory_t{});
    }
    ret = do_cpt_swapstate(
            &serializer, TRUE, headerContents->eSwapCoords, observablesHistory->swapHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    std::vector<gmx_file_position_t> outputfiles;
    ret = do_cpt_files(&serializer, TRUE, &outputfiles, nullptr, headerContents->file_version);
    if (ret)
    {
        cp_error();
    }
    do_cpt_mdmodules(headerContents->file_version, &serializer, mdModulesNotifiers, nullptr);
    if (headerContents->file_version >= CheckPointVersion::ModularSimulator)
    {
        modularSimulatorCheckpointData->deserialize(&serializer);
    }
    ret = do_cpt_footer(&serializer, headerContents->file_version);
    if (ret)
    {
        cp_error();
    }
}


void load_checkpoint(const std::filesystem::path&   fn,
                     t_fileio*                      logfio,
                     const t_commrec*               cr,
                     const ivec                     dd_nc,
                     t_inputrec*                    ir,
                     t_state*                       state,
                     ObservablesHistory*            observablesHistory,
                     gmx_bool                       reproducibilityRequested,
                     const gmx::MDModulesNotifiers& mdModulesNotifiers,
                     gmx::ReadCheckpointDataHolder* modularSimulatorCheckpointData,
                     bool                           useModularSimulator)
{
    CheckpointHeaderContents headerContents;
    if (cr->isSimulationMainRank())
    {
        /* Read the state from the checkpoint file */
        read_checkpoint(fn,
                        logfio,
                        cr,
                        dd_nc,
                        ir->eI,
                        &(ir->fepvals->init_fep_state),
                        &headerContents,
                        state,
                        observablesHistory,
                        reproducibilityRequested,
                        mdModulesNotifiers,
                        modularSimulatorCheckpointData,
                        useModularSimulator);
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(headerContents.step), &headerContents.step, cr->mpiDefaultCommunicator.comm());
        gmx::MDModulesCheckpointReadingBroadcast broadcastCheckPointData = {
            cr->mpiDefaultCommunicator.comm(), PAR(cr)
        };
        mdModulesNotifiers.checkpointingNotifier_.notify(broadcastCheckPointData);
    }
    ir->bContinuation = TRUE;
    if (ir->nsteps >= 0)
    {
        // TODO Should the following condition be <=? Currently if you
        // pass a checkpoint written by an normal completion to a restart,
        // mdrun will read all input, does some work but no steps, and
        // write successful output. But perhaps that is not desirable.
        // Note that we do not intend to support the use of mdrun
        // -nsteps to circumvent this condition.
        if (ir->nsteps + ir->init_step < headerContents.step)
        {
            char        buf[STEPSTRSIZE];
            std::string message =
                    gmx::formatString("The input requested %s steps, ", gmx_step_str(ir->nsteps, buf));
            if (ir->init_step > 0)
            {
                message += gmx::formatString("starting from step %s, ", gmx_step_str(ir->init_step, buf));
            }
            message += gmx::formatString(
                    "however the checkpoint "
                    "file has already reached step %s. The simulation will not "
                    "proceed, because either your simulation is already complete, "
                    "or your combination of input files don't match.",
                    gmx_step_str(headerContents.step, buf));
            gmx_fatal(FARGS, "%s", message.c_str());
        }
        ir->nsteps += ir->init_step - headerContents.step;
    }
    ir->init_step       = headerContents.step;
    ir->simulation_part = headerContents.simulation_part + 1;
}

void read_checkpoint_part_and_step(const std::filesystem::path& filename, int* simulation_part, int64_t* step)
{
    if (filename.empty() || !gmx_fexist(filename))
    {
        *simulation_part = 0;
        *step            = 0;
        return;
    }
    XdrSerializer serializer(filename, "r");

    CheckpointHeaderContents headerContents;
    do_cpt_header(&serializer, TRUE, nullptr, &headerContents);
    *simulation_part = headerContents.simulation_part;
    *step            = headerContents.step;
}

static CheckpointHeaderContents read_checkpoint_data(XdrSerializer*                    serializer,
                                                     t_state*                          state,
                                                     std::vector<gmx_file_position_t>* outputfiles,
                                                     gmx::ReadCheckpointDataHolder* modularSimulatorCheckpointData)
{
    CheckpointHeaderContents headerContents;
    do_cpt_header(serializer, TRUE, nullptr, &headerContents);
    state->changeNumAtoms(headerContents.natoms);
    state->ngtc          = headerContents.ngtc;
    state->nnhpres       = headerContents.nnhpres;
    state->nhchainlength = headerContents.nhchainlength;
    state->setFlags(headerContents.flags_state);
    int ret = do_cpt_state(serializer, state, nullptr);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(serializer, headerContents.flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(serializer, TRUE, headerContents.flags_enh, &enerhist, nullptr);
    if (ret)
    {
        cp_error();
    }
    PullHistory pullHist = {};
    ret = doCptPullHist(serializer, TRUE, headerContents.flagsPullHistory, &pullHist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_df_hist(
            serializer, headerContents.flags_dfh, headerContents.nlambda, &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    edsamhistory_t edsamhist = {};
    ret = do_cpt_EDstate(serializer, TRUE, headerContents.nED, &edsamhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_awh(serializer,
                     TRUE,
                     headerContents.flags_awhh,
                     state->awhHistory.get(),
                     nullptr,
                     headerContents.file_version);
    if (ret)
    {
        cp_error();
    }

    swaphistory_t swaphist = {};
    ret = do_cpt_swapstate(serializer, TRUE, headerContents.eSwapCoords, &swaphist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_files(serializer, TRUE, outputfiles, nullptr, headerContents.file_version);

    if (ret)
    {
        cp_error();
    }
    gmx::MDModulesNotifiers mdModuleNotifiers;
    do_cpt_mdmodules(headerContents.file_version, serializer, mdModuleNotifiers, nullptr);
    if (headerContents.file_version >= CheckPointVersion::ModularSimulator)
    {
        // Store modular checkpoint data into modularSimulatorCheckpointData
        modularSimulatorCheckpointData->deserialize(serializer);
    }
    ret = do_cpt_footer(serializer, headerContents.file_version);
    if (ret)
    {
        cp_error();
    }
    return headerContents;
}

void read_checkpoint_trxframe(const std::filesystem::path& filename, t_trxframe* fr)
{
    XdrSerializer                    serializer(filename, "r");
    t_state                          state;
    std::vector<gmx_file_position_t> outputfiles;
    gmx::ReadCheckpointDataHolder    modularSimulatorCheckpointData;
    CheckpointHeaderContents         headerContents =
            read_checkpoint_data(&serializer, &state, &outputfiles, &modularSimulatorCheckpointData);
    if (headerContents.isModularSimulatorCheckpoint)
    {
        gmx::ModularSimulator::readCheckpointToTrxFrame(fr, &modularSimulatorCheckpointData, headerContents);
        return;
    }

    fr->natoms    = state.numAtoms();
    fr->bStep     = TRUE;
    fr->step      = int64_to_int(headerContents.step, "conversion of checkpoint to trajectory");
    fr->bTime     = TRUE;
    fr->time      = headerContents.t;
    fr->bLambda   = TRUE;
    fr->lambda    = state.lambda[FreeEnergyPerturbationCouplingType::Fep];
    fr->fep_state = state.fep_state;
    fr->bAtoms    = FALSE;
    fr->bX        = state.hasEntry(StateEntry::X);
    if (fr->bX)
    {
        fr->x = makeRvecArray(state.x, state.numAtoms());
    }
    fr->bV = state.hasEntry(StateEntry::V);
    if (fr->bV)
    {
        fr->v = makeRvecArray(state.v, state.numAtoms());
    }
    fr->bF   = FALSE;
    fr->bBox = state.hasEntry(StateEntry::Box);
    if (fr->bBox)
    {
        copy_mat(state.box, fr->box);
    }
}

void list_checkpoint(const std::filesystem::path& fn, FILE* out)
{
    int ret;

    t_state state;

    XdrSerializer            serializer(fn, "r");
    CheckpointHeaderContents headerContents;
    do_cpt_header(&serializer, TRUE, out, &headerContents);
    state.changeNumAtoms(headerContents.natoms);
    state.ngtc          = headerContents.ngtc;
    state.nnhpres       = headerContents.nnhpres;
    state.nhchainlength = headerContents.nhchainlength;
    state.setFlags(headerContents.flags_state);
    ret = do_cpt_state(&serializer, &state, out);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(&serializer, headerContents.flags_eks, &state.ekinstate, out);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(&serializer, TRUE, headerContents.flags_enh, &enerhist, out);

    if (ret == 0)
    {
        PullHistory pullHist = {};
        ret = doCptPullHist(&serializer, TRUE, headerContents.flagsPullHistory, &pullHist, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_df_hist(
                &serializer, headerContents.flags_dfh, headerContents.nlambda, &state.dfhist, out);
    }

    if (ret == 0)
    {
        edsamhistory_t edsamhist = {};
        ret = do_cpt_EDstate(&serializer, TRUE, headerContents.nED, &edsamhist, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_awh(&serializer,
                         TRUE,
                         headerContents.flags_awhh,
                         state.awhHistory.get(),
                         out,
                         headerContents.file_version);
    }

    if (ret == 0)
    {
        swaphistory_t swaphist = {};
        ret = do_cpt_swapstate(&serializer, TRUE, headerContents.eSwapCoords, &swaphist, out);
    }

    if (ret == 0)
    {
        std::vector<gmx_file_position_t> outputfiles;
        ret = do_cpt_files(&serializer, TRUE, &outputfiles, out, headerContents.file_version);
    }
    gmx::MDModulesNotifiers mdModuleNotifiers;
    do_cpt_mdmodules(headerContents.file_version, &serializer, mdModuleNotifiers, out);
    if (headerContents.file_version >= CheckPointVersion::ModularSimulator)
    {
        gmx::ReadCheckpointDataHolder modularSimulatorCheckpointData;
        modularSimulatorCheckpointData.deserialize(&serializer);
        modularSimulatorCheckpointData.dump(out);
    }

    if (ret == 0)
    {
        ret = do_cpt_footer(&serializer, headerContents.file_version);
    }

    if (ret)
    {
        cp_warning(out);
    }
}

/* This routine cannot print tons of data, since it is called before the log file is opened. */
CheckpointHeaderContents read_checkpoint_simulation_part_and_filenames(const std::filesystem::path& filename,
                                                                       std::vector<gmx_file_position_t>* outputfiles)
{
    XdrSerializer                 serializer(filename, "r");
    t_state                       state;
    gmx::ReadCheckpointDataHolder modularSimulatorCheckpointData;
    CheckpointHeaderContents      headerContents =
            read_checkpoint_data(&serializer, &state, outputfiles, &modularSimulatorCheckpointData);
    return headerContents;
}
