/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2011,2012 by the GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

/* The source code in this file should be thread-safe.
   Please keep it that way. */
#include "gmxpre.h"

#include "checkpoint.h"

#include "config.h"

#include <cerrno>
#include <cstdlib>
#include <cstring>

#include <array>
#include <memory>

#include "buildinfo.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/awh_correlation_history.h"
#include "gromacs/mdtypes/awh_history.h"
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
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/int64_to_int.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/mdmodulenotification.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/txtdump.h"

#if GMX_FAHCORE
#    include "corewrap.h"
#endif

#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

/*! \brief Enum of values that describe the contents of a cpt file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element just above cptv_Count
 * in this enumeration, and write code below that does the right thing
 * according to the value of file_version.
 */
enum cptv
{
    cptv_Unknown = 17,                  /**< Version before numbering scheme */
    cptv_RemoveBuildMachineInformation, /**< remove functionality that makes mdrun builds non-reproducible */
    cptv_ComPrevStepAsPullGroupReference, /**< Allow using COM of previous step as pull group PBC reference */
    cptv_PullAverage, /**< Added possibility to output average pull force and position */
    cptv_MdModules,   /**< Added checkpointing for MdModules */
    cptv_Count        /**< the total number of cptv versions */
};

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
 * development branch changes, due to an extension of the cptv enum above.
 * Backward compatibility for reading old run input files is maintained
 * by checking this version number against that of the file and then using
 * the correct code path. */
static const int cpt_version = cptv_Count - 1;


const char* est_names[estNR] = { "FE-lambda",
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
                                 "barostat-integral" };

enum
{
    eeksEKIN_N,
    eeksEKINH,
    eeksDEKINDL,
    eeksMVCOS,
    eeksEKINF,
    eeksEKINO,
    eeksEKINSCALEF,
    eeksEKINSCALEH,
    eeksVSCALE,
    eeksEKINTOTAL,
    eeksNR
};

static const char* eeks_names[eeksNR] = { "Ekin_n",         "Ekinh",          "dEkindlambda",
                                          "mv_cos",         "Ekinf",          "Ekinh_old",
                                          "EkinScaleF_NHC", "EkinScaleH_NHC", "Vscale_NHC",
                                          "Ekin_Total" };

enum
{
    eenhENERGY_N,
    eenhENERGY_AVER,
    eenhENERGY_SUM,
    eenhENERGY_NSUM,
    eenhENERGY_SUM_SIM,
    eenhENERGY_NSUM_SIM,
    eenhENERGY_NSTEPS,
    eenhENERGY_NSTEPS_SIM,
    eenhENERGY_DELTA_H_NN,
    eenhENERGY_DELTA_H_LIST,
    eenhENERGY_DELTA_H_STARTTIME,
    eenhENERGY_DELTA_H_STARTLAMBDA,
    eenhNR
};

enum
{
    epullhPULL_NUMCOORDINATES,
    epullhPULL_NUMGROUPS,
    epullhPULL_NUMVALUESINXSUM,
    epullhPULL_NUMVALUESINFSUM,
    epullhNR
};

enum
{
    epullcoordh_VALUE_REF_SUM,
    epullcoordh_VALUE_SUM,
    epullcoordh_DR01_SUM,
    epullcoordh_DR23_SUM,
    epullcoordh_DR45_SUM,
    epullcoordh_FSCAL_SUM,
    epullcoordh_DYNAX_SUM,
    epullcoordh_NR
};

enum
{
    epullgrouph_X_SUM,
    epullgrouph_NR
};

static const char* eenh_names[eenhNR] = { "energy_n",
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
                                          "energy_delta_h_start_lambda" };

static const char* ePullhNames[epullhNR] = { "pullhistory_numcoordinates", "pullhistory_numgroups",
                                             "pullhistory_numvaluesinxsum",
                                             "pullhistory_numvaluesinfsum" };

/* free energy history variables -- need to be preserved over checkpoint */
enum
{
    edfhBEQUIL,
    edfhNATLAMBDA,
    edfhWLHISTO,
    edfhWLDELTA,
    edfhSUMWEIGHTS,
    edfhSUMDG,
    edfhSUMMINVAR,
    edfhSUMVAR,
    edfhACCUMP,
    edfhACCUMM,
    edfhACCUMP2,
    edfhACCUMM2,
    edfhTIJ,
    edfhTIJEMP,
    edfhNR
};
/* free energy history variable names  */
static const char* edfh_names[edfhNR] = { "bEquilibrated",
                                          "N_at_state",
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
                                          "Tij_empirical" };

/* AWH biasing history variables */
enum
{
    eawhhIN_INITIAL,
    eawhhEQUILIBRATEHISTOGRAM,
    eawhhHISTSIZE,
    eawhhNPOINTS,
    eawhhCOORDPOINT,
    eawhhUMBRELLAGRIDPOINT,
    eawhhUPDATELIST,
    eawhhLOGSCALEDSAMPLEWEIGHT,
    eawhhNUMUPDATES,
    eawhhFORCECORRELATIONGRID,
    eawhhNR
};

static const char* eawhh_names[eawhhNR] = { "awh_in_initial", "awh_equilibrateHistogram",
                                            "awh_histsize",   "awh_npoints",
                                            "awh_coordpoint", "awh_umbrellaGridpoint",
                                            "awh_updatelist", "awh_logScaledSampleWeight",
                                            "awh_numupdates", "awh_forceCorrelationGrid" };

enum
{
    epullsPREVSTEPCOM,
    epullsNR
};

static const char* epull_prev_step_com_names[epullsNR] = { "Pull groups prev step COM" };


//! Higher level vector element type, only used for formatting checkpoint dumps
enum class CptElementType
{
    integer,  //!< integer
    real,     //!< float or double, not linked to precision of type real
    real3,    //!< float[3] or double[3], not linked to precision of type real
    matrix3x3 //!< float[3][3] or double[3][3], not linked to precision of type real
};

//! \brief Parts of the checkpoint state, only used for reporting
enum class StatePart
{
    microState,         //!< The microstate of the simulated system
    kineticEnergy,      //!< Kinetic energy, needed for T/P-coupling state
    energyHistory,      //!< Energy observable statistics
    freeEnergyHistory,  //!< Free-energy state and observable statistics
    accWeightHistogram, //!< Accelerated weight histogram method state
    pullState,          //!< COM of previous step.
    pullHistory         //!< Pull history statistics (sums since last written output)
};

//! \brief Return the name of a checkpoint entry based on part and part entry
static const char* entryName(StatePart part, int ecpt)
{
    switch (part)
    {
        case StatePart::microState: return est_names[ecpt];
        case StatePart::kineticEnergy: return eeks_names[ecpt];
        case StatePart::energyHistory: return eenh_names[ecpt];
        case StatePart::freeEnergyHistory: return edfh_names[ecpt];
        case StatePart::accWeightHistogram: return eawhh_names[ecpt];
        case StatePart::pullState: return epull_prev_step_com_names[ecpt];
        case StatePart::pullHistory: return ePullhNames[ecpt];
    }

    return nullptr;
}

static void cp_warning(FILE* fp)
{
    fprintf(fp, "\nWARNING: Checkpoint file is corrupted or truncated\n\n");
}

[[noreturn]] static void cp_error()
{
    gmx_fatal(FARGS, "Checkpoint file corrupted/truncated, or maybe you are out of disk space?");
}

static void do_cpt_string_err(XDR* xd, const char* desc, gmx::ArrayRef<char> s, FILE* list)
{
    char* data = s.data();
    if (xdr_string(xd, &data, s.size()) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, data);
    }
}

static int do_cpt_int(XDR* xd, const char* desc, int* i, FILE* list)
{
    if (xdr_int(xd, i) == 0)
    {
        return -1;
    }
    if (list)
    {
        fprintf(list, "%s = %d\n", desc, *i);
    }
    return 0;
}

static int do_cpt_u_chars(XDR* xd, const char* desc, int n, unsigned char* i, FILE* list)
{
    if (list)
    {
        fprintf(list, "%s = ", desc);
    }
    bool_t res = 1;
    for (int j = 0; j < n && res; j++)
    {
        res &= xdr_u_char(xd, &i[j]);
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

static void do_cpt_int_err(XDR* xd, const char* desc, int* i, FILE* list)
{
    if (do_cpt_int(xd, desc, i, list) < 0)
    {
        cp_error();
    }
}

static void do_cpt_bool_err(XDR* xd, const char* desc, bool* b, FILE* list)
{
    int i = static_cast<int>(*b);

    if (do_cpt_int(xd, desc, &i, list) < 0)
    {
        cp_error();
    }

    *b = (i != 0);
}

static void do_cpt_step_err(XDR* xd, const char* desc, int64_t* i, FILE* list)
{
    char buf[STEPSTRSIZE];

    if (xdr_int64(xd, i) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, gmx_step_str(*i, buf));
    }
}

static void do_cpt_double_err(XDR* xd, const char* desc, double* f, FILE* list)
{
    if (xdr_double(xd, f) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %f\n", desc, *f);
    }
}

static void do_cpt_real_err(XDR* xd, real* f)
{
#if GMX_DOUBLE
    bool_t res = xdr_double(xd, f);
#else
    bool_t res = xdr_float(xd, f);
#endif
    if (res == 0)
    {
        cp_error();
    }
}

static void do_cpt_n_rvecs_err(XDR* xd, const char* desc, int n, rvec f[], FILE* list)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            do_cpt_real_err(xd, &f[i][j]);
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
    static const int value = xdr_datatype_int;
};

template<>
struct xdr_type<float>
{
    static const int value = xdr_datatype_float;
};

template<>
struct xdr_type<double>
{
    static const int value = xdr_datatype_double;
};

//! \brief Returns size in byte of an xdr_datatype
static inline unsigned int sizeOfXdrType(int xdrType)
{
    switch (xdrType)
    {
        case xdr_datatype_int: return sizeof(int);
        case xdr_datatype_float: return sizeof(float);
        case xdr_datatype_double: return sizeof(double);
        default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
    }

    return 0;
}

//! \brief Returns the XDR process function for i/o of an XDR type
static inline xdrproc_t xdrProc(int xdrType)
{
    switch (xdrType)
    {
        case xdr_datatype_int: return reinterpret_cast<xdrproc_t>(xdr_int);
        case xdr_datatype_float: return reinterpret_cast<xdrproc_t>(xdr_float);
        case xdr_datatype_double: return reinterpret_cast<xdrproc_t>(xdr_double);
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
static bool_t listXdrVector(XDR* xd, StatePart part, int ecpt, int nf, int xdrType, FILE* list, CptElementType cptElementType)
{
    bool_t res = 0;

    const unsigned int elemSize = sizeOfXdrType(xdrType);
    std::vector<char>  data(nf * elemSize);
    res = xdr_vector(xd, data.data(), nf, elemSize, xdrProc(xdrType));

    if (list != nullptr)
    {
        switch (xdrType)
        {
            case xdr_datatype_int:
                pr_ivec(list, 0, entryName(part, ecpt), reinterpret_cast<const int*>(data.data()), nf, TRUE);
                break;
            case xdr_datatype_float:
#if !GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    pr_rvecs(list, 0, entryName(part, ecpt),
                             reinterpret_cast<const rvec*>(data.data()), nf / 3);
                }
                else
#endif
                {
                    /* Note: With double precision code dumping a single precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    pr_fvec(list, 0, entryName(part, ecpt),
                            reinterpret_cast<const float*>(data.data()), nf, TRUE);
                }
                break;
            case xdr_datatype_double:
#if GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    pr_rvecs(list, 0, entryName(part, ecpt),
                             reinterpret_cast<const rvec*>(data.data()), nf / 3);
                }
                else
#endif
                {
                    /* Note: With single precision code dumping a double precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    pr_dvec(list, 0, entryName(part, ecpt),
                            reinterpret_cast<const double*>(data.data()), nf, TRUE);
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
 * If nval n<0, *nptr (with v) or vector->size() is used. On read using v,
 * the value is stored in nptr
 */
template<typename T, typename AllocatorType>
static int doVectorLow(XDR*                           xd,
                       StatePart                      part,
                       int                            ecpt,
                       int                            sflags,
                       int                            nval,
                       int*                           nptr,
                       T**                            v,
                       std::vector<T, AllocatorType>* vector,
                       FILE*                          list,
                       CptElementType                 cptElementType)
{
    GMX_RELEASE_ASSERT(list != nullptr || (v != nullptr && vector == nullptr)
                               || (v == nullptr && vector != nullptr),
                       "Without list, we should have exactly one of v and vector != NULL");

    bool_t res = 0;

    int numElemInTheFile;
    if (list == nullptr)
    {
        if (nval >= 0)
        {
            GMX_RELEASE_ASSERT(nptr == nullptr, "With nval>=0 we should have nptr==NULL");
            numElemInTheFile = nval;
        }
        else
        {
            if (v != nullptr)
            {
                GMX_RELEASE_ASSERT(nptr != nullptr, "With nval<0 we should have nptr!=NULL");
                numElemInTheFile = *nptr;
            }
            else
            {
                numElemInTheFile = vector->size();
            }
        }
    }
    /* Read/write the vector element count */
    res = xdr_int(xd, &numElemInTheFile);
    if (res == 0)
    {
        return -1;
    }
    /* Read/write the element data type */
    constexpr int xdrTypeInTheCode = xdr_type<T>::value;
    int           xdrTypeInTheFile = xdrTypeInTheCode;
    res                            = xdr_int(xd, &xdrTypeInTheFile);
    if (res == 0)
    {
        return -1;
    }

    if (list == nullptr && (sflags & (1 << ecpt)))
    {
        if (nval >= 0)
        {
            if (numElemInTheFile != nval)
            {
                gmx_fatal(FARGS,
                          "Count mismatch for state entry %s, code count is %d, file count is %d\n",
                          entryName(part, ecpt), nval, numElemInTheFile);
            }
        }
        else if (nptr != nullptr)
        {
            *nptr = numElemInTheFile;
        }

        bool typesMatch = (xdrTypeInTheFile == xdrTypeInTheCode);
        if (!typesMatch)
        {
            char buf[STRLEN];
            sprintf(buf, "mismatch for state entry %s, code precision is %s, file precision is %s",
                    entryName(part, ecpt), xdr_datatype_names[xdrTypeInTheCode],
                    xdr_datatype_names[xdrTypeInTheFile]);

            /* Matching int and real should never occur, but check anyhow */
            if (xdrTypeInTheFile == xdr_datatype_int || xdrTypeInTheCode == xdr_datatype_int)
            {
                gmx_fatal(FARGS,
                          "Type %s: incompatible checkpoint formats or corrupted checkpoint file.", buf);
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
        res = xdr_vector(xd, vChar, numElemInTheFile, sizeOfXdrType(xdrTypeInTheFile),
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
        res = listXdrVector(xd, part, ecpt, numElemInTheFile, xdrTypeInTheFile, list, cptElementType);
    }

    return 0;
}

//! \brief Read/Write a std::vector, on read checks the number of elements matches \p numElements, if specified.
template<typename T>
static int
doVector(XDR* xd, StatePart part, int ecpt, int sflags, std::vector<T>* vector, FILE* list, int numElements = -1)
{
    return doVectorLow<T>(xd, part, ecpt, sflags, numElements, nullptr, nullptr, vector, list,
                          CptElementType::real);
}

//! \brief Read/Write an ArrayRef<real>.
static int doRealArrayRef(XDR* xd, StatePart part, int ecpt, int sflags, gmx::ArrayRef<real> vector, FILE* list)
{
    real* v_real = vector.data();
    return doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, vector.size(), nullptr,
                                                   &v_real, nullptr, list, CptElementType::real);
}

//! Convert from view of RVec to view of real.
static gmx::ArrayRef<real> realArrayRefFromRVecArrayRef(gmx::ArrayRef<gmx::RVec> ofRVecs)
{
    return gmx::arrayRefFromArray<real>(reinterpret_cast<real*>(ofRVecs.data()), ofRVecs.size() * DIM);
}

//! \brief Read/Write a PaddedVector whose value_type is RVec.
template<typename PaddedVectorOfRVecType>
static int
doRvecVector(XDR* xd, StatePart part, int ecpt, int sflags, PaddedVectorOfRVecType* v, int numAtoms, FILE* list)
{
    const int numReals = numAtoms * DIM;

    if (list == nullptr)
    {
        GMX_RELEASE_ASSERT(
                sflags & (1 << ecpt),
                "When not listing, the flag for the entry should be set when requesting i/o");
        GMX_RELEASE_ASSERT(v->size() == numAtoms, "v should have sufficient size for numAtoms");

        return doRealArrayRef(xd, part, ecpt, sflags, realArrayRefFromRVecArrayRef(makeArrayRef(*v)), list);
    }
    else
    {
        // Use the rebind facility to change the value_type of the
        // allocator from RVec to real.
        using realAllocator =
                typename std::allocator_traits<typename PaddedVectorOfRVecType::allocator_type>::template rebind_alloc<real>;
        return doVectorLow<real, realAllocator>(xd, part, ecpt, sflags, numReals, nullptr, nullptr,
                                                nullptr, list, CptElementType::real);
    }
}

/* This function stores n along with the reals for reading,
 * but on reading it assumes that n matches the value in the checkpoint file,
 * a fatal error is generated when this is not the case.
 */
static int do_cpte_reals(XDR* xd, StatePart part, int ecpt, int sflags, int n, real** v, FILE* list)
{
    return doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, n, nullptr, v, nullptr,
                                                   list, CptElementType::real);
}

/* This function does the same as do_cpte_reals,
 * except that on reading it ignores the passed value of *n
 * and stores the value read from the checkpoint file in *n.
 */
static int do_cpte_n_reals(XDR* xd, StatePart part, int ecpt, int sflags, int* n, real** v, FILE* list)
{
    return doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, -1, n, v, nullptr, list,
                                                   CptElementType::real);
}

static int do_cpte_real(XDR* xd, StatePart part, int ecpt, int sflags, real* r, FILE* list)
{
    return doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, 1, nullptr, &r, nullptr,
                                                   list, CptElementType::real);
}

static int do_cpte_ints(XDR* xd, StatePart part, int ecpt, int sflags, int n, int** v, FILE* list)
{
    return doVectorLow<int, std::allocator<int>>(xd, part, ecpt, sflags, n, nullptr, v, nullptr,
                                                 list, CptElementType::integer);
}

static int do_cpte_int(XDR* xd, StatePart part, int ecpt, int sflags, int* i, FILE* list)
{
    return do_cpte_ints(xd, part, ecpt, sflags, 1, &i, list);
}

static int do_cpte_bool(XDR* xd, StatePart part, int ecpt, int sflags, bool* b, FILE* list)
{
    int i   = static_cast<int>(*b);
    int ret = do_cpte_int(xd, part, ecpt, sflags, &i, list);
    *b      = (i != 0);
    return ret;
}

static int do_cpte_doubles(XDR* xd, StatePart part, int ecpt, int sflags, int n, double** v, FILE* list)
{
    return doVectorLow<double, std::allocator<double>>(xd, part, ecpt, sflags, n, nullptr, v,
                                                       nullptr, list, CptElementType::real);
}

static int do_cpte_double(XDR* xd, StatePart part, int ecpt, int sflags, double* r, FILE* list)
{
    return do_cpte_doubles(xd, part, ecpt, sflags, 1, &r, list);
}

static int do_cpte_matrix(XDR* xd, StatePart part, int ecpt, int sflags, matrix v, FILE* list)
{
    real* vr;
    int   ret;

    vr  = &(v[0][0]);
    ret = doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, DIM * DIM, nullptr, &vr,
                                                  nullptr, nullptr, CptElementType::matrix3x3);

    if (list && ret == 0)
    {
        pr_rvecs(list, 0, entryName(part, ecpt), v, DIM);
    }

    return ret;
}


static int do_cpte_nmatrix(XDR* xd, StatePart part, int ecpt, int sflags, int n, real** v, FILE* list)
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
        reti = doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, n, nullptr, &(v[i]),
                                                       nullptr, nullptr, CptElementType::matrix3x3);
        if (list && reti == 0)
        {
            sprintf(name, "%s[%d]", entryName(part, ecpt), i);
            pr_reals(list, 0, name, v[i], n);
        }
        if (reti != 0)
        {
            ret = reti;
        }
    }
    return ret;
}

static int do_cpte_matrices(XDR* xd, StatePart part, int ecpt, int sflags, int n, matrix** v, FILE* list)
{
    bool_t  res = 0;
    matrix *vp, *va = nullptr;
    real*   vr;
    int     nf, i, j, k;
    int     ret;

    nf  = n;
    res = xdr_int(xd, &nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == nullptr && nf != n)
    {
        gmx_fatal(FARGS, "Count mismatch for state entry %s, code count is %d, file count is %d\n",
                  entryName(part, ecpt), n, nf);
    }
    if (list || !(sflags & (1 << ecpt)))
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
    ret = doVectorLow<real, std::allocator<real>>(xd, part, ecpt, sflags, nf * DIM * DIM, nullptr,
                                                  &vr, nullptr, nullptr, CptElementType::matrix3x3);
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
            pr_rvecs(list, 0, entryName(part, ecpt), vp[i], DIM);
        }
    }
    if (va)
    {
        sfree(va);
    }

    return ret;
}

static void do_cpt_header(XDR* xd, gmx_bool bRead, FILE* list, CheckpointHeaderContents* contents)
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
    res = xdr_int(xd, &magic);
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
                  magic, CPT_MAGIC1);
    }
    char fhost[255];
    if (!bRead)
    {
        gmx_gethostname(fhost, 255);
    }
    do_cpt_string_err(xd, "GROMACS version", contents->version, list);
    // The following fields are no longer ever written with meaningful
    // content, but because they precede the file version, there is no
    // good way for new code to read the old and new formats, nor a
    // good way for old code to avoid giving an error while reading a
    // new format. So we read and write a field that no longer has a
    // purpose.
    do_cpt_string_err(xd, "GROMACS build time UNUSED", contents->btime_UNUSED, list);
    do_cpt_string_err(xd, "GROMACS build user UNUSED", contents->buser_UNUSED, list);
    do_cpt_string_err(xd, "GROMACS build host UNUSED", contents->bhost_UNUSED, list);
    do_cpt_string_err(xd, "generating program", contents->fprog, list);
    do_cpt_string_err(xd, "generation time", contents->ftime, list);
    contents->file_version = cpt_version;
    do_cpt_int_err(xd, "checkpoint file version", &contents->file_version, list);
    if (contents->file_version > cpt_version)
    {
        gmx_fatal(FARGS,
                  "Attempting to read a checkpoint file of version %d with code of version %d\n",
                  contents->file_version, cpt_version);
    }
    if (contents->file_version >= 13)
    {
        do_cpt_int_err(xd, "GROMACS double precision", &contents->double_prec, list);
    }
    else
    {
        contents->double_prec = -1;
    }
    if (contents->file_version >= 12)
    {
        do_cpt_string_err(xd, "generating host", fhost, list);
    }
    do_cpt_int_err(xd, "#atoms", &contents->natoms, list);
    do_cpt_int_err(xd, "#T-coupling groups", &contents->ngtc, list);
    if (contents->file_version >= 10)
    {
        do_cpt_int_err(xd, "#Nose-Hoover T-chains", &contents->nhchainlength, list);
    }
    else
    {
        contents->nhchainlength = 1;
    }
    if (contents->file_version >= 11)
    {
        do_cpt_int_err(xd, "#Nose-Hoover T-chains for barostat ", &contents->nnhpres, list);
    }
    else
    {
        contents->nnhpres = 0;
    }
    if (contents->file_version >= 14)
    {
        do_cpt_int_err(xd, "# of total lambda states ", &contents->nlambda, list);
    }
    else
    {
        contents->nlambda = 0;
    }
    do_cpt_int_err(xd, "integrator", &contents->eIntegrator, list);
    if (contents->file_version >= 3)
    {
        do_cpt_int_err(xd, "simulation part #", &contents->simulation_part, list);
    }
    else
    {
        contents->simulation_part = 1;
    }
    if (contents->file_version >= 5)
    {
        do_cpt_step_err(xd, "step", &contents->step, list);
    }
    else
    {
        int idum = 0;
        do_cpt_int_err(xd, "step", &idum, list);
        contents->step = static_cast<int64_t>(idum);
    }
    do_cpt_double_err(xd, "t", &contents->t, list);
    do_cpt_int_err(xd, "#PP-ranks", &contents->nnodes, list);
    do_cpt_int_err(xd, "dd_nc[x]", &contents->dd_nc[XX], list);
    do_cpt_int_err(xd, "dd_nc[y]", &contents->dd_nc[YY], list);
    do_cpt_int_err(xd, "dd_nc[z]", &contents->dd_nc[ZZ], list);
    do_cpt_int_err(xd, "#PME-only ranks", &contents->npme, list);
    do_cpt_int_err(xd, "state flags", &contents->flags_state, list);
    if (contents->file_version >= 4)
    {
        do_cpt_int_err(xd, "ekin data flags", &contents->flags_eks, list);
        do_cpt_int_err(xd, "energy history flags", &contents->flags_enh, list);
    }
    else
    {
        contents->flags_eks   = 0;
        contents->flags_enh   = (contents->flags_state >> (estORIRE_DTAV + 1));
        contents->flags_state = (contents->flags_state
                                 & ~((1 << (estORIRE_DTAV + 1)) | (1 << (estORIRE_DTAV + 2))
                                     | (1 << (estORIRE_DTAV + 3))));
    }
    if (contents->file_version >= 14)
    {
        do_cpt_int_err(xd, "df history flags", &contents->flags_dfh, list);
    }
    else
    {
        contents->flags_dfh = 0;
    }

    if (contents->file_version >= 15)
    {
        do_cpt_int_err(xd, "ED data sets", &contents->nED, list);
    }
    else
    {
        contents->nED = 0;
    }

    if (contents->file_version >= 16)
    {
        do_cpt_int_err(xd, "swap", &contents->eSwapCoords, list);
    }
    else
    {
        contents->eSwapCoords = eswapNO;
    }

    if (contents->file_version >= 17)
    {
        do_cpt_int_err(xd, "AWH history flags", &contents->flags_awhh, list);
    }
    else
    {
        contents->flags_awhh = 0;
    }

    if (contents->file_version >= 18)
    {
        do_cpt_int_err(xd, "pull history flags", &contents->flagsPullHistory, list);
    }
    else
    {
        contents->flagsPullHistory = 0;
    }
}

static int do_cpt_footer(XDR* xd, int file_version)
{
    bool_t res = 0;
    int    magic;

    if (file_version >= 2)
    {
        magic = CPT_MAGIC2;
        res   = xdr_int(xd, &magic);
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

static int do_cpt_state(XDR* xd, int fflags, t_state* state, FILE* list)
{
    int             ret    = 0;
    const StatePart part   = StatePart::microState;
    const int       sflags = state->flags;
    // If reading, state->natoms was probably just read, so
    // allocations need to be managed. If writing, this won't change
    // anything that matters.
    state_change_natoms(state, state->natoms);
    for (int i = 0; (i < estNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {
                case estLAMBDA:
                    ret = doRealArrayRef(
                            xd, part, i, sflags,
                            gmx::arrayRefFromArray<real>(state->lambda.data(), state->lambda.size()),
                            list);
                    break;
                case estFEPSTATE:
                    ret = do_cpte_int(xd, part, i, sflags, &state->fep_state, list);
                    break;
                case estBOX: ret = do_cpte_matrix(xd, part, i, sflags, state->box, list); break;
                case estBOX_REL:
                    ret = do_cpte_matrix(xd, part, i, sflags, state->box_rel, list);
                    break;
                case estBOXV: ret = do_cpte_matrix(xd, part, i, sflags, state->boxv, list); break;
                case estPRES_PREV:
                    ret = do_cpte_matrix(xd, part, i, sflags, state->pres_prev, list);
                    break;
                case estSVIR_PREV:
                    ret = do_cpte_matrix(xd, part, i, sflags, state->svir_prev, list);
                    break;
                case estFVIR_PREV:
                    ret = do_cpte_matrix(xd, part, i, sflags, state->fvir_prev, list);
                    break;
                case estNH_XI:
                    ret = doVector<double>(xd, part, i, sflags, &state->nosehoover_xi, list);
                    break;
                case estNH_VXI:
                    ret = doVector<double>(xd, part, i, sflags, &state->nosehoover_vxi, list);
                    break;
                case estNHPRES_XI:
                    ret = doVector<double>(xd, part, i, sflags, &state->nhpres_xi, list);
                    break;
                case estNHPRES_VXI:
                    ret = doVector<double>(xd, part, i, sflags, &state->nhpres_vxi, list);
                    break;
                case estTHERM_INT:
                    ret = doVector<double>(xd, part, i, sflags, &state->therm_integral, list);
                    break;
                case estBAROS_INT:
                    ret = do_cpte_double(xd, part, i, sflags, &state->baros_integral, list);
                    break;
                case estVETA: ret = do_cpte_real(xd, part, i, sflags, &state->veta, list); break;
                case estVOL0: ret = do_cpte_real(xd, part, i, sflags, &state->vol0, list); break;
                case estX:
                    ret = doRvecVector(xd, part, i, sflags, &state->x, state->natoms, list);
                    break;
                case estV:
                    ret = doRvecVector(xd, part, i, sflags, &state->v, state->natoms, list);
                    break;
                /* The RNG entries are no longer written,
                 * the next 4 lines are only for reading old files.
                 */
                case estLD_RNG_NOTSUPPORTED:
                    ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list);
                    break;
                case estLD_RNGI_NOTSUPPORTED:
                    ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list);
                    break;
                case estMC_RNG_NOTSUPPORTED:
                    ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list);
                    break;
                case estMC_RNGI_NOTSUPPORTED:
                    ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list);
                    break;
                case estDISRE_INITF:
                    ret = do_cpte_real(xd, part, i, sflags, &state->hist.disre_initf, list);
                    break;
                case estDISRE_RM3TAV:
                    ret = do_cpte_n_reals(xd, part, i, sflags, &state->hist.ndisrepairs,
                                          &state->hist.disre_rm3tav, list);
                    break;
                case estORIRE_INITF:
                    ret = do_cpte_real(xd, part, i, sflags, &state->hist.orire_initf, list);
                    break;
                case estORIRE_DTAV:
                    ret = do_cpte_n_reals(xd, part, i, sflags, &state->hist.norire_Dtav,
                                          &state->hist.orire_Dtav, list);
                    break;
                case estPULLCOMPREVSTEP:
                    ret = doVector<double>(xd, part, i, sflags, &state->pull_com_prev_step, list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown state entry %d\n"
                              "You are reading a checkpoint file written by different code, which "
                              "is not supported",
                              i);
            }
        }
    }
    return ret;
}

static int do_cpt_ekinstate(XDR* xd, int fflags, ekinstate_t* ekins, FILE* list)
{
    int ret = 0;

    const StatePart part = StatePart::kineticEnergy;
    for (int i = 0; (i < eeksNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {

                case eeksEKIN_N:
                    ret = do_cpte_int(xd, part, i, fflags, &ekins->ekin_n, list);
                    break;
                case eeksEKINH:
                    ret = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinh, list);
                    break;
                case eeksEKINF:
                    ret = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinf, list);
                    break;
                case eeksEKINO:
                    ret = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinh_old, list);
                    break;
                case eeksEKINTOTAL:
                    ret = do_cpte_matrix(xd, part, i, fflags, ekins->ekin_total, list);
                    break;
                case eeksEKINSCALEF:
                    ret = doVector<double>(xd, part, i, fflags, &ekins->ekinscalef_nhc, list);
                    break;
                case eeksVSCALE:
                    ret = doVector<double>(xd, part, i, fflags, &ekins->vscale_nhc, list);
                    break;
                case eeksEKINSCALEH:
                    ret = doVector<double>(xd, part, i, fflags, &ekins->ekinscaleh_nhc, list);
                    break;
                case eeksDEKINDL:
                    ret = do_cpte_real(xd, part, i, fflags, &ekins->dekindl, list);
                    break;
                case eeksMVCOS: ret = do_cpte_real(xd, part, i, fflags, &ekins->mvcos, list); break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown ekin data state entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              i);
            }
        }
    }

    return ret;
}


static int do_cpt_swapstate(XDR* xd, gmx_bool bRead, int eSwapCoords, swaphistory_t* swapstate, FILE* list)
{
    int swap_cpt_version = 2;

    if (eSwapCoords == eswapNO)
    {
        return 0;
    }

    swapstate->bFromCpt    = bRead;
    swapstate->eSwapCoords = eSwapCoords;

    do_cpt_int_err(xd, "swap checkpoint version", &swap_cpt_version, list);
    if (bRead && swap_cpt_version < 2)
    {
        gmx_fatal(FARGS,
                  "Cannot read checkpoint files that were written with old versions"
                  "of the ion/water position swapping protocol.\n");
    }

    do_cpt_int_err(xd, "swap coupling steps", &swapstate->nAverage, list);

    /* When reading, init_swapcoords has not been called yet,
     * so we have to allocate memory first. */
    do_cpt_int_err(xd, "number of ion types", &swapstate->nIonTypes, list);
    if (bRead)
    {
        snew(swapstate->ionType, swapstate->nIonTypes);
    }

    for (int ic = 0; ic < eCompNR; ic++)
    {
        for (int ii = 0; ii < swapstate->nIonTypes; ii++)
        {
            swapstateIons_t* gs = &swapstate->ionType[ii];

            if (bRead)
            {
                do_cpt_int_err(xd, "swap requested atoms", &gs->nMolReq[ic], list);
            }
            else
            {
                do_cpt_int_err(xd, "swap requested atoms p", gs->nMolReq_p[ic], list);
            }

            if (bRead)
            {
                do_cpt_int_err(xd, "swap influx net", &gs->inflow_net[ic], list);
            }
            else
            {
                do_cpt_int_err(xd, "swap influx net p", gs->inflow_net_p[ic], list);
            }

            if (bRead && (nullptr == gs->nMolPast[ic]))
            {
                snew(gs->nMolPast[ic], swapstate->nAverage);
            }

            for (int j = 0; j < swapstate->nAverage; j++)
            {
                if (bRead)
                {
                    do_cpt_int_err(xd, "swap past atom counts", &gs->nMolPast[ic][j], list);
                }
                else
                {
                    do_cpt_int_err(xd, "swap past atom counts p", &gs->nMolPast_p[ic][j], list);
                }
            }
        }
    }

    /* Ion flux per channel */
    for (int ic = 0; ic < eChanNR; ic++)
    {
        for (int ii = 0; ii < swapstate->nIonTypes; ii++)
        {
            swapstateIons_t* gs = &swapstate->ionType[ii];

            if (bRead)
            {
                do_cpt_int_err(xd, "channel flux A->B", &gs->fluxfromAtoB[ic], list);
            }
            else
            {
                do_cpt_int_err(xd, "channel flux A->B p", gs->fluxfromAtoB_p[ic], list);
            }
        }
    }

    /* Ion flux leakage */
    if (bRead)
    {
        do_cpt_int_err(xd, "flux leakage", &swapstate->fluxleak, list);
    }
    else
    {
        do_cpt_int_err(xd, "flux leakage", swapstate->fluxleak_p, list);
    }

    /* Ion history */
    for (int ii = 0; ii < swapstate->nIonTypes; ii++)
    {
        swapstateIons_t* gs = &swapstate->ionType[ii];

        do_cpt_int_err(xd, "number of ions", &gs->nMol, list);

        if (bRead)
        {
            snew(gs->channel_label, gs->nMol);
            snew(gs->comp_from, gs->nMol);
        }

        do_cpt_u_chars(xd, "channel history", gs->nMol, gs->channel_label, list);
        do_cpt_u_chars(xd, "domain history", gs->nMol, gs->comp_from, list);
    }

    /* Save the last known whole positions to checkpoint
     * file to be able to also make multimeric channels whole in PBC */
    do_cpt_int_err(xd, "Ch0 atoms", &swapstate->nat[eChan0], list);
    do_cpt_int_err(xd, "Ch1 atoms", &swapstate->nat[eChan1], list);
    if (bRead)
    {
        snew(swapstate->xc_old_whole[eChan0], swapstate->nat[eChan0]);
        snew(swapstate->xc_old_whole[eChan1], swapstate->nat[eChan1]);
        do_cpt_n_rvecs_err(xd, "Ch0 whole x", swapstate->nat[eChan0], swapstate->xc_old_whole[eChan0], list);
        do_cpt_n_rvecs_err(xd, "Ch1 whole x", swapstate->nat[eChan1], swapstate->xc_old_whole[eChan1], list);
    }
    else
    {
        do_cpt_n_rvecs_err(xd, "Ch0 whole x", swapstate->nat[eChan0],
                           *swapstate->xc_old_whole_p[eChan0], list);
        do_cpt_n_rvecs_err(xd, "Ch1 whole x", swapstate->nat[eChan1],
                           *swapstate->xc_old_whole_p[eChan1], list);
    }

    return 0;
}


static int do_cpt_enerhist(XDR* xd, gmx_bool bRead, int fflags, energyhistory_t* enerhist, FILE* list)
{
    int ret = 0;

    if (fflags == 0)
    {
        return ret;
    }

    GMX_RELEASE_ASSERT(enerhist != nullptr,
                       "With energy history, we need a valid enerhist pointer");

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
    const StatePart    part   = StatePart::energyHistory;
    for (int i = 0; (i < eenhNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {
                case eenhENERGY_N:
                    ret = do_cpte_int(xd, part, i, fflags, &energyHistoryNumEnergies, list);
                    break;
                case eenhENERGY_AVER:
                    ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_ave, list);
                    break;
                case eenhENERGY_SUM:
                    ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_sum, list);
                    break;
                case eenhENERGY_NSUM:
                    do_cpt_step_err(xd, eenh_names[i], &enerhist->nsum, list);
                    break;
                case eenhENERGY_SUM_SIM:
                    ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_sum_sim, list);
                    break;
                case eenhENERGY_NSUM_SIM:
                    do_cpt_step_err(xd, eenh_names[i], &enerhist->nsum_sim, list);
                    break;
                case eenhENERGY_NSTEPS:
                    do_cpt_step_err(xd, eenh_names[i], &enerhist->nsteps, list);
                    break;
                case eenhENERGY_NSTEPS_SIM:
                    do_cpt_step_err(xd, eenh_names[i], &enerhist->nsteps_sim, list);
                    break;
                case eenhENERGY_DELTA_H_NN:
                {
                    int numDeltaH = 0;
                    if (!bRead && deltaH != nullptr)
                    {
                        numDeltaH = deltaH->dh.size();
                    }
                    do_cpt_int_err(xd, eenh_names[i], &numDeltaH, list);
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
                case eenhENERGY_DELTA_H_LIST:
                    for (auto dh : deltaH->dh)
                    {
                        ret = doVector<real>(xd, part, i, fflags, &dh, list);
                    }
                    break;
                case eenhENERGY_DELTA_H_STARTTIME:
                    ret = do_cpte_double(xd, part, i, fflags, &(deltaH->start_time), list);
                    break;
                case eenhENERGY_DELTA_H_STARTLAMBDA:
                    ret = do_cpte_double(xd, part, i, fflags, &(deltaH->start_lambda), list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown energy history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              i);
            }
        }
    }

    if ((fflags & (1 << eenhENERGY_SUM)) && !(fflags & (1 << eenhENERGY_SUM_SIM)))
    {
        /* Assume we have an old file format and copy sum to sum_sim */
        enerhist->ener_sum_sim = enerhist->ener_sum;
    }

    if ((fflags & (1 << eenhENERGY_NSUM)) && !(fflags & (1 << eenhENERGY_NSTEPS)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps = enerhist->nsum;
    }
    if ((fflags & (1 << eenhENERGY_NSUM_SIM)) && !(fflags & (1 << eenhENERGY_NSTEPS_SIM)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps_sim = enerhist->nsum_sim;
    }

    return ret;
}

static int doCptPullCoordHist(XDR* xd, PullCoordinateHistory* pullCoordHist, const StatePart part, FILE* list)
{
    int ret   = 0;
    int flags = 0;

    flags |= ((1 << epullcoordh_VALUE_REF_SUM) | (1 << epullcoordh_VALUE_SUM)
              | (1 << epullcoordh_DR01_SUM) | (1 << epullcoordh_DR23_SUM)
              | (1 << epullcoordh_DR45_SUM) | (1 << epullcoordh_FSCAL_SUM));

    for (int i = 0; i < epullcoordh_NR && ret == 0; i++)
    {
        switch (i)
        {
            case epullcoordh_VALUE_REF_SUM:
                ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->valueRef), list);
                break;
            case epullcoordh_VALUE_SUM:
                ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->value), list);
                break;
            case epullcoordh_DR01_SUM:
                for (int j = 0; j < DIM && ret == 0; j++)
                {
                    ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->dr01[j]), list);
                }
                break;
            case epullcoordh_DR23_SUM:
                for (int j = 0; j < DIM && ret == 0; j++)
                {
                    ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->dr23[j]), list);
                }
                break;
            case epullcoordh_DR45_SUM:
                for (int j = 0; j < DIM && ret == 0; j++)
                {
                    ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->dr45[j]), list);
                }
                break;
            case epullcoordh_FSCAL_SUM:
                ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->scalarForce), list);
                break;
            case epullcoordh_DYNAX_SUM:
                for (int j = 0; j < DIM && ret == 0; j++)
                {
                    ret = do_cpte_double(xd, part, i, flags, &(pullCoordHist->dynaX[j]), list);
                }
                break;
        }
    }

    return ret;
}

static int doCptPullGroupHist(XDR* xd, PullGroupHistory* pullGroupHist, const StatePart part, FILE* list)
{
    int ret   = 0;
    int flags = 0;

    flags |= ((1 << epullgrouph_X_SUM));

    for (int i = 0; i < epullgrouph_NR; i++)
    {
        switch (i)
        {
            case epullgrouph_X_SUM:
                for (int j = 0; j < DIM && ret == 0; j++)
                {
                    ret = do_cpte_double(xd, part, i, flags, &(pullGroupHist->x[j]), list);
                }
                break;
        }
    }

    return ret;
}


static int doCptPullHist(XDR* xd, gmx_bool bRead, int fflags, PullHistory* pullHist, const StatePart part, FILE* list)
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

    for (int i = 0; (i < epullhNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {
                case epullhPULL_NUMCOORDINATES:
                    ret = do_cpte_int(xd, part, i, fflags, &pullHistoryNumCoordinates, list);
                    break;
                case epullhPULL_NUMGROUPS:
                    do_cpt_int_err(xd, eenh_names[i], &pullHistoryNumGroups, list);
                    break;
                case epullhPULL_NUMVALUESINXSUM:
                    do_cpt_int_err(xd, eenh_names[i], &pullHist->numValuesInXSum, list);
                    break;
                case epullhPULL_NUMVALUESINFSUM:
                    do_cpt_int_err(xd, eenh_names[i], &pullHist->numValuesInFSum, list);
                    break;
                default:
                    gmx_fatal(FARGS,
                              "Unknown pull history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              i);
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
            ret = doCptPullCoordHist(xd, &(pullHist->pullCoordinateSums[i]), part, list);
        }
        for (size_t i = 0; i < pullHist->pullGroupSums.size() && ret == 0; i++)
        {
            ret = doCptPullGroupHist(xd, &(pullHist->pullGroupSums[i]), part, list);
        }
    }

    return ret;
}

static int do_cpt_df_hist(XDR* xd, int fflags, int nlambda, df_history_t** dfhistPtr, FILE* list)
{
    int ret = 0;

    if (fflags == 0)
    {
        return 0;
    }

    if (*dfhistPtr == nullptr)
    {
        snew(*dfhistPtr, 1);
        (*dfhistPtr)->nlambda = nlambda;
        init_df_history(*dfhistPtr, nlambda);
    }
    df_history_t* dfhist = *dfhistPtr;

    const StatePart part = StatePart::freeEnergyHistory;
    for (int i = 0; (i < edfhNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {
                case edfhBEQUIL:
                    ret = do_cpte_bool(xd, part, i, fflags, &dfhist->bEquil, list);
                    break;
                case edfhNATLAMBDA:
                    ret = do_cpte_ints(xd, part, i, fflags, nlambda, &dfhist->n_at_lam, list);
                    break;
                case edfhWLHISTO:
                    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->wl_histo, list);
                    break;
                case edfhWLDELTA:
                    ret = do_cpte_real(xd, part, i, fflags, &dfhist->wl_delta, list);
                    break;
                case edfhSUMWEIGHTS:
                    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_weights, list);
                    break;
                case edfhSUMDG:
                    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_dg, list);
                    break;
                case edfhSUMMINVAR:
                    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_minvar, list);
                    break;
                case edfhSUMVAR:
                    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_variance, list);
                    break;
                case edfhACCUMP:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_p, list);
                    break;
                case edfhACCUMM:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_m, list);
                    break;
                case edfhACCUMP2:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_p2, list);
                    break;
                case edfhACCUMM2:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_m2, list);
                    break;
                case edfhTIJ:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->Tij, list);
                    break;
                case edfhTIJEMP:
                    ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->Tij_empirical, list);
                    break;

                default:
                    gmx_fatal(FARGS,
                              "Unknown df history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",
                              i);
            }
        }
    }

    return ret;
}


/* This function stores the last whole configuration of the reference and
 * average structure in the .cpt file
 */
static int do_cpt_EDstate(XDR* xd, gmx_bool bRead, int nED, edsamhistory_t* EDstate, FILE* list)
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
        do_cpt_int_err(xd, buf, &EDstate->nref[i], list);
        sprintf(buf, "ED%d x_ref", i + 1);
        if (bRead)
        {
            snew(EDstate->old_sref[i], EDstate->nref[i]);
            do_cpt_n_rvecs_err(xd, buf, EDstate->nref[i], EDstate->old_sref[i], list);
        }
        else
        {
            do_cpt_n_rvecs_err(xd, buf, EDstate->nref[i], EDstate->old_sref_p[i], list);
        }

        /* Average structure SAV */
        sprintf(buf, "ED%d # of atoms in average structure", i + 1);
        do_cpt_int_err(xd, buf, &EDstate->nav[i], list);
        sprintf(buf, "ED%d x_av", i + 1);
        if (bRead)
        {
            snew(EDstate->old_sav[i], EDstate->nav[i]);
            do_cpt_n_rvecs_err(xd, buf, EDstate->nav[i], EDstate->old_sav[i], list);
        }
        else
        {
            do_cpt_n_rvecs_err(xd, buf, EDstate->nav[i], EDstate->old_sav_p[i], list);
        }
    }

    return 0;
}

static int do_cpt_correlation_grid(XDR*                         xd,
                                   gmx_bool                     bRead,
                                   gmx_unused int               fflags,
                                   gmx::CorrelationGridHistory* corrGrid,
                                   FILE*                        list,
                                   int                          eawhh)
{
    int ret = 0;

    do_cpt_int_err(xd, eawhh_names[eawhh], &(corrGrid->numCorrelationTensors), list);
    do_cpt_int_err(xd, eawhh_names[eawhh], &(corrGrid->tensorSize), list);
    do_cpt_int_err(xd, eawhh_names[eawhh], &(corrGrid->blockDataListSize), list);

    if (bRead)
    {
        initCorrelationGridHistory(corrGrid, corrGrid->numCorrelationTensors, corrGrid->tensorSize,
                                   corrGrid->blockDataListSize);
    }

    for (gmx::CorrelationBlockDataHistory& blockData : corrGrid->blockDataBuffer)
    {
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.blockSumWeight), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.blockSumSquareWeight), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.blockSumWeightX), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.blockSumWeightY), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.sumOverBlocksSquareBlockWeight), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.sumOverBlocksBlockSquareWeight), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.sumOverBlocksBlockWeightBlockWeightX), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.sumOverBlocksBlockWeightBlockWeightY), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.blockLength), list);
        do_cpt_int_err(xd, eawhh_names[eawhh], &(blockData.previousBlockIndex), list);
        do_cpt_double_err(xd, eawhh_names[eawhh], &(blockData.correlationIntegral), list);
    }

    return ret;
}

static int do_cpt_awh_bias(XDR* xd, gmx_bool bRead, int fflags, gmx::AwhBiasHistory* biasHistory, FILE* list)
{
    int ret = 0;

    gmx::AwhBiasStateHistory* state = &biasHistory->state;
    for (int i = 0; (i < eawhhNR && ret == 0); i++)
    {
        if (fflags & (1 << i))
        {
            switch (i)
            {
                case eawhhIN_INITIAL:
                    do_cpt_bool_err(xd, eawhh_names[i], &state->in_initial, list);
                    break;
                case eawhhEQUILIBRATEHISTOGRAM:
                    do_cpt_bool_err(xd, eawhh_names[i], &state->equilibrateHistogram, list);
                    break;
                case eawhhHISTSIZE:
                    do_cpt_double_err(xd, eawhh_names[i], &state->histSize, list);
                    break;
                case eawhhNPOINTS:
                {
                    int numPoints;
                    if (!bRead)
                    {
                        numPoints = biasHistory->pointState.size();
                    }
                    do_cpt_int_err(xd, eawhh_names[i], &numPoints, list);
                    if (bRead)
                    {
                        biasHistory->pointState.resize(numPoints);
                    }
                }
                break;
                case eawhhCOORDPOINT:
                    for (auto& psh : biasHistory->pointState)
                    {
                        do_cpt_double_err(xd, eawhh_names[i], &psh.target, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.free_energy, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.bias, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.weightsum_iteration, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.weightsum_covering, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.weightsum_tot, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.weightsum_ref, list);
                        do_cpt_step_err(xd, eawhh_names[i], &psh.last_update_index, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.log_pmfsum, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.visits_iteration, list);
                        do_cpt_double_err(xd, eawhh_names[i], &psh.visits_tot, list);
                    }
                    break;
                case eawhhUMBRELLAGRIDPOINT:
                    do_cpt_int_err(xd, eawhh_names[i], &(state->umbrellaGridpoint), list);
                    break;
                case eawhhUPDATELIST:
                    do_cpt_int_err(xd, eawhh_names[i], &(state->origin_index_updatelist), list);
                    do_cpt_int_err(xd, eawhh_names[i], &(state->end_index_updatelist), list);
                    break;
                case eawhhLOGSCALEDSAMPLEWEIGHT:
                    do_cpt_double_err(xd, eawhh_names[i], &(state->logScaledSampleWeight), list);
                    do_cpt_double_err(xd, eawhh_names[i], &(state->maxLogScaledSampleWeight), list);
                    break;
                case eawhhNUMUPDATES:
                    do_cpt_step_err(xd, eawhh_names[i], &(state->numUpdates), list);
                    break;
                case eawhhFORCECORRELATIONGRID:
                    ret = do_cpt_correlation_grid(xd, bRead, fflags,
                                                  &biasHistory->forceCorrelationGrid, list, i);
                    break;
                default: gmx_fatal(FARGS, "Unknown awh history entry %d\n", i);
            }
        }
    }

    return ret;
}

static int do_cpt_awh(XDR* xd, gmx_bool bRead, int fflags, gmx::AwhHistory* awhHistory, FILE* list)
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
        do_cpt_int_err(xd, "awh_nbias", &numBias, list);

        if (bRead)
        {
            awhHistory->bias.resize(numBias);
        }
        for (auto& bias : awhHistory->bias)
        {
            ret = do_cpt_awh_bias(xd, bRead, fflags, &bias, list);
            if (ret)
            {
                return ret;
            }
        }
        do_cpt_double_err(xd, "awh_potential_offset", &awhHistory->potentialOffset, list);
    }
    return ret;
}

static void do_cpt_mdmodules(int                           fileVersion,
                             t_fileio*                     checkpointFileHandle,
                             const gmx::MdModulesNotifier& mdModulesNotifier)
{
    if (fileVersion >= cptv_MdModules)
    {
        gmx::FileIOXdrSerializer serializer(checkpointFileHandle);
        gmx::KeyValueTreeObject  mdModuleCheckpointParameterTree =
                gmx::deserializeKeyValueTree(&serializer);
        gmx::MdModulesCheckpointReadingDataOnMaster mdModuleCheckpointReadingDataOnMaster = {
            mdModuleCheckpointParameterTree, fileVersion
        };
        mdModulesNotifier.notifier_.notify(mdModuleCheckpointReadingDataOnMaster);
    }
}

static int do_cpt_files(XDR* xd, gmx_bool bRead, std::vector<gmx_file_position_t>* outputfiles, FILE* list, int file_version)
{
    gmx_off_t                   offset;
    gmx_off_t                   mask = 0xFFFFFFFFL;
    int                         offset_high, offset_low;
    std::array<char, CPTSTRLEN> buf;
    GMX_RELEASE_ASSERT(outputfiles, "Must have valid outputfiles");

    // Ensure that reading pre-allocates outputfiles, while writing
    // writes what is already there.
    int nfiles = outputfiles->size();
    if (do_cpt_int(xd, "number of output files", &nfiles, list) != 0)
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
            do_cpt_string_err(xd, "output filename", buf, list);
            std::copy(std::begin(buf), std::end(buf), std::begin(outputfile.filename));

            if (do_cpt_int(xd, "file_offset_high", &offset_high, list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(xd, "file_offset_low", &offset_low, list) != 0)
            {
                return -1;
            }
            outputfile.offset = (static_cast<gmx_off_t>(offset_high) << 32)
                                | (static_cast<gmx_off_t>(offset_low) & mask);
        }
        else
        {
            do_cpt_string_err(xd, "output filename", outputfile.filename, list);
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
            if (do_cpt_int(xd, "file_offset_high", &offset_high, list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(xd, "file_offset_low", &offset_low, list) != 0)
            {
                return -1;
            }
        }
        if (file_version >= 8)
        {
            if (do_cpt_int(xd, "file_checksum_size", &outputfile.checksumSize, list) != 0)
            {
                return -1;
            }
            if (do_cpt_u_chars(xd, "file_checksum", outputfile.checksum.size(),
                               outputfile.checksum.data(), list)
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

static void mpiBarrierBeforeRename(const bool applyMpiBarrierBeforeRename, MPI_Comm mpiBarrierCommunicator)
{
    if (applyMpiBarrierBeforeRename)
    {
#if GMX_MPI
        MPI_Barrier(mpiBarrierCommunicator);
#else
        GMX_RELEASE_ASSERT(false, "Should not request a barrier without MPI");
        GMX_UNUSED_VALUE(mpiBarrierCommunicator);
#endif
    }
}

void write_checkpoint(const char*                   fn,
                      gmx_bool                      bNumberAndKeep,
                      FILE*                         fplog,
                      const t_commrec*              cr,
                      ivec                          domdecCells,
                      int                           nppnodes,
                      int                           eIntegrator,
                      int                           simulation_part,
                      gmx_bool                      bExpanded,
                      int                           elamstats,
                      int64_t                       step,
                      double                        t,
                      t_state*                      state,
                      ObservablesHistory*           observablesHistory,
                      const gmx::MdModulesNotifier& mdModulesNotifier,
                      bool                          applyMpiBarrierBeforeRename,
                      MPI_Comm                      mpiBarrierCommunicator)
{
    t_fileio* fp;
    char*     fntemp; /* the temporary checkpoint file name */
    int       npmenodes;
    char      buf[1024], suffix[5 + STEPSTRSIZE], sbuf[STEPSTRSIZE];
    t_fileio* ret;

    if (DOMAINDECOMP(cr))
    {
        npmenodes = cr->npmenodes;
    }
    else
    {
        npmenodes = 0;
    }

#if !GMX_NO_RENAME
    /* make the new temporary filename */
    snew(fntemp, std::strlen(fn) + 5 + STEPSTRSIZE);
    std::strcpy(fntemp, fn);
    fntemp[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
    sprintf(suffix, "_%s%s", "step", gmx_step_str(step, sbuf));
    std::strcat(fntemp, suffix);
    std::strcat(fntemp, fn + std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
#else
    /* if we can't rename, we just overwrite the cpt file.
     * dangerous if interrupted.
     */
    snew(fntemp, std::strlen(fn));
    std::strcpy(fntemp, fn);
#endif
    std::string timebuf = gmx_format_current_time();

    if (fplog)
    {
        fprintf(fplog, "Writing checkpoint, step %s at %s\n\n", gmx_step_str(step, buf), timebuf.c_str());
    }

    /* Get offsets for open files */
    auto outputfiles = gmx_fio_get_output_file_positions();

    fp = gmx_fio_open(fntemp, "w");

    int flags_eks;
    if (state->ekinstate.bUpToDate)
    {
        flags_eks = ((1 << eeksEKIN_N) | (1 << eeksEKINH) | (1 << eeksEKINF) | (1 << eeksEKINO)
                     | (1 << eeksEKINSCALEF) | (1 << eeksEKINSCALEH) | (1 << eeksVSCALE)
                     | (1 << eeksDEKINDL) | (1 << eeksMVCOS));
    }
    else
    {
        flags_eks = 0;
    }

    energyhistory_t* enerhist  = observablesHistory->energyHistory.get();
    int              flags_enh = 0;
    if (enerhist != nullptr && (enerhist->nsum > 0 || enerhist->nsum_sim > 0))
    {
        flags_enh |= (1 << eenhENERGY_N) | (1 << eenhENERGY_NSTEPS) | (1 << eenhENERGY_NSTEPS_SIM);
        if (enerhist->nsum > 0)
        {
            flags_enh |= ((1 << eenhENERGY_AVER) | (1 << eenhENERGY_SUM) | (1 << eenhENERGY_NSUM));
        }
        if (enerhist->nsum_sim > 0)
        {
            flags_enh |= ((1 << eenhENERGY_SUM_SIM) | (1 << eenhENERGY_NSUM_SIM));
        }
        if (enerhist->deltaHForeignLambdas != nullptr)
        {
            flags_enh |= ((1 << eenhENERGY_DELTA_H_NN) | (1 << eenhENERGY_DELTA_H_LIST)
                          | (1 << eenhENERGY_DELTA_H_STARTTIME) | (1 << eenhENERGY_DELTA_H_STARTLAMBDA));
        }
    }

    PullHistory* pullHist         = observablesHistory->pullHistory.get();
    int          flagsPullHistory = 0;
    if (pullHist != nullptr && (pullHist->numValuesInXSum > 0 || pullHist->numValuesInFSum > 0))
    {
        flagsPullHistory |= (1 << epullhPULL_NUMCOORDINATES);
        flagsPullHistory |= ((1 << epullhPULL_NUMGROUPS) | (1 << epullhPULL_NUMVALUESINXSUM)
                             | (1 << epullhPULL_NUMVALUESINFSUM));
    }

    int flags_dfh;
    if (bExpanded)
    {
        flags_dfh = ((1 << edfhBEQUIL) | (1 << edfhNATLAMBDA) | (1 << edfhSUMWEIGHTS)
                     | (1 << edfhSUMDG) | (1 << edfhTIJ) | (1 << edfhTIJEMP));
        if (EWL(elamstats))
        {
            flags_dfh |= ((1 << edfhWLDELTA) | (1 << edfhWLHISTO));
        }
        if ((elamstats == elamstatsMINVAR) || (elamstats == elamstatsBARKER)
            || (elamstats == elamstatsMETROPOLIS))
        {
            flags_dfh |= ((1 << edfhACCUMP) | (1 << edfhACCUMM) | (1 << edfhACCUMP2)
                          | (1 << edfhACCUMM2) | (1 << edfhSUMMINVAR) | (1 << edfhSUMVAR));
        }
    }
    else
    {
        flags_dfh = 0;
    }

    int flags_awhh = 0;
    if (state->awhHistory != nullptr && !state->awhHistory->bias.empty())
    {
        flags_awhh |= ((1 << eawhhIN_INITIAL) | (1 << eawhhEQUILIBRATEHISTOGRAM) | (1 << eawhhHISTSIZE)
                       | (1 << eawhhNPOINTS) | (1 << eawhhCOORDPOINT) | (1 << eawhhUMBRELLAGRIDPOINT)
                       | (1 << eawhhUPDATELIST) | (1 << eawhhLOGSCALEDSAMPLEWEIGHT)
                       | (1 << eawhhNUMUPDATES) | (1 << eawhhFORCECORRELATIONGRID));
    }

    /* We can check many more things now (CPU, acceleration, etc), but
     * it is highly unlikely to have two separate builds with exactly
     * the same version, user, time, and build host!
     */

    int nlambda = (state->dfhist ? state->dfhist->nlambda : 0);

    edsamhistory_t* edsamhist = observablesHistory->edsamHistory.get();
    int             nED       = (edsamhist ? edsamhist->nED : 0);

    swaphistory_t* swaphist    = observablesHistory->swapHistory.get();
    int            eSwapCoords = (swaphist ? swaphist->eSwapCoords : eswapNO);

    CheckpointHeaderContents headerContents = { 0,
                                                { 0 },
                                                { 0 },
                                                { 0 },
                                                { 0 },
                                                GMX_DOUBLE,
                                                { 0 },
                                                { 0 },
                                                eIntegrator,
                                                simulation_part,
                                                step,
                                                t,
                                                nppnodes,
                                                { 0 },
                                                npmenodes,
                                                state->natoms,
                                                state->ngtc,
                                                state->nnhpres,
                                                state->nhchainlength,
                                                nlambda,
                                                state->flags,
                                                flags_eks,
                                                flags_enh,
                                                flagsPullHistory,
                                                flags_dfh,
                                                flags_awhh,
                                                nED,
                                                eSwapCoords };
    std::strcpy(headerContents.version, gmx_version());
    std::strcpy(headerContents.fprog, gmx::getProgramContext().fullBinaryPath());
    std::strcpy(headerContents.ftime, timebuf.c_str());
    if (DOMAINDECOMP(cr))
    {
        copy_ivec(domdecCells, headerContents.dd_nc);
    }

    do_cpt_header(gmx_fio_getxdr(fp), FALSE, nullptr, &headerContents);

    if ((do_cpt_state(gmx_fio_getxdr(fp), state->flags, state, nullptr) < 0)
        || (do_cpt_ekinstate(gmx_fio_getxdr(fp), flags_eks, &state->ekinstate, nullptr) < 0)
        || (do_cpt_enerhist(gmx_fio_getxdr(fp), FALSE, flags_enh, enerhist, nullptr) < 0)
        || (doCptPullHist(gmx_fio_getxdr(fp), FALSE, flagsPullHistory, pullHist, StatePart::pullHistory, nullptr)
            < 0)
        || (do_cpt_df_hist(gmx_fio_getxdr(fp), flags_dfh, nlambda, &state->dfhist, nullptr) < 0)
        || (do_cpt_EDstate(gmx_fio_getxdr(fp), FALSE, nED, edsamhist, nullptr) < 0)
        || (do_cpt_awh(gmx_fio_getxdr(fp), FALSE, flags_awhh, state->awhHistory.get(), nullptr) < 0)
        || (do_cpt_swapstate(gmx_fio_getxdr(fp), FALSE, eSwapCoords, swaphist, nullptr) < 0)
        || (do_cpt_files(gmx_fio_getxdr(fp), FALSE, &outputfiles, nullptr, headerContents.file_version) < 0))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    // Checkpointing MdModules
    {
        gmx::KeyValueTreeBuilder          builder;
        gmx::MdModulesWriteCheckpointData mdModulesWriteCheckpoint = { builder.rootObject(),
                                                                       headerContents.file_version };
        mdModulesNotifier.notifier_.notify(mdModulesWriteCheckpoint);
        auto                     tree = builder.build();
        gmx::FileIOXdrSerializer serializer(fp);
        gmx::serializeKeyValueTree(tree, &serializer);
    }

    do_cpt_footer(gmx_fio_getxdr(fp), headerContents.file_version);

    /* we really, REALLY, want to make sure to physically write the checkpoint,
       and all the files it depends on, out to disk. Because we've
       opened the checkpoint with gmx_fio_open(), it's in our list
       of open files.  */
    ret = gmx_fio_all_output_fsync();

    if (ret)
    {
        char buf[STRLEN];
        sprintf(buf, "Cannot fsync '%s'; maybe you are out of disk space?", gmx_fio_getname(ret));

        if (getenv(GMX_IGNORE_FSYNC_FAILURE_ENV) == nullptr)
        {
            gmx_file(buf);
        }
        else
        {
            gmx_warning("%s", buf);
        }
    }

    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    /* we don't move the checkpoint if the user specified they didn't want it,
       or if the fsyncs failed */
#if !GMX_NO_RENAME
    if (!bNumberAndKeep && !ret)
    {
        if (gmx_fexist(fn))
        {
            /* Rename the previous checkpoint file */
            mpiBarrierBeforeRename(applyMpiBarrierBeforeRename, mpiBarrierCommunicator);

            std::strcpy(buf, fn);
            buf[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
            std::strcat(buf, "_prev");
            std::strcat(buf, fn + std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
            if (!GMX_FAHCORE)
            {
                /* we copy here so that if something goes wrong between now and
                 * the rename below, there's always a state.cpt.
                 * If renames are atomic (such as in POSIX systems),
                 * this copying should be unneccesary.
                 */
                gmx_file_copy(fn, buf, FALSE);
                /* We don't really care if this fails:
                 * there's already a new checkpoint.
                 */
            }
            else
            {
                gmx_file_rename(fn, buf);
            }
        }

        /* Rename the checkpoint file from the temporary to the final name */
        mpiBarrierBeforeRename(applyMpiBarrierBeforeRename, mpiBarrierCommunicator);

        if (gmx_file_rename(fntemp, fn) != 0)
        {
            gmx_file("Cannot rename checkpoint file; maybe you are out of disk space?");
        }
    }
#endif /* GMX_NO_RENAME */

    sfree(fntemp);

#if GMX_FAHCORE
    /*code for alternate checkpointing scheme.  moved from top of loop over
       steps */
    fcRequestCheckPoint();
    if (fcCheckPointParallel(cr->nodeid, NULL, 0) == 0)
    {
        gmx_fatal(3, __FILE__, __LINE__, "Checkpoint error on step %d\n", step);
    }
#endif /* end GMX_FAHCORE block */
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
                "single/double precision will lead to precision or performance loss.\n";
        if (fplog)
        {
            fprintf(fplog, "%s\n", msg_precision_difference);
        }
    }

    gmx_bool mm = (versionDiffers || precisionDiffers);

    if (reproducibilityRequested)
    {
        check_string(fplog, "Program name", gmx::getProgramContext().fullBinaryPath(),
                     headerContents.fprog, &mm);

        check_int(fplog, "#ranks", cr->nnodes, headerContents.nnodes, &mm);
    }

    if (cr->nnodes > 1 && reproducibilityRequested)
    {
        check_int(fplog, "#PME-ranks", cr->npmenodes, headerContents.npme, &mm);

        int npp = cr->nnodes;
        if (cr->npmenodes >= 0)
        {
            npp -= cr->npmenodes;
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

static void read_checkpoint(const char*                   fn,
                            t_fileio*                     logfio,
                            const t_commrec*              cr,
                            const ivec                    dd_nc,
                            int                           eIntegrator,
                            int*                          init_fep_state,
                            CheckpointHeaderContents*     headerContents,
                            t_state*                      state,
                            ObservablesHistory*           observablesHistory,
                            gmx_bool                      reproducibilityRequested,
                            const gmx::MdModulesNotifier& mdModulesNotifier)
{
    t_fileio* fp;
    char      buf[STEPSTRSIZE];
    int       ret;

    fp = gmx_fio_open(fn, "r");
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, nullptr, headerContents);

    // If we are appending, then we don't want write to the open log
    // file because we still need to compute a checksum for it. In
    // that case, the filehandle will be nullptr. Otherwise, we report
    // to the new log file about the checkpoint file that we are
    // reading from.
    FILE* fplog = gmx_fio_getfp(logfio);
    if (fplog)
    {
        fprintf(fplog, "\n");
        fprintf(fplog, "Reading checkpoint file %s\n", fn);
        fprintf(fplog, "  file generated by:     %s\n", headerContents->fprog);
        fprintf(fplog, "  file generated at:     %s\n", headerContents->ftime);
        fprintf(fplog, "  GROMACS double prec.:  %d\n", headerContents->double_prec);
        fprintf(fplog, "  simulation part #:     %d\n", headerContents->simulation_part);
        fprintf(fplog, "  step:                  %s\n", gmx_step_str(headerContents->step, buf));
        fprintf(fplog, "  time:                  %f\n", headerContents->t);
        fprintf(fplog, "\n");
    }

    if (headerContents->natoms != state->natoms)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d atoms, while the current system consists "
                  "of %d atoms",
                  headerContents->natoms, state->natoms);
    }
    if (headerContents->ngtc != state->ngtc)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d T-coupling groups, while the current "
                  "system consists of %d T-coupling groups",
                  headerContents->ngtc, state->ngtc);
    }
    if (headerContents->nnhpres != state->nnhpres)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system of %d NH-pressure-coupling variables, while the "
                  "current system consists of %d NH-pressure-coupling variables",
                  headerContents->nnhpres, state->nnhpres);
    }

    int nlambdaHistory = (state->dfhist ? state->dfhist->nlambda : 0);
    if (headerContents->nlambda != nlambdaHistory)
    {
        gmx_fatal(FARGS,
                  "Checkpoint file is for a system with %d lambda states, while the current system "
                  "consists of %d lambda states",
                  headerContents->nlambda, nlambdaHistory);
    }

    init_gtc_state(state, state->ngtc, state->nnhpres,
                   headerContents->nhchainlength); /* need to keep this here to keep the tpr format working */
    /* write over whatever was read; we use the number of Nose-Hoover chains from the checkpoint */

    if (headerContents->eIntegrator != eIntegrator)
    {
        gmx_fatal(FARGS,
                  "Cannot change integrator during a checkpoint restart. Perhaps you should make a "
                  "new .tpr with grompp -f new.mdp -t %s",
                  fn);
    }

    if (headerContents->flags_state != state->flags)
    {
        gmx_fatal(FARGS,
                  "Cannot change a simulation algorithm during a checkpoint restart. Perhaps you "
                  "should make a new .tpr with grompp -f new.mdp -t %s",
                  fn);
    }

    if (MASTER(cr))
    {
        check_match(fplog, cr, dd_nc, *headerContents, reproducibilityRequested);
    }

    ret             = do_cpt_state(gmx_fio_getxdr(fp), headerContents->flags_state, state, nullptr);
    *init_fep_state = state->fep_state; /* there should be a better way to do this than setting it
                                           here. Investigate for 5.0. */
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), headerContents->flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }
    state->ekinstate.hasReadEkinState = (((headerContents->flags_eks & (1 << eeksEKINH)) != 0)
                                         || ((headerContents->flags_eks & (1 << eeksEKINF)) != 0)
                                         || ((headerContents->flags_eks & (1 << eeksEKINO)) != 0)
                                         || (((headerContents->flags_eks & (1 << eeksEKINSCALEF))
                                              | (headerContents->flags_eks & (1 << eeksEKINSCALEH))
                                              | (headerContents->flags_eks & (1 << eeksVSCALE)))
                                             != 0));

    if (headerContents->flags_enh && observablesHistory->energyHistory == nullptr)
    {
        observablesHistory->energyHistory = std::make_unique<energyhistory_t>();
    }
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE, headerContents->flags_enh,
                          observablesHistory->energyHistory.get(), nullptr);
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
        ret = doCptPullHist(gmx_fio_getxdr(fp), TRUE, headerContents->flagsPullHistory,
                            observablesHistory->pullHistory.get(), StatePart::pullHistory, nullptr);
        if (ret)
        {
            cp_error();
        }
    }

    if (headerContents->file_version < 6)
    {
        gmx_fatal(FARGS,
                  "Continuing from checkpoint files written before GROMACS 4.5 is not supported");
    }

    ret = do_cpt_df_hist(gmx_fio_getxdr(fp), headerContents->flags_dfh, headerContents->nlambda,
                         &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->nED > 0 && observablesHistory->edsamHistory == nullptr)
    {
        observablesHistory->edsamHistory = std::make_unique<edsamhistory_t>(edsamhistory_t{});
    }
    ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, headerContents->nED,
                         observablesHistory->edsamHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->flags_awhh != 0 && state->awhHistory == nullptr)
    {
        state->awhHistory = std::make_shared<gmx::AwhHistory>();
    }
    ret = do_cpt_awh(gmx_fio_getxdr(fp), TRUE, headerContents->flags_awhh, state->awhHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    if (headerContents->eSwapCoords != eswapNO && observablesHistory->swapHistory == nullptr)
    {
        observablesHistory->swapHistory = std::make_unique<swaphistory_t>(swaphistory_t{});
    }
    ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, headerContents->eSwapCoords,
                           observablesHistory->swapHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    std::vector<gmx_file_position_t> outputfiles;
    ret = do_cpt_files(gmx_fio_getxdr(fp), TRUE, &outputfiles, nullptr, headerContents->file_version);
    if (ret)
    {
        cp_error();
    }
    do_cpt_mdmodules(headerContents->file_version, fp, mdModulesNotifier);
    ret = do_cpt_footer(gmx_fio_getxdr(fp), headerContents->file_version);
    if (ret)
    {
        cp_error();
    }
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}


void load_checkpoint(const char*                   fn,
                     t_fileio*                     logfio,
                     const t_commrec*              cr,
                     const ivec                    dd_nc,
                     t_inputrec*                   ir,
                     t_state*                      state,
                     ObservablesHistory*           observablesHistory,
                     gmx_bool                      reproducibilityRequested,
                     const gmx::MdModulesNotifier& mdModulesNotifier)
{
    CheckpointHeaderContents headerContents;
    if (SIMMASTER(cr))
    {
        /* Read the state from the checkpoint file */
        read_checkpoint(fn, logfio, cr, dd_nc, ir->eI, &(ir->fepvals->init_fep_state), &headerContents,
                        state, observablesHistory, reproducibilityRequested, mdModulesNotifier);
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(headerContents.step), &headerContents.step, cr);
        gmx::MdModulesCheckpointReadingBroadcast broadcastCheckPointData = { *cr, headerContents.file_version };
        mdModulesNotifier.notifier_.notify(broadcastCheckPointData);
    }
    ir->bContinuation = TRUE;
    // TODO Should the following condition be <=? Currently if you
    // pass a checkpoint written by an normal completion to a restart,
    // mdrun will read all input, does some work but no steps, and
    // write successful output. But perhaps that is not desirable.
    if ((ir->nsteps >= 0) && (ir->nsteps < headerContents.step))
    {
        // Note that we do not intend to support the use of mdrun
        // -nsteps to circumvent this condition.
        char nstepsString[STEPSTRSIZE], stepString[STEPSTRSIZE];
        gmx_step_str(ir->nsteps, nstepsString);
        gmx_step_str(headerContents.step, stepString);
        gmx_fatal(FARGS,
                  "The input requested %s steps, however the checkpoint "
                  "file has already reached step %s. The simulation will not "
                  "proceed, because either your simulation is already complete, "
                  "or your combination of input files don't match.",
                  nstepsString, stepString);
    }
    if (ir->nsteps >= 0)
    {
        ir->nsteps += ir->init_step - headerContents.step;
    }
    ir->init_step       = headerContents.step;
    ir->simulation_part = headerContents.simulation_part + 1;
}

void read_checkpoint_part_and_step(const char* filename, int* simulation_part, int64_t* step)
{
    t_fileio* fp;

    if (filename == nullptr || !gmx_fexist(filename) || ((fp = gmx_fio_open(filename, "r")) == nullptr))
    {
        *simulation_part = 0;
        *step            = 0;
        return;
    }

    CheckpointHeaderContents headerContents;
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, nullptr, &headerContents);
    gmx_fio_close(fp);
    *simulation_part = headerContents.simulation_part;
    *step            = headerContents.step;
}

static CheckpointHeaderContents read_checkpoint_data(t_fileio*                         fp,
                                                     t_state*                          state,
                                                     std::vector<gmx_file_position_t>* outputfiles)
{
    CheckpointHeaderContents headerContents;
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, nullptr, &headerContents);
    state->natoms        = headerContents.natoms;
    state->ngtc          = headerContents.ngtc;
    state->nnhpres       = headerContents.nnhpres;
    state->nhchainlength = headerContents.nhchainlength;
    state->flags         = headerContents.flags_state;
    int ret              = do_cpt_state(gmx_fio_getxdr(fp), state->flags, state, nullptr);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), headerContents.flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE, headerContents.flags_enh, &enerhist, nullptr);
    if (ret)
    {
        cp_error();
    }
    PullHistory pullHist = {};
    ret = doCptPullHist(gmx_fio_getxdr(fp), TRUE, headerContents.flagsPullHistory, &pullHist,
                        StatePart::pullHistory, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_df_hist(gmx_fio_getxdr(fp), headerContents.flags_dfh, headerContents.nlambda,
                         &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    edsamhistory_t edsamhist = {};
    ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, headerContents.nED, &edsamhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_awh(gmx_fio_getxdr(fp), TRUE, headerContents.flags_awhh, state->awhHistory.get(), nullptr);
    if (ret)
    {
        cp_error();
    }

    swaphistory_t swaphist = {};
    ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, headerContents.eSwapCoords, &swaphist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_files(gmx_fio_getxdr(fp), TRUE, outputfiles, nullptr, headerContents.file_version);

    if (ret)
    {
        cp_error();
    }
    gmx::MdModulesNotifier mdModuleNotifier;
    do_cpt_mdmodules(headerContents.file_version, fp, mdModuleNotifier);
    ret = do_cpt_footer(gmx_fio_getxdr(fp), headerContents.file_version);
    if (ret)
    {
        cp_error();
    }
    return headerContents;
}

void read_checkpoint_trxframe(t_fileio* fp, t_trxframe* fr)
{
    t_state                          state;
    std::vector<gmx_file_position_t> outputfiles;
    CheckpointHeaderContents headerContents = read_checkpoint_data(fp, &state, &outputfiles);

    fr->natoms    = state.natoms;
    fr->bStep     = TRUE;
    fr->step      = int64_to_int(headerContents.step, "conversion of checkpoint to trajectory");
    fr->bTime     = TRUE;
    fr->time      = headerContents.t;
    fr->bLambda   = TRUE;
    fr->lambda    = state.lambda[efptFEP];
    fr->fep_state = state.fep_state;
    fr->bAtoms    = FALSE;
    fr->bX        = ((state.flags & (1 << estX)) != 0);
    if (fr->bX)
    {
        fr->x = makeRvecArray(state.x, state.natoms);
    }
    fr->bV = ((state.flags & (1 << estV)) != 0);
    if (fr->bV)
    {
        fr->v = makeRvecArray(state.v, state.natoms);
    }
    fr->bF   = FALSE;
    fr->bBox = ((state.flags & (1 << estBOX)) != 0);
    if (fr->bBox)
    {
        copy_mat(state.box, fr->box);
    }
}

void list_checkpoint(const char* fn, FILE* out)
{
    t_fileio* fp;
    int       ret;

    t_state state;

    fp = gmx_fio_open(fn, "r");
    CheckpointHeaderContents headerContents;
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, out, &headerContents);
    state.natoms        = headerContents.natoms;
    state.ngtc          = headerContents.ngtc;
    state.nnhpres       = headerContents.nnhpres;
    state.nhchainlength = headerContents.nhchainlength;
    state.flags         = headerContents.flags_state;
    ret                 = do_cpt_state(gmx_fio_getxdr(fp), state.flags, &state, out);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), headerContents.flags_eks, &state.ekinstate, out);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE, headerContents.flags_enh, &enerhist, out);

    if (ret == 0)
    {
        PullHistory pullHist = {};
        ret = doCptPullHist(gmx_fio_getxdr(fp), TRUE, headerContents.flagsPullHistory, &pullHist,
                            StatePart::pullHistory, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_df_hist(gmx_fio_getxdr(fp), headerContents.flags_dfh, headerContents.nlambda,
                             &state.dfhist, out);
    }

    if (ret == 0)
    {
        edsamhistory_t edsamhist = {};
        ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, headerContents.nED, &edsamhist, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_awh(gmx_fio_getxdr(fp), TRUE, headerContents.flags_awhh, state.awhHistory.get(), out);
    }

    if (ret == 0)
    {
        swaphistory_t swaphist = {};
        ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, headerContents.eSwapCoords, &swaphist, out);
    }

    if (ret == 0)
    {
        std::vector<gmx_file_position_t> outputfiles;
        ret = do_cpt_files(gmx_fio_getxdr(fp), TRUE, &outputfiles, out, headerContents.file_version);
    }

    if (ret == 0)
    {
        ret = do_cpt_footer(gmx_fio_getxdr(fp), headerContents.file_version);
    }

    if (ret)
    {
        cp_warning(out);
    }
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}

/* This routine cannot print tons of data, since it is called before the log file is opened. */
CheckpointHeaderContents read_checkpoint_simulation_part_and_filenames(t_fileio* fp,
                                                                       std::vector<gmx_file_position_t>* outputfiles)
{
    t_state                  state;
    CheckpointHeaderContents headerContents = read_checkpoint_data(fp, &state, outputfiles);
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
    return headerContents;
}
