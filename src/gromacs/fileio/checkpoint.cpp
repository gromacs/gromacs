/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <fcntl.h>
#if GMX_NATIVE_WINDOWS
#include <io.h>
#include <sys/locking.h>
#endif

#include "buildinfo.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/int64_to_int.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/txtdump.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819
#define CPTSTRLEN 1024

/* cpt_version should normally only be changed
 * when the header or footer format changes.
 * The state data format itself is backward and forward compatible.
 * But old code can not read a new entry that is present in the file
 * (but can read a new format when new entries are not present).
 */
static const int cpt_version = 16;


const char *est_names[estNR] =
{
    "FE-lambda",
    "box", "box-rel", "box-v", "pres_prev",
    "nosehoover-xi", "thermostat-integral",
    "x", "v", "sdx-unsupported", "CGp", "LD-rng-unsupported", "LD-rng-i-unsupported",
    "disre_initf", "disre_rm3tav",
    "orire_initf", "orire_Dtav",
    "svir_prev", "nosehoover-vxi", "v_eta", "vol0", "nhpres_xi", "nhpres_vxi", "fvir_prev", "fep_state", "MC-rng-unsupported", "MC-rng-i-unsupported"
};

enum {
    eeksEKIN_N, eeksEKINH, eeksDEKINDL, eeksMVCOS, eeksEKINF, eeksEKINO, eeksEKINSCALEF, eeksEKINSCALEH, eeksVSCALE, eeksEKINTOTAL, eeksNR
};

const char *eeks_names[eeksNR] =
{
    "Ekin_n", "Ekinh", "dEkindlambda", "mv_cos",
    "Ekinf", "Ekinh_old", "EkinScaleF_NHC", "EkinScaleH_NHC", "Vscale_NHC", "Ekin_Total"
};

enum {
    eenhENERGY_N, eenhENERGY_AVER, eenhENERGY_SUM, eenhENERGY_NSUM,
    eenhENERGY_SUM_SIM, eenhENERGY_NSUM_SIM,
    eenhENERGY_NSTEPS, eenhENERGY_NSTEPS_SIM,
    eenhENERGY_DELTA_H_NN,
    eenhENERGY_DELTA_H_LIST,
    eenhENERGY_DELTA_H_STARTTIME,
    eenhENERGY_DELTA_H_STARTLAMBDA,
    eenhNR
};

const char *eenh_names[eenhNR] =
{
    "energy_n", "energy_aver", "energy_sum", "energy_nsum",
    "energy_sum_sim", "energy_nsum_sim",
    "energy_nsteps", "energy_nsteps_sim",
    "energy_delta_h_nn",
    "energy_delta_h_list",
    "energy_delta_h_start_time",
    "energy_delta_h_start_lambda"
};

/* free energy history variables -- need to be preserved over checkpoint */
enum {
    edfhBEQUIL, edfhNATLAMBDA, edfhWLHISTO, edfhWLDELTA, edfhSUMWEIGHTS, edfhSUMDG, edfhSUMMINVAR, edfhSUMVAR,
    edfhACCUMP, edfhACCUMM, edfhACCUMP2, edfhACCUMM2, edfhTIJ, edfhTIJEMP, edfhNR
};
/* free energy history variable names  */
const char *edfh_names[edfhNR] =
{
    "bEquilibrated", "N_at_state", "Wang-Landau Histogram", "Wang-Landau Delta", "Weights", "Free Energies", "minvar", "variance",
    "accumulated_plus", "accumulated_minus", "accumulated_plus_2",  "accumulated_minus_2", "Tij", "Tij_empirical"
};

//! Higher level vector element type, only used for formatting checkpoint dumps
enum class CptElementType
{
    integer,   //!< integer
    real,      //!< float or double, not linked to precision of type real
    real3,     //!< float[3] or double[3], not linked to precision of type real
    matrix3x3  //!< float[3][3] or double[3][3], not linked to precision of type real
};

//! \brief Parts of the checkpoint state, only used for reporting
enum class StatePart
{
    microState,       //!< The microstate of the simulated system
    kineticEnergy,    //!< Kinetic energy, needed for T/P-coupling state
    energyHistory,    //!< Energy observable statistics
    freeEnergyHistory //!< Free-energy state and observable statistics
};

//! \brief Return the name of a checkpoint entry based on part and part entry
static const char *entryName(StatePart part, int ecpt)
{
    switch (part)
    {
        case StatePart::microState:        return est_names [ecpt];
        case StatePart::kineticEnergy:     return eeks_names[ecpt];
        case StatePart::energyHistory:     return eenh_names[ecpt];
        case StatePart::freeEnergyHistory: return edfh_names[ecpt];
    }

    return nullptr;
}

static void cp_warning(FILE *fp)
{
    fprintf(fp, "\nWARNING: Checkpoint file is corrupted or truncated\n\n");
}

static void cp_error()
{
    gmx_fatal(FARGS, "Checkpoint file corrupted/truncated, or maybe you are out of disk space?");
}

static void do_cpt_string_err(XDR *xd, gmx_bool bRead, const char *desc, char **s, FILE *list)
{
    if (bRead)
    {
        snew(*s, CPTSTRLEN);
    }
    if (xdr_string(xd, s, CPTSTRLEN) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, *s);
        sfree(*s);
    }
}

static int do_cpt_int(XDR *xd, const char *desc, int *i, FILE *list)
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

static int do_cpt_u_chars(XDR *xd, const char *desc, int n, unsigned char *i, FILE *list)
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

static void do_cpt_int_err(XDR *xd, const char *desc, int *i, FILE *list)
{
    if (do_cpt_int(xd, desc, i, list) < 0)
    {
        cp_error();
    }
}

static void do_cpt_step_err(XDR *xd, const char *desc, gmx_int64_t *i, FILE *list)
{
    char   buf[STEPSTRSIZE];

    if (xdr_int64(xd, i) == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list, "%s = %s\n", desc, gmx_step_str(*i, buf));
    }
}

static void do_cpt_double_err(XDR *xd, const char *desc, double *f, FILE *list)
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

static void do_cpt_real_err(XDR *xd, real *f)
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

static void do_cpt_n_rvecs_err(XDR *xd, const char *desc, int n, rvec f[], FILE *list)
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

template <typename T>
struct xdr_type
{
};

template <>
struct xdr_type<int>
{
    // cppcheck-suppress unusedStructMember
    static const int value = xdr_datatype_int;
};

template <>
struct xdr_type<float>
{
    // cppcheck-suppress unusedStructMember
    static const int value = xdr_datatype_float;
};

template <>
struct xdr_type<double>
{
    // cppcheck-suppress unusedStructMember
    static const int value = xdr_datatype_double;
};

//! \brief Returns size in byte of an xdr_datatype
static inline unsigned int sizeOfXdrType(int xdrType)
{
    switch (xdrType)
    {
        case xdr_datatype_int:
            return sizeof(int);
            break;
        case xdr_datatype_float:
            return sizeof(float);
            break;
        case xdr_datatype_double:
            return sizeof(double);
            break;
        default: GMX_RELEASE_ASSERT(false, "XDR data type not implemented");
    }

    return 0;
}

//! \brief Returns the XDR process function for i/o of an XDR type
static inline xdrproc_t xdrProc(int xdrType)
{
    switch (xdrType)
    {
        case xdr_datatype_int:
            return reinterpret_cast<xdrproc_t>(xdr_int);
            break;
        case xdr_datatype_float:
            return reinterpret_cast<xdrproc_t>(xdr_float);
            break;
        case xdr_datatype_double:
            return reinterpret_cast<xdrproc_t>(xdr_double);
            break;
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
static bool_t listXdrVector(XDR *xd, StatePart part, int ecpt, int nf, int xdrType,
                            FILE *list, CptElementType cptElementType)
{
    bool_t             res = 0;

    const unsigned int elemSize = sizeOfXdrType(xdrType);
    std::vector<char>  data(nf*elemSize);
    res = xdr_vector(xd, data.data(), nf, elemSize, xdrProc(xdrType));

    if (list != nullptr)
    {
        switch (xdrType)
        {
            case xdr_datatype_int:
                pr_ivec(list, 0, entryName(part, ecpt), reinterpret_cast<const int *>(data.data()), nf, TRUE);
                break;
            case xdr_datatype_float:
#if !GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    // cppcheck-suppress invalidPointerCast
                    pr_rvecs(list, 0, entryName(part, ecpt), reinterpret_cast<const rvec *>(data.data()), nf/3);
                }
                else
#endif
                {
                    /* Note: With double precision code dumping a single precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    // cppcheck-suppress invalidPointerCast
                    pr_fvec(list, 0, entryName(part, ecpt), reinterpret_cast<const float *>(data.data()), nf, TRUE);
                }
                break;
            case xdr_datatype_double:
#if GMX_DOUBLE
                if (cptElementType == CptElementType::real3)
                {
                    // cppcheck-suppress invalidPointerCast
                    pr_rvecs(list, 0, entryName(part, ecpt), reinterpret_cast<const rvec *>(data.data()), nf/3);
                }
                else
#endif
                {
                    /* Note: With single precision code dumping a double precision rvec will produce float iso rvec print, but that's a minor annoyance */
                    // cppcheck-suppress invalidPointerCast
                    pr_dvec(list, 0, entryName(part, ecpt), reinterpret_cast<const double *>(data.data()), nf, TRUE);
                }
                break;
            default: GMX_RELEASE_ASSERT(false, "Data type not implemented for listing");
        }
    }

    return res;
}

//! \brief Convert a double array, typed char*, to float
static void convertArrayRealPrecision(const char *c, float *v, int n)
{
    // cppcheck-suppress invalidPointerCast
    const double *d = reinterpret_cast<const double *>(c);
    for (int i = 0; i < n; i++)
    {
        v[i] = static_cast<float>(d[i]);
    }
}

//! \brief Convert a float array, typed char*, to double
static void convertArrayRealPrecision(const char *c, double *v, int n)
{
    // cppcheck-suppress invalidPointerCast
    const float *f = reinterpret_cast<const float *>(c);
    for (int i = 0; i < n; i++)
    {
        v[i] = static_cast<double>(f[i]);
    }
}

//! \brief Generate an error for trying to convert to integer
static void convertArrayRealPrecision(const char gmx_unused *c, int gmx_unused *v, int gmx_unused n)
{
    GMX_RELEASE_ASSERT(false, "We only expect type mismatches between float and double, not integer");
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
 * When not listing, we use either v or vector, depending on which is !=NULL.
 * If nval >= 0, nval is used; on read this should match the passed value.
 * If nval n<0, *nptr (with v) or vector->size() is used. On read using v,
 * the value is stored in nptr
 */
template<typename T>
static int doVectorLow(XDR *xd, StatePart part, int ecpt, int sflags,
                       int nval, int *nptr,
                       T **v, std::vector<T> *vector,
                       FILE *list, CptElementType cptElementType)
{
    GMX_RELEASE_ASSERT(list != nullptr || (v != nullptr && vector == nullptr) || (v == nullptr && vector != nullptr), "Without list, we should have exactly one of v and vector != NULL");

    bool_t res = 0;

    int    numElemInTheFile;
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
                // cppcheck-suppress nullPointer
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
    gmx_constexpr int xdrTypeInTheCode = xdr_type<T>::value;
    int               xdrTypeInTheFile = xdrTypeInTheCode;
    res = xdr_int(xd, &xdrTypeInTheFile);
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
                gmx_fatal(FARGS, "Count mismatch for state entry %s, code count is %d, file count is %d\n", entryName(part, ecpt), nval, numElemInTheFile);
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
                    entryName(part, ecpt),
                    xdr_datatype_names[xdrTypeInTheCode],
                    xdr_datatype_names[xdrTypeInTheFile]);

            /* Matching int and real should never occur, but check anyhow */
            if (xdrTypeInTheFile == xdr_datatype_int ||
                xdrTypeInTheCode == xdr_datatype_int)
            {
                gmx_fatal(FARGS, "Type %s: incompatible checkpoint formats or corrupted checkpoint file.", buf);
            }
            fprintf(stderr, "Precision %s\n", buf);
        }

        T *vp;
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
             * PaddedRVecVector has a size of numElemInThefile and we
             * don't want to lose that padding here.
             */
            if (vector->size() < static_cast<unsigned int>(numElemInTheFile))
            {
                vector->resize(numElemInTheFile);
            }
            vp = vector->data();
        }

        char *vChar;
        if (typesMatch)
        {
            vChar = reinterpret_cast<char *>(vp);
        }
        else
        {
            snew(vChar, numElemInTheFile*sizeOfXdrType(xdrTypeInTheFile));
        }
        res = xdr_vector(xd, vChar,
                         numElemInTheFile, sizeOfXdrType(xdrTypeInTheFile),
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
        res = listXdrVector(xd, part, ecpt, numElemInTheFile, xdrTypeInTheFile,
                            list, cptElementType);
    }

    return 0;
}

//! \brief Read/Write a std::vector.
template <typename T>
static int doVector(XDR *xd, StatePart part, int ecpt, int sflags,
                    std::vector<T> *vector, FILE *list)
{
    return doVectorLow<T>(xd, part, ecpt, sflags, -1, nullptr, nullptr, vector, list, CptElementType::real);
}

//! \brief Read/Write an ArrayRef<real>.
static int doRealArrayRef(XDR *xd, StatePart part, int ecpt, int sflags,
                          gmx::ArrayRef<real> vector, FILE *list)
{
    real *v_real = vector.data();
    return doVectorLow<real>(xd, part, ecpt, sflags, vector.size(), nullptr, &v_real, nullptr, list, CptElementType::real);
}

//! \brief Read/Write a PaddedRVecVector.
static int doPaddedRvecVector(XDR *xd, StatePart part, int ecpt, int sflags,
                              gmx::PaddedRVecVector *vector, FILE *list)
{
    real *v_real;

    if (list == nullptr && (sflags & (1 << ecpt)))
    {
        v_real = vector->data()->as_vec();
    }
    else
    {
        v_real = nullptr;
    }
    // The current invariant of a PaddedRVecVector is that its size is
    // one larger than necessary to store the data. Make sure that we
    // read/write only the valid data, and don't leak to the outside
    // world that currently we find it convenient internally to
    // allocate one extra element.
    gmx::ArrayRef<real> ref(v_real, v_real + (vector->size()-1) * DIM);

    return doRealArrayRef(xd, part, ecpt, sflags, ref, list);
}

/* This function stores n along with the reals for reading,
 * but on reading it assumes that n matches the value in the checkpoint file,
 * a fatal error is generated when this is not the case.
 */
static int do_cpte_reals(XDR *xd, StatePart part, int ecpt, int sflags,
                         int n, real **v, FILE *list)
{
    return doVectorLow<real>(xd, part, ecpt, sflags, n, nullptr, v, nullptr, list, CptElementType::real);
}

/* This function does the same as do_cpte_reals,
 * except that on reading it ignores the passed value of *n
 * and stores the value read from the checkpoint file in *n.
 */
static int do_cpte_n_reals(XDR *xd, StatePart part, int ecpt, int sflags,
                           int *n, real **v, FILE *list)
{
    return doVectorLow<real>(xd, part, ecpt, sflags, -1, n, v, nullptr, list, CptElementType::real);
}

static int do_cpte_real(XDR *xd, StatePart part, int ecpt, int sflags,
                        real *r, FILE *list)
{
    return doVectorLow<real>(xd, part, ecpt, sflags, 1, nullptr, &r, nullptr, list, CptElementType::real);
}

static int do_cpte_ints(XDR *xd, StatePart part, int ecpt, int sflags,
                        int n, int **v, FILE *list)
{
    return doVectorLow<int>(xd, part, ecpt, sflags, n, nullptr, v, nullptr, list, CptElementType::integer);
}

static int do_cpte_int(XDR *xd, StatePart part, int ecpt, int sflags,
                       int *i, FILE *list)
{
    return do_cpte_ints(xd, part, ecpt, sflags, 1, &i, list);
}

static int do_cpte_doubles(XDR *xd, StatePart part, int ecpt, int sflags,
                           int n, double **v, FILE *list)
{
    return doVectorLow<double>(xd, part, ecpt, sflags, n, nullptr, v, nullptr, list, CptElementType::real);
}

static int do_cpte_double(XDR *xd, StatePart part, int ecpt, int sflags,
                          double *r, FILE *list)
{
    return do_cpte_doubles(xd, part, ecpt, sflags, 1, &r, list);
}

static int do_cpte_matrix(XDR *xd, StatePart part, int ecpt, int sflags,
                          matrix v, FILE *list)
{
    real *vr;
    int   ret;

    vr  = &(v[0][0]);
    ret = doVectorLow<real>(xd, part, ecpt, sflags,
                            DIM*DIM, nullptr, &vr, nullptr, nullptr, CptElementType::matrix3x3);

    if (list && ret == 0)
    {
        pr_rvecs(list, 0, entryName(part, ecpt), v, DIM);
    }

    return ret;
}


static int do_cpte_nmatrix(XDR *xd, StatePart part, int ecpt, int sflags,
                           int n, real **v, FILE *list)
{
    int   i;
    int   ret, reti;
    char  name[CPTSTRLEN];

    ret = 0;
    if (v == nullptr)
    {
        snew(v, n);
    }
    for (i = 0; i < n; i++)
    {
        reti = doVectorLow<real>(xd, part, ecpt, sflags, n, nullptr, &(v[i]), nullptr, nullptr, CptElementType::matrix3x3);
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

static int do_cpte_matrices(XDR *xd, StatePart part, int ecpt, int sflags,
                            int n, matrix **v, FILE *list)
{
    bool_t  res = 0;
    matrix *vp, *va = nullptr;
    real   *vr;
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
        gmx_fatal(FARGS, "Count mismatch for state entry %s, code count is %d, file count is %d\n", entryName(part, ecpt), n, nf);
    }
    if (list || !(sflags & (1<<ecpt)))
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
    snew(vr, nf*DIM*DIM);
    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                vr[(i*DIM+j)*DIM+k] = vp[i][j][k];
            }
        }
    }
    ret = doVectorLow<real>(xd, part, ecpt, sflags,
                            nf*DIM*DIM, nullptr, &vr, nullptr, nullptr,
                            CptElementType::matrix3x3);
    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                vp[i][j][k] = vr[(i*DIM+j)*DIM+k];
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

static void do_cpt_header(XDR *xd, gmx_bool bRead, int *file_version,
                          char **version, char **btime, char **buser, char **bhost,
                          int *double_prec,
                          char **fprog, char **ftime,
                          int *eIntegrator, int *simulation_part,
                          gmx_int64_t *step, double *t,
                          int *nnodes, int *dd_nc, int *npme,
                          int *natoms, int *ngtc, int *nnhpres, int *nhchainlength,
                          int *nlambda, int *flags_state,
                          int *flags_eks, int *flags_enh, int *flags_dfh,
                          int *nED, int *eSwapCoords,
                          FILE *list)
{
    bool_t res = 0;
    int    magic;
    int    idum = 0;
    char  *fhost;

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
        gmx_fatal(FARGS, "The checkpoint file is empty/corrupted, or maybe you are out of disk space?");
    }
    if (magic != CPT_MAGIC1)
    {
        gmx_fatal(FARGS, "Start of file magic number mismatch, checkpoint file has %d, should be %d\n"
                  "The checkpoint file is corrupted or not a checkpoint file",
                  magic, CPT_MAGIC1);
    }
    if (!bRead)
    {
        snew(fhost, 255);
        gmx_gethostname(fhost, 255);
    }
    do_cpt_string_err(xd, bRead, "GROMACS version", version, list);
    do_cpt_string_err(xd, bRead, "GROMACS build time", btime, list);
    do_cpt_string_err(xd, bRead, "GROMACS build user", buser, list);
    do_cpt_string_err(xd, bRead, "GROMACS build host", bhost, list);
    do_cpt_string_err(xd, bRead, "generating program", fprog, list);
    do_cpt_string_err(xd, bRead, "generation time", ftime, list);
    *file_version = cpt_version;
    do_cpt_int_err(xd, "checkpoint file version", file_version, list);
    if (*file_version > cpt_version)
    {
        gmx_fatal(FARGS, "Attempting to read a checkpoint file of version %d with code of version %d\n", *file_version, cpt_version);
    }
    if (*file_version >= 13)
    {
        do_cpt_int_err(xd, "GROMACS double precision", double_prec, list);
    }
    else
    {
        *double_prec = -1;
    }
    if (*file_version >= 12)
    {
        do_cpt_string_err(xd, bRead, "generating host", &fhost, list);
        if (list == nullptr)
        {
            sfree(fhost);
        }
    }
    do_cpt_int_err(xd, "#atoms", natoms, list);
    do_cpt_int_err(xd, "#T-coupling groups", ngtc, list);
    if (*file_version >= 10)
    {
        do_cpt_int_err(xd, "#Nose-Hoover T-chains", nhchainlength, list);
    }
    else
    {
        *nhchainlength = 1;
    }
    if (*file_version >= 11)
    {
        do_cpt_int_err(xd, "#Nose-Hoover T-chains for barostat ", nnhpres, list);
    }
    else
    {
        *nnhpres = 0;
    }
    if (*file_version >= 14)
    {
        do_cpt_int_err(xd, "# of total lambda states ", nlambda, list);
    }
    else
    {
        *nlambda = 0;
    }
    do_cpt_int_err(xd, "integrator", eIntegrator, list);
    if (*file_version >= 3)
    {
        do_cpt_int_err(xd, "simulation part #", simulation_part, list);
    }
    else
    {
        *simulation_part = 1;
    }
    if (*file_version >= 5)
    {
        do_cpt_step_err(xd, "step", step, list);
    }
    else
    {
        do_cpt_int_err(xd, "step", &idum, list);
        *step = idum;
    }
    do_cpt_double_err(xd, "t", t, list);
    do_cpt_int_err(xd, "#PP-ranks", nnodes, list);
    idum = 1;
    do_cpt_int_err(xd, "dd_nc[x]", dd_nc ? &(dd_nc[0]) : &idum, list);
    do_cpt_int_err(xd, "dd_nc[y]", dd_nc ? &(dd_nc[1]) : &idum, list);
    do_cpt_int_err(xd, "dd_nc[z]", dd_nc ? &(dd_nc[2]) : &idum, list);
    do_cpt_int_err(xd, "#PME-only ranks", npme, list);
    do_cpt_int_err(xd, "state flags", flags_state, list);
    if (*file_version >= 4)
    {
        do_cpt_int_err(xd, "ekin data flags", flags_eks, list);
        do_cpt_int_err(xd, "energy history flags", flags_enh, list);
    }
    else
    {
        *flags_eks   = 0;
        *flags_enh   = (*flags_state >> (estORIRE_DTAV+1));
        *flags_state = (*flags_state & ~((1<<(estORIRE_DTAV+1)) |
                                         (1<<(estORIRE_DTAV+2)) |
                                         (1<<(estORIRE_DTAV+3))));
    }
    if (*file_version >= 14)
    {
        do_cpt_int_err(xd, "df history flags", flags_dfh, list);
    }
    else
    {
        *flags_dfh = 0;
    }

    if (*file_version >= 15)
    {
        do_cpt_int_err(xd, "ED data sets", nED, list);
    }
    else
    {
        *nED = 0;
    }
    if (*file_version >= 16)
    {
        do_cpt_int_err(xd, "swap", eSwapCoords, list);
    }
    else
    {
        *eSwapCoords = eswapNO;
    }
}

static int do_cpt_footer(XDR *xd, int file_version)
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

static int do_cpt_state(XDR *xd,
                        int fflags, t_state *state,
                        FILE *list)
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
        if (fflags & (1<<i))
        {
            switch (i)
            {
                case estLAMBDA:  ret      = doRealArrayRef(xd, part, i, sflags, gmx::arrayRefFromArray<real>(state->lambda.data(), state->lambda.size()), list); break;
                case estFEPSTATE: ret     = do_cpte_int (xd, part, i, sflags, &state->fep_state, list); break;
                case estBOX:     ret      = do_cpte_matrix(xd, part, i, sflags, state->box, list); break;
                case estBOX_REL: ret      = do_cpte_matrix(xd, part, i, sflags, state->box_rel, list); break;
                case estBOXV:    ret      = do_cpte_matrix(xd, part, i, sflags, state->boxv, list); break;
                case estPRES_PREV: ret    = do_cpte_matrix(xd, part, i, sflags, state->pres_prev, list); break;
                case estSVIR_PREV:  ret   = do_cpte_matrix(xd, part, i, sflags, state->svir_prev, list); break;
                case estFVIR_PREV:  ret   = do_cpte_matrix(xd, part, i, sflags, state->fvir_prev, list); break;
                case estNH_XI:   ret      = doVector<double>(xd, part, i, sflags, &state->nosehoover_xi, list); break;
                case estNH_VXI:  ret      = doVector<double>(xd, part, i, sflags, &state->nosehoover_vxi, list); break;
                case estNHPRES_XI:   ret  = doVector<double>(xd, part, i, sflags, &state->nhpres_xi, list); break;
                case estNHPRES_VXI:  ret  = doVector<double>(xd, part, i, sflags, &state->nhpres_vxi, list); break;
                case estTC_INT:  ret      = doVector<double>(xd, part, i, sflags, &state->therm_integral, list); break;
                case estVETA:    ret      = do_cpte_real(xd, part, i, sflags, &state->veta, list); break;
                case estVOL0:    ret      = do_cpte_real(xd, part, i, sflags, &state->vol0, list); break;
                case estX:       ret      = doPaddedRvecVector(xd, part, i, sflags, &state->x, list); break;
                case estV:       ret      = doPaddedRvecVector(xd, part, i, sflags, &state->v, list); break;
                /* The RNG entries are no longer written,
                 * the next 4 lines are only for reading old files.
                 */
                case estLD_RNG_NOTSUPPORTED:  ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list); break;
                case estLD_RNGI_NOTSUPPORTED: ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list); break;
                case estMC_RNG_NOTSUPPORTED:  ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list); break;
                case estMC_RNGI_NOTSUPPORTED: ret = do_cpte_ints(xd, part, i, sflags, 0, nullptr, list); break;
                case estDISRE_INITF:  ret         = do_cpte_real (xd, part, i, sflags, &state->hist.disre_initf, list); break;
                case estDISRE_RM3TAV: ret         = do_cpte_n_reals(xd, part, i, sflags, &state->hist.ndisrepairs, &state->hist.disre_rm3tav, list); break;
                case estORIRE_INITF:  ret         = do_cpte_real (xd, part, i, sflags, &state->hist.orire_initf, list); break;
                case estORIRE_DTAV:   ret         = do_cpte_n_reals(xd, part, i, sflags, &state->hist.norire_Dtav, &state->hist.orire_Dtav, list); break;
                default:
                    gmx_fatal(FARGS, "Unknown state entry %d\n"
                              "You are reading a checkpoint file written by different code, which is not supported", i);
            }
        }
    }

    return ret;
}

static int do_cpt_ekinstate(XDR *xd, int fflags, ekinstate_t *ekins,
                            FILE *list)
{
    int             ret  = 0;

    const StatePart part = StatePart::kineticEnergy;
    for (int i = 0; (i < eeksNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {

                case eeksEKIN_N:     ret = do_cpte_int(xd, part, i, fflags, &ekins->ekin_n, list); break;
                case eeksEKINH:     ret  = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinh, list); break;
                case eeksEKINF:      ret = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinf, list); break;
                case eeksEKINO:      ret = do_cpte_matrices(xd, part, i, fflags, ekins->ekin_n, &ekins->ekinh_old, list); break;
                case eeksEKINTOTAL:  ret = do_cpte_matrix(xd, part, i, fflags, ekins->ekin_total, list); break;
                case eeksEKINSCALEF: ret = doVector<double>(xd, part, i, fflags, &ekins->ekinscalef_nhc, list); break;
                case eeksVSCALE:     ret = doVector<double>(xd, part, i, fflags, &ekins->vscale_nhc, list); break;
                case eeksEKINSCALEH: ret = doVector<double>(xd, part, i, fflags, &ekins->ekinscaleh_nhc, list); break;
                case eeksDEKINDL:   ret  = do_cpte_real(xd, part, i, fflags, &ekins->dekindl, list); break;
                case eeksMVCOS:      ret = do_cpte_real(xd, part, i, fflags, &ekins->mvcos, list); break;
                default:
                    gmx_fatal(FARGS, "Unknown ekin data state entry %d\n"
                              "You are probably reading a new checkpoint file with old code", i);
            }
        }
    }

    return ret;
}


static int do_cpt_swapstate(XDR *xd, gmx_bool bRead,
                            int eSwapCoords, swapstate_t **swapstatePtr, FILE *list)
{
    int swap_cpt_version = 2;

    if (eSwapCoords == eswapNO)
    {
        return 0;
    }

    if (*swapstatePtr == nullptr)
    {
        snew(*swapstatePtr, 1);
    }
    swapstate_t *swapstate = *swapstatePtr;
    swapstate->bFromCpt    = bRead;
    swapstate->eSwapCoords = eSwapCoords;

    do_cpt_int_err(xd, "swap checkpoint version", &swap_cpt_version, list);
    if (bRead && swap_cpt_version < 2)
    {
        gmx_fatal(FARGS, "Cannot read checkpoint files that were written with old versions"
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
            swapstateIons_t *gs = &swapstate->ionType[ii];

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

            if (bRead && (nullptr == gs->nMolPast[ic]) )
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
            swapstateIons_t *gs = &swapstate->ionType[ii];

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
        swapstateIons_t *gs = &swapstate->ionType[ii];

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
        do_cpt_n_rvecs_err(xd, "Ch0 whole x", swapstate->nat[eChan0], *swapstate->xc_old_whole_p[eChan0], list);
        do_cpt_n_rvecs_err(xd, "Ch1 whole x", swapstate->nat[eChan1], *swapstate->xc_old_whole_p[eChan1], list);
    }

    return 0;
}


static int do_cpt_enerhist(XDR *xd, gmx_bool bRead,
                           int fflags, energyhistory_t *enerhist,
                           FILE *list)
{
    int ret = 0;

    /* This is stored/read for backward compatibility */
    int  energyHistoryNumEnergies = 0;
    if (bRead)
    {
        enerhist->nsteps     = 0;
        enerhist->nsum       = 0;
        enerhist->nsteps_sim = 0;
        enerhist->nsum_sim   = 0;
    }
    else
    {
        energyHistoryNumEnergies = enerhist->ener_sum_sim.size();
    }

    delta_h_history_t *deltaH = enerhist->deltaHForeignLambdas.get();
    const StatePart    part   = StatePart::energyHistory;
    for (int i = 0; (i < eenhNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {
                case eenhENERGY_N:     ret = do_cpte_int(xd, part, i, fflags, &energyHistoryNumEnergies, list); break;
                case eenhENERGY_AVER:  ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_ave, list); break;
                case eenhENERGY_SUM:   ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_sum, list); break;
                case eenhENERGY_NSUM:  do_cpt_step_err(xd, eenh_names[i], &enerhist->nsum, list); break;
                case eenhENERGY_SUM_SIM: ret = doVector<double>(xd, part, i, fflags, &enerhist->ener_sum_sim, list); break;
                case eenhENERGY_NSUM_SIM:   do_cpt_step_err(xd, eenh_names[i], &enerhist->nsum_sim, list); break;
                case eenhENERGY_NSTEPS:     do_cpt_step_err(xd, eenh_names[i], &enerhist->nsteps, list); break;
                case eenhENERGY_NSTEPS_SIM: do_cpt_step_err(xd, eenh_names[i], &enerhist->nsteps_sim, list); break;
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
                            enerhist->deltaHForeignLambdas.reset(new delta_h_history_t);
                            deltaH = enerhist->deltaHForeignLambdas.get();
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
                    ret = do_cpte_double(xd, part, i, fflags, &(deltaH->start_time), list); break;
                case eenhENERGY_DELTA_H_STARTLAMBDA:
                    ret = do_cpte_double(xd, part, i, fflags, &(deltaH->start_lambda), list); break;
                default:
                    gmx_fatal(FARGS, "Unknown energy history entry %d\n"
                              "You are probably reading a new checkpoint file with old code", i);
            }
        }
    }

    if ((fflags & (1<<eenhENERGY_SUM)) && !(fflags & (1<<eenhENERGY_SUM_SIM)))
    {
        /* Assume we have an old file format and copy sum to sum_sim */
        enerhist->ener_sum_sim = enerhist->ener_sum;
    }

    if ( (fflags & (1<<eenhENERGY_NSUM)) &&
         !(fflags & (1<<eenhENERGY_NSTEPS)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps = enerhist->nsum;
    }
    if ( (fflags & (1<<eenhENERGY_NSUM_SIM)) &&
         !(fflags & (1<<eenhENERGY_NSTEPS_SIM)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps_sim = enerhist->nsum_sim;
    }

    return ret;
}

static int do_cpt_df_hist(XDR *xd, int fflags, int nlambda, df_history_t **dfhistPtr, FILE *list)
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
    df_history_t    *dfhist = *dfhistPtr;

    const StatePart  part   = StatePart::freeEnergyHistory;
    for (int i = 0; (i < edfhNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {
                case edfhBEQUIL:       ret = do_cpte_int(xd, part, i, fflags, &dfhist->bEquil, list); break;
                case edfhNATLAMBDA:    ret = do_cpte_ints(xd, part, i, fflags, nlambda, &dfhist->n_at_lam, list); break;
                case edfhWLHISTO:      ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->wl_histo, list); break;
                case edfhWLDELTA:      ret = do_cpte_real(xd, part, i, fflags, &dfhist->wl_delta, list); break;
                case edfhSUMWEIGHTS:   ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_weights, list); break;
                case edfhSUMDG:        ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_dg, list); break;
                case edfhSUMMINVAR:    ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_minvar, list); break;
                case edfhSUMVAR:       ret = do_cpte_reals(xd, part, i, fflags, nlambda, &dfhist->sum_variance, list); break;
                case edfhACCUMP:       ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_p, list); break;
                case edfhACCUMM:       ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_m, list); break;
                case edfhACCUMP2:      ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_p2, list); break;
                case edfhACCUMM2:      ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->accum_m2, list); break;
                case edfhTIJ:          ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->Tij, list); break;
                case edfhTIJEMP:       ret = do_cpte_nmatrix(xd, part, i, fflags, nlambda, dfhist->Tij_empirical, list); break;

                default:
                    gmx_fatal(FARGS, "Unknown df history entry %d\n"
                              "You are probably reading a new checkpoint file with old code", i);
            }
        }
    }

    return ret;
}


/* This function stores the last whole configuration of the reference and
 * average structure in the .cpt file
 */
static int do_cpt_EDstate(XDR *xd, gmx_bool bRead,
                          int nED, edsamstate_t **EDstatePtr, FILE *list)
{
    if (nED == 0)
    {
        return 0;
    }

    if (*EDstatePtr == nullptr)
    {
        snew(*EDstatePtr, 1);
    }
    edsamstate_t *EDstate = *EDstatePtr;

    EDstate->bFromCpt     = bRead;
    EDstate->nED          = nED;

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
        sprintf(buf, "ED%d # of atoms in reference structure", i+1);
        do_cpt_int_err(xd, buf, &EDstate->nref[i], list);
        sprintf(buf, "ED%d x_ref", i+1);
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
        sprintf(buf, "ED%d # of atoms in average structure", i+1);
        do_cpt_int_err(xd, buf, &EDstate->nav[i], list);
        sprintf(buf, "ED%d x_av", i+1);
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


static int do_cpt_files(XDR *xd, gmx_bool bRead,
                        gmx_file_position_t **p_outputfiles, int *nfiles,
                        FILE *list, int file_version)
{
    int                  i;
    gmx_off_t            offset;
    gmx_off_t            mask = 0xFFFFFFFFL;
    int                  offset_high, offset_low;
    char                *buf;
    gmx_file_position_t *outputfiles;

    if (do_cpt_int(xd, "number of output files", nfiles, list) != 0)
    {
        return -1;
    }

    if (bRead)
    {
        snew(*p_outputfiles, *nfiles);
    }

    outputfiles = *p_outputfiles;

    for (i = 0; i < *nfiles; i++)
    {
        /* 64-bit XDR numbers are not portable, so it is stored as separate high/low fractions */
        if (bRead)
        {
            do_cpt_string_err(xd, bRead, "output filename", &buf, list);
            std::strncpy(outputfiles[i].filename, buf, CPTSTRLEN-1);
            if (list == nullptr)
            {
                sfree(buf);
            }

            if (do_cpt_int(xd, "file_offset_high", &offset_high, list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(xd, "file_offset_low", &offset_low, list) != 0)
            {
                return -1;
            }
            outputfiles[i].offset = (static_cast<gmx_off_t>(offset_high) << 32 ) | ( static_cast<gmx_off_t>(offset_low) & mask );
        }
        else
        {
            buf = outputfiles[i].filename;
            do_cpt_string_err(xd, bRead, "output filename", &buf, list);
            /* writing */
            offset      = outputfiles[i].offset;
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
            if (do_cpt_int(xd, "file_checksum_size", &(outputfiles[i].chksum_size),
                           list) != 0)
            {
                return -1;
            }
            if (do_cpt_u_chars(xd, "file_checksum", 16, outputfiles[i].chksum, list) != 0)
            {
                return -1;
            }
        }
        else
        {
            outputfiles[i].chksum_size = -1;
        }
    }
    return 0;
}


void write_checkpoint(const char *fn, gmx_bool bNumberAndKeep,
                      FILE *fplog, t_commrec *cr,
                      ivec domdecCells, int nppnodes,
                      int eIntegrator, int simulation_part,
                      gmx_bool bExpanded, int elamstats,
                      gmx_int64_t step, double t,
                      t_state *state, energyhistory_t *enerhist)
{
    t_fileio            *fp;
    int                  file_version;
    char                *version;
    char                *btime;
    char                *buser;
    char                *bhost;
    int                  double_prec;
    char                *fprog;
    char                *fntemp; /* the temporary checkpoint file name */
    char                 timebuf[STRLEN];
    int                  npmenodes;
    char                 buf[1024], suffix[5+STEPSTRSIZE], sbuf[STEPSTRSIZE];
    gmx_file_position_t *outputfiles;
    int                  noutputfiles;
    char                *ftime;
    int                  flags_eks, flags_enh, flags_dfh;
    t_fileio            *ret;

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
    snew(fntemp, std::strlen(fn)+5+STEPSTRSIZE);
    std::strcpy(fntemp, fn);
    fntemp[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
    sprintf(suffix, "_%s%s", "step", gmx_step_str(step, sbuf));
    std::strcat(fntemp, suffix);
    std::strcat(fntemp, fn+std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
#else
    /* if we can't rename, we just overwrite the cpt file.
     * dangerous if interrupted.
     */
    snew(fntemp, std::strlen(fn));
    std::strcpy(fntemp, fn);
#endif
    gmx_format_current_time(timebuf, STRLEN);

    if (fplog)
    {
        fprintf(fplog, "Writing checkpoint, step %s at %s\n\n",
                gmx_step_str(step, buf), timebuf);
    }

    /* Get offsets for open files */
    gmx_fio_get_output_file_positions(&outputfiles, &noutputfiles);

    fp = gmx_fio_open(fntemp, "w");

    if (state->ekinstate.bUpToDate)
    {
        flags_eks =
            ((1<<eeksEKIN_N) | (1<<eeksEKINH) | (1<<eeksEKINF) |
             (1<<eeksEKINO) | (1<<eeksEKINSCALEF) | (1<<eeksEKINSCALEH) |
             (1<<eeksVSCALE) | (1<<eeksDEKINDL) | (1<<eeksMVCOS));
    }
    else
    {
        flags_eks = 0;
    }

    flags_enh = 0;
    if (enerhist->nsum > 0 || enerhist->nsum_sim > 0)
    {
        flags_enh |= (1<<eenhENERGY_N) | (1<<eenhENERGY_NSTEPS) | (1<<eenhENERGY_NSTEPS_SIM);
        if (enerhist->nsum > 0)
        {
            flags_enh |= ((1<<eenhENERGY_AVER) | (1<<eenhENERGY_SUM) |
                          (1<<eenhENERGY_NSUM));
        }
        if (enerhist->nsum_sim > 0)
        {
            flags_enh |= ((1<<eenhENERGY_SUM_SIM) | (1<<eenhENERGY_NSUM_SIM));
        }
        if (enerhist->deltaHForeignLambdas != nullptr)
        {
            flags_enh |= ( (1<< eenhENERGY_DELTA_H_NN) |
                           (1<< eenhENERGY_DELTA_H_LIST) |
                           (1<< eenhENERGY_DELTA_H_STARTTIME) |
                           (1<< eenhENERGY_DELTA_H_STARTLAMBDA) );
        }
    }

    if (bExpanded)
    {
        flags_dfh = ((1<<edfhBEQUIL) | (1<<edfhNATLAMBDA) | (1<<edfhSUMWEIGHTS) |  (1<<edfhSUMDG)  |
                     (1<<edfhTIJ) | (1<<edfhTIJEMP));
        if (EWL(elamstats))
        {
            flags_dfh |= ((1<<edfhWLDELTA) | (1<<edfhWLHISTO));
        }
        if ((elamstats == elamstatsMINVAR) || (elamstats == elamstatsBARKER) || (elamstats == elamstatsMETROPOLIS))
        {
            flags_dfh |= ((1<<edfhACCUMP) | (1<<edfhACCUMM) | (1<<edfhACCUMP2) | (1<<edfhACCUMM2)
                          | (1<<edfhSUMMINVAR) | (1<<edfhSUMVAR));
        }
    }
    else
    {
        flags_dfh = 0;
    }

    /* We can check many more things now (CPU, acceleration, etc), but
     * it is highly unlikely to have two separate builds with exactly
     * the same version, user, time, and build host!
     */

    version = gmx_strdup(gmx_version());
    btime   = gmx_strdup(BUILD_TIME);
    buser   = gmx_strdup(BUILD_USER);
    bhost   = gmx_strdup(BUILD_HOST);

    double_prec = GMX_DOUBLE;
    fprog       = gmx_strdup(gmx::getProgramContext().fullBinaryPath());

    ftime   = &(timebuf[0]);

    int nlambda     = (state->dfhist ? state->dfhist->nlambda : 0);
    int nED         = (state->edsamstate ? state->edsamstate->nED : 0);
    int eSwapCoords = (state->swapstate ? state->swapstate->eSwapCoords : eswapNO);

    do_cpt_header(gmx_fio_getxdr(fp), FALSE, &file_version,
                  &version, &btime, &buser, &bhost, &double_prec, &fprog, &ftime,
                  &eIntegrator, &simulation_part, &step, &t, &nppnodes,
                  DOMAINDECOMP(cr) ? domdecCells : nullptr, &npmenodes,
                  &state->natoms, &state->ngtc, &state->nnhpres,
                  &state->nhchainlength, &nlambda, &state->flags, &flags_eks, &flags_enh, &flags_dfh,
                  &nED, &eSwapCoords,
                  nullptr);

    sfree(version);
    sfree(btime);
    sfree(buser);
    sfree(bhost);
    sfree(fprog);

    if ((do_cpt_state(gmx_fio_getxdr(fp), state->flags, state, nullptr) < 0)        ||
        (do_cpt_ekinstate(gmx_fio_getxdr(fp), flags_eks, &state->ekinstate, nullptr) < 0) ||
        (do_cpt_enerhist(gmx_fio_getxdr(fp), FALSE, flags_enh, enerhist, nullptr) < 0)  ||
        (do_cpt_df_hist(gmx_fio_getxdr(fp), flags_dfh, nlambda, &state->dfhist, nullptr) < 0)  ||
        (do_cpt_EDstate(gmx_fio_getxdr(fp), FALSE, nED, &state->edsamstate, nullptr) < 0)      ||
        (do_cpt_swapstate(gmx_fio_getxdr(fp), FALSE, eSwapCoords, &state->swapstate, nullptr) < 0) ||
        (do_cpt_files(gmx_fio_getxdr(fp), FALSE, &outputfiles, &noutputfiles, nullptr,
                      file_version) < 0))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    do_cpt_footer(gmx_fio_getxdr(fp), file_version);

    /* we really, REALLY, want to make sure to physically write the checkpoint,
       and all the files it depends on, out to disk. Because we've
       opened the checkpoint with gmx_fio_open(), it's in our list
       of open files.  */
    ret = gmx_fio_all_output_fsync();

    if (ret)
    {
        char buf[STRLEN];
        sprintf(buf,
                "Cannot fsync '%s'; maybe you are out of disk space?",
                gmx_fio_getname(ret));

        if (getenv(GMX_IGNORE_FSYNC_FAILURE_ENV) == nullptr)
        {
            gmx_file(buf);
        }
        else
        {
            gmx_warning(buf);
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
            std::strcpy(buf, fn);
            buf[std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
            std::strcat(buf, "_prev");
            std::strcat(buf, fn+std::strlen(fn) - std::strlen(ftp2ext(fn2ftp(fn))) - 1);
#ifndef GMX_FAHCORE
            /* we copy here so that if something goes wrong between now and
             * the rename below, there's always a state.cpt.
             * If renames are atomic (such as in POSIX systems),
             * this copying should be unneccesary.
             */
            gmx_file_copy(fn, buf, FALSE);
            /* We don't really care if this fails:
             * there's already a new checkpoint.
             */
#else
            gmx_file_rename(fn, buf);
#endif
        }
        if (gmx_file_rename(fntemp, fn) != 0)
        {
            gmx_file("Cannot rename checkpoint file; maybe you are out of disk space?");
        }
    }
#endif  /* GMX_NO_RENAME */

    sfree(outputfiles);
    sfree(fntemp);

#ifdef GMX_FAHCORE
    /*code for alternate checkpointing scheme.  moved from top of loop over
       steps */
    fcRequestCheckPoint();
    if (fcCheckPointParallel( cr->nodeid, NULL, 0) == 0)
    {
        gmx_fatal( 3, __FILE__, __LINE__, "Checkpoint error on step %d\n", step );
    }
#endif /* end GMX_FAHCORE block */
}

static void print_flag_mismatch(FILE *fplog, int sflags, int fflags)
{
    int i;

    fprintf(fplog, "\nState entry mismatch between the simulation and the checkpoint file\n");
    fprintf(fplog, "Entries which are not present in the checkpoint file will not be updated\n");
    fprintf(fplog, "  %24s    %11s    %11s\n", "", "simulation", "checkpoint");
    for (i = 0; i < estNR; i++)
    {
        if ((sflags & (1<<i)) || (fflags & (1<<i)))
        {
            fprintf(fplog, "  %24s    %11s    %11s\n",
                    est_names[i],
                    (sflags & (1<<i)) ? "  present  " : "not present",
                    (fflags & (1<<i)) ? "  present  " : "not present");
        }
    }
}

static void check_int(FILE *fplog, const char *type, int p, int f, gmx_bool *mm)
{
    FILE *fp = fplog ? fplog : stderr;

    if (p != f)
    {
        fprintf(fp, "  %s mismatch,\n", type);
        fprintf(fp, "    current program: %d\n", p);
        fprintf(fp, "    checkpoint file: %d\n", f);
        fprintf(fp, "\n");
        *mm = TRUE;
    }
}

static void check_string(FILE *fplog, const char *type, const char *p,
                         const char *f, gmx_bool *mm)
{
    FILE *fp = fplog ? fplog : stderr;

    if (std::strcmp(p, f) != 0)
    {
        fprintf(fp, "  %s mismatch,\n", type);
        fprintf(fp, "    current program: %s\n", p);
        fprintf(fp, "    checkpoint file: %s\n", f);
        fprintf(fp, "\n");
        *mm = TRUE;
    }
}

static void check_match(FILE *fplog,
                        char *version,
                        char *btime, char *buser, char *bhost, int double_prec,
                        char *fprog,
                        const t_commrec *cr, int npp_f, int npme_f,
                        ivec dd_nc, ivec dd_nc_f,
                        gmx_bool reproducibilityRequested)
{
    /* Note that this check_string on the version will also print a message
     * when only the minor version differs. But we only print a warning
     * message further down with reproducibilityRequested=TRUE.
     */
    gmx_bool versionDiffers = FALSE;
    check_string(fplog, "Version", gmx_version(), version, &versionDiffers);

    gmx_bool precisionDiffers = FALSE;
    check_int   (fplog, "Double prec.", GMX_DOUBLE, double_prec, &precisionDiffers);
    if (precisionDiffers)
    {
        const char msg_precision_difference[] =
            "You are continuing a simulation with a different precision. Not matching\n"
            "single/double precision will lead to precision or performance loss.\n";
        fprintf(stderr, "%s\n", msg_precision_difference);
        if (fplog)
        {
            fprintf(fplog, "%s\n", msg_precision_difference);
        }
    }

    gmx_bool mm = (versionDiffers || precisionDiffers);

    if (reproducibilityRequested)
    {
        check_string(fplog, "Build time", BUILD_TIME, btime, &mm);
        check_string(fplog, "Build user", BUILD_USER, buser, &mm);
        check_string(fplog, "Build host", BUILD_HOST, bhost, &mm);
        check_string(fplog, "Program name", gmx::getProgramContext().fullBinaryPath(), fprog, &mm);

        check_int   (fplog, "#ranks", cr->nnodes, npp_f+npme_f, &mm);
    }

    if (cr->nnodes > 1 && reproducibilityRequested)
    {
        check_int (fplog, "#PME-ranks", cr->npmenodes, npme_f, &mm);

        int npp = cr->nnodes;
        if (cr->npmenodes >= 0)
        {
            npp -= cr->npmenodes;
        }
        if (npp == npp_f)
        {
            check_int (fplog, "#DD-cells[x]", dd_nc[XX], dd_nc_f[XX], &mm);
            check_int (fplog, "#DD-cells[y]", dd_nc[YY], dd_nc_f[YY], &mm);
            check_int (fplog, "#DD-cells[z]", dd_nc[ZZ], dd_nc_f[ZZ], &mm);
        }
    }

    if (mm)
    {
        /* Gromacs should be able to continue from checkpoints between
         * different patch level versions, but we do not guarantee
         * compatibility between different major/minor versions - check this.
         */
        int        gmx_major;
        int        cpt_major;
        sscanf(gmx_version(), "%5d", &gmx_major);
        int        ret                 = sscanf(version, "%5d", &cpt_major);
        gmx_bool   majorVersionDiffers = (ret < 1 || gmx_major != cpt_major);

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

        const char msg_logdetails[] =
            "See the log file for details.\n";

        if (majorVersionDiffers)
        {
            fprintf(stderr, "%s%s\n", msg_major_version_difference, fplog ? msg_logdetails : "");

            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_major_version_difference);
            }
        }
        else if (reproducibilityRequested)
        {
            /* Major & minor versions match at least, but something is different. */
            fprintf(stderr, "%s%s\n", msg_mismatch_notice, fplog ? msg_logdetails : "");
            if (fplog)
            {
                fprintf(fplog, "%s\n", msg_mismatch_notice);
            }
        }
    }
}

static void read_checkpoint(const char *fn, FILE **pfplog,
                            const t_commrec *cr,
                            ivec dd_nc, int *npme,
                            int eIntegrator, int *init_fep_state, gmx_int64_t *step, double *t,
                            t_state *state, gmx_bool *bReadEkin,
                            energyhistory_t *enerhist,
                            int *simulation_part,
                            gmx_bool bAppendOutputFiles, gmx_bool bForceAppend,
                            gmx_bool reproducibilityRequested)
{
    t_fileio            *fp;
    int                  i, j, rc;
    int                  file_version;
    char                *version, *btime, *buser, *bhost, *fprog, *ftime;
    int                  double_prec;
    char                 buf[STEPSTRSIZE];
    int                  eIntegrator_f, nppnodes_f, npmenodes_f;
    ivec                 dd_nc_f;
    int                  natoms, ngtc, nnhpres, nhchainlength, nlambda, fflags, flags_eks, flags_enh, flags_dfh;
    int                  nED, eSwapCoords;
    int                  d;
    int                  ret;
    gmx_file_position_t *outputfiles;
    int                  nfiles;
    t_fileio            *chksum_file;
    FILE               * fplog = *pfplog;
    unsigned char        digest[16];
#if !defined __native_client__ && !GMX_NATIVE_WINDOWS
    struct flock         fl; /* don't initialize here: the struct order is OS
                                dependent! */
#endif

    const char *int_warn =
        "WARNING: The checkpoint file was generated with integrator %s,\n"
        "         while the simulation uses integrator %s\n\n";

#if !defined __native_client__ && !GMX_NATIVE_WINDOWS
    fl.l_type   = F_WRLCK;
    fl.l_whence = SEEK_SET;
    fl.l_start  = 0;
    fl.l_len    = 0;
    fl.l_pid    = 0;
#endif

    fp = gmx_fio_open(fn, "r");
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, &file_version,
                  &version, &btime, &buser, &bhost, &double_prec, &fprog, &ftime,
                  &eIntegrator_f, simulation_part, step, t,
                  &nppnodes_f, dd_nc_f, &npmenodes_f,
                  &natoms, &ngtc, &nnhpres, &nhchainlength, &nlambda,
                  &fflags, &flags_eks, &flags_enh, &flags_dfh,
                  &nED, &eSwapCoords, nullptr);

    if (bAppendOutputFiles &&
        file_version >= 13 && double_prec != GMX_DOUBLE)
    {
        gmx_fatal(FARGS, "Output file appending requested, but the code and checkpoint file precision (single/double) don't match");
    }

    if (cr == nullptr || MASTER(cr))
    {
        fprintf(stderr, "\nReading checkpoint file %s generated: %s\n\n",
                fn, ftime);
    }

    /* This will not be written if we do appending, since fplog is still NULL then */
    if (fplog)
    {
        fprintf(fplog, "\n");
        fprintf(fplog, "Reading checkpoint file %s\n", fn);
        fprintf(fplog, "  file generated by:     %s\n", fprog);
        fprintf(fplog, "  file generated at:     %s\n", ftime);
        fprintf(fplog, "  GROMACS build time:    %s\n", btime);
        fprintf(fplog, "  GROMACS build user:    %s\n", buser);
        fprintf(fplog, "  GROMACS build host:    %s\n", bhost);
        fprintf(fplog, "  GROMACS double prec.:  %d\n", double_prec);
        fprintf(fplog, "  simulation part #:     %d\n", *simulation_part);
        fprintf(fplog, "  step:                  %s\n", gmx_step_str(*step, buf));
        fprintf(fplog, "  time:                  %f\n", *t);
        fprintf(fplog, "\n");
    }

    if (natoms != state->natoms)
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d atoms, while the current system consists of %d atoms", natoms, state->natoms);
    }
    if (ngtc != state->ngtc)
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d T-coupling groups, while the current system consists of %d T-coupling groups", ngtc, state->ngtc);
    }
    if (nnhpres != state->nnhpres)
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system of %d NH-pressure-coupling variables, while the current system consists of %d NH-pressure-coupling variables", nnhpres, state->nnhpres);
    }

    int nlambdaHistory = (state->dfhist ? state->dfhist->nlambda : 0);
    if (nlambda != nlambdaHistory)
    {
        gmx_fatal(FARGS, "Checkpoint file is for a system with %d lambda states, while the current system consists of %d lambda states", nlambda, nlambdaHistory);
    }

    init_gtc_state(state, state->ngtc, state->nnhpres, nhchainlength); /* need to keep this here to keep the tpr format working */
    /* write over whatever was read; we use the number of Nose-Hoover chains from the checkpoint */

    if (eIntegrator_f != eIntegrator)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, int_warn, EI(eIntegrator_f), EI(eIntegrator));
        }
        if (bAppendOutputFiles)
        {
            gmx_fatal(FARGS,
                      "Output file appending requested, but input/checkpoint integrators do not match.\n"
                      "Stopping the run to prevent you from ruining all your data...\n"
                      "If you _really_ know what you are doing, try with the -noappend option.\n");
        }
        if (fplog)
        {
            fprintf(fplog, int_warn, EI(eIntegrator_f), EI(eIntegrator));
        }
    }

    if (!PAR(cr))
    {
        *npme = 0;
    }
    else if (cr->nnodes == nppnodes_f + npmenodes_f)
    {
        if (*npme < 0)
        {
            *npme = npmenodes_f;
        }
        int nppnodes = cr->nnodes - *npme;
        if (nppnodes == nppnodes_f)
        {
            for (d = 0; d < DIM; d++)
            {
                if (dd_nc[d] == 0)
                {
                    dd_nc[d] = dd_nc_f[d];
                }
            }
        }
    }

    if (fflags != state->flags)
    {

        if (MASTER(cr))
        {
            if (bAppendOutputFiles)
            {
                gmx_fatal(FARGS,
                          "Output file appending requested, but input and checkpoint states are not identical.\n"
                          "Stopping the run to prevent you from ruining all your data...\n"
                          "You can try with the -noappend option, and get more info in the log file.\n");
            }

            if (getenv("GMX_ALLOW_CPT_MISMATCH") == nullptr)
            {
                gmx_fatal(FARGS, "You seem to have switched ensemble, integrator, T and/or P-coupling algorithm between the cpt and tpr file. The recommended way of doing this is passing the cpt file to grompp (with option -t) instead of to mdrun. If you know what you are doing, you can override this error by setting the env.var. GMX_ALLOW_CPT_MISMATCH");
            }
            else
            {
                fprintf(stderr,
                        "WARNING: The checkpoint state entries do not match the simulation,\n"
                        "         see the log file for details\n\n");
            }
        }

        if (fplog)
        {
            print_flag_mismatch(fplog, state->flags, fflags);
        }
    }
    else
    {
        if (MASTER(cr))
        {
            check_match(fplog, version, btime, buser, bhost, double_prec, fprog,
                        cr, nppnodes_f, npmenodes_f, dd_nc, dd_nc_f,
                        reproducibilityRequested);
        }
    }
    ret             = do_cpt_state(gmx_fio_getxdr(fp), fflags, state, nullptr);
    *init_fep_state = state->fep_state;  /* there should be a better way to do this than setting it here.
                                            Investigate for 5.0. */
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }
    *bReadEkin = ((flags_eks & (1<<eeksEKINH)) || (flags_eks & (1<<eeksEKINF)) || (flags_eks & (1<<eeksEKINO)) ||
                  ((flags_eks & (1<<eeksEKINSCALEF)) | (flags_eks & (1<<eeksEKINSCALEH)) | (flags_eks & (1<<eeksVSCALE))));

    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE,
                          flags_enh, enerhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    if (file_version < 6)
    {
        const char *warn = "Reading checkpoint file in old format, assuming that the run that generated this file started at step 0, if this is not the case the averages stored in the energy file will be incorrect.";

        fprintf(stderr, "\nWARNING: %s\n\n", warn);
        if (fplog)
        {
            fprintf(fplog, "\nWARNING: %s\n\n", warn);
        }
        enerhist->nsum     = *step;
        enerhist->nsum_sim = *step;
    }

    ret = do_cpt_df_hist(gmx_fio_getxdr(fp), flags_dfh, nlambda, &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, nED, &state->edsamstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, eSwapCoords, &state->swapstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_files(gmx_fio_getxdr(fp), TRUE, &outputfiles, &nfiles, nullptr, file_version);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_footer(gmx_fio_getxdr(fp), file_version);
    if (ret)
    {
        cp_error();
    }
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }

    sfree(fprog);
    sfree(ftime);
    sfree(btime);
    sfree(buser);
    sfree(bhost);

    /* If the user wants to append to output files,
     * we use the file pointer positions of the output files stored
     * in the checkpoint file and truncate the files such that any frames
     * written after the checkpoint time are removed.
     * All files are md5sum checked such that we can be sure that
     * we do not truncate other (maybe imprortant) files.
     */
    if (bAppendOutputFiles)
    {
        if (fn2ftp(outputfiles[0].filename) != efLOG)
        {
            /* make sure first file is log file so that it is OK to use it for
             * locking
             */
            gmx_fatal(FARGS, "The first output file should always be the log "
                      "file but instead is: %s. Cannot do appending because of this condition.", outputfiles[0].filename);
        }
        for (i = 0; i < nfiles; i++)
        {
            if (outputfiles[i].offset < 0)
            {
                gmx_fatal(FARGS, "The original run wrote a file called '%s' which "
                          "is larger than 2 GB, but mdrun did not support large file"
                          " offsets. Can not append. Run mdrun with -noappend",
                          outputfiles[i].filename);
            }
#ifdef GMX_FAHCORE
            chksum_file = gmx_fio_open(outputfiles[i].filename, "a");

#else
            chksum_file = gmx_fio_open(outputfiles[i].filename, "r+");

            /* lock log file */
            if (i == 0)
            {
                /* Note that there are systems where the lock operation
                 * will succeed, but a second process can also lock the file.
                 * We should probably try to detect this.
                 */
#if defined __native_client__
                errno = ENOSYS;
                if (1)

#elif GMX_NATIVE_WINDOWS
                if (_locking(fileno(gmx_fio_getfp(chksum_file)), _LK_NBLCK, LONG_MAX) == -1)
#else
                if (fcntl(fileno(gmx_fio_getfp(chksum_file)), F_SETLK, &fl) == -1)
#endif
                {
                    if (errno == ENOSYS)
                    {
                        if (!bForceAppend)
                        {
                            gmx_fatal(FARGS, "File locking is not supported on this system. Use -noappend or specify -append explicitly to append anyhow.");
                        }
                        else
                        {
                            fprintf(stderr, "\nNOTE: File locking is not supported on this system, will not lock %s\n\n", outputfiles[i].filename);
                            if (fplog)
                            {
                                fprintf(fplog, "\nNOTE: File locking not supported on this system, will not lock %s\n\n", outputfiles[i].filename);
                            }
                        }
                    }
                    else if (errno == EACCES || errno == EAGAIN)
                    {
                        gmx_fatal(FARGS, "Failed to lock: %s. Already running "
                                  "simulation?", outputfiles[i].filename);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Failed to lock: %s. %s.",
                                  outputfiles[i].filename, std::strerror(errno));
                    }
                }
            }

            /* compute md5 chksum */
            if (outputfiles[i].chksum_size != -1)
            {
                if (gmx_fio_get_file_md5(chksum_file, outputfiles[i].offset,
                                         digest) != outputfiles[i].chksum_size) /*at the end of the call the file position is at the end of the file*/
                {
                    gmx_fatal(FARGS, "Can't read %d bytes of '%s' to compute checksum. The file has been replaced or its contents have been modified. Cannot do appending because of this condition.",
                              outputfiles[i].chksum_size,
                              outputfiles[i].filename);
                }
            }
            if (i == 0)  /*log file needs to be seeked in case we need to truncate (other files are truncated below)*/
            {
                if (gmx_fio_seek(chksum_file, outputfiles[i].offset))
                {
                    gmx_fatal(FARGS, "Seek error! Failed to truncate log-file: %s.", std::strerror(errno));
                }
            }
#endif

            if (i == 0) /*open log file here - so that lock is never lifted
                           after chksum is calculated */
            {
                *pfplog = gmx_fio_getfp(chksum_file);
            }
            else
            {
                gmx_fio_close(chksum_file);
            }
#ifndef GMX_FAHCORE
            /* compare md5 chksum */
            if (outputfiles[i].chksum_size != -1 &&
                memcmp(digest, outputfiles[i].chksum, 16) != 0)
            {
                if (debug)
                {
                    fprintf(debug, "chksum for %s: ", outputfiles[i].filename);
                    for (j = 0; j < 16; j++)
                    {
                        fprintf(debug, "%02x", digest[j]);
                    }
                    fprintf(debug, "\n");
                }
                gmx_fatal(FARGS, "Checksum wrong for '%s'. The file has been replaced or its contents have been modified. Cannot do appending because of this condition.",
                          outputfiles[i].filename);
            }
#endif


            if (i != 0) /*log file is already seeked to correct position */
            {
#if !GMX_NATIVE_WINDOWS || !defined(GMX_FAHCORE)
                /* For FAHCORE, we do this elsewhere*/
                rc = gmx_truncate(outputfiles[i].filename, outputfiles[i].offset);
                if (rc != 0)
                {
                    gmx_fatal(FARGS, "Truncation of file %s failed. Cannot do appending because of this failure.", outputfiles[i].filename);
                }
#endif
            }
        }
    }

    sfree(outputfiles);
}


void load_checkpoint(const char *fn, FILE **fplog,
                     const t_commrec *cr, ivec dd_nc, int *npme,
                     t_inputrec *ir, t_state *state,
                     gmx_bool *bReadEkin,
                     energyhistory_t *enerhist,
                     gmx_bool bAppend, gmx_bool bForceAppend,
                     gmx_bool reproducibilityRequested)
{
    gmx_int64_t     step;
    double          t;

    if (SIMMASTER(cr))
    {
        /* Read the state from the checkpoint file */
        read_checkpoint(fn, fplog,
                        cr, dd_nc, npme,
                        ir->eI, &(ir->fepvals->init_fep_state), &step, &t,
                        state, bReadEkin, enerhist,
                        &ir->simulation_part, bAppend, bForceAppend,
                        reproducibilityRequested);
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(*npme), npme, cr);
        gmx_bcast(DIM*sizeof(dd_nc[0]), dd_nc, cr);
        gmx_bcast(sizeof(step), &step, cr);
        gmx_bcast(sizeof(*bReadEkin), bReadEkin, cr);
    }
    ir->bContinuation    = TRUE;
    if (ir->nsteps >= 0)
    {
        ir->nsteps          += ir->init_step - step;
    }
    ir->init_step        = step;
    ir->simulation_part += 1;
}

void read_checkpoint_part_and_step(const char  *filename,
                                   int         *simulation_part,
                                   gmx_int64_t *step)
{
    int       file_version;
    char     *version, *btime, *buser, *bhost, *fprog, *ftime;
    int       double_prec;
    int       eIntegrator;
    int       nppnodes, npme;
    ivec      dd_nc;
    int       nlambda;
    int       flags_eks, flags_enh, flags_dfh;
    double    t;
    t_state   state;
    int       nED, eSwapCoords;
    t_fileio *fp;

    if (filename == nullptr ||
        !gmx_fexist(filename) ||
        (!(fp = gmx_fio_open(filename, "r"))))
    {
        *simulation_part = 0;
        *step            = 0;
        return;
    }

    /* Not calling initializing state before use is nasty, but all we
       do is read into its member variables and throw the struct away
       again immediately. */

    do_cpt_header(gmx_fio_getxdr(fp), TRUE, &file_version,
                  &version, &btime, &buser, &bhost, &double_prec, &fprog, &ftime,
                  &eIntegrator, simulation_part, step, &t, &nppnodes, dd_nc, &npme,
                  &state.natoms, &state.ngtc, &state.nnhpres, &state.nhchainlength,
                  &nlambda, &state.flags, &flags_eks, &flags_enh, &flags_dfh,
                  &nED, &eSwapCoords, nullptr);

    gmx_fio_close(fp);
}

static void read_checkpoint_data(t_fileio *fp, int *simulation_part,
                                 gmx_int64_t *step, double *t, t_state *state,
                                 int *nfiles, gmx_file_position_t **outputfiles)
{
    int                  file_version;
    char                *version, *btime, *buser, *bhost, *fprog, *ftime;
    int                  double_prec;
    int                  eIntegrator;
    int                  nppnodes, npme;
    ivec                 dd_nc;
    int                  nlambda;
    int                  flags_eks, flags_enh, flags_dfh;
    int                  nED, eSwapCoords;
    int                  nfiles_loc;
    gmx_file_position_t *files_loc = nullptr;
    int                  ret;

    do_cpt_header(gmx_fio_getxdr(fp), TRUE, &file_version,
                  &version, &btime, &buser, &bhost, &double_prec, &fprog, &ftime,
                  &eIntegrator, simulation_part, step, t, &nppnodes, dd_nc, &npme,
                  &state->natoms, &state->ngtc, &state->nnhpres, &state->nhchainlength,
                  &nlambda, &state->flags, &flags_eks, &flags_enh, &flags_dfh,
                  &nED, &eSwapCoords, nullptr);
    ret =
        do_cpt_state(gmx_fio_getxdr(fp), state->flags, state, nullptr);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), flags_eks, &state->ekinstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE,
                          flags_enh, &enerhist, nullptr);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_df_hist(gmx_fio_getxdr(fp), flags_dfh, nlambda, &state->dfhist, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, nED, &state->edsamstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, eSwapCoords, &state->swapstate, nullptr);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_files(gmx_fio_getxdr(fp), TRUE,
                       &files_loc,
                       &nfiles_loc,
                       nullptr, file_version);
    if (outputfiles != nullptr)
    {
        *outputfiles = files_loc;
    }
    else
    {
        sfree(files_loc);
    }
    if (nfiles != nullptr)
    {
        *nfiles = nfiles_loc;
    }

    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_footer(gmx_fio_getxdr(fp), file_version);
    if (ret)
    {
        cp_error();
    }

    sfree(fprog);
    sfree(ftime);
    sfree(btime);
    sfree(buser);
    sfree(bhost);
}

void
read_checkpoint_state(const char *fn, int *simulation_part,
                      gmx_int64_t *step, double *t, t_state *state)
{
    t_fileio *fp;

    fp = gmx_fio_open(fn, "r");
    read_checkpoint_data(fp, simulation_part, step, t, state, nullptr, nullptr);
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}

void read_checkpoint_trxframe(t_fileio *fp, t_trxframe *fr)
{
    t_state         state;
    int             simulation_part;
    gmx_int64_t     step;
    double          t;

    read_checkpoint_data(fp, &simulation_part, &step, &t, &state, nullptr, nullptr);

    fr->natoms  = state.natoms;
    fr->bTitle  = FALSE;
    fr->bStep   = TRUE;
    fr->step    = gmx_int64_to_int(step,
                                   "conversion of checkpoint to trajectory");
    fr->bTime      = TRUE;
    fr->time       = t;
    fr->bLambda    = TRUE;
    fr->lambda     = state.lambda[efptFEP];
    fr->fep_state  = state.fep_state;
    fr->bAtoms     = FALSE;
    fr->bX         = (state.flags & (1<<estX));
    if (fr->bX)
    {
        fr->x   = getRvecArrayFromPaddedRVecVector(&state.x, state.natoms);
    }
    fr->bV      = (state.flags & (1<<estV));
    if (fr->bV)
    {
        fr->v   = getRvecArrayFromPaddedRVecVector(&state.v, state.natoms);
    }
    fr->bF      = FALSE;
    fr->bBox    = (state.flags & (1<<estBOX));
    if (fr->bBox)
    {
        copy_mat(state.box, fr->box);
    }
}

void list_checkpoint(const char *fn, FILE *out)
{
    t_fileio            *fp;
    int                  file_version;
    char                *version, *btime, *buser, *bhost, *fprog, *ftime;
    int                  double_prec;
    int                  eIntegrator, simulation_part, nppnodes, npme;
    gmx_int64_t          step;
    double               t;
    ivec                 dd_nc;
    int                  nlambda;
    int                  flags_eks, flags_enh, flags_dfh;
    int                  nED, eSwapCoords;
    int                  ret;
    gmx_file_position_t *outputfiles;
    int                  nfiles;

    t_state              state;

    fp = gmx_fio_open(fn, "r");
    do_cpt_header(gmx_fio_getxdr(fp), TRUE, &file_version,
                  &version, &btime, &buser, &bhost, &double_prec, &fprog, &ftime,
                  &eIntegrator, &simulation_part, &step, &t, &nppnodes, dd_nc, &npme,
                  &state.natoms, &state.ngtc, &state.nnhpres, &state.nhchainlength,
                  &nlambda, &state.flags,
                  &flags_eks, &flags_enh, &flags_dfh, &nED, &eSwapCoords,
                  out);
    ret = do_cpt_state(gmx_fio_getxdr(fp), state.flags, &state, out);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp), flags_eks, &state.ekinstate, out);
    if (ret)
    {
        cp_error();
    }

    energyhistory_t enerhist;
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp), TRUE,
                          flags_enh, &enerhist, out);

    if (ret == 0)
    {
        ret = do_cpt_df_hist(gmx_fio_getxdr(fp),
                             flags_dfh, nlambda, &state.dfhist, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_EDstate(gmx_fio_getxdr(fp), TRUE, nED, &state.edsamstate, out);
    }

    if (ret == 0)
    {
        ret = do_cpt_swapstate(gmx_fio_getxdr(fp), TRUE, eSwapCoords, &state.swapstate, out);
    }

    if (ret == 0)
    {
        do_cpt_files(gmx_fio_getxdr(fp), TRUE, &outputfiles, &nfiles, out, file_version);
    }

    if (ret == 0)
    {
        ret = do_cpt_footer(gmx_fio_getxdr(fp), file_version);
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
void
read_checkpoint_simulation_part_and_filenames(t_fileio             *fp,
                                              int                  *simulation_part,
                                              int                  *nfiles,
                                              gmx_file_position_t **outputfiles)
{
    gmx_int64_t step = 0;
    double      t;
    t_state     state;

    read_checkpoint_data(fp, simulation_part, &step, &t, &state,
                         nfiles, outputfiles);
    if (gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of disk space?");
    }
}
