/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#include "gmxpre.h"

#include "update.h"

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <memory>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

struct gmx_sd_const_t
{
    double em = 0;
};

struct gmx_sd_sigma_t
{
    real V = 0;
};

struct gmx_stochd_t
{
    /* BD stuff */
    std::vector<real> bd_rf;
    /* SD stuff */
    std::vector<gmx_sd_const_t> sdc;
    std::vector<gmx_sd_sigma_t> sdsig;
    /* andersen temperature control stuff */
    std::vector<bool> randomize_group;
    std::vector<real> boltzfac;

    explicit gmx_stochd_t(const t_inputrec* ir);
};

//! pImpled implementation for Update
class Update::Impl
{
public:
    //! Constructor
    Impl(const t_inputrec* ir, BoxDeformation* boxDeformation);
    //! Destructor
    ~Impl() = default;
    //! stochastic dynamics struct
    std::unique_ptr<gmx_stochd_t> sd;
    //! xprime for constraint algorithms
    PaddedVector<RVec> xp;
    //! Box deformation handler (or nullptr if inactive).
    BoxDeformation* deform = nullptr;
};

Update::Update(const t_inputrec* ir, BoxDeformation* boxDeformation) :
    impl_(new Impl(ir, boxDeformation)){};

Update::~Update(){};

gmx_stochd_t* Update::sd() const
{
    return impl_->sd.get();
}

PaddedVector<RVec>* Update::xp()
{
    return &impl_->xp;
}

BoxDeformation* Update::deform() const
{
    return impl_->deform;
}

/*! \brief Sets the velocities of virtual sites to zero */
static void clearVsiteVelocities(int start, int nrend, const unsigned short* particleType, rvec* gmx_restrict v)
{
    for (int a = start; a < nrend; a++)
    {
        if (particleType[a] == eptVSite)
        {
            clear_rvec(v[a]);
        }
    }
}

/*! \brief Sets the number of different temperature coupling values */
enum class NumTempScaleValues
{
    single,  //!< Single T-scaling value (either one group or all values =1)
    multiple //!< Multiple T-scaling values, need to use T-group indices
};

/*! \brief Sets if to apply no or diagonal Parrinello-Rahman pressure scaling
 *
 * Note that this enum is only used in updateMDLeapfrogSimple(), which does
 * not handle fully anistropic Parrinello-Rahman scaling, so we only have
 * options \p no and \p diagonal here and no anistropic option.
 */
enum class ApplyParrinelloRahmanVScaling
{
    no,      //!< Do not apply velocity scaling (not a PR-coupling run or step)
    diagonal //!< Apply velocity scaling using a diagonal matrix
};

/*! \brief Integrate using leap-frog with T-scaling and optionally diagonal Parrinello-Rahman p-coupling
 *
 * \tparam       numTempScaleValues     The number of different T-couple values
 * \tparam       applyPRVScaling        Apply Parrinello-Rahman velocity scaling
 * \param[in]    start                  Index of first atom to update
 * \param[in]    nrend                  Last atom to update: \p nrend - 1
 * \param[in]    dt                     The time step
 * \param[in]    dtPressureCouple       Time step for pressure coupling
 * \param[in]    invMassPerDim          1/mass per atom and dimension
 * \param[in]    tcstat                 Temperature coupling information
 * \param[in]    cTC                    T-coupling group index per atom
 * \param[in]    pRVScaleMatrixDiagonal Parrinello-Rahman v-scale matrix diagonal
 * \param[in]    x                      Input coordinates
 * \param[out]   xprime                 Updated coordinates
 * \param[inout] v                      Velocities
 * \param[in]    f                      Forces
 *
 * We expect this template to get good SIMD acceleration by most compilers,
 * unlike the more complex general template.
 * Note that we might get even better SIMD acceleration when we introduce
 * aligned (and padded) memory, possibly with some hints for the compilers.
 */
template<NumTempScaleValues numTempScaleValues, ApplyParrinelloRahmanVScaling applyPRVScaling>
static void updateMDLeapfrogSimple(int         start,
                                   int         nrend,
                                   real        dt,
                                   real        dtPressureCouple,
                                   const rvec* gmx_restrict          invMassPerDim,
                                   gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                   const unsigned short*             cTC,
                                   const rvec                        pRVScaleMatrixDiagonal,
                                   const rvec* gmx_restrict x,
                                   rvec* gmx_restrict xprime,
                                   rvec* gmx_restrict v,
                                   const rvec* gmx_restrict f)
{
    real lambdaGroup;

    if (numTempScaleValues == NumTempScaleValues::single)
    {
        lambdaGroup = tcstat[0].lambda;
    }

    for (int a = start; a < nrend; a++)
    {
        if (numTempScaleValues == NumTempScaleValues::multiple)
        {
            lambdaGroup = tcstat[cTC[a]].lambda;
        }

        for (int d = 0; d < DIM; d++)
        {
            /* Note that using rvec invMassPerDim results in more efficient
             * SIMD code, but this increases the cache pressure.
             * For large systems with PME on the CPU this slows down the
             * (then already slow) update by 20%. If all data remains in cache,
             * using rvec is much faster.
             */
            real vNew = lambdaGroup * v[a][d] + f[a][d] * invMassPerDim[a][d] * dt;

            if (applyPRVScaling == ApplyParrinelloRahmanVScaling::diagonal)
            {
                vNew -= dtPressureCouple * pRVScaleMatrixDiagonal[d] * v[a][d];
            }
            v[a][d]      = vNew;
            xprime[a][d] = x[a][d] + vNew * dt;
        }
    }
}

#if GMX_SIMD && GMX_SIMD_HAVE_REAL
#    define GMX_HAVE_SIMD_UPDATE 1
#else
#    define GMX_HAVE_SIMD_UPDATE 0
#endif

#if GMX_HAVE_SIMD_UPDATE

/*! \brief Load (aligned) the contents of GMX_SIMD_REAL_WIDTH rvec elements sequentially into 3 SIMD registers
 *
 * The loaded output is:
 * \p r0: { r[index][XX], r[index][YY], ... }
 * \p r1: { ... }
 * \p r2: { ..., r[index+GMX_SIMD_REAL_WIDTH-1][YY], r[index+GMX_SIMD_REAL_WIDTH-1][ZZ] }
 *
 * \param[in]  r      Real to an rvec array, has to be aligned to SIMD register width
 * \param[in]  index  Index of the first rvec triplet of reals to load
 * \param[out] r0     Pointer to first SIMD register
 * \param[out] r1     Pointer to second SIMD register
 * \param[out] r2     Pointer to third SIMD register
 */
static inline void simdLoadRvecs(const rvec* r, int index, SimdReal* r0, SimdReal* r1, SimdReal* r2)
{
    const real* realPtr = r[index];

    GMX_ASSERT(isSimdAligned(realPtr), "Pointer should be SIMD aligned");

    *r0 = simdLoad(realPtr + 0 * GMX_SIMD_REAL_WIDTH);
    *r1 = simdLoad(realPtr + 1 * GMX_SIMD_REAL_WIDTH);
    *r2 = simdLoad(realPtr + 2 * GMX_SIMD_REAL_WIDTH);
}

/*! \brief Store (aligned) 3 SIMD registers sequentially to GMX_SIMD_REAL_WIDTH rvec elements
 *
 * The stored output is:
 * \p r[index] = { { r0[0], r0[1], ... }
 * ...
 * \p r[index+GMX_SIMD_REAL_WIDTH-1] =  { ... , r2[GMX_SIMD_REAL_WIDTH-2], r2[GMX_SIMD_REAL_WIDTH-1] }
 *
 * \param[out] r      Pointer to an rvec array, has to be aligned to SIMD register width
 * \param[in]  index  Index of the first rvec triplet of reals to store to
 * \param[in]  r0     First SIMD register
 * \param[in]  r1     Second SIMD register
 * \param[in]  r2     Third SIMD register
 */
static inline void simdStoreRvecs(rvec* r, int index, SimdReal r0, SimdReal r1, SimdReal r2)
{
    real* realPtr = r[index];

    GMX_ASSERT(isSimdAligned(realPtr), "Pointer should be SIMD aligned");

    store(realPtr + 0 * GMX_SIMD_REAL_WIDTH, r0);
    store(realPtr + 1 * GMX_SIMD_REAL_WIDTH, r1);
    store(realPtr + 2 * GMX_SIMD_REAL_WIDTH, r2);
}

/*! \brief Integrate using leap-frog with single group T-scaling and SIMD
 *
 * \param[in]    start                  Index of first atom to update
 * \param[in]    nrend                  Last atom to update: \p nrend - 1
 * \param[in]    dt                     The time step
 * \param[in]    invMass                1/mass per atom
 * \param[in]    tcstat                 Temperature coupling information
 * \param[in]    x                      Input coordinates
 * \param[out]   xprime                 Updated coordinates
 * \param[inout] v                      Velocities
 * \param[in]    f                      Forces
 */
static void updateMDLeapfrogSimpleSimd(int         start,
                                       int         nrend,
                                       real        dt,
                                       const real* gmx_restrict          invMass,
                                       gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                       const rvec* gmx_restrict x,
                                       rvec* gmx_restrict xprime,
                                       rvec* gmx_restrict v,
                                       const rvec* gmx_restrict f)
{
    SimdReal timestep(dt);
    SimdReal lambdaSystem(tcstat[0].lambda);

    /* We declare variables here, since code is often slower when declaring them inside the loop */

    /* Note: We should implement a proper PaddedVector, so we don't need this check */
    GMX_ASSERT(isSimdAligned(invMass), "invMass should be aligned");

    for (int a = start; a < nrend; a += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal invMass0, invMass1, invMass2;
        expandScalarsToTriplets(simdLoad(invMass + a), &invMass0, &invMass1, &invMass2);

        SimdReal v0, v1, v2;
        SimdReal f0, f1, f2;
        simdLoadRvecs(v, a, &v0, &v1, &v2);
        simdLoadRvecs(f, a, &f0, &f1, &f2);

        v0 = fma(f0 * invMass0, timestep, lambdaSystem * v0);
        v1 = fma(f1 * invMass1, timestep, lambdaSystem * v1);
        v2 = fma(f2 * invMass2, timestep, lambdaSystem * v2);

        simdStoreRvecs(v, a, v0, v1, v2);

        SimdReal x0, x1, x2;
        simdLoadRvecs(x, a, &x0, &x1, &x2);

        SimdReal xprime0 = fma(v0, timestep, x0);
        SimdReal xprime1 = fma(v1, timestep, x1);
        SimdReal xprime2 = fma(v2, timestep, x2);

        simdStoreRvecs(xprime, a, xprime0, xprime1, xprime2);
    }
}

#endif // GMX_HAVE_SIMD_UPDATE

/*! \brief Sets the NEMD acceleration type */
enum class AccelerationType
{
    none,
    group,
    cosine
};

/*! \brief Integrate using leap-frog with support for everything
 *
 * \tparam       accelerationType       Type of NEMD acceleration
 * \param[in]    start                  Index of first atom to update
 * \param[in]    nrend                  Last atom to update: \p nrend - 1
 * \param[in]    doNoseHoover           If to apply Nose-Hoover T-coupling
 * \param[in]    dt                     The time step
 * \param[in]    dtPressureCouple       Time step for pressure coupling, is 0 when no pressure
 * coupling should be applied at this step \param[in]    ir                     The input parameter
 * record \param[in]    md                     Atom properties \param[in]    ekind Kinetic energy
 * data \param[in]    box                    The box dimensions \param[in]    x Input coordinates \param[out]
 * xprime                 Updated coordinates \param[inout] v                      Velocities \param[in]
 * f                      Forces \param[in]    nh_vxi                 Nose-Hoover velocity scaling
 * factors \param[in]    M                      Parrinello-Rahman scaling matrix
 */
template<AccelerationType accelerationType>
static void updateMDLeapfrogGeneral(int                   start,
                                    int                   nrend,
                                    bool                  doNoseHoover,
                                    real                  dt,
                                    real                  dtPressureCouple,
                                    const t_inputrec*     ir,
                                    const t_mdatoms*      md,
                                    const gmx_ekindata_t* ekind,
                                    const matrix          box,
                                    const rvec* gmx_restrict x,
                                    rvec* gmx_restrict xprime,
                                    rvec* gmx_restrict v,
                                    const rvec* gmx_restrict f,
                                    const double* gmx_restrict nh_vxi,
                                    const matrix               M)
{
    /* This is a version of the leap-frog integrator that supports
     * all combinations of T-coupling, P-coupling and NEMD.
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */

    gmx::ArrayRef<const t_grp_tcstat> tcstat  = ekind->tcstat;
    gmx::ArrayRef<const t_grp_acc>    grpstat = ekind->grpstat;
    const unsigned short*             cTC     = md->cTC;
    const unsigned short*             cACC    = md->cACC;
    const rvec*                       accel   = ir->opts.acc;

    const rvec* gmx_restrict invMassPerDim = md->invMassPerDim;

    /* Initialize group values, changed later when multiple groups are used */
    int  ga       = 0;
    int  gt       = 0;
    real factorNH = 0;

    real omega_Z = 2 * static_cast<real>(M_PI) / box[ZZ][ZZ];

    for (int n = start; n < nrend; n++)
    {
        if (cTC)
        {
            gt = cTC[n];
        }
        real lg = tcstat[gt].lambda;

        rvec vRel;
        real cosineZ, vCosine;
#ifdef __INTEL_COMPILER
#    pragma warning(disable : 280)
#endif
        switch (accelerationType)
        {
            case AccelerationType::none: copy_rvec(v[n], vRel); break;
            case AccelerationType::group:
                if (cACC)
                {
                    ga = cACC[n];
                }
                /* Avoid scaling the group velocity */
                rvec_sub(v[n], grpstat[ga].u, vRel);
                break;
            case AccelerationType::cosine:
                cosineZ = std::cos(x[n][ZZ] * omega_Z);
                vCosine = cosineZ * ekind->cosacc.vcos;
                /* Avoid scaling the cosine profile velocity */
                copy_rvec(v[n], vRel);
                vRel[XX] -= vCosine;
                break;
        }

        if (doNoseHoover)
        {
            /* Here we account for multiple time stepping, by increasing
             * the Nose-Hoover correction by nsttcouple
             */
            factorNH = 0.5 * ir->nsttcouple * dt * nh_vxi[gt];
        }

        for (int d = 0; d < DIM; d++)
        {
            real vNew = (lg * vRel[d]
                         + (f[n][d] * invMassPerDim[n][d] * dt - factorNH * vRel[d]
                            - dtPressureCouple * iprod(M[d], vRel)))
                        / (1 + factorNH);
            switch (accelerationType)
            {
                case AccelerationType::none: break;
                case AccelerationType::group:
                    /* Add back the mean velocity and apply acceleration */
                    vNew += grpstat[ga].u[d] + accel[ga][d] * dt;
                    break;
                case AccelerationType::cosine:
                    if (d == XX)
                    {
                        /* Add back the mean velocity and apply acceleration */
                        vNew += vCosine + cosineZ * ekind->cosacc.cos_accel * dt;
                    }
                    break;
            }
            v[n][d]      = vNew;
            xprime[n][d] = x[n][d] + vNew * dt;
        }
    }
}

/*! \brief Handles the Leap-frog MD x and v integration */
static void do_update_md(int                   start,
                         int                   nrend,
                         int64_t               step,
                         real                  dt,
                         const t_inputrec*     ir,
                         const t_mdatoms*      md,
                         const gmx_ekindata_t* ekind,
                         const matrix          box,
                         const rvec* gmx_restrict x,
                         rvec* gmx_restrict xprime,
                         rvec* gmx_restrict v,
                         const rvec* gmx_restrict f,
                         const double* gmx_restrict nh_vxi,
                         const matrix               M)
{
    GMX_ASSERT(nrend == start || xprime != x,
               "For SIMD optimization certain compilers need to have xprime != x");
    GMX_ASSERT(ir->eI == eiMD,
               "Leap-frog integrator was called while another integrator was requested");

    /* Note: Berendsen pressure scaling is handled after do_update_md() */
    bool doTempCouple = (ir->etc != etcNO && do_per_step(step + ir->nsttcouple - 1, ir->nsttcouple));
    bool doNoseHoover = (ir->etc == etcNOSEHOOVER && doTempCouple);
    bool doParrinelloRahman =
            (ir->epc == epcPARRINELLORAHMAN && do_per_step(step + ir->nstpcouple - 1, ir->nstpcouple));
    bool doPROffDiagonal = (doParrinelloRahman && (M[YY][XX] != 0 || M[ZZ][XX] != 0 || M[ZZ][YY] != 0));

    real dtPressureCouple = (doParrinelloRahman ? ir->nstpcouple * dt : 0);

    /* NEMD (also cosine) acceleration is applied in updateMDLeapFrogGeneral */
    bool doAcceleration = (ekind->bNEMD || ekind->cosacc.cos_accel != 0);

    if (doNoseHoover || doPROffDiagonal || doAcceleration)
    {
        matrix stepM;
        if (!doParrinelloRahman)
        {
            /* We should not apply PR scaling at this step */
            clear_mat(stepM);
        }
        else
        {
            copy_mat(M, stepM);
        }

        if (!doAcceleration)
        {
            updateMDLeapfrogGeneral<AccelerationType::none>(start, nrend, doNoseHoover, dt,
                                                            dtPressureCouple, ir, md, ekind, box, x,
                                                            xprime, v, f, nh_vxi, stepM);
        }
        else if (ekind->bNEMD)
        {
            updateMDLeapfrogGeneral<AccelerationType::group>(start, nrend, doNoseHoover, dt,
                                                             dtPressureCouple, ir, md, ekind, box,
                                                             x, xprime, v, f, nh_vxi, stepM);
        }
        else
        {
            updateMDLeapfrogGeneral<AccelerationType::cosine>(start, nrend, doNoseHoover, dt,
                                                              dtPressureCouple, ir, md, ekind, box,
                                                              x, xprime, v, f, nh_vxi, stepM);
        }
    }
    else
    {
        /* Use a simple and thus more efficient integration loop. */
        /* The simple loop does not check for particle type (so it can
         * be vectorized), which means we need to clear the velocities
         * of virtual sites in advance, when present. Note that vsite
         * velocities are computed after update and constraints from
         * their displacement.
         */
        if (md->haveVsites)
        {
            /* Note: The overhead of this loop is completely neligible */
            clearVsiteVelocities(start, nrend, md->ptype, v);
        }

        /* We determine if we have a single T-coupling lambda value for all
         * atoms. That allows for better SIMD acceleration in the template.
         * If we do not do temperature coupling (in the run or this step),
         * all scaling values are 1, so we effectively have a single value.
         */
        bool haveSingleTempScaleValue = (!doTempCouple || ekind->ngtc == 1);

        /* Extract some pointers needed by all cases */
        const unsigned short*             cTC           = md->cTC;
        gmx::ArrayRef<const t_grp_tcstat> tcstat        = ekind->tcstat;
        const rvec*                       invMassPerDim = md->invMassPerDim;

        if (doParrinelloRahman)
        {
            GMX_ASSERT(!doPROffDiagonal,
                       "updateMDLeapfrogSimple only support diagonal Parrinello-Rahman scaling "
                       "matrices");

            rvec diagM;
            for (int d = 0; d < DIM; d++)
            {
                diagM[d] = M[d][d];
            }

            if (haveSingleTempScaleValue)
            {
                updateMDLeapfrogSimple<NumTempScaleValues::single, ApplyParrinelloRahmanVScaling::diagonal>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, diagM, x,
                        xprime, v, f);
            }
            else
            {
                updateMDLeapfrogSimple<NumTempScaleValues::multiple, ApplyParrinelloRahmanVScaling::diagonal>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, diagM, x,
                        xprime, v, f);
            }
        }
        else
        {
            if (haveSingleTempScaleValue)
            {
                /* Note that modern compilers are pretty good at vectorizing
                 * updateMDLeapfrogSimple(). But the SIMD version will still
                 * be faster because invMass lowers the cache pressure
                 * compared to invMassPerDim.
                 */
#if GMX_HAVE_SIMD_UPDATE
                /* Check if we can use invmass instead of invMassPerDim */
                if (!md->havePartiallyFrozenAtoms)
                {
                    updateMDLeapfrogSimpleSimd(start, nrend, dt, md->invmass, tcstat, x, xprime, v, f);
                }
                else
#endif
                {
                    updateMDLeapfrogSimple<NumTempScaleValues::single, ApplyParrinelloRahmanVScaling::no>(
                            start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, nullptr,
                            x, xprime, v, f);
                }
            }
            else
            {
                updateMDLeapfrogSimple<NumTempScaleValues::multiple, ApplyParrinelloRahmanVScaling::no>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, nullptr, x,
                        xprime, v, f);
            }
        }
    }
}

static void do_update_vv_vel(int                  start,
                             int                  nrend,
                             real                 dt,
                             const rvec           accel[],
                             const ivec           nFreeze[],
                             const real           invmass[],
                             const unsigned short ptype[],
                             const unsigned short cFREEZE[],
                             const unsigned short cACC[],
                             rvec                 v[],
                             const rvec           f[],
                             gmx_bool             bExtended,
                             real                 veta,
                             real                 alpha)
{
    int  gf = 0, ga = 0;
    int  n, d;
    real g, mv1, mv2;

    if (bExtended)
    {
        g   = 0.25 * dt * veta * alpha;
        mv1 = std::exp(-g);
        mv2 = gmx::series_sinhx(g);
    }
    else
    {
        mv1 = 1.0;
        mv2 = 1.0;
    }
    for (n = start; n < nrend; n++)
    {
        real w_dt = invmass[n] * dt;
        if (cFREEZE)
        {
            gf = cFREEZE[n];
        }
        if (cACC)
        {
            ga = cACC[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
            {
                v[n][d] = mv1 * (mv1 * v[n][d] + 0.5 * (w_dt * mv2 * f[n][d])) + 0.5 * accel[ga][d] * dt;
            }
            else
            {
                v[n][d] = 0.0;
            }
        }
    }
} /* do_update_vv_vel */

static void do_update_vv_pos(int                  start,
                             int                  nrend,
                             real                 dt,
                             const ivec           nFreeze[],
                             const unsigned short ptype[],
                             const unsigned short cFREEZE[],
                             const rvec           x[],
                             rvec                 xprime[],
                             const rvec           v[],
                             gmx_bool             bExtended,
                             real                 veta)
{
    int  gf = 0;
    int  n, d;
    real g, mr1, mr2;

    /* Would it make more sense if Parrinello-Rahman was put here? */
    if (bExtended)
    {
        g   = 0.5 * dt * veta;
        mr1 = std::exp(g);
        mr2 = gmx::series_sinhx(g);
    }
    else
    {
        mr1 = 1.0;
        mr2 = 1.0;
    }

    for (n = start; n < nrend; n++)
    {

        if (cFREEZE)
        {
            gf = cFREEZE[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
            {
                xprime[n][d] = mr1 * (mr1 * x[n][d] + mr2 * dt * v[n][d]);
            }
            else
            {
                xprime[n][d] = x[n][d];
            }
        }
    }
} /* do_update_vv_pos */

gmx_stochd_t::gmx_stochd_t(const t_inputrec* ir)
{
    const t_grpopts* opts = &ir->opts;
    const int        ngtc = opts->ngtc;

    if (ir->eI == eiBD)
    {
        bd_rf.resize(ngtc);
    }
    else if (EI_SD(ir->eI))
    {
        sdc.resize(ngtc);
        sdsig.resize(ngtc);

        for (int gt = 0; gt < ngtc; gt++)
        {
            if (opts->tau_t[gt] > 0)
            {
                sdc[gt].em = std::exp(-ir->delta_t / opts->tau_t[gt]);
            }
            else
            {
                /* No friction and noise on this group */
                sdc[gt].em = 1;
            }
        }
    }
    else if (ETC_ANDERSEN(ir->etc))
    {
        randomize_group.resize(ngtc);
        boltzfac.resize(ngtc);

        /* for now, assume that all groups, if randomized, are randomized at the same rate, i.e. tau_t is the same. */
        /* since constraint groups don't necessarily match up with temperature groups! This is checked in readir.c */

        for (int gt = 0; gt < ngtc; gt++)
        {
            real reft = std::max<real>(0, opts->ref_t[gt]);
            if ((opts->tau_t[gt] > 0)
                && (reft > 0)) /* tau_t or ref_t = 0 means that no randomization is done */
            {
                randomize_group[gt] = true;
                boltzfac[gt]        = BOLTZ * opts->ref_t[gt];
            }
            else
            {
                randomize_group[gt] = false;
            }
        }
    }
}

void update_temperature_constants(gmx_stochd_t* sd, const t_inputrec* ir)
{
    if (ir->eI == eiBD)
    {
        if (ir->bd_fric != 0)
        {
            for (int gt = 0; gt < ir->opts.ngtc; gt++)
            {
                sd->bd_rf[gt] = std::sqrt(2.0 * BOLTZ * ir->opts.ref_t[gt] / (ir->bd_fric * ir->delta_t));
            }
        }
        else
        {
            for (int gt = 0; gt < ir->opts.ngtc; gt++)
            {
                sd->bd_rf[gt] = std::sqrt(2.0 * BOLTZ * ir->opts.ref_t[gt]);
            }
        }
    }
    if (ir->eI == eiSD1)
    {
        for (int gt = 0; gt < ir->opts.ngtc; gt++)
        {
            real kT = BOLTZ * ir->opts.ref_t[gt];
            /* The mass is accounted for later, since this differs per atom */
            sd->sdsig[gt].V = std::sqrt(kT * (1 - sd->sdc[gt].em * sd->sdc[gt].em));
        }
    }
}

Update::Impl::Impl(const t_inputrec* ir, BoxDeformation* boxDeformation)
{
    sd = std::make_unique<gmx_stochd_t>(ir);
    update_temperature_constants(sd.get(), ir);
    xp.resizeWithPadding(0);
    deform = boxDeformation;
}

void Update::setNumAtoms(int nAtoms)
{

    impl_->xp.resizeWithPadding(nAtoms);
}

/*! \brief Sets the SD update type */
enum class SDUpdate : int
{
    ForcesOnly,
    FrictionAndNoiseOnly,
    Combined
};

/*! \brief SD integrator update
 *
 * Two phases are required in the general case of a constrained
 * update, the first phase from the contribution of forces, before
 * applying constraints, and then a second phase applying the friction
 * and noise, and then further constraining. For details, see
 * Goga2012.
 *
 * Without constraints, the two phases can be combined, for
 * efficiency.
 *
 * Thus three instantiations of this templated function will be made,
 * two with only one contribution, and one with both contributions. */
template<SDUpdate updateType>
static void doSDUpdateGeneral(const gmx_stochd_t&  sd,
                              int                  start,
                              int                  nrend,
                              real                 dt,
                              const rvec           accel[],
                              const ivec           nFreeze[],
                              const real           invmass[],
                              const unsigned short ptype[],
                              const unsigned short cFREEZE[],
                              const unsigned short cACC[],
                              const unsigned short cTC[],
                              const rvec           x[],
                              rvec                 xprime[],
                              rvec                 v[],
                              const rvec           f[],
                              int64_t              step,
                              int                  seed,
                              const int*           gatindex)
{
    // cTC, cACC and cFREEZE can be nullptr any time, but various
    // instantiations do not make sense with particular pointer
    // values.
    if (updateType == SDUpdate::ForcesOnly)
    {
        GMX_ASSERT(f != nullptr, "SD update with only forces requires forces");
        GMX_ASSERT(cTC == nullptr, "SD update with only forces cannot handle temperature groups");
    }
    if (updateType == SDUpdate::FrictionAndNoiseOnly)
    {
        GMX_ASSERT(f == nullptr, "SD update with only noise cannot handle forces");
        GMX_ASSERT(cACC == nullptr, "SD update with only noise cannot handle acceleration groups");
    }
    if (updateType == SDUpdate::Combined)
    {
        GMX_ASSERT(f != nullptr, "SD update with forces and noise requires forces");
    }

    // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
    gmx::ThreeFry2x64<0>                       rng(seed, gmx::RandomDomain::UpdateCoordinates);
    gmx::TabulatedNormalDistribution<real, 14> dist;

    for (int n = start; n < nrend; n++)
    {
        int globalAtomIndex = gatindex ? gatindex[n] : n;
        rng.restart(step, globalAtomIndex);
        dist.reset();

        real inverseMass = invmass[n];
        real invsqrtMass = std::sqrt(inverseMass);

        int freezeGroup       = cFREEZE ? cFREEZE[n] : 0;
        int accelerationGroup = cACC ? cACC[n] : 0;
        int temperatureGroup  = cTC ? cTC[n] : 0;

        for (int d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[freezeGroup][d])
            {
                if (updateType == SDUpdate::ForcesOnly)
                {
                    real vn = v[n][d] + (inverseMass * f[n][d] + accel[accelerationGroup][d]) * dt;
                    v[n][d] = vn;
                    // Simple position update.
                    xprime[n][d] = x[n][d] + v[n][d] * dt;
                }
                else if (updateType == SDUpdate::FrictionAndNoiseOnly)
                {
                    real vn = v[n][d];
                    v[n][d] = (vn * sd.sdc[temperatureGroup].em
                               + invsqrtMass * sd.sdsig[temperatureGroup].V * dist(rng));
                    // The previous phase already updated the
                    // positions with a full v*dt term that must
                    // now be half removed.
                    xprime[n][d] = xprime[n][d] + 0.5 * (v[n][d] - vn) * dt;
                }
                else
                {
                    real vn = v[n][d] + (inverseMass * f[n][d] + accel[accelerationGroup][d]) * dt;
                    v[n][d] = (vn * sd.sdc[temperatureGroup].em
                               + invsqrtMass * sd.sdsig[temperatureGroup].V * dist(rng));
                    // Here we include half of the friction+noise
                    // update of v into the position update.
                    xprime[n][d] = x[n][d] + 0.5 * (vn + v[n][d]) * dt;
                }
            }
            else
            {
                // When using constraints, the update is split into
                // two phases, but we only need to zero the update of
                // virtual, shell or frozen particles in at most one
                // of the phases.
                if (updateType != SDUpdate::FrictionAndNoiseOnly)
                {
                    v[n][d]      = 0.0;
                    xprime[n][d] = x[n][d];
                }
            }
        }
    }
}

static void do_update_bd(int                  start,
                         int                  nrend,
                         real                 dt,
                         const ivec           nFreeze[],
                         const real           invmass[],
                         const unsigned short ptype[],
                         const unsigned short cFREEZE[],
                         const unsigned short cTC[],
                         const rvec           x[],
                         rvec                 xprime[],
                         rvec                 v[],
                         const rvec           f[],
                         real                 friction_coefficient,
                         const real*          rf,
                         int64_t              step,
                         int                  seed,
                         const int*           gatindex)
{
    /* note -- these appear to be full step velocities . . .  */
    int  gf = 0, gt = 0;
    real vn;
    real invfr = 0;
    int  n, d;
    // Use 1 bit of internal counters to give us 2*2 64-bits values per stream
    // Each 64-bit value is enough for 4 normal distribution table numbers.
    gmx::ThreeFry2x64<0>                       rng(seed, gmx::RandomDomain::UpdateCoordinates);
    gmx::TabulatedNormalDistribution<real, 14> dist;

    if (friction_coefficient != 0)
    {
        invfr = 1.0 / friction_coefficient;
    }

    for (n = start; (n < nrend); n++)
    {
        int ng = gatindex ? gatindex[n] : n;

        rng.restart(step, ng);
        dist.reset();

        if (cFREEZE)
        {
            gf = cFREEZE[n];
        }
        if (cTC)
        {
            gt = cTC[n];
        }
        for (d = 0; (d < DIM); d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
            {
                if (friction_coefficient != 0)
                {
                    vn = invfr * f[n][d] + rf[gt] * dist(rng);
                }
                else
                {
                    /* NOTE: invmass = 2/(mass*friction_constant*dt) */
                    vn = 0.5 * invmass[n] * f[n][d] * dt
                         + std::sqrt(0.5 * invmass[n]) * rf[gt] * dist(rng);
                }

                v[n][d]      = vn;
                xprime[n][d] = x[n][d] + vn * dt;
            }
            else
            {
                v[n][d]      = 0.0;
                xprime[n][d] = x[n][d];
            }
        }
    }
}

static void calc_ke_part_normal(ArrayRef<const RVec> v,
                                const t_grpopts*     opts,
                                const t_mdatoms*     md,
                                gmx_ekindata_t*      ekind,
                                t_nrnb*              nrnb,
                                gmx_bool             bEkinAveVel)
{
    int                         g;
    gmx::ArrayRef<t_grp_tcstat> tcstat  = ekind->tcstat;
    gmx::ArrayRef<t_grp_acc>    grpstat = ekind->grpstat;

    /* three main: VV with AveVel, vv with AveEkin, leap with AveEkin.  Leap with AveVel is also
       an option, but not supported now.
       bEkinAveVel: If TRUE, we sum into ekin, if FALSE, into ekinh.
     */

    /* group velocities are calculated in update_ekindata and
     * accumulated in acumulate_groups.
     * Now the partial global and groups ekin.
     */
    for (g = 0; (g < opts->ngtc); g++)
    {
        copy_mat(tcstat[g].ekinh, tcstat[g].ekinh_old);
        if (bEkinAveVel)
        {
            clear_mat(tcstat[g].ekinf);
            tcstat[g].ekinscalef_nhc = 1.0; /* need to clear this -- logic is complicated! */
        }
        else
        {
            clear_mat(tcstat[g].ekinh);
        }
    }
    ekind->dekindl_old = ekind->dekindl;
    int nthread        = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        // This OpenMP only loops over arrays and does not call any functions
        // or memory allocation. It should not be able to throw, so for now
        // we do not need a try/catch wrapper.
        int     start_t, end_t, n;
        int     ga, gt;
        rvec    v_corrt;
        real    hm;
        int     d, m;
        matrix* ekin_sum;
        real*   dekindl_sum;

        start_t = ((thread + 0) * md->homenr) / nthread;
        end_t   = ((thread + 1) * md->homenr) / nthread;

        ekin_sum    = ekind->ekin_work[thread];
        dekindl_sum = ekind->dekindl_work[thread];

        for (gt = 0; gt < opts->ngtc; gt++)
        {
            clear_mat(ekin_sum[gt]);
        }
        *dekindl_sum = 0.0;

        ga = 0;
        gt = 0;
        for (n = start_t; n < end_t; n++)
        {
            if (md->cACC)
            {
                ga = md->cACC[n];
            }
            if (md->cTC)
            {
                gt = md->cTC[n];
            }
            hm = 0.5 * md->massT[n];

            for (d = 0; (d < DIM); d++)
            {
                v_corrt[d] = v[n][d] - grpstat[ga].u[d];
            }
            for (d = 0; (d < DIM); d++)
            {
                for (m = 0; (m < DIM); m++)
                {
                    /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
                    ekin_sum[gt][m][d] += hm * v_corrt[m] * v_corrt[d];
                }
            }
            if (md->nMassPerturbed && md->bPerturbed[n])
            {
                *dekindl_sum += 0.5 * (md->massB[n] - md->massA[n]) * iprod(v_corrt, v_corrt);
            }
        }
    }

    ekind->dekindl = 0;
    for (int thread = 0; thread < nthread; thread++)
    {
        for (g = 0; g < opts->ngtc; g++)
        {
            if (bEkinAveVel)
            {
                m_add(tcstat[g].ekinf, ekind->ekin_work[thread][g], tcstat[g].ekinf);
            }
            else
            {
                m_add(tcstat[g].ekinh, ekind->ekin_work[thread][g], tcstat[g].ekinh);
            }
        }

        ekind->dekindl += *ekind->dekindl_work[thread];
    }

    inc_nrnb(nrnb, eNR_EKIN, md->homenr);
}

static void calc_ke_part_visc(const matrix         box,
                              ArrayRef<const RVec> x,
                              ArrayRef<const RVec> v,
                              const t_grpopts*     opts,
                              const t_mdatoms*     md,
                              gmx_ekindata_t*      ekind,
                              t_nrnb*              nrnb,
                              gmx_bool             bEkinAveVel)
{
    int                         start = 0, homenr = md->homenr;
    int                         g, d, n, m, gt = 0;
    rvec                        v_corrt;
    real                        hm;
    gmx::ArrayRef<t_grp_tcstat> tcstat = ekind->tcstat;
    t_cos_acc*                  cosacc = &(ekind->cosacc);
    real                        dekindl;
    real                        fac, cosz;
    double                      mvcos;

    for (g = 0; g < opts->ngtc; g++)
    {
        copy_mat(ekind->tcstat[g].ekinh, ekind->tcstat[g].ekinh_old);
        clear_mat(ekind->tcstat[g].ekinh);
    }
    ekind->dekindl_old = ekind->dekindl;

    fac     = 2 * M_PI / box[ZZ][ZZ];
    mvcos   = 0;
    dekindl = 0;
    for (n = start; n < start + homenr; n++)
    {
        if (md->cTC)
        {
            gt = md->cTC[n];
        }
        hm = 0.5 * md->massT[n];

        /* Note that the times of x and v differ by half a step */
        /* MRS -- would have to be changed for VV */
        cosz = std::cos(fac * x[n][ZZ]);
        /* Calculate the amplitude of the new velocity profile */
        mvcos += 2 * cosz * md->massT[n] * v[n][XX];

        copy_rvec(v[n], v_corrt);
        /* Subtract the profile for the kinetic energy */
        v_corrt[XX] -= cosz * cosacc->vcos;
        for (d = 0; (d < DIM); d++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
                if (bEkinAveVel)
                {
                    tcstat[gt].ekinf[m][d] += hm * v_corrt[m] * v_corrt[d];
                }
                else
                {
                    tcstat[gt].ekinh[m][d] += hm * v_corrt[m] * v_corrt[d];
                }
            }
        }
        if (md->nPerturbed && md->bPerturbed[n])
        {
            /* The minus sign here might be confusing.
             * The kinetic contribution from dH/dl doesn't come from
             * d m(l)/2 v^2 / dl, but rather from d p^2/2m(l) / dl,
             * where p are the momenta. The difference is only a minus sign.
             */
            dekindl -= 0.5 * (md->massB[n] - md->massA[n]) * iprod(v_corrt, v_corrt);
        }
    }
    ekind->dekindl = dekindl;
    cosacc->mvcos  = mvcos;

    inc_nrnb(nrnb, eNR_EKIN, homenr);
}

void calc_ke_part(ArrayRef<const RVec> x,
                  ArrayRef<const RVec> v,
                  const matrix         box,
                  const t_grpopts*     opts,
                  const t_mdatoms*     md,
                  gmx_ekindata_t*      ekind,
                  t_nrnb*              nrnb,
                  gmx_bool             bEkinAveVel)
{
    if (ekind->cosacc.cos_accel == 0)
    {
        calc_ke_part_normal(v, opts, md, ekind, nrnb, bEkinAveVel);
    }
    else
    {
        calc_ke_part_visc(box, x, v, opts, md, ekind, nrnb, bEkinAveVel);
    }
}

extern void init_ekinstate(ekinstate_t* ekinstate, const t_inputrec* ir)
{
    ekinstate->ekin_n = ir->opts.ngtc;
    snew(ekinstate->ekinh, ekinstate->ekin_n);
    snew(ekinstate->ekinf, ekinstate->ekin_n);
    snew(ekinstate->ekinh_old, ekinstate->ekin_n);
    ekinstate->ekinscalef_nhc.resize(ekinstate->ekin_n);
    ekinstate->ekinscaleh_nhc.resize(ekinstate->ekin_n);
    ekinstate->vscale_nhc.resize(ekinstate->ekin_n);
    ekinstate->dekindl          = 0;
    ekinstate->mvcos            = 0;
    ekinstate->hasReadEkinState = false;
}

void update_ekinstate(ekinstate_t* ekinstate, const gmx_ekindata_t* ekind)
{
    int i;

    for (i = 0; i < ekinstate->ekin_n; i++)
    {
        copy_mat(ekind->tcstat[i].ekinh, ekinstate->ekinh[i]);
        copy_mat(ekind->tcstat[i].ekinf, ekinstate->ekinf[i]);
        copy_mat(ekind->tcstat[i].ekinh_old, ekinstate->ekinh_old[i]);
        ekinstate->ekinscalef_nhc[i] = ekind->tcstat[i].ekinscalef_nhc;
        ekinstate->ekinscaleh_nhc[i] = ekind->tcstat[i].ekinscaleh_nhc;
        ekinstate->vscale_nhc[i]     = ekind->tcstat[i].vscale_nhc;
    }

    copy_mat(ekind->ekin, ekinstate->ekin_total);
    ekinstate->dekindl = ekind->dekindl;
    ekinstate->mvcos   = ekind->cosacc.mvcos;
}

void restore_ekinstate_from_state(const t_commrec* cr, gmx_ekindata_t* ekind, const ekinstate_t* ekinstate)
{
    int i, n;

    if (MASTER(cr))
    {
        for (i = 0; i < ekinstate->ekin_n; i++)
        {
            copy_mat(ekinstate->ekinh[i], ekind->tcstat[i].ekinh);
            copy_mat(ekinstate->ekinf[i], ekind->tcstat[i].ekinf);
            copy_mat(ekinstate->ekinh_old[i], ekind->tcstat[i].ekinh_old);
            ekind->tcstat[i].ekinscalef_nhc = ekinstate->ekinscalef_nhc[i];
            ekind->tcstat[i].ekinscaleh_nhc = ekinstate->ekinscaleh_nhc[i];
            ekind->tcstat[i].vscale_nhc     = ekinstate->vscale_nhc[i];
        }

        copy_mat(ekinstate->ekin_total, ekind->ekin);

        ekind->dekindl      = ekinstate->dekindl;
        ekind->cosacc.mvcos = ekinstate->mvcos;
        n                   = ekinstate->ekin_n;
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(n), &n, cr->mpi_comm_mygroup);
        for (i = 0; i < n; i++)
        {
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinh[0][0]), ekind->tcstat[i].ekinh[0],
                      cr->mpi_comm_mygroup);
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinf[0][0]), ekind->tcstat[i].ekinf[0],
                      cr->mpi_comm_mygroup);
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinh_old[0][0]),
                      ekind->tcstat[i].ekinh_old[0], cr->mpi_comm_mygroup);

            gmx_bcast(sizeof(ekind->tcstat[i].ekinscalef_nhc), &(ekind->tcstat[i].ekinscalef_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].ekinscaleh_nhc), &(ekind->tcstat[i].ekinscaleh_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].vscale_nhc), &(ekind->tcstat[i].vscale_nhc),
                      cr->mpi_comm_mygroup);
        }
        gmx_bcast(DIM * DIM * sizeof(ekind->ekin[0][0]), ekind->ekin[0], cr->mpi_comm_mygroup);

        gmx_bcast(sizeof(ekind->dekindl), &ekind->dekindl, cr->mpi_comm_mygroup);
        gmx_bcast(sizeof(ekind->cosacc.mvcos), &ekind->cosacc.mvcos, cr->mpi_comm_mygroup);
    }
}

void update_tcouple(int64_t           step,
                    const t_inputrec* inputrec,
                    t_state*          state,
                    gmx_ekindata_t*   ekind,
                    const t_extmass*  MassQ,
                    const t_mdatoms*  md)

{
    // This condition was explicitly checked in previous version, but should have never been satisfied
    GMX_ASSERT(!(EI_VV(inputrec->eI)
                 && (inputrecNvtTrotter(inputrec) || inputrecNptTrotter(inputrec)
                     || inputrecNphTrotter(inputrec))),
               "Temperature coupling was requested with velocity verlet and trotter");

    bool doTemperatureCoupling = false;

    // For VV temperature coupling parameters are updated on the current
    // step, for the others - one step before.
    if (inputrec->etc == etcNO)
    {
        doTemperatureCoupling = false;
    }
    else if (EI_VV(inputrec->eI))
    {
        doTemperatureCoupling = do_per_step(step, inputrec->nsttcouple);
    }
    else
    {
        doTemperatureCoupling = do_per_step(step + inputrec->nsttcouple - 1, inputrec->nsttcouple);
    }

    if (doTemperatureCoupling)
    {
        real dttc = inputrec->nsttcouple * inputrec->delta_t;

        // TODO: berendsen_tcoupl(...), nosehoover_tcoupl(...) and vrescale_tcoupl(...) update
        //      temperature coupling parameters, which should be reflected in the name of these
        //      subroutines
        switch (inputrec->etc)
        {
            case etcNO: break;
            case etcBERENDSEN:
                berendsen_tcoupl(inputrec, ekind, dttc, state->therm_integral);
                break;
            case etcNOSEHOOVER:
                nosehoover_tcoupl(&(inputrec->opts), ekind, dttc, state->nosehoover_xi.data(),
                                  state->nosehoover_vxi.data(), MassQ);
                break;
            case etcVRESCALE:
                vrescale_tcoupl(inputrec, step, ekind, dttc, state->therm_integral.data());
                break;
        }
        /* rescale in place here */
        if (EI_VV(inputrec->eI))
        {
            rescale_velocities(ekind, md, 0, md->homenr, state->v.rvec_array());
        }
    }
    else
    {
        // Set the T scaling lambda to 1 to have no scaling
        // TODO: Do we have to do it on every non-t-couple step?
        for (int i = 0; (i < inputrec->opts.ngtc); i++)
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

void getThreadAtomRange(int numThreads, int threadIndex, int numAtoms, int* startAtom, int* endAtom)
{
#if GMX_HAVE_SIMD_UPDATE
    constexpr int blockSize = GMX_SIMD_REAL_WIDTH;
#else
    constexpr int blockSize = 1;
#endif
    int numBlocks = (numAtoms + blockSize - 1) / blockSize;

    *startAtom = ((numBlocks * threadIndex) / numThreads) * blockSize;
    *endAtom   = ((numBlocks * (threadIndex + 1)) / numThreads) * blockSize;
    if (threadIndex == numThreads - 1)
    {
        *endAtom = numAtoms;
    }
}

void update_pcouple_before_coordinates(FILE*             fplog,
                                       int64_t           step,
                                       const t_inputrec* inputrec,
                                       t_state*          state,
                                       matrix            parrinellorahmanMu,
                                       matrix            M,
                                       gmx_bool          bInitStep)
{
    /* Berendsen P-coupling is completely handled after the coordinate update.
     * Trotter P-coupling is handled by separate calls to trotter_update().
     */
    if (inputrec->epc == epcPARRINELLORAHMAN
        && do_per_step(step + inputrec->nstpcouple - 1, inputrec->nstpcouple))
    {
        real dtpc = inputrec->nstpcouple * inputrec->delta_t;

        parrinellorahman_pcoupl(fplog, step, inputrec, dtpc, state->pres_prev, state->box,
                                state->box_rel, state->boxv, M, parrinellorahmanMu, bInitStep);
    }
}

void constrain_velocities(int64_t step,
                          real* dvdlambda, /* the contribution to be added to the bonded interactions */
                          t_state*          state,
                          tensor            vir_part,
                          gmx::Constraints* constr,
                          gmx_bool          bCalcVir,
                          bool              do_log,
                          bool              do_ene)
{
    if (!constr)
    {
        return;
    }

    /*
     *  Steps (7C, 8C)
     *  APPLY CONSTRAINTS:
     *  BLOCK SHAKE
     */

    {
        tensor vir_con;

        /* clear out constraints before applying */
        clear_mat(vir_part);

        /* Constrain the coordinates upd->xp */
        constr->apply(do_log, do_ene, step, 1, 1.0, state->x.arrayRefWithPadding(),
                      state->v.arrayRefWithPadding(), state->v.arrayRefWithPadding().unpaddedArrayRef(),
                      state->box, state->lambda[efptBONDED], dvdlambda, ArrayRefWithPadding<RVec>(),
                      bCalcVir ? &vir_con : nullptr, ConstraintVariable::Velocities);

        if (bCalcVir)
        {
            m_add(vir_part, vir_con, vir_part);
        }
    }
}

void constrain_coordinates(int64_t step,
                           real* dvdlambda, /* the contribution to be added to the bonded interactions */
                           t_state*          state,
                           tensor            vir_part,
                           Update*           upd,
                           gmx::Constraints* constr,
                           gmx_bool          bCalcVir,
                           bool              do_log,
                           bool              do_ene)
{
    if (!constr)
    {
        return;
    }

    {
        tensor vir_con;

        /* clear out constraints before applying */
        clear_mat(vir_part);

        /* Constrain the coordinates upd->xp */
        constr->apply(do_log, do_ene, step, 1, 1.0, state->x.arrayRefWithPadding(),
                      upd->xp()->arrayRefWithPadding(), ArrayRef<RVec>(), state->box,
                      state->lambda[efptBONDED], dvdlambda, state->v.arrayRefWithPadding(),
                      bCalcVir ? &vir_con : nullptr, ConstraintVariable::Positions);

        if (bCalcVir)
        {
            m_add(vir_part, vir_con, vir_part);
        }
    }
}

void update_sd_second_half(int64_t step,
                           real* dvdlambda, /* the contribution to be added to the bonded interactions */
                           const t_inputrec* inputrec, /* input record and box stuff	*/
                           const t_mdatoms*  md,
                           t_state*          state,
                           const t_commrec*  cr,
                           t_nrnb*           nrnb,
                           gmx_wallcycle_t   wcycle,
                           Update*           upd,
                           gmx::Constraints* constr,
                           bool              do_log,
                           bool              do_ene)
{
    if (!constr)
    {
        return;
    }
    if (inputrec->eI == eiSD1)
    {
        int homenr = md->homenr;

        /* Cast delta_t from double to real to make the integrators faster.
         * The only reason for having delta_t double is to get accurate values
         * for t=delta_t*step when step is larger than float precision.
         * For integration dt the accuracy of real suffices, since with
         * integral += dt*integrand the increment is nearly always (much) smaller
         * than the integral (and the integrand has real precision).
         */
        real dt = inputrec->delta_t;

        wallcycle_start(wcycle, ewcUPDATE);

        int nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)
        for (int th = 0; th < nth; th++)
        {
            try
            {
                int start_th, end_th;
                getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

                doSDUpdateGeneral<SDUpdate::FrictionAndNoiseOnly>(
                        *upd->sd(), start_th, end_th, dt, inputrec->opts.acc, inputrec->opts.nFreeze,
                        md->invmass, md->ptype, md->cFREEZE, nullptr, md->cTC, state->x.rvec_array(),
                        upd->xp()->rvec_array(), state->v.rvec_array(), nullptr, step, inputrec->ld_seed,
                        DOMAINDECOMP(cr) ? cr->dd->globalAtomIndices.data() : nullptr);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        inc_nrnb(nrnb, eNR_UPDATE, homenr);
        wallcycle_stop(wcycle, ewcUPDATE);

        /* Constrain the coordinates upd->xp for half a time step */
        constr->apply(do_log, do_ene, step, 1, 0.5, state->x.arrayRefWithPadding(),
                      upd->xp()->arrayRefWithPadding(), ArrayRef<RVec>(), state->box,
                      state->lambda[efptBONDED], dvdlambda, state->v.arrayRefWithPadding(), nullptr,
                      ConstraintVariable::Positions);
    }
}

void finish_update(const t_inputrec*       inputrec, /* input record and box stuff	*/
                   const t_mdatoms*        md,
                   t_state*                state,
                   gmx_wallcycle_t         wcycle,
                   Update*                 upd,
                   const gmx::Constraints* constr)
{
    /* NOTE: Currently we always integrate to a temporary buffer and
     * then copy the results back here.
     */

    wallcycle_start_nocount(wcycle, ewcUPDATE);

    const int homenr = md->homenr;
    auto      xp     = makeConstArrayRef(*upd->xp()).subArray(0, homenr);
    auto      x      = makeArrayRef(state->x).subArray(0, homenr);

    if (md->havePartiallyFrozenAtoms && constr != nullptr)
    {
        /* We have atoms that are frozen along some, but not all dimensions,
         * then constraints will have moved them also along the frozen dimensions.
         * To freeze such degrees of freedom we do not copy them back here.
         */
        const ivec* nFreeze = inputrec->opts.nFreeze;

        for (int i = 0; i < homenr; i++)
        {
            const int g = md->cFREEZE[i];

            for (int d = 0; d < DIM; d++)
            {
                if (nFreeze[g][d] == 0)
                {
                    x[i][d] = xp[i][d];
                }
            }
        }
    }
    else
    {
        /* We have no frozen atoms or fully frozen atoms which have not
         * been moved by the update, so we can simply copy all coordinates.
         */
        int gmx_unused nth = gmx_omp_nthreads_get(emntUpdate);
#pragma omp parallel for num_threads(nth) schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            // Trivial statement, does not throw
            x[i] = xp[i];
        }
    }

    wallcycle_stop(wcycle, ewcUPDATE);
}

void update_pcouple_after_coordinates(FILE*             fplog,
                                      int64_t           step,
                                      const t_inputrec* inputrec,
                                      const t_mdatoms*  md,
                                      const matrix      pressure,
                                      const matrix      forceVirial,
                                      const matrix      constraintVirial,
                                      matrix            pressureCouplingMu,
                                      t_state*          state,
                                      t_nrnb*           nrnb,
                                      Update*           upd,
                                      const bool        scaleCoordinates)
{
    int start  = 0;
    int homenr = md->homenr;

    /* Cast to real for faster code, no loss in precision (see comment above) */
    real dt = inputrec->delta_t;


    /* now update boxes */
    switch (inputrec->epc)
    {
        case (epcNO): break;
        case (epcBERENDSEN):
            if (do_per_step(step, inputrec->nstpcouple))
            {
                real dtpc = inputrec->nstpcouple * dt;
                berendsen_pcoupl(fplog, step, inputrec, dtpc, pressure, state->box, forceVirial,
                                 constraintVirial, pressureCouplingMu, &state->baros_integral);
                berendsen_pscale(inputrec, pressureCouplingMu, state->box, state->box_rel, start,
                                 homenr, state->x.rvec_array(), md->cFREEZE, nrnb, scaleCoordinates);
            }
            break;
        case (epcPARRINELLORAHMAN):
            if (do_per_step(step + inputrec->nstpcouple - 1, inputrec->nstpcouple))
            {
                /* The box velocities were updated in do_pr_pcoupl,
                 * but we dont change the box vectors until we get here
                 * since we need to be able to shift/unshift above.
                 */
                real dtpc = inputrec->nstpcouple * dt;
                for (int i = 0; i < DIM; i++)
                {
                    for (int m = 0; m <= i; m++)
                    {
                        state->box[i][m] += dtpc * state->boxv[i][m];
                    }
                }
                preserve_box_shape(inputrec, state->box_rel, state->box);

                /* Scale the coordinates */
                if (scaleCoordinates)
                {
                    auto x = state->x.rvec_array();
                    for (int n = start; n < start + homenr; n++)
                    {
                        tmvmul_ur0(pressureCouplingMu, x[n], x[n]);
                    }
                }
            }
            break;
        case (epcMTTK):
            switch (inputrec->epct)
            {
                case (epctISOTROPIC):
                    /* DIM * eta = ln V.  so DIM*eta_new = DIM*eta_old + DIM*dt*veta =>
                       ln V_new = ln V_old + 3*dt*veta => V_new = V_old*exp(3*dt*veta) =>
                       Side length scales as exp(veta*dt) */

                    msmul(state->box, std::exp(state->veta * dt), state->box);

                    /* Relate veta to boxv.  veta = d(eta)/dT = (1/DIM)*1/V dV/dT.
                       o               If we assume isotropic scaling, and box length scaling
                       factor L, then V = L^DIM (det(M)).  So dV/dt = DIM
                       L^(DIM-1) dL/dt det(M), and veta = (1/L) dL/dt.  The
                       determinant of B is L^DIM det(M), and the determinant
                       of dB/dt is (dL/dT)^DIM det (M).  veta will be
                       (det(dB/dT)/det(B))^(1/3).  Then since M =
                       B_new*(vol_new)^(1/3), dB/dT_new = (veta_new)*B(new). */

                    msmul(state->box, state->veta, state->boxv);
                    break;
                default: break;
            }
            break;
        default: break;
    }

    if (upd->deform())
    {
        auto localX = makeArrayRef(state->x).subArray(start, homenr);
        upd->deform()->apply(localX, state->box, step);
    }
}

void update_coords(int64_t           step,
                   const t_inputrec* inputrec, /* input record and box stuff	*/
                   const t_mdatoms*  md,
                   t_state*          state,
                   gmx::ArrayRefWithPadding<const gmx::RVec> f,
                   const t_fcdata*                           fcd,
                   const gmx_ekindata_t*                     ekind,
                   const matrix                              M,
                   Update*                                   upd,
                   int                                       UpdatePart,
                   const t_commrec* cr, /* these shouldn't be here -- need to think about it */
                   const gmx::Constraints* constr)
{
    gmx_bool bDoConstr = (nullptr != constr);

    /* Running the velocity half does nothing except for velocity verlet */
    if ((UpdatePart == etrtVELOCITY1 || UpdatePart == etrtVELOCITY2) && !EI_VV(inputrec->eI))
    {
        gmx_incons("update_coords called for velocity without VV integrator");
    }

    int homenr = md->homenr;

    /* Cast to real for faster code, no loss in precision (see comment above) */
    real dt = inputrec->delta_t;

    /* We need to update the NMR restraint history when time averaging is used */
    if (state->flags & (1 << estDISRE_RM3TAV))
    {
        update_disres_history(fcd, &state->hist);
    }
    if (state->flags & (1 << estORIRE_DTAV))
    {
        update_orires_history(fcd, &state->hist);
    }

    /* ############# START The update of velocities and positions ######### */
    int nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            const rvec* x_rvec  = state->x.rvec_array();
            rvec*       xp_rvec = upd->xp()->rvec_array();
            rvec*       v_rvec  = state->v.rvec_array();
            const rvec* f_rvec  = as_rvec_array(f.unpaddedArrayRef().data());

            switch (inputrec->eI)
            {
                case (eiMD):
                    do_update_md(start_th, end_th, step, dt, inputrec, md, ekind, state->box,
                                 x_rvec, xp_rvec, v_rvec, f_rvec, state->nosehoover_vxi.data(), M);
                    break;
                case (eiSD1):
                    if (bDoConstr)
                    {
                        // With constraints, the SD update is done in 2 parts
                        doSDUpdateGeneral<SDUpdate::ForcesOnly>(
                                *upd->sd(), start_th, end_th, dt, inputrec->opts.acc, inputrec->opts.nFreeze,
                                md->invmass, md->ptype, md->cFREEZE, md->cACC, nullptr, x_rvec,
                                xp_rvec, v_rvec, f_rvec, step, inputrec->ld_seed, nullptr);
                    }
                    else
                    {
                        doSDUpdateGeneral<SDUpdate::Combined>(
                                *upd->sd(), start_th, end_th, dt, inputrec->opts.acc,
                                inputrec->opts.nFreeze, md->invmass, md->ptype, md->cFREEZE, md->cACC,
                                md->cTC, x_rvec, xp_rvec, v_rvec, f_rvec, step, inputrec->ld_seed,
                                DOMAINDECOMP(cr) ? cr->dd->globalAtomIndices.data() : nullptr);
                    }
                    break;
                case (eiBD):
                    do_update_bd(start_th, end_th, dt, inputrec->opts.nFreeze, md->invmass,
                                 md->ptype, md->cFREEZE, md->cTC, x_rvec, xp_rvec, v_rvec, f_rvec,
                                 inputrec->bd_fric, upd->sd()->bd_rf.data(), step, inputrec->ld_seed,
                                 DOMAINDECOMP(cr) ? cr->dd->globalAtomIndices.data() : nullptr);
                    break;
                case (eiVV):
                case (eiVVAK):
                {
                    gmx_bool bExtended = (inputrec->etc == etcNOSEHOOVER || inputrec->epc == epcPARRINELLORAHMAN
                                          || inputrec->epc == epcMTTK);

                    /* assuming barostat coupled to group 0 */
                    real alpha = 1.0 + DIM / static_cast<real>(inputrec->opts.nrdf[0]);
                    switch (UpdatePart)
                    {
                        case etrtVELOCITY1:
                        case etrtVELOCITY2:
                            do_update_vv_vel(start_th, end_th, dt, inputrec->opts.acc,
                                             inputrec->opts.nFreeze, md->invmass, md->ptype, md->cFREEZE,
                                             md->cACC, v_rvec, f_rvec, bExtended, state->veta, alpha);
                            break;
                        case etrtPOSITION:
                            do_update_vv_pos(start_th, end_th, dt, inputrec->opts.nFreeze, md->ptype,
                                             md->cFREEZE, x_rvec, xp_rvec, v_rvec, bExtended, state->veta);
                            break;
                    }
                    break;
                }
                default: gmx_fatal(FARGS, "Don't know how to update coordinates");
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

extern gmx_bool update_randomize_velocities(const t_inputrec*        ir,
                                            int64_t                  step,
                                            const t_commrec*         cr,
                                            const t_mdatoms*         md,
                                            gmx::ArrayRef<gmx::RVec> v,
                                            const Update*            upd,
                                            const gmx::Constraints*  constr)
{

    real rate = (ir->delta_t) / ir->opts.tau_t[0];

    if (ir->etc == etcANDERSEN && constr != nullptr)
    {
        /* Currently, Andersen thermostat does not support constrained
           systems. Functionality exists in the andersen_tcoupl
           function in GROMACS 4.5.7 to allow this combination. That
           code could be ported to the current random-number
           generation approach, but has not yet been done because of
           lack of time and resources. */
        gmx_fatal(FARGS,
                  "Normal Andersen is currently not supported with constraints, use massive "
                  "Andersen instead");
    }

    /* proceed with andersen if 1) it's fixed probability per
       particle andersen or 2) it's massive andersen and it's tau_t/dt */
    if ((ir->etc == etcANDERSEN) || do_per_step(step, roundToInt(1.0 / rate)))
    {
        andersen_tcoupl(ir, step, cr, md, v, rate, upd->sd()->randomize_group, upd->sd()->boltzfac);
        return TRUE;
    }
    return FALSE;
}
