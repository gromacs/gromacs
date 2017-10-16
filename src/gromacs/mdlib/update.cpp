/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <math.h>
#include <stdio.h>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/*For debugging, start at v(-dt/2) for velolcity verlet -- uncomment next line */
/*#define STARTFROMDT2*/

typedef struct {
    double em;
} gmx_sd_const_t;

typedef struct {
    real V;
} gmx_sd_sigma_t;

typedef struct {
    /* BD stuff */
    real           *bd_rf;
    /* SD stuff */
    gmx_sd_const_t *sdc;
    gmx_sd_sigma_t *sdsig;
    /* andersen temperature control stuff */
    gmx_bool       *randomize_group;
    real           *boltzfac;
} gmx_stochd_t;

struct gmx_update_t
{
    gmx_stochd_t     *sd;
    /* xprime for constraint algorithms */
    PaddedRVecVector  xp;

    /* Variables for the deform algorithm */
    gmx_int64_t       deformref_step;
    matrix            deformref_box;
};

static bool isTemperatureCouplingStep(gmx_int64_t step, const t_inputrec *ir)
{
    /* We should only couple after a step where energies were determined (for leapfrog versions)
       or the step energies are determined, for velocity verlet versions */
    int offset;
    if (EI_VV(ir->eI))
    {
        offset = 0;
    }
    else
    {
        offset = 1;
    }
    return ir->etc != etcNO &&
           (ir->nsttcouple == 1 ||
            do_per_step(step + ir->nsttcouple - offset, ir->nsttcouple));
}

static bool isPressureCouplingStep(gmx_int64_t step, const t_inputrec *ir)
{
    GMX_ASSERT(ir->epc != epcMTTK, "MTTK pressure coupling is not handled here");

    int offset;
    if (ir->epc == epcBERENDSEN)
    {
        offset = 0;
    }
    else
    {
        offset = 1;
    }
    /* We should only couple after a step where pressures were determined */
    return ir->epc != etcNO &&
           (ir->nstpcouple == 1 ||
            do_per_step(step + ir->nstpcouple - offset, ir->nstpcouple));
}

/*! \brief Sets the number of different temperature coupling values */
enum class NumTempScaleValues
{
    single,   //!< Single T-scaling value (either one group or all values =1)
    multiple  //!< Multiple T-scaling values, need to use T-group indices
};

/*! \brief Sets if to apply no or diagonal Parrinello-Rahman pressure scaling
 *
 * Note that this enum is only used in updateMDLeapfrogSimple(), which does
 * not handle fully anistropic Parrinello-Rahman scaling, so we only have
 * options \p no and \p diagonal here and no anistropic option.
 */
enum class ApplyParrinelloRahmanVScaling
{
    no,       //!< Do not apply velocity scaling (not a PR-coupling run or step)
    diagonal  //!< Apply velocity scaling using a diagonal matrix
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
template<NumTempScaleValues            numTempScaleValues,
         ApplyParrinelloRahmanVScaling applyPRVScaling>
static void
updateMDLeapfrogSimple(int                       start,
                       int                       nrend,
                       real                      dt,
                       real                      dtPressureCouple,
                       const rvec * gmx_restrict invMassPerDim,
                       const t_grp_tcstat      * tcstat,
                       const unsigned short    * cTC,
                       const rvec                pRVScaleMatrixDiagonal,
                       const rvec * gmx_restrict x,
                       rvec       * gmx_restrict xprime,
                       rvec       * gmx_restrict v,
                       const rvec * gmx_restrict f)
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
            real vNew = lambdaGroup*v[a][d] + f[a][d]*invMassPerDim[a][d]*dt;

            if (applyPRVScaling == ApplyParrinelloRahmanVScaling::diagonal)
            {
                vNew -= dtPressureCouple*pRVScaleMatrixDiagonal[d]*v[a][d];
            }
            v[a][d]      = vNew;
            xprime[a][d] = x[a][d] + vNew*dt;
        }
    }
}

#if GMX_SIMD && GMX_SIMD_HAVE_REAL
#define GMX_HAVE_SIMD_UPDATE 1
#else
#define GMX_HAVE_SIMD_UPDATE 0
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
static inline void simdLoadRvecs(const rvec *r,
                                 int         index,
                                 SimdReal   *r0,
                                 SimdReal   *r1,
                                 SimdReal   *r2)
{
    const real *realPtr = r[index];

    GMX_ASSERT(isSimdAligned(realPtr), "Pointer should be SIMD aligned");

    *r0 = simdLoad(realPtr + 0*GMX_SIMD_REAL_WIDTH);
    *r1 = simdLoad(realPtr + 1*GMX_SIMD_REAL_WIDTH);
    *r2 = simdLoad(realPtr + 2*GMX_SIMD_REAL_WIDTH);
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
static inline void simdStoreRvecs(rvec     *r,
                                  int       index,
                                  SimdReal  r0,
                                  SimdReal  r1,
                                  SimdReal  r2)
{
    real *realPtr = r[index];

    GMX_ASSERT(isSimdAligned(realPtr), "Pointer should be SIMD aligned");

    store(realPtr + 0*GMX_SIMD_REAL_WIDTH, r0);
    store(realPtr + 1*GMX_SIMD_REAL_WIDTH, r1);
    store(realPtr + 2*GMX_SIMD_REAL_WIDTH, r2);
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
static void
updateMDLeapfrogSimpleSimd(int                       start,
                           int                       nrend,
                           real                      dt,
                           const real * gmx_restrict invMass,
                           const t_grp_tcstat      * tcstat,
                           const rvec * gmx_restrict x,
                           rvec       * gmx_restrict xprime,
                           rvec       * gmx_restrict v,
                           const rvec * gmx_restrict f)
{
    SimdReal                  timestep(dt);
    SimdReal                  lambdaSystem(tcstat[0].lambda);

    /* We declare variables here, since code is often slower when declaring them inside the loop */

    /* Note: We should implement a proper PaddedVector, so we don't need this check */
    GMX_ASSERT(isSimdAligned(invMass), "invMass should be aligned");

    for (int a = start; a < nrend; a += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal invMass0, invMass1, invMass2;
        expandScalarsToTriplets(simdLoad(invMass + a),
                                &invMass0, &invMass1, &invMass2);

        SimdReal v0, v1, v2;
        SimdReal f0, f1, f2;
        simdLoadRvecs(v, a, &v0, &v1, &v2);
        simdLoadRvecs(f, a, &f0, &f1, &f2);

        v0 = fma(f0*invMass0, timestep, lambdaSystem*v0);
        v1 = fma(f1*invMass1, timestep, lambdaSystem*v1);
        v2 = fma(f2*invMass2, timestep, lambdaSystem*v2);

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
    none, group, cosine
};

/*! \brief Integrate using leap-frog with support for everything
 *
 * \tparam       accelerationType       Type of NEMD acceleration
 * \param[in]    start                  Index of first atom to update
 * \param[in]    nrend                  Last atom to update: \p nrend - 1
 * \param[in]    doNoseHoover           If to apply Nose-Hoover T-coupling
 * \param[in]    dt                     The time step
 * \param[in]    dtPressureCouple       Time step for pressure coupling, is 0 when no pressure coupling should be applied at this step
 * \param[in]    ir                     The input parameter record
 * \param[in]    md                     Atom properties
 * \param[in]    ekind                  Kinetic energy data
 * \param[in]    box                    The box dimensions
 * \param[in]    x                      Input coordinates
 * \param[out]   xprime                 Updated coordinates
 * \param[inout] v                      Velocities
 * \param[in]    f                      Forces
 * \param[in]    nh_vxi                 Nose-Hoover velocity scaling factors
 * \param[in]    M                      Parrinello-Rahman scaling matrix
 */
template<AccelerationType accelerationType>
static void
updateMDLeapfrogGeneral(int                         start,
                        int                         nrend,
                        bool                        doNoseHoover,
                        real                        dt,
                        real                        dtPressureCouple,
                        const t_inputrec          * ir,
                        const t_mdatoms           * md,
                        const gmx_ekindata_t      * ekind,
                        const matrix                box,
                        const rvec   * gmx_restrict x,
                        rvec         * gmx_restrict xprime,
                        rvec         * gmx_restrict v,
                        const rvec   * gmx_restrict f,
                        const double * gmx_restrict nh_vxi,
                        const matrix                M)
{
    /* This is a version of the leap-frog integrator that supports
     * all combinations of T-coupling, P-coupling and NEMD.
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */

    const unsigned short    * cTC           = md->cTC;
    const t_grp_tcstat      * tcstat        = ekind->tcstat;

    const unsigned short    * cACC          = md->cACC;
    const rvec              * accel         = ir->opts.acc;
    const t_grp_acc         * grpstat       = ekind->grpstat;

    const rvec * gmx_restrict invMassPerDim = md->invMassPerDim;

    /* Initialize group values, changed later when multiple groups are used */
    int  ga       = 0;
    int  gt       = 0;
    real factorNH = 0;

    for (int n = start; n < nrend; n++)
    {
        if (cTC)
        {
            gt  = cTC[n];
        }
        real lg = tcstat[gt].lambda;

        rvec vRel;
        real cosineZ, vCosine;
#ifdef __INTEL_COMPILER
#pragma warning( disable : 280 )
#endif
        switch (accelerationType)
        {
            case AccelerationType::none:
                copy_rvec(v[n], vRel);
                break;
            case AccelerationType::group:
                if (cACC)
                {
                    ga = cACC[n];
                }
                /* Avoid scaling the group velocity */
                rvec_sub(v[n], grpstat[ga].u, vRel);
                break;
            case AccelerationType::cosine:
                cosineZ = std::cos(x[n][ZZ]*static_cast<real>(M_PI)/box[ZZ][ZZ]);
                vCosine = cosineZ*ekind->cosacc.vcos;
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
            factorNH = 0.5*ir->nsttcouple*dt*nh_vxi[gt];
        }

        for (int d = 0; d < DIM; d++)
        {
            real vNew =
                (lg*vRel[d] + (f[n][d]*invMassPerDim[n][d]*dt - factorNH*vRel[d]
                               - dtPressureCouple*iprod(M[d], vRel)))/(1 + factorNH);
            switch (accelerationType)
            {
                case AccelerationType::none:
                    break;
                case AccelerationType::group:
                    /* Add back the mean velocity and apply acceleration */
                    vNew += grpstat[ga].u[d] + accel[ga][d]*dt;
                    break;
                case AccelerationType::cosine:
                    if (d == XX)
                    {
                        /* Add back the mean velocity and apply acceleration */
                        vNew += vCosine + cosineZ*ekind->cosacc.cos_accel*dt;
                    }
                    break;
            }
            v[n][d]       = vNew;
            xprime[n][d]  = x[n][d] + vNew*dt;
        }
    }
}

/*! \brief Handles the Leap-frog MD x and v integration */
static void do_update_md(int                         start,
                         int                         nrend,
                         gmx_int64_t                 step,
                         real                        dt,
                         const t_inputrec          * ir,
                         const t_mdatoms           * md,
                         const gmx_ekindata_t      * ekind,
                         const matrix                box,
                         const rvec   * gmx_restrict x,
                         rvec         * gmx_restrict xprime,
                         rvec         * gmx_restrict v,
                         const rvec   * gmx_restrict f,
                         const double * gmx_restrict nh_vxi,
                         const matrix                M)
{
    GMX_ASSERT(nrend == start || xprime != x, "For SIMD optimization certain compilers need to have xprime != x");

    /* Note: Berendsen pressure scaling is handled after do_update_md() */
    bool doTempCouple       = isTemperatureCouplingStep(step, ir);
    bool doNoseHoover       = (ir->etc == etcNOSEHOOVER && doTempCouple);
    bool doParrinelloRahman = (ir->epc == epcPARRINELLORAHMAN && isPressureCouplingStep(step, ir));
    bool doPROffDiagonal    = (doParrinelloRahman && (M[YY][XX] != 0 || M[ZZ][XX] != 0 || M[ZZ][YY] != 0));

    real dtPressureCouple   = (doParrinelloRahman ? ir->nstpcouple*dt : 0);

    /* NEMD (also cosine) acceleration is applied in updateMDLeapFrogGeneral */
    bool doAcceleration     = (ekind->bNEMD || ekind->cosacc.cos_accel != 0);

    if (doNoseHoover || doPROffDiagonal || doAcceleration)
    {
        if (!doAcceleration)
        {
            updateMDLeapfrogGeneral<AccelerationType::none>
                (start, nrend, doNoseHoover, dt, dtPressureCouple,
                ir, md, ekind, box, x, xprime, v, f, nh_vxi, M);
        }
        else if (ekind->bNEMD)
        {
            updateMDLeapfrogGeneral<AccelerationType::group>
                (start, nrend, doNoseHoover, dt, dtPressureCouple,
                ir, md, ekind, box, x, xprime, v, f, nh_vxi, M);
        }
        else
        {
            updateMDLeapfrogGeneral<AccelerationType::cosine>
                (start, nrend, doNoseHoover, dt, dtPressureCouple,
                ir, md, ekind, box, x, xprime, v, f, nh_vxi, M);
        }
    }
    else
    {
        /* We determine if we have a single T-coupling lambda value for all
         * atoms. That allows for better SIMD acceleration in the template.
         * If we do not do temperature coupling (in the run or this step),
         * all scaling values are 1, so we effectively have a single value.
         */
        bool haveSingleTempScaleValue = (!doTempCouple || ekind->ngtc == 1);

        /* Extract some pointers needed by all cases */
        const unsigned short *cTC           = md->cTC;
        const t_grp_tcstat   *tcstat        = ekind->tcstat;
        const rvec           *invMassPerDim = md->invMassPerDim;

        if (doParrinelloRahman)
        {
            GMX_ASSERT(!doPROffDiagonal, "updateMDLeapfrogSimple only support diagonal Parrinello-Rahman scaling matrices");

            rvec diagM;
            for (int d = 0; d < DIM; d++)
            {
                diagM[d] = M[d][d];
            }

            if (haveSingleTempScaleValue)
            {
                updateMDLeapfrogSimple
                <NumTempScaleValues::single,
                 ApplyParrinelloRahmanVScaling::diagonal>
                    (start, nrend, dt, dtPressureCouple,
                    invMassPerDim, tcstat, cTC, diagM, x, xprime, v, f);
            }
            else
            {
                updateMDLeapfrogSimple
                <NumTempScaleValues::multiple,
                 ApplyParrinelloRahmanVScaling::diagonal>
                    (start, nrend, dt, dtPressureCouple,
                    invMassPerDim, tcstat, cTC, diagM, x, xprime, v, f);
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
                    updateMDLeapfrogSimpleSimd
                        (start, nrend, dt, md->invmass, tcstat, x, xprime, v, f);
                }
                else
#endif
                {
                    updateMDLeapfrogSimple
                    <NumTempScaleValues::single,
                     ApplyParrinelloRahmanVScaling::no>
                        (start, nrend, dt, dtPressureCouple,
                        invMassPerDim, tcstat, cTC, nullptr, x, xprime, v, f);
                }
            }
            else
            {
                updateMDLeapfrogSimple
                <NumTempScaleValues::multiple,
                 ApplyParrinelloRahmanVScaling::no>
                    (start, nrend, dt, dtPressureCouple,
                    invMassPerDim, tcstat, cTC, nullptr, x, xprime, v, f);
            }
        }
    }
}

static void do_update_vv_vel(int start, int nrend, real dt,
                             rvec accel[], ivec nFreeze[], real invmass[],
                             unsigned short ptype[], unsigned short cFREEZE[],
                             unsigned short cACC[], rvec v[], const rvec f[],
                             gmx_bool bExtended, real veta, real alpha)
{
    int    gf = 0, ga = 0;
    int    n, d;
    real   g, mv1, mv2;

    if (bExtended)
    {
        g        = 0.25*dt*veta*alpha;
        mv1      = std::exp(-g);
        mv2      = gmx::series_sinhx(g);
    }
    else
    {
        mv1      = 1.0;
        mv2      = 1.0;
    }
    for (n = start; n < nrend; n++)
    {
        real w_dt = invmass[n]*dt;
        if (cFREEZE)
        {
            gf   = cFREEZE[n];
        }
        if (cACC)
        {
            ga   = cACC[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
            {
                v[n][d]             = mv1*(mv1*v[n][d] + 0.5*(w_dt*mv2*f[n][d]))+0.5*accel[ga][d]*dt;
            }
            else
            {
                v[n][d]        = 0.0;
            }
        }
    }
} /* do_update_vv_vel */

static void do_update_vv_pos(int start, int nrend, real dt,
                             ivec nFreeze[],
                             unsigned short ptype[], unsigned short cFREEZE[],
                             const rvec x[], rvec xprime[], rvec v[],
                             gmx_bool bExtended, real veta)
{
    int    gf = 0;
    int    n, d;
    real   g, mr1, mr2;

    /* Would it make more sense if Parrinello-Rahman was put here? */
    if (bExtended)
    {
        g        = 0.5*dt*veta;
        mr1      = std::exp(g);
        mr2      = gmx::series_sinhx(g);
    }
    else
    {
        mr1      = 1.0;
        mr2      = 1.0;
    }

    for (n = start; n < nrend; n++)
    {

        if (cFREEZE)
        {
            gf   = cFREEZE[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
            {
                xprime[n][d]   = mr1*(mr1*x[n][d]+mr2*dt*v[n][d]);
            }
            else
            {
                xprime[n][d]   = x[n][d];
            }
        }
    }
} /* do_update_vv_pos */

static gmx_stochd_t *init_stochd(const t_inputrec *ir)
{
    gmx_stochd_t   *sd;

    snew(sd, 1);

    const t_grpopts *opts = &ir->opts;
    int              ngtc = opts->ngtc;

    if (ir->eI == eiBD)
    {
        snew(sd->bd_rf, ngtc);
    }
    else if (EI_SD(ir->eI))
    {
        snew(sd->sdc, ngtc);
        snew(sd->sdsig, ngtc);

        gmx_sd_const_t *sdc = sd->sdc;

        for (int gt = 0; gt < ngtc; gt++)
        {
            if (opts->tau_t[gt] > 0)
            {
                sdc[gt].em  = std::exp(-ir->delta_t/opts->tau_t[gt]);
            }
            else
            {
                /* No friction and noise on this group */
                sdc[gt].em  = 1;
            }
        }
    }
    else if (ETC_ANDERSEN(ir->etc))
    {
        snew(sd->randomize_group, ngtc);
        snew(sd->boltzfac, ngtc);

        /* for now, assume that all groups, if randomized, are randomized at the same rate, i.e. tau_t is the same. */
        /* since constraint groups don't necessarily match up with temperature groups! This is checked in readir.c */

        for (int gt = 0; gt < ngtc; gt++)
        {
            real reft = std::max<real>(0, opts->ref_t[gt]);
            if ((opts->tau_t[gt] > 0) && (reft > 0))  /* tau_t or ref_t = 0 means that no randomization is done */
            {
                sd->randomize_group[gt] = TRUE;
                sd->boltzfac[gt]        = BOLTZ*opts->ref_t[gt];
            }
            else
            {
                sd->randomize_group[gt] = FALSE;
            }
        }
    }

    return sd;
}

void update_temperature_constants(gmx_update_t *upd, const t_inputrec *ir)
{
    if (ir->eI == eiBD)
    {
        if (ir->bd_fric != 0)
        {
            for (int gt = 0; gt < ir->opts.ngtc; gt++)
            {
                upd->sd->bd_rf[gt] = std::sqrt(2.0*BOLTZ*ir->opts.ref_t[gt]/(ir->bd_fric*ir->delta_t));
            }
        }
        else
        {
            for (int gt = 0; gt < ir->opts.ngtc; gt++)
            {
                upd->sd->bd_rf[gt] = std::sqrt(2.0*BOLTZ*ir->opts.ref_t[gt]);
            }
        }
    }
    if (ir->eI == eiSD1)
    {
        for (int gt = 0; gt < ir->opts.ngtc; gt++)
        {
            real kT = BOLTZ*ir->opts.ref_t[gt];
            /* The mass is accounted for later, since this differs per atom */
            upd->sd->sdsig[gt].V  = std::sqrt(kT*(1 - upd->sd->sdc[gt].em*upd->sd->sdc[gt].em));
        }
    }
}

gmx_update_t *init_update(const t_inputrec *ir)
{
    gmx_update_t *upd = new(gmx_update_t);

    if (ir->eI == eiBD || EI_SD(ir->eI) || ir->etc == etcVRESCALE || ETC_ANDERSEN(ir->etc))
    {
        upd->sd    = init_stochd(ir);
    }

    update_temperature_constants(upd, ir);

    upd->xp.resize(0);

    return upd;
}

void update_realloc(gmx_update_t *upd, int natoms)
{
    GMX_ASSERT(upd, "upd must be allocated before its fields can be reallocated");

    upd->xp.resize(gmx::paddedRVecVectorSize(natoms));
}

static void do_update_sd1(gmx_stochd_t *sd,
                          int start, int nrend, real dt,
                          rvec accel[], ivec nFreeze[],
                          real invmass[], unsigned short ptype[],
                          unsigned short cFREEZE[], unsigned short cACC[],
                          unsigned short cTC[],
                          const rvec x[], rvec xprime[], rvec v[], const rvec f[],
                          gmx_bool bDoConstr,
                          gmx_bool bFirstHalfConstr,
                          gmx_int64_t step, int seed, int* gatindex)
{
    gmx_sd_const_t *sdc;
    gmx_sd_sigma_t *sig;
    int             gf = 0, ga = 0, gt = 0;
    real            ism;
    int             n, d;

    // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
    gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
    gmx::TabulatedNormalDistribution<real, 14> dist;

    sdc = sd->sdc;
    sig = sd->sdsig;

    if (!bDoConstr)
    {
        for (n = start; n < nrend; n++)
        {
            int  ng = gatindex ? gatindex[n] : n;

            rng.restart(step, ng);
            dist.reset();

            ism = std::sqrt(invmass[n]);

            if (cFREEZE)
            {
                gf  = cFREEZE[n];
            }
            if (cACC)
            {
                ga  = cACC[n];
            }
            if (cTC)
            {
                gt  = cTC[n];
            }

            for (d = 0; d < DIM; d++)
            {
                if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
                {
                    real sd_V, vn;

                    sd_V         = ism*sig[gt].V*dist(rng);
                    vn           = v[n][d] + (invmass[n]*f[n][d] + accel[ga][d])*dt;
                    v[n][d]      = vn*sdc[gt].em + sd_V;
                    /* Here we include half of the friction+noise
                     * update of v into the integration of x.
                     */
                    xprime[n][d] = x[n][d] + 0.5*(vn + v[n][d])*dt;
                }
                else
                {
                    v[n][d]      = 0.0;
                    xprime[n][d] = x[n][d];
                }
            }
        }
    }
    else
    {
        /* We do have constraints */
        if (bFirstHalfConstr)
        {
            /* First update without friction and noise */
            real im;

            for (n = start; n < nrend; n++)
            {
                im = invmass[n];

                if (cFREEZE)
                {
                    gf  = cFREEZE[n];
                }
                if (cACC)
                {
                    ga  = cACC[n];
                }

                for (d = 0; d < DIM; d++)
                {
                    if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
                    {
                        v[n][d]      = v[n][d] + (im*f[n][d] + accel[ga][d])*dt;
                        xprime[n][d] = x[n][d] +  v[n][d]*dt;
                    }
                    else
                    {
                        v[n][d]      = 0.0;
                        xprime[n][d] = x[n][d];
                    }
                }
            }
        }
        else
        {
            /* Update friction and noise only */
            for (n = start; n < nrend; n++)
            {
                int  ng = gatindex ? gatindex[n] : n;

                rng.restart(step, ng);
                dist.reset();

                ism = std::sqrt(invmass[n]);

                if (cFREEZE)
                {
                    gf  = cFREEZE[n];
                }
                if (cTC)
                {
                    gt  = cTC[n];
                }

                for (d = 0; d < DIM; d++)
                {
                    if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
                    {
                        real sd_V, vn;

                        sd_V         = ism*sig[gt].V*dist(rng);
                        vn           = v[n][d];
                        v[n][d]      = vn*sdc[gt].em + sd_V;
                        /* Add the friction and noise contribution only */
                        xprime[n][d] = xprime[n][d] + 0.5*(v[n][d] - vn)*dt;
                    }
                }
            }
        }
    }
}

static void do_update_bd(int start, int nrend, real dt,
                         ivec nFreeze[],
                         real invmass[], unsigned short ptype[],
                         unsigned short cFREEZE[], unsigned short cTC[],
                         const rvec x[], rvec xprime[], rvec v[],
                         const rvec f[], real friction_coefficient,
                         real *rf, gmx_int64_t step, int seed,
                         int* gatindex)
{
    /* note -- these appear to be full step velocities . . .  */
    int    gf = 0, gt = 0;
    real   vn;
    real   invfr = 0;
    int    n, d;
    // Use 1 bit of internal counters to give us 2*2 64-bits values per stream
    // Each 64-bit value is enough for 4 normal distribution table numbers.
    gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
    gmx::TabulatedNormalDistribution<real, 14> dist;

    if (friction_coefficient != 0)
    {
        invfr = 1.0/friction_coefficient;
    }

    for (n = start; (n < nrend); n++)
    {
        int  ng  = gatindex ? gatindex[n] : n;

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
                    vn = invfr*f[n][d] + rf[gt]*dist(rng);
                }
                else
                {
                    /* NOTE: invmass = 2/(mass*friction_constant*dt) */
                    vn = 0.5*invmass[n]*f[n][d]*dt
                        + std::sqrt(0.5*invmass[n])*rf[gt]*dist(rng);
                }

                v[n][d]      = vn;
                xprime[n][d] = x[n][d]+vn*dt;
            }
            else
            {
                v[n][d]      = 0.0;
                xprime[n][d] = x[n][d];
            }
        }
    }
}

static void dump_it_all(FILE gmx_unused *fp, const char gmx_unused *title,
                        int gmx_unused natoms,
                        const gmx_unused PaddedRVecVector *x,
                        const gmx_unused PaddedRVecVector *xp,
                        const gmx_unused PaddedRVecVector *v,
                        const gmx_unused PaddedRVecVector *f)
{
#ifdef DEBUG
    if (fp)
    {
        fprintf(fp, "%s\n", title);
        pr_rvecs(fp, 0, "x", as_rvec_array(x->data()), natoms);
        pr_rvecs(fp, 0, "xp", as_rvec_array(xp->data()), natoms);
        pr_rvecs(fp, 0, "v", as_rvec_array(v->data()), natoms);
        if (f != NULL)
        {
            pr_rvecs(fp, 0, "f", as_rvec_array(f->data()), natoms);
        }
    }
#endif
}

static void calc_ke_part_normal(rvec v[], t_grpopts *opts, t_mdatoms *md,
                                gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel)
{
    int           g;
    t_grp_tcstat *tcstat  = ekind->tcstat;
    t_grp_acc    *grpstat = ekind->grpstat;
    int           nthread, thread;

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
            tcstat[g].ekinscalef_nhc = 1.0;   /* need to clear this -- logic is complicated! */
        }
        else
        {
            clear_mat(tcstat[g].ekinh);
        }
    }
    ekind->dekindl_old = ekind->dekindl;
    nthread            = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        // This OpenMP only loops over arrays and does not call any functions
        // or memory allocation. It should not be able to throw, so for now
        // we do not need a try/catch wrapper.
        int     start_t, end_t, n;
        int     ga, gt;
        rvec    v_corrt;
        real    hm;
        int     d, m;
        matrix *ekin_sum;
        real   *dekindl_sum;

        start_t = ((thread+0)*md->homenr)/nthread;
        end_t   = ((thread+1)*md->homenr)/nthread;

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
            hm   = 0.5*md->massT[n];

            for (d = 0; (d < DIM); d++)
            {
                v_corrt[d]  = v[n][d]  - grpstat[ga].u[d];
            }
            for (d = 0; (d < DIM); d++)
            {
                for (m = 0; (m < DIM); m++)
                {
                    /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
                    ekin_sum[gt][m][d] += hm*v_corrt[m]*v_corrt[d];
                }
            }
            if (md->nMassPerturbed && md->bPerturbed[n])
            {
                *dekindl_sum +=
                    0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt, v_corrt);
            }
        }
    }

    ekind->dekindl = 0;
    for (thread = 0; thread < nthread; thread++)
    {
        for (g = 0; g < opts->ngtc; g++)
        {
            if (bEkinAveVel)
            {
                m_add(tcstat[g].ekinf, ekind->ekin_work[thread][g],
                      tcstat[g].ekinf);
            }
            else
            {
                m_add(tcstat[g].ekinh, ekind->ekin_work[thread][g],
                      tcstat[g].ekinh);
            }
        }

        ekind->dekindl += *ekind->dekindl_work[thread];
    }

    inc_nrnb(nrnb, eNR_EKIN, md->homenr);
}

static void calc_ke_part_visc(matrix box, rvec x[], rvec v[],
                              t_grpopts *opts, t_mdatoms *md,
                              gmx_ekindata_t *ekind,
                              t_nrnb *nrnb, gmx_bool bEkinAveVel)
{
    int           start = 0, homenr = md->homenr;
    int           g, d, n, m, gt = 0;
    rvec          v_corrt;
    real          hm;
    t_grp_tcstat *tcstat = ekind->tcstat;
    t_cos_acc    *cosacc = &(ekind->cosacc);
    real          dekindl;
    real          fac, cosz;
    double        mvcos;

    for (g = 0; g < opts->ngtc; g++)
    {
        copy_mat(ekind->tcstat[g].ekinh, ekind->tcstat[g].ekinh_old);
        clear_mat(ekind->tcstat[g].ekinh);
    }
    ekind->dekindl_old = ekind->dekindl;

    fac     = 2*M_PI/box[ZZ][ZZ];
    mvcos   = 0;
    dekindl = 0;
    for (n = start; n < start+homenr; n++)
    {
        if (md->cTC)
        {
            gt = md->cTC[n];
        }
        hm   = 0.5*md->massT[n];

        /* Note that the times of x and v differ by half a step */
        /* MRS -- would have to be changed for VV */
        cosz         = std::cos(fac*x[n][ZZ]);
        /* Calculate the amplitude of the new velocity profile */
        mvcos       += 2*cosz*md->massT[n]*v[n][XX];

        copy_rvec(v[n], v_corrt);
        /* Subtract the profile for the kinetic energy */
        v_corrt[XX] -= cosz*cosacc->vcos;
        for (d = 0; (d < DIM); d++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* if we're computing a full step velocity, v_corrt[d] has v(t).  Otherwise, v(t+dt/2) */
                if (bEkinAveVel)
                {
                    tcstat[gt].ekinf[m][d] += hm*v_corrt[m]*v_corrt[d];
                }
                else
                {
                    tcstat[gt].ekinh[m][d] += hm*v_corrt[m]*v_corrt[d];
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
            dekindl -= 0.5*(md->massB[n] - md->massA[n])*iprod(v_corrt, v_corrt);
        }
    }
    ekind->dekindl = dekindl;
    cosacc->mvcos  = mvcos;

    inc_nrnb(nrnb, eNR_EKIN, homenr);
}

void calc_ke_part(t_state *state, t_grpopts *opts, t_mdatoms *md,
                  gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel)
{
    if (ekind->cosacc.cos_accel == 0)
    {
        calc_ke_part_normal(as_rvec_array(state->v.data()), opts, md, ekind, nrnb, bEkinAveVel);
    }
    else
    {
        calc_ke_part_visc(state->box, as_rvec_array(state->x.data()), as_rvec_array(state->v.data()), opts, md, ekind, nrnb, bEkinAveVel);
    }
}

extern void init_ekinstate(ekinstate_t *ekinstate, const t_inputrec *ir)
{
    ekinstate->ekin_n = ir->opts.ngtc;
    snew(ekinstate->ekinh, ekinstate->ekin_n);
    snew(ekinstate->ekinf, ekinstate->ekin_n);
    snew(ekinstate->ekinh_old, ekinstate->ekin_n);
    ekinstate->ekinscalef_nhc.resize(ekinstate->ekin_n);
    ekinstate->ekinscaleh_nhc.resize(ekinstate->ekin_n);
    ekinstate->vscale_nhc.resize(ekinstate->ekin_n);
    ekinstate->dekindl = 0;
    ekinstate->mvcos   = 0;
}

void update_ekinstate(ekinstate_t *ekinstate, gmx_ekindata_t *ekind)
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

void restore_ekinstate_from_state(t_commrec *cr,
                                  gmx_ekindata_t *ekind, const ekinstate_t *ekinstate)
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
        gmx_bcast(sizeof(n), &n, cr);
        for (i = 0; i < n; i++)
        {
            gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinh[0][0]),
                      ekind->tcstat[i].ekinh[0], cr);
            gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinf[0][0]),
                      ekind->tcstat[i].ekinf[0], cr);
            gmx_bcast(DIM*DIM*sizeof(ekind->tcstat[i].ekinh_old[0][0]),
                      ekind->tcstat[i].ekinh_old[0], cr);

            gmx_bcast(sizeof(ekind->tcstat[i].ekinscalef_nhc),
                      &(ekind->tcstat[i].ekinscalef_nhc), cr);
            gmx_bcast(sizeof(ekind->tcstat[i].ekinscaleh_nhc),
                      &(ekind->tcstat[i].ekinscaleh_nhc), cr);
            gmx_bcast(sizeof(ekind->tcstat[i].vscale_nhc),
                      &(ekind->tcstat[i].vscale_nhc), cr);
        }
        gmx_bcast(DIM*DIM*sizeof(ekind->ekin[0][0]),
                  ekind->ekin[0], cr);

        gmx_bcast(sizeof(ekind->dekindl), &ekind->dekindl, cr);
        gmx_bcast(sizeof(ekind->cosacc.mvcos), &ekind->cosacc.mvcos, cr);
    }
}

void set_deform_reference_box(gmx_update_t *upd, gmx_int64_t step, matrix box)
{
    upd->deformref_step = step;
    copy_mat(box, upd->deformref_box);
}

static void deform(gmx_update_t *upd,
                   int start, int homenr, rvec x[], matrix box,
                   const t_inputrec *ir, gmx_int64_t step)
{
    matrix bnew, invbox, mu;
    real   elapsed_time;
    int    i, j;

    elapsed_time = (step + 1 - upd->deformref_step)*ir->delta_t;
    copy_mat(box, bnew);
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            if (ir->deform[i][j] != 0)
            {
                bnew[i][j] =
                    upd->deformref_box[i][j] + elapsed_time*ir->deform[i][j];
            }
        }
    }
    /* We correct the off-diagonal elements,
     * which can grow indefinitely during shearing,
     * so the shifts do not get messed up.
     */
    for (i = 1; i < DIM; i++)
    {
        for (j = i-1; j >= 0; j--)
        {
            while (bnew[i][j] - box[i][j] > 0.5*bnew[j][j])
            {
                rvec_dec(bnew[i], bnew[j]);
            }
            while (bnew[i][j] - box[i][j] < -0.5*bnew[j][j])
            {
                rvec_inc(bnew[i], bnew[j]);
            }
        }
    }
    gmx::invertBoxMatrix(box, invbox);
    copy_mat(bnew, box);
    mmul_ur0(box, invbox, mu);

    for (i = start; i < start+homenr; i++)
    {
        x[i][XX] = mu[XX][XX]*x[i][XX]+mu[YY][XX]*x[i][YY]+mu[ZZ][XX]*x[i][ZZ];
        x[i][YY] = mu[YY][YY]*x[i][YY]+mu[ZZ][YY]*x[i][ZZ];
        x[i][ZZ] = mu[ZZ][ZZ]*x[i][ZZ];
    }
}

void update_tcouple(gmx_int64_t       step,
                    t_inputrec       *inputrec,
                    t_state          *state,
                    gmx_ekindata_t   *ekind,
                    t_extmass        *MassQ,
                    t_mdatoms        *md)

{
    bool doTemperatureCoupling = false;

    /* if using vv with trotter decomposition methods, we do this elsewhere in the code */
    if (!(EI_VV(inputrec->eI) &&
          (inputrecNvtTrotter(inputrec) || inputrecNptTrotter(inputrec) || inputrecNphTrotter(inputrec))))
    {
        doTemperatureCoupling = isTemperatureCouplingStep(step, inputrec);
    }

    if (doTemperatureCoupling)
    {
        real dttc = inputrec->nsttcouple*inputrec->delta_t;

        switch (inputrec->etc)
        {
            case etcNO:
                break;
            case etcBERENDSEN:
                berendsen_tcoupl(inputrec, ekind, dttc, state->therm_integral);
                break;
            case etcNOSEHOOVER:
                nosehoover_tcoupl(&(inputrec->opts), ekind, dttc,
                                  state->nosehoover_xi.data(), state->nosehoover_vxi.data(), MassQ);
                break;
            case etcVRESCALE:
                vrescale_tcoupl(inputrec, step, ekind, dttc,
                                state->therm_integral.data());
                break;
        }
        /* rescale in place here */
        if (EI_VV(inputrec->eI))
        {
            rescale_velocities(ekind, md, 0, md->homenr, as_rvec_array(state->v.data()));
        }
    }
    else
    {
        /* Set the T scaling lambda to 1 to have no scaling */
        for (int i = 0; (i < inputrec->opts.ngtc); i++)
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

/*! \brief Computes the atom range for a thread to operate on, ensured SIMD aligned ranges
 *
 * \param[in]  numThreads   The number of threads to divide atoms over
 * \param[in]  threadIndex  The thread to get the range for
 * \param[in]  numAtoms     The total number of atoms (on this rank)
 * \param[out] startAtom    The start of the atom range
 * \param[out] endAtom      The end of the atom range, note that this is in general not a multiple of the SIMD width
 */
static void getThreadAtomRange(int numThreads, int threadIndex, int numAtoms,
                               int *startAtom, int *endAtom)
{
#if GMX_HAVE_SIMD_UPDATE
    constexpr int blockSize = GMX_SIMD_REAL_WIDTH;
#else
    constexpr int blockSize = 1;
#endif
    int           numBlocks = (numAtoms + blockSize - 1)/blockSize;

    *startAtom    = ((numBlocks* threadIndex     )/numThreads)*blockSize;
    *endAtom      = ((numBlocks*(threadIndex + 1))/numThreads)*blockSize;
    if (threadIndex == numThreads - 1)
    {
        *endAtom  = numAtoms;
    }
}

void update_pcouple_before_coordinates(FILE             *fplog,
                                       gmx_int64_t       step,
                                       const t_inputrec *inputrec,
                                       t_state          *state,
                                       matrix            parrinellorahmanMu,
                                       matrix            M,
                                       gmx_bool          bInitStep)
{
    /* Berendsen P-coupling is completely handled after the coordinate update.
     * Trotter P-coupling is handled by separate calls to trotter_update().
     */
    if (inputrec->epc == epcPARRINELLORAHMAN &&
        isPressureCouplingStep(step, inputrec))
    {
        real dtpc = inputrec->nstpcouple*inputrec->delta_t;

        parrinellorahman_pcoupl(fplog, step, inputrec, dtpc, state->pres_prev,
                                state->box, state->box_rel, state->boxv,
                                M, parrinellorahmanMu, bInitStep);
    }
}

void update_constraints(FILE             *fplog,
                        gmx_int64_t       step,
                        real             *dvdlambda, /* the contribution to be added to the bonded interactions */
                        t_inputrec       *inputrec,  /* input record and box stuff	*/
                        t_mdatoms        *md,
                        t_state          *state,
                        gmx_bool          bMolPBC,
                        t_graph          *graph,
                        PaddedRVecVector *force,     /* forces on home particles */
                        t_idef           *idef,
                        tensor            vir_part,
                        t_commrec        *cr,
                        t_nrnb           *nrnb,
                        gmx_wallcycle_t   wcycle,
                        gmx_update_t     *upd,
                        gmx_constr_t      constr,
                        gmx_bool          bFirstHalf,
                        gmx_bool          bCalcVir)
{
    gmx_bool             bLastStep, bLog = FALSE, bEner = FALSE, bDoConstr = FALSE;
    tensor               vir_con;
    int                  nth, th;

    if (constr)
    {
        bDoConstr = TRUE;
    }
    if (bFirstHalf && !EI_VV(inputrec->eI))
    {
        bDoConstr = FALSE;
    }

    /* for now, SD update is here -- though it really seems like it
       should be reformulated as a velocity verlet method, since it has two parts */

    int  homenr = md->homenr;

    /* Cast delta_t from double to real to make the integrators faster.
     * The only reason for having delta_t double is to get accurate values
     * for t=delta_t*step when step is larger than float precision.
     * For integration dt the accuracy of real suffices, since with
     * integral += dt*integrand the increment is nearly always (much) smaller
     * than the integral (and the integrand has real precision).
     */
    real dt     = inputrec->delta_t;

    /*
     *  Steps (7C, 8C)
     *  APPLY CONSTRAINTS:
     *  BLOCK SHAKE

     * When doing PR pressure coupling we have to constrain the
     * bonds in each iteration. If we are only using Nose-Hoover tcoupling
     * it is enough to do this once though, since the relative velocities
     * after this will be normal to the bond vector
     */

    if (bDoConstr)
    {
        /* clear out constraints before applying */
        clear_mat(vir_part);

        bLastStep = (step == inputrec->init_step+inputrec->nsteps);
        bLog      = (do_per_step(step, inputrec->nstlog) || bLastStep || (step < 0));
        bEner     = (do_per_step(step, inputrec->nstenergy) || bLastStep);
        /* Constrain the coordinates upd->xp */
        wallcycle_start(wcycle, ewcCONSTR);
        if (EI_VV(inputrec->eI) && bFirstHalf)
        {
            constrain(nullptr, bLog, bEner, constr, idef,
                      inputrec, cr, step, 1, 1.0, md,
                      as_rvec_array(state->x.data()), as_rvec_array(state->v.data()), as_rvec_array(state->v.data()),
                      bMolPBC, state->box,
                      state->lambda[efptBONDED], dvdlambda,
                      nullptr, bCalcVir ? &vir_con : nullptr, nrnb, econqVeloc);
        }
        else
        {
            constrain(nullptr, bLog, bEner, constr, idef,
                      inputrec, cr, step, 1, 1.0, md,
                      as_rvec_array(state->x.data()), as_rvec_array(upd->xp.data()), nullptr,
                      bMolPBC, state->box,
                      state->lambda[efptBONDED], dvdlambda,
                      as_rvec_array(state->v.data()), bCalcVir ? &vir_con : nullptr, nrnb, econqCoord);
        }
        wallcycle_stop(wcycle, ewcCONSTR);

        where();

        dump_it_all(fplog, "After Shake",
                    state->natoms, &state->x, &upd->xp, &state->v, force);

        if (bCalcVir)
        {
            m_add(vir_part, vir_con, vir_part);
            if (debug)
            {
                pr_rvecs(debug, 0, "constraint virial", vir_part, DIM);
            }
        }
    }

    where();

    if (inputrec->eI == eiSD1 && bDoConstr && !bFirstHalf)
    {
        wallcycle_start(wcycle, ewcUPDATE);

        nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)
        for (th = 0; th < nth; th++)
        {
            try
            {
                int start_th, end_th;
                getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

                /* The second part of the SD integration */
                do_update_sd1(upd->sd,
                              start_th, end_th, dt,
                              inputrec->opts.acc, inputrec->opts.nFreeze,
                              md->invmass, md->ptype,
                              md->cFREEZE, md->cACC, md->cTC,
                              as_rvec_array(state->x.data()), as_rvec_array(upd->xp.data()), as_rvec_array(state->v.data()), as_rvec_array(force->data()),
                              bDoConstr, FALSE,
                              step, inputrec->ld_seed,
                              DOMAINDECOMP(cr) ? cr->dd->gatindex : nullptr);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        }
        inc_nrnb(nrnb, eNR_UPDATE, homenr);
        wallcycle_stop(wcycle, ewcUPDATE);

        if (bDoConstr)
        {
            /* Constrain the coordinates upd->xp for half a time step */
            wallcycle_start(wcycle, ewcCONSTR);

            constrain(nullptr, bLog, bEner, constr, idef,
                      inputrec, cr, step, 1, 0.5, md,
                      as_rvec_array(state->x.data()), as_rvec_array(upd->xp.data()), nullptr,
                      bMolPBC, state->box,
                      state->lambda[efptBONDED], dvdlambda,
                      as_rvec_array(state->v.data()), nullptr, nrnb, econqCoord);

            wallcycle_stop(wcycle, ewcCONSTR);
        }
    }

    /* We must always unshift after updating coordinates; if we did not shake
       x was shifted in do_force */

    if (!(bFirstHalf)) /* in the first half of vv, no shift. */
    {
        /* NOTE This part of the update actually does not belong with
         * the constraints, since we also call it without constraints.
         * But currently we always integrate to a temporary buffer and
         * then copy the results back here.
         */
        wallcycle_start_nocount(wcycle, ewcUPDATE);

        if (md->cFREEZE != nullptr && constr != nullptr)
        {
            /* If we have atoms that are frozen along some, but not all
             * dimensions, the constraints will have moved them also along
             * the frozen dimensions. To freeze such degrees of freedom
             * we copy them back here to later copy them forward. It would
             * be more elegant and slightly more efficient to copies zero
             * times instead of twice, but the graph case below prevents this.
             */
            const ivec *nFreeze                     = inputrec->opts.nFreeze;
            bool        partialFreezeAndConstraints = false;
            for (int g = 0; g < inputrec->opts.ngfrz; g++)
            {
                int numFreezeDim = nFreeze[g][XX] + nFreeze[g][YY] + nFreeze[g][ZZ];
                if (numFreezeDim > 0 && numFreezeDim < 3)
                {
                    partialFreezeAndConstraints = true;
                }
            }
            if (partialFreezeAndConstraints)
            {
                for (int i = 0; i < homenr; i++)
                {
                    int g = md->cFREEZE[i];

                    for (int d = 0; d < DIM; d++)
                    {
                        if (nFreeze[g][d])
                        {
                            upd->xp[i][d] = state->x[i][d];
                        }
                    }
                }
            }
        }

        if (graph && (graph->nnodes > 0))
        {
            unshift_x(graph, state->box, as_rvec_array(state->x.data()), as_rvec_array(upd->xp.data()));
            if (TRICLINIC(state->box))
            {
                inc_nrnb(nrnb, eNR_SHIFTX, 2*graph->nnodes);
            }
            else
            {
                inc_nrnb(nrnb, eNR_SHIFTX, graph->nnodes);
            }
        }
        else
        {
            /* The copy is performance sensitive, so use a bare pointer */
            rvec *xp = as_rvec_array(upd->xp.data());
#ifndef __clang_analyzer__
            // cppcheck-suppress unreadVariable
            nth = gmx_omp_nthreads_get(emntUpdate);
#endif
#pragma omp parallel for num_threads(nth) schedule(static)
            for (int i = 0; i < homenr; i++)
            {
                // Trivial statement, does not throw
                copy_rvec(xp[i], state->x[i]);
            }
        }
        wallcycle_stop(wcycle, ewcUPDATE);

        dump_it_all(fplog, "After unshift",
                    state->natoms, &state->x, &upd->xp, &state->v, force);
    }
/* ############# END the update of velocities and positions ######### */
}

void update_pcouple_after_coordinates(FILE             *fplog,
                                      gmx_int64_t       step,
                                      const t_inputrec *inputrec,
                                      const t_mdatoms  *md,
                                      const matrix      pressure,
                                      const matrix      forceVirial,
                                      const matrix      constraintVirial,
                                      const matrix      parrinellorahmanMu,
                                      t_state          *state,
                                      t_nrnb           *nrnb,
                                      gmx_update_t     *upd)
{
    int  start  = 0;
    int  homenr = md->homenr;

    /* Cast to real for faster code, no loss in precision (see comment above) */
    real dt     = inputrec->delta_t;

    where();

    /* now update boxes */
    switch (inputrec->epc)
    {
        case (epcNO):
            break;
        case (epcBERENDSEN):
            if (isPressureCouplingStep(step, inputrec))
            {
                real   dtpc = inputrec->nstpcouple*dt;
                matrix mu;
                berendsen_pcoupl(fplog, step, inputrec, dtpc,
                                 pressure, state->box,
                                 forceVirial, constraintVirial,
                                 mu, &state->baros_integral);
                berendsen_pscale(inputrec, mu, state->box, state->box_rel,
                                 start, homenr, as_rvec_array(state->x.data()),
                                 md->cFREEZE, nrnb);
            }
            break;
        case (epcPARRINELLORAHMAN):
            if (isPressureCouplingStep(step, inputrec))
            {
                /* The box velocities were updated in do_pr_pcoupl,
                 * but we dont change the box vectors until we get here
                 * since we need to be able to shift/unshift above.
                 */
                real dtpc = inputrec->nstpcouple*dt;
                for (int i = 0; i < DIM; i++)
                {
                    for (int m = 0; m <= i; m++)
                    {
                        state->box[i][m] += dtpc*state->boxv[i][m];
                    }
                }
                preserve_box_shape(inputrec, state->box_rel, state->box);

                /* Scale the coordinates */
                for (int n = start; n < start + homenr; n++)
                {
                    tmvmul_ur0(parrinellorahmanMu, state->x[n], state->x[n]);
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

                    msmul(state->box, std::exp(state->veta*dt), state->box);

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
                default:
                    break;
            }
            break;
        default:
            break;
    }

    if (inputrecDeform(inputrec))
    {
        deform(upd, start, homenr, as_rvec_array(state->x.data()), state->box, inputrec, step);
    }
    where();
    dump_it_all(fplog, "After update",
                state->natoms, &state->x, &upd->xp, &state->v, nullptr);
}

void update_coords(FILE             *fplog,
                   gmx_int64_t       step,
                   t_inputrec       *inputrec,  /* input record and box stuff	*/
                   t_mdatoms        *md,
                   t_state          *state,
                   PaddedRVecVector *f,    /* forces on home particles */
                   t_fcdata         *fcd,
                   gmx_ekindata_t   *ekind,
                   matrix            M,
                   gmx_update_t     *upd,
                   int               UpdatePart,
                   t_commrec        *cr, /* these shouldn't be here -- need to think about it */
                   gmx_constr_t      constr)
{
    gmx_bool bDoConstr = (nullptr != constr);

    /* Running the velocity half does nothing except for velocity verlet */
    if ((UpdatePart == etrtVELOCITY1 || UpdatePart == etrtVELOCITY2) &&
        !EI_VV(inputrec->eI))
    {
        gmx_incons("update_coords called for velocity without VV integrator");
    }

    int  homenr = md->homenr;

    /* Cast to real for faster code, no loss in precision (see comment above) */
    real dt     = inputrec->delta_t;

    /* We need to update the NMR restraint history when time averaging is used */
    if (state->flags & (1<<estDISRE_RM3TAV))
    {
        update_disres_history(fcd, &state->hist);
    }
    if (state->flags & (1<<estORIRE_DTAV))
    {
        update_orires_history(fcd, &state->hist);
    }

    /* ############# START The update of velocities and positions ######### */
    where();
    dump_it_all(fplog, "Before update",
                state->natoms, &state->x, &upd->xp, &state->v, f);

    int nth = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

#ifndef NDEBUG
            /* Strictly speaking, we would only need this check with SIMD
             * and for the actual SIMD width. But since the code currently
             * always adds padding for GMX_REAL_MAX_SIMD_WIDTH, we check that.
             */
            size_t homenrSimdPadded = ((homenr + GMX_REAL_MAX_SIMD_WIDTH - 1)/GMX_REAL_MAX_SIMD_WIDTH)*GMX_REAL_MAX_SIMD_WIDTH;
            GMX_ASSERT(state->x.size() >= homenrSimdPadded, "state->x needs to be padded for SIMD access");
            GMX_ASSERT(upd->xp.size()  >= homenrSimdPadded, "upd->xp needs to be padded for SIMD access");
            GMX_ASSERT(state->v.size() >= homenrSimdPadded, "state->v needs to be padded for SIMD access");
            GMX_ASSERT(f->size()       >= homenrSimdPadded, "f needs to be padded for SIMD access");
#endif

            const rvec *x_rvec  = as_rvec_array(state->x.data());
            rvec       *xp_rvec = as_rvec_array(upd->xp.data());
            rvec       *v_rvec  = as_rvec_array(state->v.data());
            const rvec *f_rvec  = as_rvec_array(f->data());

            switch (inputrec->eI)
            {
                case (eiMD):
                    do_update_md(start_th, end_th, step, dt,
                                 inputrec, md, ekind, state->box,
                                 x_rvec, xp_rvec, v_rvec, f_rvec,
                                 state->nosehoover_vxi.data(), M);
                    break;
                case (eiSD1):
                    /* With constraints, the SD1 update is done in 2 parts */
                    do_update_sd1(upd->sd,
                                  start_th, end_th, dt,
                                  inputrec->opts.acc, inputrec->opts.nFreeze,
                                  md->invmass, md->ptype,
                                  md->cFREEZE, md->cACC, md->cTC,
                                  x_rvec, xp_rvec, v_rvec, f_rvec,
                                  bDoConstr, TRUE,
                                  step, inputrec->ld_seed, DOMAINDECOMP(cr) ? cr->dd->gatindex : nullptr);
                    break;
                case (eiBD):
                    do_update_bd(start_th, end_th, dt,
                                 inputrec->opts.nFreeze, md->invmass, md->ptype,
                                 md->cFREEZE, md->cTC,
                                 x_rvec, xp_rvec, v_rvec, f_rvec,
                                 inputrec->bd_fric,
                                 upd->sd->bd_rf,
                                 step, inputrec->ld_seed, DOMAINDECOMP(cr) ? cr->dd->gatindex : nullptr);
                    break;
                case (eiVV):
                case (eiVVAK):
                {
                    gmx_bool bExtended = (inputrec->etc == etcNOSEHOOVER ||
                                          inputrec->epc == epcPARRINELLORAHMAN ||
                                          inputrec->epc == epcMTTK);

                    /* assuming barostat coupled to group 0 */
                    real alpha = 1.0 + DIM/static_cast<real>(inputrec->opts.nrdf[0]);
                    switch (UpdatePart)
                    {
                        case etrtVELOCITY1:
                        case etrtVELOCITY2:
                            do_update_vv_vel(start_th, end_th, dt,
                                             inputrec->opts.acc, inputrec->opts.nFreeze,
                                             md->invmass, md->ptype,
                                             md->cFREEZE, md->cACC,
                                             v_rvec, f_rvec,
                                             bExtended, state->veta, alpha);
                            break;
                        case etrtPOSITION:
                            do_update_vv_pos(start_th, end_th, dt,
                                             inputrec->opts.nFreeze,
                                             md->ptype, md->cFREEZE,
                                             x_rvec, xp_rvec, v_rvec,
                                             bExtended, state->veta);
                            break;
                    }
                    break;
                }
                default:
                    gmx_fatal(FARGS, "Don't know how to update coordinates");
                    break;
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

}


void correct_ekin(FILE *log, int start, int end, rvec v[], rvec vcm, real mass[],
                  real tmass, tensor ekin)
{
    /*
     * This is a debugging routine. It should not be called for production code
     *
     * The kinetic energy should calculated according to:
     *   Ekin = 1/2 m (v-vcm)^2
     * However the correction is not always applied, since vcm may not be
     * known in time and we compute
     *   Ekin' = 1/2 m v^2 instead
     * This can be corrected afterwards by computing
     *   Ekin = Ekin' + 1/2 m ( -2 v vcm + vcm^2)
     * or in hsorthand:
     *   Ekin = Ekin' - m v vcm + 1/2 m vcm^2
     */
    int    i, j, k;
    real   m, tm;
    rvec   hvcm, mv;
    tensor dekin;

    /* Local particles */
    clear_rvec(mv);

    /* Processor dependent part. */
    tm = 0;
    for (i = start; (i < end); i++)
    {
        m      = mass[i];
        tm    += m;
        for (j = 0; (j < DIM); j++)
        {
            mv[j] += m*v[i][j];
        }
    }
    /* Shortcut */
    svmul(1/tmass, vcm, vcm);
    svmul(0.5, vcm, hvcm);
    clear_mat(dekin);
    for (j = 0; (j < DIM); j++)
    {
        for (k = 0; (k < DIM); k++)
        {
            dekin[j][k] += vcm[k]*(tm*hvcm[j]-mv[j]);
        }
    }
    pr_rvecs(log, 0, "dekin", dekin, DIM);
    pr_rvecs(log, 0, " ekin", ekin, DIM);
    fprintf(log, "dekin = %g, ekin = %g  vcm = (%8.4f %8.4f %8.4f)\n",
            trace(dekin), trace(ekin), vcm[XX], vcm[YY], vcm[ZZ]);
    fprintf(log, "mv = (%8.4f %8.4f %8.4f)\n",
            mv[XX], mv[YY], mv[ZZ]);
}

extern gmx_bool update_randomize_velocities(t_inputrec *ir, gmx_int64_t step, const t_commrec *cr,
                                            t_mdatoms *md, t_state *state, gmx_update_t *upd, gmx_constr_t constr)
{

    real rate = (ir->delta_t)/ir->opts.tau_t[0];

    if (ir->etc == etcANDERSEN && constr != nullptr)
    {
        /* Currently, Andersen thermostat does not support constrained
           systems. Functionality exists in the andersen_tcoupl
           function in GROMACS 4.5.7 to allow this combination. That
           code could be ported to the current random-number
           generation approach, but has not yet been done because of
           lack of time and resources. */
        gmx_fatal(FARGS, "Normal Andersen is currently not supported with constraints, use massive Andersen instead");
    }

    /* proceed with andersen if 1) it's fixed probability per
       particle andersen or 2) it's massive andersen and it's tau_t/dt */
    if ((ir->etc == etcANDERSEN) || do_per_step(step, static_cast<int>(1.0/rate + 0.5)))
    {
        andersen_tcoupl(ir, step, cr, md, state, rate,
                        upd->sd->randomize_group, upd->sd->boltzfac);
        return TRUE;
    }
    return FALSE;
}
