/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
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

    explicit gmx_stochd_t(const t_inputrec& inputRecord);
};

//! pImpled implementation for Update
class Update::Impl
{
public:
    //! Constructor
    Impl(const t_inputrec& inputRecord, BoxDeformation* boxDeformation);
    //! Destructor
    ~Impl() = default;

    void update_coords(const t_inputrec&                                inputRecord,
                       int64_t                                          step,
                       int                                              homenr,
                       bool                                             havePartiallyFrozenAtoms,
                       gmx::ArrayRef<const ParticleType>                ptype,
                       gmx::ArrayRef<const real>                        invMass,
                       gmx::ArrayRef<const rvec>                        invMassPerDim,
                       t_state*                                         state,
                       const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                       t_fcdata*                                        fcdata,
                       const gmx_ekindata_t*                            ekind,
                       const matrix                                     M,
                       int                                              UpdatePart,
                       const t_commrec*                                 cr,
                       bool                                             haveConstraints);

    void finish_update(const t_inputrec&                   inputRecord,
                       bool                                havePartiallyFrozenAtoms,
                       int                                 homenr,
                       gmx::ArrayRef<const unsigned short> cFREEZE,
                       t_state*                            state,
                       gmx_wallcycle*                      wcycle,
                       bool                                haveConstraints);

    void update_sd_second_half(const t_inputrec&                 inputRecord,
                               int64_t                           step,
                               real*                             dvdlambda,
                               int                               homenr,
                               gmx::ArrayRef<const ParticleType> ptype,
                               gmx::ArrayRef<const real>         invMass,
                               t_state*                          state,
                               const t_commrec*                  cr,
                               t_nrnb*                           nrnb,
                               gmx_wallcycle*                    wcycle,
                               gmx::Constraints*                 constr,
                               bool                              do_log,
                               bool                              do_ene);

    void update_for_constraint_virial(const t_inputrec&         inputRecord,
                                      int                       homenr,
                                      bool                      havePartiallyFrozenAtoms,
                                      gmx::ArrayRef<const real> invmass,
                                      gmx::ArrayRef<const rvec> invMassPerDim,
                                      const t_state&            state,
                                      const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                      const gmx_ekindata_t&                            ekind);

    void update_temperature_constants(const t_inputrec& inputRecord);

    const std::vector<bool>& getAndersenRandomizeGroup() const { return sd_.randomize_group; }

    const std::vector<real>& getBoltzmanFactor() const { return sd_.boltzfac; }

    PaddedVector<RVec>* xp() { return &xp_; }

    BoxDeformation* deform() const { return deform_; }

    //! Group index for freezing
    gmx::ArrayRef<const unsigned short> cFREEZE_;
    //! Group index for temperature coupling
    gmx::ArrayRef<const unsigned short> cTC_;
    //! Group index for accleration groups
    gmx::ArrayRef<const unsigned short> cAcceleration_;

private:
    //! stochastic dynamics struct
    gmx_stochd_t sd_;
    //! xprime for constraint algorithms
    PaddedVector<RVec> xp_;
    //! Box deformation handler (or nullptr if inactive).
    BoxDeformation* deform_ = nullptr;
};

Update::Update(const t_inputrec& inputRecord, BoxDeformation* boxDeformation) :
    impl_(new Impl(inputRecord, boxDeformation)){};

Update::~Update(){};

const std::vector<bool>& Update::getAndersenRandomizeGroup() const
{
    return impl_->getAndersenRandomizeGroup();
}

const std::vector<real>& Update::getBoltzmanFactor() const
{
    return impl_->getBoltzmanFactor();
}

PaddedVector<RVec>* Update::xp()
{
    return impl_->xp();
}

BoxDeformation* Update::deform() const
{
    return impl_->deform();
}

void Update::update_coords(const t_inputrec&                 inputRecord,
                           int64_t                           step,
                           const int                         homenr,
                           const bool                        havePartiallyFrozenAtoms,
                           gmx::ArrayRef<const ParticleType> ptype,
                           gmx::ArrayRef<const real>         invMass,
                           gmx::ArrayRef<const rvec>         invMassPerDim,
                           t_state*                          state,
                           const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                           t_fcdata*                                        fcdata,
                           const gmx_ekindata_t*                            ekind,
                           const matrix                                     M,
                           int                                              updatePart,
                           const t_commrec*                                 cr,
                           const bool                                       haveConstraints)
{
    return impl_->update_coords(inputRecord,
                                step,
                                homenr,
                                havePartiallyFrozenAtoms,
                                ptype,
                                invMass,
                                invMassPerDim,
                                state,
                                f,
                                fcdata,
                                ekind,
                                M,
                                updatePart,
                                cr,
                                haveConstraints);
}

void Update::finish_update(const t_inputrec& inputRecord,
                           const bool        havePartiallyFrozenAtoms,
                           const int         homenr,
                           t_state*          state,
                           gmx_wallcycle*    wcycle,
                           const bool        haveConstraints)
{
    return impl_->finish_update(
            inputRecord, havePartiallyFrozenAtoms, homenr, impl_->cFREEZE_, state, wcycle, haveConstraints);
}

void Update::update_sd_second_half(const t_inputrec&                 inputRecord,
                                   int64_t                           step,
                                   real*                             dvdlambda,
                                   const int                         homenr,
                                   gmx::ArrayRef<const ParticleType> ptype,
                                   gmx::ArrayRef<const real>         invMass,
                                   t_state*                          state,
                                   const t_commrec*                  cr,
                                   t_nrnb*                           nrnb,
                                   gmx_wallcycle*                    wcycle,
                                   gmx::Constraints*                 constr,
                                   bool                              do_log,
                                   bool                              do_ene)
{
    return impl_->update_sd_second_half(
            inputRecord, step, dvdlambda, homenr, ptype, invMass, state, cr, nrnb, wcycle, constr, do_log, do_ene);
}

void Update::update_for_constraint_virial(const t_inputrec&         inputRecord,
                                          const int                 homenr,
                                          const bool                havePartiallyFrozenAtoms,
                                          gmx::ArrayRef<const real> invmass,
                                          gmx::ArrayRef<const rvec> invMassPerDim,
                                          const t_state&            state,
                                          const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                          const gmx_ekindata_t&                            ekind)
{
    return impl_->update_for_constraint_virial(
            inputRecord, homenr, havePartiallyFrozenAtoms, invmass, invMassPerDim, state, f, ekind);
}

void Update::update_temperature_constants(const t_inputrec& inputRecord)
{
    return impl_->update_temperature_constants(inputRecord);
}

/*! \brief Sets whether we store the updated velocities */
enum class StoreUpdatedVelocities
{
    yes, //!< Store the updated velocities
    no   //!< Do not store the updated velocities
};

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
 * \tparam       storeUpdatedVelocities Tells whether we should store the updated velocities
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
 * \param[inout] v                      Velocities, type either rvec* or const rvec*
 * \param[in]    f                      Forces
 *
 * We expect this template to get good SIMD acceleration by most compilers,
 * unlike the more complex general template.
 * Note that we might get even better SIMD acceleration when we introduce
 * aligned (and padded) memory, possibly with some hints for the compilers.
 */
template<StoreUpdatedVelocities storeUpdatedVelocities, NumTempScaleValues numTempScaleValues, ApplyParrinelloRahmanVScaling applyPRVScaling, typename VelocityType>
static std::enable_if_t<std::is_same<VelocityType, rvec>::value || std::is_same<VelocityType, const rvec>::value, void>
updateMDLeapfrogSimple(int                                 start,
                       int                                 nrend,
                       real                                dt,
                       real                                dtPressureCouple,
                       gmx::ArrayRef<const rvec>           invMassPerDim,
                       gmx::ArrayRef<const t_grp_tcstat>   tcstat,
                       gmx::ArrayRef<const unsigned short> cTC,
                       const rvec                          pRVScaleMatrixDiagonal,
                       const rvec* gmx_restrict            x,
                       rvec* gmx_restrict                  xprime,
                       VelocityType* gmx_restrict          v,
                       const rvec* gmx_restrict            f)
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
            if constexpr (storeUpdatedVelocities == StoreUpdatedVelocities::yes)
            {
                v[a][d] = vNew;
            }
            // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
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
 * \tparam       storeUpdatedVelocities Tells whether we should store the updated velocities
 * \param[in]    start                  Index of first atom to update
 * \param[in]    nrend                  Last atom to update: \p nrend - 1
 * \param[in]    dt                     The time step
 * \param[in]    invMass                1/mass per atom
 * \param[in]    tcstat                 Temperature coupling information
 * \param[in]    x                      Input coordinates
 * \param[out]   xprime                 Updated coordinates
 * \param[inout] v                      Velocities, type either rvec* or const rvec*
 * \param[in]    f                      Forces
 */
template<StoreUpdatedVelocities storeUpdatedVelocities, typename VelocityType>
static std::enable_if_t<std::is_same<VelocityType, rvec>::value || std::is_same<VelocityType, const rvec>::value, void>
updateMDLeapfrogSimpleSimd(int                               start,
                           int                               nrend,
                           real                              dt,
                           gmx::ArrayRef<const real>         invMass,
                           gmx::ArrayRef<const t_grp_tcstat> tcstat,
                           const rvec* gmx_restrict          x,
                           rvec* gmx_restrict                xprime,
                           VelocityType* gmx_restrict        v,
                           const rvec* gmx_restrict          f)
{
    SimdReal timestep(dt);
    SimdReal lambdaSystem(tcstat[0].lambda);

    /* We declare variables here, since code is often slower when declaring them inside the loop */

    /* Note: We should implement a proper PaddedVector, so we don't need this check */
    GMX_ASSERT(isSimdAligned(invMass.data()), "invMass should be aligned");

    for (int a = start; a < nrend; a += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal invMass0, invMass1, invMass2;
        expandScalarsToTriplets(simdLoad(invMass.data() + a), &invMass0, &invMass1, &invMass2);

        SimdReal v0, v1, v2;
        SimdReal f0, f1, f2;
        simdLoadRvecs(v, a, &v0, &v1, &v2);
        simdLoadRvecs(f, a, &f0, &f1, &f2);

        v0 = fma(f0 * invMass0, timestep, lambdaSystem * v0);
        v1 = fma(f1 * invMass1, timestep, lambdaSystem * v1);
        v2 = fma(f2 * invMass2, timestep, lambdaSystem * v2);

        if constexpr (storeUpdatedVelocities == StoreUpdatedVelocities::yes)
        {
            simdStoreRvecs(v, a, v0, v1, v2);
        }

        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
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

/*! \brief Integrate using leap-frog with support for everything.
 *
 * \tparam        accelerationType  Type of NEMD acceleration.
 * \param[in]     start             Index of first atom to update.
 * \param[in]     nrend             Last atom to update: \p nrend - 1.
 * \param[in]     doNoseHoover      If to apply Nose-Hoover T-coupling.
 * \param[in]     dt                The time step.
 * \param[in]     dtPressureCouple  Time step for pressure coupling, is 0 when no pressure
 *                                  coupling should be applied at this step.
 * \param[in]     cTC               Temperature coupling group indices
 * \param[in]     cAcceleration     Acceleration group indices
 * \param[in]     acceleration      Acceleration per group.
 * \param[in]     invMassPerDim     Inverse mass per dimension
 * \param[in]     ekind             Kinetic energy data.
 * \param[in]     box               The box dimensions.
 * \param[in]     x                 Input coordinates.
 * \param[out]    xprime            Updated coordinates.
 * \param[inout]  v                 Velocities.
 * \param[in]     f                 Forces.
 * \param[in]     nh_vxi            Nose-Hoover velocity scaling factors.
 * \param[in]     nsttcouple        Frequency of the temperature coupling steps.
 * \param[in]     M                 Parrinello-Rahman scaling matrix.
 */
template<AccelerationType accelerationType>
static void updateMDLeapfrogGeneral(int                                 start,
                                    int                                 nrend,
                                    bool                                doNoseHoover,
                                    real                                dt,
                                    real                                dtPressureCouple,
                                    gmx::ArrayRef<const unsigned short> cTC,
                                    gmx::ArrayRef<const unsigned short> cAcceleration,
                                    const rvec* gmx_restrict            acceleration,
                                    gmx::ArrayRef<const rvec>           invMassPerDim,
                                    const gmx_ekindata_t*               ekind,
                                    const matrix                        box,
                                    const rvec* gmx_restrict            x,
                                    rvec* gmx_restrict                  xprime,
                                    rvec* gmx_restrict                  v,
                                    const rvec* gmx_restrict            f,
                                    const double* gmx_restrict          nh_vxi,
                                    const int                           nsttcouple,
                                    const matrix                        M)
{
    /* This is a version of the leap-frog integrator that supports
     * all combinations of T-coupling, P-coupling and NEMD.
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */

    gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind->tcstat;

    /* Initialize group values, changed later when multiple groups are used */
    int ga = 0;
    int gt = 0;

    real omega_Z = 2 * static_cast<real>(M_PI) / box[ZZ][ZZ];

    for (int n = start; n < nrend; n++)
    {
        if (!cTC.empty())
        {
            gt = cTC[n];
        }
        real lg = tcstat[gt].lambda;

        rvec vRel;
        real cosineZ, vCosine;
        switch (accelerationType)
        {
            case AccelerationType::none: copy_rvec(v[n], vRel); break;
            case AccelerationType::group:
                if (!cAcceleration.empty())
                {
                    ga = cAcceleration[n];
                }
                /* With constant acceleration we do scale the velocity of the accelerated groups */
                copy_rvec(v[n], vRel);
                break;
            case AccelerationType::cosine:
                cosineZ = std::cos(x[n][ZZ] * omega_Z);
                vCosine = cosineZ * ekind->cosacc.vcos;
                /* Avoid scaling the cosine profile velocity */
                copy_rvec(v[n], vRel);
                vRel[XX] -= vCosine;
                break;
        }

        real factorNH = 0.0;
        if (doNoseHoover)
        {
            /* Here we account for multiple time stepping, by increasing
             * the Nose-Hoover correction by nsttcouple
             * TODO: This can be pre-computed.
             */
            factorNH = 0.5 * nsttcouple * dt * nh_vxi[gt];
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
                    /* Apply the constant acceleration */
                    vNew += acceleration[ga][d] * dt;
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
static void do_update_md(int                                  start,
                         int                                  nrend,
                         real                                 dt,
                         int64_t                              step,
                         const rvec* gmx_restrict             x,
                         rvec* gmx_restrict                   xprime,
                         rvec* gmx_restrict                   v,
                         const rvec* gmx_restrict             f,
                         const TemperatureCoupling            etc,
                         const PressureCoupling               epc,
                         const int                            nsttcouple,
                         const int                            nstpcouple,
                         gmx::ArrayRef<const unsigned short>  cTC,
                         const bool                           useConstantAcceleration,
                         gmx::ArrayRef<const unsigned short>  cAcceleration,
                         const rvec*                          acceleration,
                         gmx::ArrayRef<const real> gmx_unused invmass,
                         gmx::ArrayRef<const rvec>            invMassPerDim,
                         const gmx_ekindata_t*                ekind,
                         const matrix                         box,
                         const double* gmx_restrict           nh_vxi,
                         const matrix                         M,
                         bool gmx_unused                      havePartiallyFrozenAtoms)
{
    GMX_ASSERT(nrend == start || xprime != x,
               "For SIMD optimization certain compilers need to have xprime != x");

    /* Note: Berendsen pressure scaling is handled after do_update_md() */
    bool doTempCouple =
            (etc != TemperatureCoupling::No && do_per_step(step + nsttcouple - 1, nsttcouple));
    bool doNoseHoover       = (etc == TemperatureCoupling::NoseHoover && doTempCouple);
    bool doParrinelloRahman = (epc == PressureCoupling::ParrinelloRahman
                               && do_per_step(step + nstpcouple - 1, nstpcouple));
    bool doPROffDiagonal = (doParrinelloRahman && (M[YY][XX] != 0 || M[ZZ][XX] != 0 || M[ZZ][YY] != 0));

    real dtPressureCouple = (doParrinelloRahman ? nstpcouple * dt : 0);

    /* NEMD (also cosine) acceleration is applied in updateMDLeapFrogGeneral */
    const bool doAcceleration = (useConstantAcceleration || ekind->cosacc.cos_accel != 0);

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
            updateMDLeapfrogGeneral<AccelerationType::none>(start,
                                                            nrend,
                                                            doNoseHoover,
                                                            dt,
                                                            dtPressureCouple,
                                                            cTC,
                                                            cAcceleration,
                                                            acceleration,
                                                            invMassPerDim,
                                                            ekind,
                                                            box,
                                                            x,
                                                            xprime,
                                                            v,
                                                            f,
                                                            nh_vxi,
                                                            nsttcouple,
                                                            stepM);
        }
        else if (useConstantAcceleration)
        {
            updateMDLeapfrogGeneral<AccelerationType::group>(start,
                                                             nrend,
                                                             doNoseHoover,
                                                             dt,
                                                             dtPressureCouple,
                                                             cTC,
                                                             cAcceleration,
                                                             acceleration,
                                                             invMassPerDim,
                                                             ekind,
                                                             box,
                                                             x,
                                                             xprime,
                                                             v,
                                                             f,
                                                             nh_vxi,
                                                             nsttcouple,
                                                             stepM);
        }
        else
        {
            updateMDLeapfrogGeneral<AccelerationType::cosine>(start,
                                                              nrend,
                                                              doNoseHoover,
                                                              dt,
                                                              dtPressureCouple,
                                                              cTC,
                                                              cAcceleration,
                                                              acceleration,
                                                              invMassPerDim,
                                                              ekind,
                                                              box,
                                                              x,
                                                              xprime,
                                                              v,
                                                              f,
                                                              nh_vxi,
                                                              nsttcouple,
                                                              stepM);
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
        gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind->tcstat;

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
                updateMDLeapfrogSimple<StoreUpdatedVelocities::yes, NumTempScaleValues::single, ApplyParrinelloRahmanVScaling::diagonal>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, diagM, x, xprime, v, f);
            }
            else
            {
                updateMDLeapfrogSimple<StoreUpdatedVelocities::yes, NumTempScaleValues::multiple, ApplyParrinelloRahmanVScaling::diagonal>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, diagM, x, xprime, v, f);
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
                if (!havePartiallyFrozenAtoms)
                {
                    updateMDLeapfrogSimpleSimd<StoreUpdatedVelocities::yes>(
                            start, nrend, dt, invmass, tcstat, x, xprime, v, f);
                }
                else
#endif
                {
                    updateMDLeapfrogSimple<StoreUpdatedVelocities::yes, NumTempScaleValues::single, ApplyParrinelloRahmanVScaling::no>(
                            start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, nullptr, x, xprime, v, f);
                }
            }
            else
            {
                updateMDLeapfrogSimple<StoreUpdatedVelocities::yes, NumTempScaleValues::multiple, ApplyParrinelloRahmanVScaling::no>(
                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, nullptr, x, xprime, v, f);
            }
        }
    }
}
/*! \brief Handles the Leap-frog MD x and v integration */
static void doUpdateMDDoNotUpdateVelocities(int                      start,
                                            int                      nrend,
                                            real                     dt,
                                            const rvec* gmx_restrict x,
                                            rvec* gmx_restrict       xprime,
                                            const rvec* gmx_restrict v,
                                            const rvec* gmx_restrict f,
                                            bool gmx_unused          havePartiallyFrozenAtoms,
                                            gmx::ArrayRef<const real> gmx_unused invmass,
                                            gmx::ArrayRef<const rvec>            invMassPerDim,
                                            const gmx_ekindata_t&                ekind)
{
    GMX_ASSERT(nrend == start || xprime != x,
               "For SIMD optimization certain compilers need to have xprime != x");

    gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind.tcstat;

    /* Check if we can use invmass instead of invMassPerDim */
#if GMX_HAVE_SIMD_UPDATE
    if (!havePartiallyFrozenAtoms)
    {
        updateMDLeapfrogSimpleSimd<StoreUpdatedVelocities::no>(
                start, nrend, dt, invmass, tcstat, x, xprime, v, f);
    }
    else
#endif
    {
        updateMDLeapfrogSimple<StoreUpdatedVelocities::no, NumTempScaleValues::single, ApplyParrinelloRahmanVScaling::no>(
                start, nrend, dt, dt, invMassPerDim, tcstat, gmx::ArrayRef<const unsigned short>(), nullptr, x, xprime, v, f);
    }
}

static void do_update_vv_vel(int                                 start,
                             int                                 nrend,
                             real                                dt,
                             gmx::ArrayRef<const ivec>           nFreeze,
                             gmx::ArrayRef<const unsigned short> cAcceleration,
                             const rvec*                         acceleration,
                             gmx::ArrayRef<const real>           invmass,
                             gmx::ArrayRef<const ParticleType>   ptype,
                             gmx::ArrayRef<const unsigned short> cFREEZE,
                             rvec                                v[],
                             const rvec                          f[],
                             gmx_bool                            bExtended,
                             real                                veta,
                             real                                alpha)
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
        if (!cFREEZE.empty())
        {
            gf = cFREEZE[n];
        }
        if (!cAcceleration.empty())
        {
            ga = cAcceleration[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != ParticleType::Shell) && !nFreeze[gf][d])
            {
                v[n][d] = mv1 * (mv1 * v[n][d] + 0.5 * (w_dt * mv2 * f[n][d]))
                          + 0.5 * acceleration[ga][d] * dt;
            }
            else
            {
                v[n][d] = 0.0;
            }
        }
    }
} /* do_update_vv_vel */

static void do_update_vv_pos(int                                 start,
                             int                                 nrend,
                             real                                dt,
                             gmx::ArrayRef<const ivec>           nFreeze,
                             gmx::ArrayRef<const ParticleType>   ptype,
                             gmx::ArrayRef<const unsigned short> cFREEZE,
                             const rvec                          x[],
                             rvec                                xprime[],
                             const rvec                          v[],
                             gmx_bool                            bExtended,
                             real                                veta)
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

        if (!cFREEZE.empty())
        {
            gf = cFREEZE[n];
        }

        for (d = 0; d < DIM; d++)
        {
            if ((ptype[n] != ParticleType::Shell) && !nFreeze[gf][d])
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

gmx_stochd_t::gmx_stochd_t(const t_inputrec& inputRecord)
{
    const t_grpopts* opts = &inputRecord.opts;
    const int        ngtc = opts->ngtc;

    if (inputRecord.eI == IntegrationAlgorithm::BD)
    {
        bd_rf.resize(ngtc);
    }
    else if (EI_SD(inputRecord.eI))
    {
        sdc.resize(ngtc);
        sdsig.resize(ngtc);

        for (int gt = 0; gt < ngtc; gt++)
        {
            if (opts->tau_t[gt] > 0)
            {
                sdc[gt].em = std::exp(-inputRecord.delta_t / opts->tau_t[gt]);
            }
            else
            {
                /* No friction and noise on this group */
                sdc[gt].em = 1;
            }
        }
    }
    else if (ETC_ANDERSEN(inputRecord.etc))
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
                boltzfac[gt]        = gmx::c_boltz * opts->ref_t[gt];
            }
            else
            {
                randomize_group[gt] = false;
            }
        }
    }
}

void Update::Impl::update_temperature_constants(const t_inputrec& inputRecord)
{
    if (inputRecord.eI == IntegrationAlgorithm::BD)
    {
        if (inputRecord.bd_fric != 0)
        {
            for (int gt = 0; gt < inputRecord.opts.ngtc; gt++)
            {
                sd_.bd_rf[gt] = std::sqrt(2.0 * gmx::c_boltz * inputRecord.opts.ref_t[gt]
                                          / (inputRecord.bd_fric * inputRecord.delta_t));
            }
        }
        else
        {
            for (int gt = 0; gt < inputRecord.opts.ngtc; gt++)
            {
                sd_.bd_rf[gt] = std::sqrt(2.0 * gmx::c_boltz * inputRecord.opts.ref_t[gt]);
            }
        }
    }
    if (inputRecord.eI == IntegrationAlgorithm::SD1)
    {
        for (int gt = 0; gt < inputRecord.opts.ngtc; gt++)
        {
            real kT = gmx::c_boltz * inputRecord.opts.ref_t[gt];
            /* The mass is accounted for later, since this differs per atom */
            sd_.sdsig[gt].V = std::sqrt(kT * (1 - sd_.sdc[gt].em * sd_.sdc[gt].em));
        }
    }
}

Update::Impl::Impl(const t_inputrec& inputRecord, BoxDeformation* boxDeformation) :
    sd_(inputRecord), deform_(boxDeformation)
{
    update_temperature_constants(inputRecord);
    xp_.resizeWithPadding(0);
}

void Update::updateAfterPartition(int                                 numAtoms,
                                  gmx::ArrayRef<const unsigned short> cFREEZE,
                                  gmx::ArrayRef<const unsigned short> cTC,
                                  gmx::ArrayRef<const unsigned short> cAcceleration)
{

    impl_->xp()->resizeWithPadding(numAtoms);
    impl_->cFREEZE_       = cFREEZE;
    impl_->cTC_           = cTC;
    impl_->cAcceleration_ = cAcceleration;
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
static void doSDUpdateGeneral(const gmx_stochd_t&                 sd,
                              int                                 start,
                              int                                 nrend,
                              real                                dt,
                              gmx::ArrayRef<const ivec>           nFreeze,
                              gmx::ArrayRef<const real>           invmass,
                              gmx::ArrayRef<const ParticleType>   ptype,
                              gmx::ArrayRef<const unsigned short> cFREEZE,
                              gmx::ArrayRef<const unsigned short> cTC,
                              gmx::ArrayRef<const unsigned short> cAcceleration,
                              const rvec*                         acceleration,
                              const rvec                          x[],
                              rvec                                xprime[],
                              rvec                                v[],
                              const rvec                          f[],
                              int64_t                             step,
                              int                                 seed,
                              const int*                          gatindex)
{
    // cTC, cACC and cFREEZE can be nullptr any time, but various
    // instantiations do not make sense with particular pointer
    // values.
    if (updateType == SDUpdate::ForcesOnly)
    {
        GMX_ASSERT(f != nullptr, "SD update with only forces requires forces");
        GMX_ASSERT(cTC.empty(), "SD update with only forces cannot handle temperature groups");
    }
    if (updateType == SDUpdate::FrictionAndNoiseOnly)
    {
        GMX_ASSERT(f == nullptr, "SD update with only noise cannot handle forces");
        GMX_ASSERT(cAcceleration.empty(),
                   "SD update with only noise cannot handle acceleration groups");
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

        int freezeGroup       = !cFREEZE.empty() ? cFREEZE[n] : 0;
        int accelerationGroup = !cAcceleration.empty() ? cAcceleration[n] : 0;
        int temperatureGroup  = !cTC.empty() ? cTC[n] : 0;

        for (int d = 0; d < DIM; d++)
        {
            if ((ptype[n] != ParticleType::Shell) && !nFreeze[freezeGroup][d])
            {
                if (updateType == SDUpdate::ForcesOnly)
                {
                    real vn = v[n][d] + (inverseMass * f[n][d] + acceleration[accelerationGroup][d]) * dt;
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
                    real vn = v[n][d] + (inverseMass * f[n][d] + acceleration[accelerationGroup][d]) * dt;
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

static void do_update_sd(int                                 start,
                         int                                 nrend,
                         real                                dt,
                         int64_t                             step,
                         const rvec* gmx_restrict            x,
                         rvec* gmx_restrict                  xprime,
                         rvec* gmx_restrict                  v,
                         const rvec* gmx_restrict            f,
                         gmx::ArrayRef<const ivec>           nFreeze,
                         gmx::ArrayRef<const real>           invmass,
                         gmx::ArrayRef<const ParticleType>   ptype,
                         gmx::ArrayRef<const unsigned short> cFREEZE,
                         gmx::ArrayRef<const unsigned short> cTC,
                         gmx::ArrayRef<const unsigned short> cAcceleration,
                         const rvec*                         acceleration,
                         int                                 seed,
                         const t_commrec*                    cr,
                         const gmx_stochd_t&                 sd,
                         bool                                haveConstraints)
{
    if (haveConstraints)
    {
        // With constraints, the SD update is done in 2 parts
        doSDUpdateGeneral<SDUpdate::ForcesOnly>(sd,
                                                start,
                                                nrend,
                                                dt,
                                                nFreeze,
                                                invmass,
                                                ptype,
                                                cFREEZE,
                                                gmx::ArrayRef<const unsigned short>(),
                                                cAcceleration,
                                                acceleration,
                                                x,
                                                xprime,
                                                v,
                                                f,
                                                step,
                                                seed,
                                                nullptr);
    }
    else
    {
        doSDUpdateGeneral<SDUpdate::Combined>(
                sd,
                start,
                nrend,
                dt,
                nFreeze,
                invmass,
                ptype,
                cFREEZE,
                cTC,
                cAcceleration,
                acceleration,
                x,
                xprime,
                v,
                f,
                step,
                seed,
                haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr);
    }
}

static void do_update_bd(int                                 start,
                         int                                 nrend,
                         real                                dt,
                         int64_t                             step,
                         const rvec* gmx_restrict            x,
                         rvec* gmx_restrict                  xprime,
                         rvec* gmx_restrict                  v,
                         const rvec* gmx_restrict            f,
                         gmx::ArrayRef<const ivec>           nFreeze,
                         gmx::ArrayRef<const real>           invmass,
                         gmx::ArrayRef<const ParticleType>   ptype,
                         gmx::ArrayRef<const unsigned short> cFREEZE,
                         gmx::ArrayRef<const unsigned short> cTC,
                         real                                friction_coefficient,
                         const real*                         rf,
                         int                                 seed,
                         const int*                          gatindex)
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

        if (!cFREEZE.empty())
        {
            gf = cFREEZE[n];
        }
        if (!cTC.empty())
        {
            gt = cTC[n];
        }
        for (d = 0; (d < DIM); d++)
        {
            if ((ptype[n] != ParticleType::Shell) && !nFreeze[gf][d])
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
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinh[0][0]),
                      ekind->tcstat[i].ekinh[0],
                      cr->mpi_comm_mygroup);
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinf[0][0]),
                      ekind->tcstat[i].ekinf[0],
                      cr->mpi_comm_mygroup);
            gmx_bcast(DIM * DIM * sizeof(ekind->tcstat[i].ekinh_old[0][0]),
                      ekind->tcstat[i].ekinh_old[0],
                      cr->mpi_comm_mygroup);

            gmx_bcast(sizeof(ekind->tcstat[i].ekinscalef_nhc),
                      &(ekind->tcstat[i].ekinscalef_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].ekinscaleh_nhc),
                      &(ekind->tcstat[i].ekinscaleh_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].vscale_nhc), &(ekind->tcstat[i].vscale_nhc), cr->mpi_comm_mygroup);
        }
        gmx_bcast(DIM * DIM * sizeof(ekind->ekin[0][0]), ekind->ekin[0], cr->mpi_comm_mygroup);

        gmx_bcast(sizeof(ekind->dekindl), &ekind->dekindl, cr->mpi_comm_mygroup);
        gmx_bcast(sizeof(ekind->cosacc.mvcos), &ekind->cosacc.mvcos, cr->mpi_comm_mygroup);
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

void Update::Impl::update_sd_second_half(const t_inputrec&                 inputRecord,
                                         int64_t                           step,
                                         real*                             dvdlambda,
                                         int                               homenr,
                                         gmx::ArrayRef<const ParticleType> ptype,
                                         gmx::ArrayRef<const real>         invMass,
                                         t_state*                          state,
                                         const t_commrec*                  cr,
                                         t_nrnb*                           nrnb,
                                         gmx_wallcycle*                    wcycle,
                                         gmx::Constraints*                 constr,
                                         bool                              do_log,
                                         bool                              do_ene)
{
    if (!constr)
    {
        return;
    }
    if (inputRecord.eI == IntegrationAlgorithm::SD1)
    {

        /* Cast delta_t from double to real to make the integrators faster.
         * The only reason for having delta_t double is to get accurate values
         * for t=delta_t*step when step is larger than float precision.
         * For integration dt the accuracy of real suffices, since with
         * integral += dt*integrand the increment is nearly always (much) smaller
         * than the integral (and the integrand has real precision).
         */
        real dt = inputRecord.delta_t;

        wallcycle_start(wcycle, WallCycleCounter::Update);

        int nth = gmx_omp_nthreads_get(ModuleMultiThread::Update);

#pragma omp parallel for num_threads(nth) schedule(static)
        for (int th = 0; th < nth; th++)
        {
            try
            {
                int start_th, end_th;
                getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

                doSDUpdateGeneral<SDUpdate::FrictionAndNoiseOnly>(
                        sd_,
                        start_th,
                        end_th,
                        dt,
                        gmx::arrayRefFromArray(inputRecord.opts.nFreeze, inputRecord.opts.ngfrz),
                        invMass,
                        ptype,
                        cFREEZE_,
                        cTC_,
                        cAcceleration_,
                        inputRecord.opts.acceleration,
                        state->x.rvec_array(),
                        xp_.rvec_array(),
                        state->v.rvec_array(),
                        nullptr,
                        step,
                        inputRecord.ld_seed,
                        haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        inc_nrnb(nrnb, eNR_UPDATE, homenr);
        wallcycle_stop(wcycle, WallCycleCounter::Update);

        /* Constrain the coordinates upd->xp for half a time step */
        bool computeVirial = false;
        constr->apply(do_log,
                      do_ene,
                      step,
                      1,
                      0.5,
                      state->x.arrayRefWithPadding(),
                      xp_.arrayRefWithPadding(),
                      ArrayRef<RVec>(),
                      state->box,
                      state->lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Bonded)],
                      dvdlambda,
                      state->v.arrayRefWithPadding(),
                      computeVirial,
                      nullptr,
                      ConstraintVariable::Positions);
    }
}

void Update::Impl::finish_update(const t_inputrec&                   inputRecord,
                                 const bool                          havePartiallyFrozenAtoms,
                                 const int                           homenr,
                                 gmx::ArrayRef<const unsigned short> cFREEZE,
                                 t_state*                            state,
                                 gmx_wallcycle*                      wcycle,
                                 const bool                          haveConstraints)
{
    /* NOTE: Currently we always integrate to a temporary buffer and
     * then copy the results back here.
     */

    wallcycle_start_nocount(wcycle, WallCycleCounter::Update);

    auto xp = makeConstArrayRef(xp_).subArray(0, homenr);
    auto x  = makeArrayRef(state->x).subArray(0, homenr);

    if (havePartiallyFrozenAtoms && haveConstraints)
    {
        /* We have atoms that are frozen along some, but not all dimensions,
         * then constraints will have moved them also along the frozen dimensions.
         * To freeze such degrees of freedom we do not copy them back here.
         */
        const ivec* nFreeze = inputRecord.opts.nFreeze;

        for (int i = 0; i < homenr; i++)
        {
            const int g = cFREEZE[i];

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
        int gmx_unused nth = gmx_omp_nthreads_get(ModuleMultiThread::Update);
#pragma omp parallel for num_threads(nth) schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            // Trivial statement, does not throw
            x[i] = xp[i];
        }
    }

    wallcycle_stop(wcycle, WallCycleCounter::Update);
}

void Update::Impl::update_coords(const t_inputrec&                 inputRecord,
                                 int64_t                           step,
                                 int                               homenr,
                                 bool                              havePartiallyFrozenAtoms,
                                 gmx::ArrayRef<const ParticleType> ptype,
                                 gmx::ArrayRef<const real>         invMass,
                                 gmx::ArrayRef<const rvec>         invMassPerDim,
                                 t_state*                          state,
                                 const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                 t_fcdata*                                        fcdata,
                                 const gmx_ekindata_t*                            ekind,
                                 const matrix                                     M,
                                 int                                              updatePart,
                                 const t_commrec*                                 cr,
                                 const bool                                       haveConstraints)
{
    /* Running the velocity half does nothing except for velocity verlet */
    if ((updatePart == etrtVELOCITY1 || updatePart == etrtVELOCITY2) && !EI_VV(inputRecord.eI))
    {
        gmx_incons("update_coords called for velocity without VV integrator");
    }

    /* Cast to real for faster code, no loss in precision (see comment above) */
    real dt = inputRecord.delta_t;

    /* We need to update the NMR restraint history when time averaging is used */
    if (state->flags & enumValueToBitMask(StateEntry::DisreRm3Tav))
    {
        update_disres_history(*fcdata->disres, &state->hist);
    }
    if (state->flags & enumValueToBitMask(StateEntry::OrireDtav))
    {
        GMX_ASSERT(fcdata, "Need valid fcdata");
        fcdata->orires->updateHistory();
    }

    /* ############# START The update of velocities and positions ######### */
    int nth = gmx_omp_nthreads_get(ModuleMultiThread::Update);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            const rvec* x_rvec  = state->x.rvec_array();
            rvec*       xp_rvec = xp_.rvec_array();
            rvec*       v_rvec  = state->v.rvec_array();
            const rvec* f_rvec  = as_rvec_array(f.unpaddedConstArrayRef().data());

            switch (inputRecord.eI)
            {
                case (IntegrationAlgorithm::MD):
                    do_update_md(start_th,
                                 end_th,
                                 dt,
                                 step,
                                 x_rvec,
                                 xp_rvec,
                                 v_rvec,
                                 f_rvec,
                                 inputRecord.etc,
                                 inputRecord.epc,
                                 inputRecord.nsttcouple,
                                 inputRecord.nstpcouple,
                                 cTC_,
                                 inputRecord.useConstantAcceleration,
                                 cAcceleration_,
                                 inputRecord.opts.acceleration,
                                 invMass,
                                 invMassPerDim,
                                 ekind,
                                 state->box,
                                 state->nosehoover_vxi.data(),
                                 M,
                                 havePartiallyFrozenAtoms);
                    break;
                case (IntegrationAlgorithm::SD1):
                    do_update_sd(start_th,
                                 end_th,
                                 dt,
                                 step,
                                 x_rvec,
                                 xp_rvec,
                                 v_rvec,
                                 f_rvec,
                                 gmx::arrayRefFromArray(inputRecord.opts.nFreeze, inputRecord.opts.ngfrz),
                                 invMass,
                                 ptype,
                                 cFREEZE_,
                                 cTC_,
                                 cAcceleration_,
                                 inputRecord.opts.acceleration,
                                 inputRecord.ld_seed,
                                 cr,
                                 sd_,
                                 haveConstraints);
                    break;
                case (IntegrationAlgorithm::BD):
                    do_update_bd(start_th,
                                 end_th,
                                 dt,
                                 step,
                                 x_rvec,
                                 xp_rvec,
                                 v_rvec,
                                 f_rvec,
                                 gmx::arrayRefFromArray(inputRecord.opts.nFreeze, inputRecord.opts.ngfrz),
                                 invMass,
                                 ptype,
                                 cFREEZE_,
                                 cTC_,
                                 inputRecord.bd_fric,
                                 sd_.bd_rf.data(),
                                 inputRecord.ld_seed,
                                 haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr);
                    break;
                case (IntegrationAlgorithm::VV):
                case (IntegrationAlgorithm::VVAK):
                {
                    gmx_bool bExtended = (inputRecord.etc == TemperatureCoupling::NoseHoover
                                          || inputRecord.epc == PressureCoupling::ParrinelloRahman
                                          || inputRecord.epc == PressureCoupling::Mttk);

                    /* assuming barostat coupled to group 0 */
                    real alpha = 1.0 + DIM / static_cast<real>(inputRecord.opts.nrdf[0]);
                    switch (updatePart)
                    {
                        case etrtVELOCITY1:
                        case etrtVELOCITY2:
                            do_update_vv_vel(start_th,
                                             end_th,
                                             dt,
                                             gmx::arrayRefFromArray(inputRecord.opts.nFreeze,
                                                                    inputRecord.opts.ngfrz),
                                             cAcceleration_,
                                             inputRecord.opts.acceleration,
                                             invMass,
                                             ptype,
                                             cFREEZE_,
                                             v_rvec,
                                             f_rvec,
                                             bExtended,
                                             state->veta,
                                             alpha);
                            break;
                        case etrtPOSITION:
                            do_update_vv_pos(start_th,
                                             end_th,
                                             dt,
                                             gmx::arrayRefFromArray(inputRecord.opts.nFreeze,
                                                                    inputRecord.opts.ngfrz),
                                             ptype,
                                             cFREEZE_,
                                             x_rvec,
                                             xp_rvec,
                                             v_rvec,
                                             bExtended,
                                             state->veta);
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

void Update::Impl::update_for_constraint_virial(const t_inputrec&         inputRecord,
                                                int                       homenr,
                                                bool                      havePartiallyFrozenAtoms,
                                                gmx::ArrayRef<const real> invmass,
                                                gmx::ArrayRef<const rvec> invMassPerDim,
                                                const t_state&            state,
                                                const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                                const gmx_ekindata_t& ekind)
{
    GMX_ASSERT(inputRecord.eI == IntegrationAlgorithm::MD || inputRecord.eI == IntegrationAlgorithm::SD1,
               "Only leap-frog is supported here");

    // Cast to real for faster code, no loss in precision
    const real dt = inputRecord.delta_t;

    const int nth = gmx_omp_nthreads_get(ModuleMultiThread::Update);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            int start_th, end_th;
            getThreadAtomRange(nth, th, homenr, &start_th, &end_th);

            const rvec* x_rvec  = state.x.rvec_array();
            rvec*       xp_rvec = xp_.rvec_array();
            const rvec* v_rvec  = state.v.rvec_array();
            const rvec* f_rvec  = as_rvec_array(f.unpaddedConstArrayRef().data());

            doUpdateMDDoNotUpdateVelocities(
                    start_th, end_th, dt, x_rvec, xp_rvec, v_rvec, f_rvec, havePartiallyFrozenAtoms, invmass, invMassPerDim, ekind);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}
