/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "update.h"

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <type_traits>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/arrayrefwithpadding.h"
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
#include "gromacs/random/seed.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/template_mp.h"

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

/*! \brief Sets the NEMD acceleration type */
enum class AccelerationType
{
    None,
    Group,
    Cosine,
    BoxDeformation,
    Count
};

//! pImpled implementation for Update
class Update::Impl
{
public:
    //! Constructor
    Impl(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind, BoxDeformation* boxDeformation);
    //! Destructor
    ~Impl() = default;

    void update_coords(const t_inputrec&                                inputRecord,
                       int64_t                                          step,
                       int                                              homenr,
                       bool                                             havePartiallyFrozenAtoms,
                       gmx::ArrayRef<const ParticleType>                ptype,
                       gmx::ArrayRef<const real>                        invMass,
                       gmx::ArrayRef<const gmx::RVec>                   invMassPerDim,
                       t_state*                                         state,
                       const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                       t_fcdata*                                        fcdata,
                       const gmx_ekindata_t*                            ekind,
                       const Matrix3x3&                                 parrinelloRahmanM,
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

    void update_for_constraint_virial(const t_inputrec&              inputRecord,
                                      int                            homenr,
                                      bool                           havePartiallyFrozenAtoms,
                                      gmx::ArrayRef<const real>      invmass,
                                      gmx::ArrayRef<const gmx::RVec> invMassPerDim,
                                      const t_state&                 state,
                                      const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                      const gmx_ekindata_t&                            ekind);

    void update_temperature_constants(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind);

    const std::vector<bool>& getAndersenRandomizeGroup() const { return sd_.randomize_group; }

    const std::vector<real>& getBoltzmanFactor() const { return sd_.boltzfac; }

    PaddedVector<RVec>* xp() { return &xp_; }

    BoxDeformation* deform() const { return deform_; }

    //! Group index for freezing
    gmx::ArrayRef<const unsigned short> cFREEZE_;
    //! Group index for temperature coupling
    gmx::ArrayRef<const unsigned short> cTC_;
    //! Group index for acceleration groups
    gmx::ArrayRef<const unsigned short> cAcceleration_;

private:
    //! Type of acceleration used in the simulation
    AccelerationType accelerationType_;
    //! stochastic dynamics struct
    gmx_stochd_t sd_;
    //! xprime for constraint algorithms
    PaddedVector<RVec> xp_;
    //! Box deformation handler (or nullptr if inactive).
    BoxDeformation* deform_ = nullptr;
};

Update::Update(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind, BoxDeformation* boxDeformation) :
    impl_(new Impl(inputRecord, ekind, boxDeformation)){};

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
                           gmx::ArrayRef<const gmx::RVec>    invMassPerDim,
                           t_state*                          state,
                           const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                           t_fcdata*                                        fcdata,
                           const gmx_ekindata_t*                            ekind,
                           const Matrix3x3&                                 parrinelloRahmanM,
                           int                                              updatePart,
                           const t_commrec*                                 cr,
                           const bool                                       haveConstraints)
{
    impl_->update_coords(inputRecord,
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
                         parrinelloRahmanM,
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
    impl_->finish_update(
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
    impl_->update_sd_second_half(
            inputRecord, step, dvdlambda, homenr, ptype, invMass, state, cr, nrnb, wcycle, constr, do_log, do_ene);
}

void Update::update_for_constraint_virial(const t_inputrec&              inputRecord,
                                          const int                      homenr,
                                          const bool                     havePartiallyFrozenAtoms,
                                          gmx::ArrayRef<const real>      invmass,
                                          gmx::ArrayRef<const gmx::RVec> invMassPerDim,
                                          const t_state&                 state,
                                          const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                          const gmx_ekindata_t&                            ekind)
{
    impl_->update_for_constraint_virial(
            inputRecord, homenr, havePartiallyFrozenAtoms, invmass, invMassPerDim, state, f, ekind);
}

void Update::update_temperature_constants(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind)
{
    impl_->update_temperature_constants(inputRecord, ekind);
}

/*! \brief Sets whether we store the updated velocities */
enum class StoreUpdatedVelocities
{
    Yes, //!< Store the updated velocities
    No,  //!< Do not store the updated velocities
    Count
};

/*! \brief Sets the number of different temperature coupling values */
enum class NumTempScaleValues
{
    None,     //!< No temperature scaling
    Single,   //!< Single T-scaling value (either one group or all values =1)
    Multiple, //!< Multiple T-scaling values, need to use T-group indices
    Count
};

//! Describes the properties of the Parrinello-Rahman pressure scaling matrix
enum class ParrinelloRahmanVelocityScaling
{
    No,          //!< Do not apply velocity scaling (not a PR-coupling run or step)
    Diagonal,    //!< Apply velocity scaling using a diagonal matrix
    Anisotropic, //!< Apply velocity scaling using a matrix with off-diagonal elements
    Count
};

/*! \brief Integrate using leap-frog with T-scaling and optionally diagonal Parrinello-Rahman p-coupling
 *
 * \tparam       storeUpdatedVelocities Tells whether we should store the updated velocities
 * \tparam       numTempScaleValues     The number of different T-couple values
 * \tparam       applyPRVScaling        Apply Parrinello-Rahman velocity scaling
 * \tparam       parrinelloRahmanVelocityScaling  The properties of the Parrinello-Rahman velocity scaling matrix
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
template<StoreUpdatedVelocities storeUpdatedVelocities, NumTempScaleValues numTempScaleValues, ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling, typename VelocityType>
static std::enable_if_t<std::is_same<VelocityType, rvec>::value || std::is_same<VelocityType, const rvec>::value, void>
updateMDLeapfrogSimple(int                                 start,
                       int                                 nrend,
                       real                                dt,
                       real                                dtPressureCouple,
                       gmx::ArrayRef<const gmx::RVec>      invMassPerDim,
                       gmx::ArrayRef<const t_grp_tcstat>   tcstat,
                       gmx::ArrayRef<const unsigned short> cTC,
                       const gmx::RVec                     pRVScaleMatrixDiagonal,
                       const rvec* gmx_restrict            x,
                       rvec* gmx_restrict                  xprime,
                       VelocityType* gmx_restrict          v,
                       const rvec* gmx_restrict            f)
{
    real lambdaGroup;

    if (numTempScaleValues == NumTempScaleValues::None)
    {
        lambdaGroup = 1.0_real;
    }
    else if (numTempScaleValues == NumTempScaleValues::Single)
    {
        lambdaGroup = tcstat[0].lambda;
    }

    for (int a = start; a < nrend; a++)
    {
        if (numTempScaleValues == NumTempScaleValues::Multiple)
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

            if (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
            {
                vNew -= dtPressureCouple * pRVScaleMatrixDiagonal[d] * v[a][d];
            }
            if constexpr (storeUpdatedVelocities == StoreUpdatedVelocities::Yes)
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
using UpdateSimdReal   = SimdReal;
using UpdateSimdTraits = gmx::internal::SimdTraits<UpdateSimdReal>;
#else
#    define GMX_HAVE_SIMD_UPDATE 0
using UpdateSimdReal = real;
struct UpdateSimdTraits
{
    using type                 = UpdateSimdReal;
    static constexpr int width = 1;
};
#endif

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
static inline void simdLoadRvecs(const rvec* r, int index, UpdateSimdReal* r0, UpdateSimdReal* r1, UpdateSimdReal* r2)
{
    const real* realPtr = r[index];

    *r0 = load<UpdateSimdReal>(realPtr + 0 * UpdateSimdTraits::width);
    *r1 = load<UpdateSimdReal>(realPtr + 1 * UpdateSimdTraits::width);
    *r2 = load<UpdateSimdReal>(realPtr + 2 * UpdateSimdTraits::width);
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
static inline void simdStoreRvecs(rvec* r, int index, UpdateSimdReal r0, UpdateSimdReal r1, UpdateSimdReal r2)
{
    real* realPtr = r[index];

    store(realPtr + 0 * UpdateSimdTraits::width, r0);
    store(realPtr + 1 * UpdateSimdTraits::width, r1);
    store(realPtr + 2 * UpdateSimdTraits::width, r2);
}

/*! \brief Integrate using leap-frog with single group T-scaling and SIMD
 *
 * \tparam       storeUpdatedVelocities Tells whether we should store the updated velocities
 * \tparam       numTempScaleValues     The number of different T-couple values
 * \tparam       VelocityType           Either rvec or const rvec according to whether we store velocities
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
template<StoreUpdatedVelocities storeUpdatedVelocities, NumTempScaleValues numTempScaleValues, typename VelocityType>
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
    UpdateSimdReal timestep(dt);
    UpdateSimdReal lambdaSystem(tcstat[0].lambda);

    for (int a = start; a < nrend; a += UpdateSimdTraits::width)
    {
        UpdateSimdReal invMass0, invMass1, invMass2;
        expandScalarsToTriplets(load<UpdateSimdReal>(invMass.data() + a), &invMass0, &invMass1, &invMass2);

        UpdateSimdReal v0, v1, v2;
        UpdateSimdReal f0, f1, f2;
        simdLoadRvecs(v, a, &v0, &v1, &v2);
        simdLoadRvecs(f, a, &f0, &f1, &f2);

        if constexpr (numTempScaleValues == NumTempScaleValues::None)
        {
            v0 = gmx::fma(f0 * invMass0, timestep, v0);
            v1 = gmx::fma(f1 * invMass1, timestep, v1);
            v2 = gmx::fma(f2 * invMass2, timestep, v2);
        }
        else
        {
            static_assert(
                    numTempScaleValues == NumTempScaleValues::Single,
                    "Invalid multiple temperature-scaling values (ie. groups) for SIMD leapfrog");
            v0 = gmx::fma(f0 * invMass0, timestep, lambdaSystem * v0);
            v1 = gmx::fma(f1 * invMass1, timestep, lambdaSystem * v1);
            v2 = gmx::fma(f2 * invMass2, timestep, lambdaSystem * v2);
        }
        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        if constexpr (storeUpdatedVelocities == StoreUpdatedVelocities::Yes)
        {
            simdStoreRvecs(v, a, v0, v1, v2);
        }

        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        UpdateSimdReal x0, x1, x2;
        simdLoadRvecs(x, a, &x0, &x1, &x2);

        UpdateSimdReal xprime0 = gmx::fma(v0, timestep, x0);
        UpdateSimdReal xprime1 = gmx::fma(v1, timestep, x1);
        UpdateSimdReal xprime2 = gmx::fma(v2, timestep, x2);

        simdStoreRvecs(xprime, a, xprime0, xprime1, xprime2);
    }
}

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
 * \param[in]     boxDeformation    Velocity of box components in units of nm/ps
 * \param[in]     invMassPerDim     Inverse mass per dimension
 * \param[in]     ekind             Kinetic energy data.
 * \param[in]     box               The box dimensions.
 * \param[in]     x                 Input coordinates.
 * \param[out]    xprime            Updated coordinates.
 * \param[inout]  v                 Velocities.
 * \param[in]     f                 Forces.
 * \param[in]     nh_vxi            Nose-Hoover velocity scaling factors.
 * \param[in]     nsttcouple        Frequency of the temperature coupling steps.
 * \param[in]     parrinelloRahmanM Parrinello-Rahman scaling matrix.
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
                                    const matrix                        boxDeformation,
                                    gmx::ArrayRef<const gmx::RVec>      invMassPerDim,
                                    const gmx_ekindata_t*               ekind,
                                    const matrix                        box,
                                    const rvec* gmx_restrict            x,
                                    rvec* gmx_restrict                  xprime,
                                    rvec* gmx_restrict                  v,
                                    const rvec* gmx_restrict            f,
                                    const double* gmx_restrict          nh_vxi,
                                    const int                           nsttcouple,
                                    const Matrix3x3&                    parrinelloRahmanM)
{
    /* This is a version of the leap-frog integrator that supports
     * all combinations of T-coupling, P-coupling and NEMD.
     * Nose-Hoover uses the reversible leap-frog integrator from
     * Holian et al. Phys Rev E 52(3) : 2338, 1995
     */

    gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind->tcstat;

    // Matrix that multiplied with the velocity gives the flow
    matrix deformFlowMatrix;
    // The velocity of the system COM after subtracting the flow profile
    rvec systemVelocity;
    if (accelerationType == AccelerationType::BoxDeformation)
    {
        setBoxDeformationFlowMatrix(boxDeformation, box, deformFlowMatrix);

        const SystemMomentum& systemMomentum = ekind->systemMomenta->momentumHalfStep;
        for (int d = 0; d < DIM; d++)
        {
            systemVelocity[d] = systemMomentum.momentum[d] / systemMomentum.mass;
        }
    }

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
        rvec vFlow;
        real cosineZ, vCosine;
        copy_rvec(v[n], vRel);
        switch (accelerationType)
        {
            case AccelerationType::None: break;
            case AccelerationType::Group:
                if (!cAcceleration.empty())
                {
                    ga = cAcceleration[n];
                }
                /* With constant acceleration we do scale the velocity of the accelerated groups */
                break;
            case AccelerationType::Cosine:
                cosineZ = std::cos(x[n][ZZ] * omega_Z);
                vCosine = cosineZ * ekind->cosacc.vcos;
                /* Avoid scaling the cosine profile velocity */
                vRel[XX] -= vCosine;
                break;
            case AccelerationType::BoxDeformation:
                for (int d = 0; d < DIM; d++)
                {
                    vFlow[d] = iprod(x[n], deformFlowMatrix[d]) - systemVelocity[d];
                    vRel[d] -= vFlow[d];
                }
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

        const RVec parrinelloRahmanScaledVelocity =
                dtPressureCouple * multiplyVectorByMatrix(parrinelloRahmanM, vRel);
        for (int d = 0; d < DIM; d++)
        {
            real vNew = (lg * vRel[d]
                         + (f[n][d] * invMassPerDim[n][d] * dt - factorNH * vRel[d]
                            - parrinelloRahmanScaledVelocity[d]))
                        / (1 + factorNH);
            switch (accelerationType)
            {
                case AccelerationType::None: break;
                case AccelerationType::Group:
                    /* Apply the constant acceleration */
                    vNew += acceleration[ga][d] * dt;
                    break;
                case AccelerationType::Cosine:
                    if (d == XX)
                    {
                        /* Add back the mean velocity and apply acceleration */
                        vNew += vCosine + cosineZ * ekind->cosacc.cos_accel * dt;
                    }
                    break;
                case AccelerationType::BoxDeformation: vNew += vFlow[d]; break;
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
                         const AccelerationType               simulationAccelerationType,
                         gmx::ArrayRef<const unsigned short>  cAcceleration,
                         const rvec*                          acceleration,
                         const matrix                         boxDeformation,
                         gmx::ArrayRef<const real> gmx_unused invmass,
                         gmx::ArrayRef<const gmx::RVec>       invMassPerDim,
                         const gmx_ekindata_t*                ekind,
                         const matrix                         box,
                         const double* gmx_restrict           nh_vxi,
                         const Matrix3x3&                     parrinelloRahmanM,
                         bool gmx_unused                      havePartiallyFrozenAtoms)
{
    GMX_ASSERT(nrend == start || xprime != x,
               "For SIMD optimization certain compilers need to have xprime != x");

    /* Note: Berendsen pressure scaling is handled after do_update_md() */
    const bool doTempCouple =
            (etc != TemperatureCoupling::No && do_per_step(step + nsttcouple - 1, nsttcouple));
    const bool doNoseHoover = (etc == TemperatureCoupling::NoseHoover && doTempCouple);
    const bool doParrinelloRahmanThisStep = (epc == PressureCoupling::ParrinelloRahman
                                             && do_per_step(step + nstpcouple - 1, nstpcouple));
    ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling =
            (doParrinelloRahmanThisStep ? ((parrinelloRahmanM(YY, XX) != 0 || parrinelloRahmanM(ZZ, XX) != 0
                                            || parrinelloRahmanM(ZZ, YY) != 0)
                                                   ? ParrinelloRahmanVelocityScaling::Anisotropic
                                                   : ParrinelloRahmanVelocityScaling::Diagonal)
                                        : ParrinelloRahmanVelocityScaling::No);

    real dtPressureCouple =
            ((parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No) ? nstpcouple * dt : 0);

    /* The forces for constant and cosine acceleration are applied in update, so we need to use
     * the general update function as that is the only one that implements acceleration.
     * With box deformation, the temperature coupling should not scale the flow profile,
     * so we need to subtract that which is only done in the general update function.
     * At non Tcouple steps we should do a normal update to avoid accessing invalid
     * flow data.
     */
    const AccelerationType stepAccelerationType =
            (simulationAccelerationType == AccelerationType::BoxDeformation && !doTempCouple
                     ? AccelerationType::None
                     : simulationAccelerationType);

    if (doNoseHoover || (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Anisotropic)
        || stepAccelerationType != AccelerationType::None)
    {
        // If there's no Parrinello-Rahman scaling this step, we need to pass a zero matrix instead
        Matrix3x3        zero = { { 0._real } };
        const Matrix3x3& parrinelloRahmanMToUseThisStep =
                parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No ? parrinelloRahmanM
                                                                                       : zero;

        dispatchTemplatedFunction(
                [=](auto stepAccelerationType) {
                    return updateMDLeapfrogGeneral<stepAccelerationType>(start,
                                                                         nrend,
                                                                         doNoseHoover,
                                                                         dt,
                                                                         dtPressureCouple,
                                                                         cTC,
                                                                         cAcceleration,
                                                                         acceleration,
                                                                         boxDeformation,
                                                                         invMassPerDim,
                                                                         ekind,
                                                                         box,
                                                                         x,
                                                                         xprime,
                                                                         v,
                                                                         f,
                                                                         nh_vxi,
                                                                         nsttcouple,
                                                                         parrinelloRahmanMToUseThisStep);
                },
                stepAccelerationType);
    }
    else
    {
        /* We determine if we have a single T-coupling lambda value for all
         * atoms. That allows for better SIMD acceleration in the template.
         * If we do not do temperature coupling (in the run or this step),
         * all scaling values are 1, so we effectively have a single value.
         */
        const int          ntcg = ekind->numTemperatureCouplingGroups();
        NumTempScaleValues numTempScaleValues =
                (!doTempCouple || ntcg == 0)
                        ? NumTempScaleValues::None
                        : (ntcg == 1 ? NumTempScaleValues::Single : NumTempScaleValues::Multiple);

        /* Extract some pointers needed by all cases */
        gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind->tcstat;

        if (parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No)
        {
            GMX_ASSERT((parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal),
                       "updateMDLeapfrogSimple only supports diagonal Parrinello-Rahman scaling "
                       "matrices when scaling is active");
            RVec diagM = {
                parrinelloRahmanM(XX, XX),
                parrinelloRahmanM(YY, YY),
                parrinelloRahmanM(ZZ, ZZ),
            };

            dispatchTemplatedFunction(
                    [=](auto numTempScaleValues, auto parrinelloRahmanVelocityScaling) {
                        return updateMDLeapfrogSimple<StoreUpdatedVelocities::Yes, numTempScaleValues, parrinelloRahmanVelocityScaling>(
                                start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, diagM, x, xprime, v, f);
                    },
                    numTempScaleValues,
                    parrinelloRahmanVelocityScaling);
        }
        else
        {
            if (GMX_HAVE_SIMD_UPDATE && numTempScaleValues != NumTempScaleValues::Multiple
                && !havePartiallyFrozenAtoms)
            {
                // Because no atoms are partially frozen, we can use
                // invmass instead of invMassPerDim.
                if (numTempScaleValues == NumTempScaleValues::Single)
                {
                    updateMDLeapfrogSimpleSimd<StoreUpdatedVelocities::Yes, NumTempScaleValues::Single>(
                            start, nrend, dt, invmass, tcstat, x, xprime, v, f);
                }
                else
                {
                    GMX_ASSERT(numTempScaleValues == NumTempScaleValues::None,
                               "Must not have temperature scaling values here");
                    updateMDLeapfrogSimpleSimd<StoreUpdatedVelocities::Yes, NumTempScaleValues::None>(
                            start, nrend, dt, invmass, tcstat, x, xprime, v, f);
                }
            }
            else
            {
                dispatchTemplatedFunction(
                        [=](auto numTempScaleValues) {
                            /* Note that modern compilers are pretty good at vectorizing
                             * updateMDLeapfrogSimple(). But the SIMD version will still
                             * be faster because invMass lowers the cache pressure
                             * compared to invMassPerDim.
                             */
                            {
                                updateMDLeapfrogSimple<StoreUpdatedVelocities::Yes, numTempScaleValues, ParrinelloRahmanVelocityScaling::No>(
                                        start, nrend, dt, dtPressureCouple, invMassPerDim, tcstat, cTC, {}, x, xprime, v, f);
                            }
                        },
                        numTempScaleValues);
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
                                            gmx::ArrayRef<const gmx::RVec>       invMassPerDim,
                                            const gmx_ekindata_t&                ekind)
{
    GMX_ASSERT(nrend == start || xprime != x,
               "For SIMD optimization certain compilers need to have xprime != x");

    gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind.tcstat;

    /* Check if we can use invmass instead of invMassPerDim */
    if (GMX_HAVE_SIMD_UPDATE && !havePartiallyFrozenAtoms)
    {
        updateMDLeapfrogSimpleSimd<StoreUpdatedVelocities::No, NumTempScaleValues::Single>(
                start, nrend, dt, invmass, tcstat, x, xprime, v, f);
    }
    else
    {
        updateMDLeapfrogSimple<StoreUpdatedVelocities::No, NumTempScaleValues::Single, ParrinelloRahmanVelocityScaling::No>(
                start, nrend, dt, dt, invMassPerDim, tcstat, gmx::ArrayRef<const unsigned short>(), {}, x, xprime, v, f);
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

void Update::Impl::update_temperature_constants(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind)
{
    const int ntcg = ekind.numTemperatureCouplingGroups();

    if (inputRecord.eI == IntegrationAlgorithm::BD)
    {
        if (inputRecord.bd_fric != 0)
        {
            for (int gt = 0; gt < ntcg; gt++)
            {
                sd_.bd_rf[gt] = std::sqrt(2.0 * gmx::c_boltz * ekind.currentReferenceTemperature(gt)
                                          / (inputRecord.bd_fric * inputRecord.delta_t));
            }
        }
        else
        {
            for (int gt = 0; gt < ntcg; gt++)
            {
                sd_.bd_rf[gt] = std::sqrt(2.0 * gmx::c_boltz * ekind.currentReferenceTemperature(gt));
            }
        }
    }
    if (inputRecord.eI == IntegrationAlgorithm::SD1)
    {
        for (int gt = 0; gt < ntcg; gt++)
        {
            real kT = gmx::c_boltz * ekind.currentReferenceTemperature(gt);
            /* The mass is accounted for later, since this differs per atom */
            sd_.sdsig[gt].V = std::sqrt(kT * (1 - sd_.sdc[gt].em * sd_.sdc[gt].em));
        }
    }
}

Update::Impl::Impl(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind, BoxDeformation* boxDeformation) :
    accelerationType_(inputRecord.useConstantAcceleration
                              ? AccelerationType::Group
                              : ((inputRecord.cos_accel != 0) ? AccelerationType::Cosine
                                                              : (ir_haveBoxDeformation(inputRecord)
                                                                         ? AccelerationType::BoxDeformation
                                                                         : AccelerationType::None))),
    sd_(inputRecord),
    deform_(boxDeformation)
{
    update_temperature_constants(inputRecord, ekind);
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
    Combined,
    Count
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
                              const int*                          gatindex,
                              real                                dtPressureCouple,
                              const Matrix3x3&                    parrinelloRahmanM)
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

        RVec parrinelloRahmanScaledVelocity;
        if (updateType != SDUpdate::FrictionAndNoiseOnly)
        {
            parrinelloRahmanScaledVelocity =
                    dtPressureCouple * multiplyVectorByMatrix(parrinelloRahmanM, v[n]);
        }
        for (int d = 0; d < DIM; d++)
        {
            if ((ptype[n] != ParticleType::Shell) && !nFreeze[freezeGroup][d])
            {
                if (updateType == SDUpdate::ForcesOnly)
                {
                    real vn = v[n][d] + (inverseMass * f[n][d] + acceleration[accelerationGroup][d]) * dt
                              - parrinelloRahmanScaledVelocity[d];
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
                    real vn = v[n][d] + (inverseMass * f[n][d] + acceleration[accelerationGroup][d]) * dt
                              - parrinelloRahmanScaledVelocity[d];
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
                         bool                                haveConstraints,
                         const PressureCoupling              pressureCoupling,
                         const int                           nstpcouple,
                         const Matrix3x3&                    parrinelloRahmanM)
{
    const bool doParrinelloRahmanThisStep = (pressureCoupling == PressureCoupling::ParrinelloRahman
                                             && do_per_step(step + nstpcouple - 1, nstpcouple));
    ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling =
            (doParrinelloRahmanThisStep ? ParrinelloRahmanVelocityScaling::Anisotropic
                                        : ParrinelloRahmanVelocityScaling::No);
    // If there's no Parrinello-Rahman scaling this step, we need to pass a zero matrix instead
    Matrix3x3        zero = { { 0._real } };
    const Matrix3x3& parrinelloRahmanMToUseThisStep =
            parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No ? parrinelloRahmanM
                                                                                   : zero;
    const real dtPressureCouple =
            ((parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No) ? nstpcouple * dt : 0);

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
                                                nullptr,
                                                dtPressureCouple,
                                                parrinelloRahmanMToUseThisStep);
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
                (cr != nullptr && haveDDAtomOrdering(*cr)) ? cr->dd->globalAtomIndices.data() : nullptr,
                dtPressureCouple,
                parrinelloRahmanMToUseThisStep);
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

void update_ekinstate(ekinstate_t* ekinstate, const gmx_ekindata_t* ekind, const bool sumEkin, const t_commrec* cr)
{
    /* Note that it might seem like we are storing the current kinetic energy at time t+dt/2 here
     * as we are operating on ekinh and dekindl. But this function is called before the kinetic
     * energy is computed. Thus we are actually operating on the kinetic energy at time t-dt/2.
     * The kinetic energy computation called right after this will copy the values stored here
     * to ekind->ekinh_old and ekind->dekindl_old.
     */

    const bool reduceEkin = (sumEkin && havePPDomainDecomposition(cr));

    if (reduceEkin)
    {
        /* The kinetic energy terms from t-dt/2 have not been reduced yet.
         * We need to checkpoint reduced values. We reduce here using double
         * precision so we get binary identical reduced results compared with
         * the reduction in compte_globals() which also uses double precision.
         */
        const int           ntcg = ekind->numTemperatureCouplingGroups();
        std::vector<double> buffer(ntcg * 2 * DIM * DIM + 1);
        int                 bufIndex = 0;
        for (int g = 0; g < ntcg; g++)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    buffer[bufIndex++] = ekind->tcstat[g].ekinh[i][j];
                }
            }
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    buffer[bufIndex++] = ekind->tcstat[g].ekinf[i][j];
                }
            }
        }
        buffer[bufIndex++] = ekind->dekindl;

        gmx_sumd(bufIndex, buffer.data(), cr);

        // Extract to ekinstate on the main rank only
        if (MAIN(cr))
        {
            bufIndex = 0;
            for (int g = 0; g < ekinstate->ekin_n; g++)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        ekinstate->ekinh[g][i][j] = buffer[bufIndex++];
                    }
                }
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        ekinstate->ekinf[g][i][j] = buffer[bufIndex++];
                    }
                }
            }
            ekinstate->dekindl = buffer[bufIndex++];
        }
    }

    if (MAIN(cr))
    {
        if (!reduceEkin)
        {
            for (int g = 0; g < ekinstate->ekin_n; g++)
            {
                copy_mat(ekind->tcstat[g].ekinh, ekinstate->ekinh[g]);
                copy_mat(ekind->tcstat[g].ekinf, ekinstate->ekinf[g]);
            }
            ekinstate->dekindl = ekind->dekindl;
        }

        /* These terms are likely not part of the state at all and can be removed
         * as they are (re)computed when restarting from a checkpoint.
         */
        for (int g = 0; g < ekinstate->ekin_n; g++)
        {
            ekinstate->ekinscalef_nhc[g] = ekind->tcstat[g].ekinscalef_nhc;
            ekinstate->ekinscaleh_nhc[g] = ekind->tcstat[g].ekinscaleh_nhc;
            ekinstate->vscale_nhc[g]     = ekind->tcstat[g].vscale_nhc;
        }
        ekinstate->mvcos = ekind->cosacc.mvcos;
    }
}

void restore_ekinstate_from_state(const t_commrec* cr, gmx_ekindata_t* ekind, const ekinstate_t* ekinstate)
{
    int i, n;

    if (MAIN(cr))
    {
        for (i = 0; i < ekinstate->ekin_n; i++)
        {
            copy_mat(ekinstate->ekinh[i], ekind->tcstat[i].ekinh);
            copy_mat(ekinstate->ekinf[i], ekind->tcstat[i].ekinf);
            ekind->tcstat[i].ekinscalef_nhc = ekinstate->ekinscalef_nhc[i];
            ekind->tcstat[i].ekinscaleh_nhc = ekinstate->ekinscaleh_nhc[i];
            ekind->tcstat[i].vscale_nhc     = ekinstate->vscale_nhc[i];
        }

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

            gmx_bcast(sizeof(ekind->tcstat[i].ekinscalef_nhc),
                      &(ekind->tcstat[i].ekinscalef_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].ekinscaleh_nhc),
                      &(ekind->tcstat[i].ekinscaleh_nhc),
                      cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(ekind->tcstat[i].vscale_nhc), &(ekind->tcstat[i].vscale_nhc), cr->mpi_comm_mygroup);
        }

        gmx_bcast(sizeof(ekind->dekindl), &ekind->dekindl, cr->mpi_comm_mygroup);
        gmx_bcast(sizeof(ekind->cosacc.mvcos), &ekind->cosacc.mvcos, cr->mpi_comm_mygroup);
    }
}

void getThreadAtomRange(int numThreads, int threadIndex, int numAtoms, int* startAtom, int* endAtom)
{
    constexpr int blockSize = UpdateSimdTraits::width;
    const int     numBlocks = divideRoundUp(numAtoms, blockSize);

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

        Matrix3x3 parrinelloRahmanM{ { 0._real } };
        real      dtPressureCouple = 0;

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
                        haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr,
                        dtPressureCouple,
                        parrinelloRahmanM);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        inc_nrnb(nrnb, eNR_UPDATE, homenr);
        wallcycle_stop(wcycle, WallCycleCounter::Update);

        /* Constrain the coordinates upd->xp for half a time step */
        bool computeVirial = false;
        constr->apply(do_log || do_ene,
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
                                 gmx::ArrayRef<const gmx::RVec>    invMassPerDim,
                                 t_state*                          state,
                                 const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                 t_fcdata*                                        fcdata,
                                 const gmx_ekindata_t*                            ekind,
                                 const Matrix3x3&                                 parrinelloRahmanM,
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
    if (state->hasEntry(StateEntry::DisreRm3Tav))
    {
        update_disres_history(*fcdata->disres, &state->hist);
    }
    if (state->hasEntry(StateEntry::OrireDtav))
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
                                 inputRecord.pressureCouplingOptions.epc,
                                 inputRecord.nsttcouple,
                                 inputRecord.pressureCouplingOptions.nstpcouple,
                                 cTC_,
                                 accelerationType_,
                                 cAcceleration_,
                                 inputRecord.opts.acceleration,
                                 inputRecord.deform,
                                 invMass,
                                 invMassPerDim,
                                 ekind,
                                 state->box,
                                 state->nosehoover_vxi.data(),
                                 parrinelloRahmanM,
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
                                 haveConstraints,
                                 inputRecord.pressureCouplingOptions.epc,
                                 inputRecord.pressureCouplingOptions.nstpcouple,
                                 parrinelloRahmanM);
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
                    gmx_bool bExtended =
                            (inputRecord.etc == TemperatureCoupling::NoseHoover
                             || inputRecord.pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
                             || inputRecord.pressureCouplingOptions.epc == PressureCoupling::Mttk);

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
                                                gmx::ArrayRef<const gmx::RVec> invMassPerDim,
                                                const t_state&                 state,
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
