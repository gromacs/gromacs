/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

/*! \libinternal
 * \defgroup module_awh Accelerated weight histogram (AWH) method
 * \ingroup group_mdrun
 * \brief
 * Implements the "accelerated weight histogram" sampling method.
 *
 * This class provides the interface between the AWH module and
 * other modules using it. Currently AWH can only act on COM pull
 * reaction coordinates, but this can easily be extended to other
 * types of reaction coordinates.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 */

/*! \libinternal \file
 *
 * \brief
 * Declares the Awh class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_awh
 */

#ifndef GMX_AWH_H
#define GMX_AWH_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_multisim_t;
struct gmx_wallcycle;
struct pull_t;
class t_state;
struct t_commrec;
struct t_inputrec;
struct t_mdatoms;

namespace gmx
{

struct AwhHistory;
struct AwhParams;
class Bias;
class BiasCoupledToSystem;

/*! \libinternal \brief The accelerated weight histogram method (AWH).
 *
 * AWH calculates the free energy along an order parameter of the system.
 * Free energy barriers are overcome by adaptively tuning a bias along the parameter
 * such that the biased distribution along the parameter converges toward a chosen target distribution.
 * The fundamental equation governing the tuning is: log(target) = bias - free energy, where
 * the bias and free energy are initially unknown. Typically the target distribution is simply
 * chosen uniform, such that the bias completely flattens the free energy landscape.
 *
 * This class implements AWH for the case when the order parameter corresponds to a reaction coordinate,
 * here referred to as coordinate for short, i.e. a function of the system configuration.
 * The bias is coupled to the system by a bias potential: either in the form of an harmonic ("umbrella") potential
 * Monte-Carlo (MC) "jumping" around the current coordinate value, or as a smooth average of the umbrellas.
 *
 * The AWH module is organizes as follows:
 * The Awh class is the interface between the outside and inside of the module.
 * The Awh class contains one or more BiasCoupledToSystem objects.
 * The BiasCoupledToSystem class takes care of the reaction coordinate input
 * and force output for the single Bias object it containts.
 * The Bias class is a container and wrapper for a object BiasState + helpers.
 * All computation takes place in the BiasState object and its sub-classes.
 * The Bias class also contains a BiasWriter object that takes care of i/o.
 * Additionally, there are, currently quite messy, files that care of
 * parameter reading and checkpointing. These should be rewritten when we
 * have proper, general C++ modules that can take of these tasks.
 *
 * The basic use of Awh in mdrun consists of 2 method calls:
 * Call the constructor Awh() after the pull module has been initialized.
 * Call applyBiasForcesAndUpdateBias() at every MD step after the pull
 * potential calculation function has been called.
 *
 * In grompp the pull potential provider should be registerd using
 * registerAwhWithPull() so grompp can check for unregistered potentials.
 *
 * The main tasks of AWH are:
 * - calculate and set the bias force given the current coordinate value.
 * - after accumulating a number of coordinate samples, update the free energy estimate and the bias.
 *
 * AWH currently relies on the pull code for the first task. Pull provides AWH with updated coordinate values
 * and distributes the bias force that AWH calculates to the atoms making up the coordinate. This
 * also means that there are some order dependencies where pull functions need to be called before AWH
 * functions (see below).
 *
 * The implementation is quite general. There can be multiple independent AWH biases coupled to the system
 * simultaneously. This makes sense if the system is made up of several fairly independent parts,
 * like monomers in a protein. Each bias acts on exactly one, possibly multidimensional, coordinate.
 * Each coordinate dimension maps to exactly one pull coordinate. Thus, an n-dimensional
 * biased coordinate is defined by a set of n pull coordinates. Periodicity is taken care of for coordinate
 * dimensions that require it (ddihedral angles). For increased parallelism, there is the option of
 * having multiple communicating simulations sharing all samples. All simulations would then share a single
 * bias and free energy estimate. Alternatively, one may partition the sampling domain into smaller
 * subdomains with some overlap and have multiple independent simulations sample each subdomain.
 *
 */
class Awh
{
    public:
        /*! \brief Construct an AWH at the start of a simulation.
         *
         * AWH will here also register itself with the pull struct as the
         * potential provider for the pull coordinates given as AWH coordinates
         * in the user input. This allows AWH to later apply the bias force to
         * these coordinate in \ref Awh::applyBiasForcesAndUpdateBias.
         *
         * \param[in,out] fplog           General output file, normally md.log, can be nullptr.
         * \param[in]     ir              General input parameters.
         * \param[in]     cr              Struct for communication, can be nullptr.
         * \param[in]     awhParams       AWH input parameters.
         * \param[in,out] pull_work       Pull struct which AWH will register the bias into, has to be initialized.
         */
        Awh(FILE              *fplog,
            const t_inputrec  &ir,
            const t_commrec   *cr,
            const AwhParams   &awhParams,
            pull_t            *pull_work);

        /*! \brief Destructor. */
        ~Awh();

        /*! \brief Peform an AWH update, to be called every MD step.
         *
         * An update has two tasks: apply the bias force and improve
         * the bias and the free energy estimate that AWH keeps internally.
         *
         * For the first task, AWH retrieves the pull coordinate values from the pull struct.
         * With these, the bias potential and forces are calculated.
         * The bias force together with the atom forces and virial
         * are passed on to pull which applies the bias force to the atoms.
         * This is done at every step.
         *
         * Secondly, coordinate values are regularly sampled and kept by AWH.
         * Convergence of the bias and free energy estimate is achieved by
         * updating the AWH bias state after a certain number of samples has been collected.
         *
         * \note Requires that pull_potential from pull.h has been called first
         * since AWH needs the current coordinate values.
         *
         * \param[in,out] pull_work    Pull working struct.
         * \param[in]   mdatoms        Atom properties.
         * \param[in] ePBC             Type of periodic boundary conditions.
         * \param[in] box              Box vectors.
         * \param[in,out] force        Forces.
         * \param[in,out] virial       The virial
         * \param[in] ms               Struct for multi-simulation communication.
         * \param[in]     t            Time.
         * \param[in]     step         Time step.
         * \param[in,out] wallcycle    Wallcycle counter.
         * \param[in,out] fplog        General output file, normally md.log.
         * \returns the potential energy for the bias.
         */
        real applyBiasForcesAndUpdateBias(struct pull_t          *pull_work,
                                          int                     ePBC,
                                          const t_mdatoms        &mdatoms,
                                          const matrix            box,
                                          rvec                   *force,
                                          tensor                  virial,
                                          const gmx_multisim_t   *ms,
                                          double                  t,
                                          gmx_int64_t             step,
                                          gmx_wallcycle          *wallcycle,
                                          FILE                   *fplog);

        /*! \brief
         * Update the AWH history in preparation for writing to checkpoint file.
         *
         * \param[in,out] awhHistory  AWH history to set.
         */
        void updateHistory(AwhHistory *awhHistory) const;

        /*! \brief
         * Allocate and initialize an AWH history with the given AWH state.
         *
         * This function will be called at the start of a new simulation.
         * Note that only constant data will be initialized here.
         * History data is set by \ref Awh::updateHistory.
         *
         * \param[in,out] awhHistory  AWH history to initialize.
         */
        void initHistoryFromState(AwhHistory *awhHistory) const;

        /*! \brief Restore the AWH state from the given history.
         *
         * Should be called with a valid t_commrec on all ranks.
         * Should pass a valid awhHistory on the master rank.
         *
         * \param[in] awhHistory  AWH history to restore from.
         * \param[in] cr          Struct for communication.
         */
        void restoreStateFromHistory(const AwhHistory *awhHistory,
                                     const t_commrec  *cr);

        /*! \brief Returns string "AWH" for registering AWH as an external potential provider with the pull module.
         */
        static const char *externalPotentialString();

        /*! \brief Register the AWH biased coordinates with pull.
         *
         * This function is public because it needs to be called at
         * the pre-processing stage.
         * Pull requires all external potentials to register themselves
         * before the end of pre-processing and before the first MD step.
         *
         * \param[in]     awhParams  The AWH parameters.
         * \param[in,out] pull_work  Pull struct which AWH will register the bias into.
         */
        static void registerAwhWithPull(const AwhParams &awhParams,
                                        struct pull_t   *pull_work);

    private:
        std::vector<BiasCoupledToSystem> biasCoupledToSystem_; /**< AWH biases and definitions of their coupling to the system. */
        const gmx_int64_t                seed_;                /**< Random seed for MC jumping with umbrella type bias potential. */
        double                           potentialOffset_;     /**< The offset of the bias potential which changes due to bias updates. */
};

}      // namespace gmx

#endif /* GMX_AWH_H */
