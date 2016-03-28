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
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_multisim_t;
struct gmx_wallcycle;
struct pull_work_t;
struct pull_t;
class t_state;
struct t_commrec;
struct t_enxframe;
struct t_inputrec;
struct t_mdatoms;

namespace gmx
{

struct AwhHistory;
struct AwhParams;
class Bias;
struct BiasCoupledToSystem;
class ForceWithVirial;

/*! \libinternal
 * \brief Coupling of the accelerated weight histogram method (AWH) with the system.
 *
 * AWH calculates the free energy along order parameters of the system.
 * Free energy barriers are overcome by adaptively tuning a bias potential along
 * the order parameter such that the biased distribution along the parameter
 * converges toward a chosen target distribution.
 *
 * The Awh class takes care of the coupling between the system and the AWH
 * bias(es). The Awh class contains one or more BiasCoupledToSystem objects.
 * The BiasCoupledToSystem class takes care of the reaction coordinate input
 * and force output for the single Bias object it containts.
 *
 * \todo Update parameter reading and checkpointing, when general C++ framework is ready.
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
         * \param[in,out] fplog             General output file, normally md.log, can be nullptr.
         * \param[in]     inputRecord       General input parameters (as set up by grompp).
         * \param[in]     commRecord        Struct for communication, can be nullptr.
         * \param[in]     awhParams         AWH input parameters, consistent with the relevant parts of \p inputRecord (as set up by grompp).
         * \param[in]     biasInitFilename  Name of file to read PMF and target from.
         * \param[in,out] pull_work         Pointer to a pull struct which AWH will couple to, has to be initialized, is assumed not to change during the lifetime of the Awh object.
         */
        Awh(FILE              *fplog,
            const t_inputrec  &inputRecord,
            const t_commrec   *commRecord,
            const AwhParams   &awhParams,
            const std::string &biasInitFilename,
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
         * since AWH needs the current coordinate values (the pull code checks
         * for this).
         *
         * \param[in]     mdatoms          Atom properties.
         * \param[in]     ePBC             Type of periodic boundary conditions.
         * \param[in]     box              Box vectors.
         * \param[in,out] forceWithVirial  Force and virial buffers, should cover at least the local atoms.
         * \param[in]     t                Time.
         * \param[in]     step             The current MD step.
         * \param[in,out] wallcycle        Wallcycle counter, can be nullptr.
         * \param[in,out] fplog            General output file, normally md.log, can be nullptr.
         * \returns the potential energy for the bias.
         */
        real applyBiasForcesAndUpdateBias(int                     ePBC,
                                          const t_mdatoms        &mdatoms,
                                          const matrix            box,
                                          gmx::ForceWithVirial   *forceWithVirial,
                                          double                  t,
                                          gmx_int64_t             step,
                                          gmx_wallcycle          *wallcycle,
                                          FILE                   *fplog);

        /*! \brief
         * Update the AWH history in preparation for writing to checkpoint file.
         *
         * Should be called at least on the master rank at checkpoint steps.
         *
         * Should be called with a valid \p awhHistory (is checked).
         *
         * \param[in,out] awhHistory  AWH history to set.
         */
        void updateHistory(AwhHistory *awhHistory) const;

        /*! \brief
         * Allocate and initialize an AWH history with the given AWH state.
         *
         * This function should be called at the start of a new simulation
         * at least on the master rank.
         * Note that only constant data will be initialized here.
         * History data is set by \ref Awh::updateHistory.
         *
         * \returns a shared pointer to the AWH history on the rank that does I/O, nullptr otherwise.
         */
        std::shared_ptr<AwhHistory> initHistoryFromState() const;

        /*! \brief Restore the AWH state from the given history.
         *
         * Should be called on all ranks (for internal MPI broadcast).
         * Should pass a point to an AwhHistory on the master rank that
         * is compatible with the AWH setup in this simulation. Will throw
         * an exception if it is not compatible.
         *
         * \param[in] awhHistory  AWH history to restore from.
         */
        void restoreStateFromHistory(const AwhHistory *awhHistory);

        /*! \brief Fills the AWH data block of an energy frame with data at certain steps.
         *
         * \param[in]     step  The current MD step.
         * \param[in,out] fr    Energy data frame.
         */
        void writeToEnergyFrame(gmx_int64_t  step,
                                t_enxframe  *fr) const;

        /*! \brief Returns string "AWH" for registering AWH as an external potential provider with the pull module.
         */
        static const char *externalPotentialString();

        /*! \brief Register the AWH biased coordinates with pull.
         *
         * This function is public because it needs to be called by grompp
         * (and is otherwise only called by Awh()).
         * Pull requires all external potentials to register themselves
         * before the end of pre-processing and before the first MD step.
         * If this has not happened, pull with throw an error.
         *
         * \param[in]     awhParams  The AWH parameters.
         * \param[in,out] pull_work  Pull struct which AWH will register the bias into.
         */
        static void registerAwhWithPull(const AwhParams &awhParams,
                                        pull_t          *pull_work);

    private:
        /*! \brief Returns whether we need to write output at the current step.
         *
         * \param[in]     step             The current MD step.
         */
        bool isOutputStep(gmx_int64_t step) const;

    private:
        std::vector<BiasCoupledToSystem> biasCoupledToSystem_; /**< AWH biases and definitions of their coupling to the system. */
        const gmx_int64_t                seed_;                /**< Random seed for MC jumping with umbrella type bias potential. */
        const int                        nstout_;              /**< Interval in steps for writing to energy file. */
        const t_commrec                 *commRecord_;          /**< Pointer to the communication record. */
        pull_t                          *pull_;                /**< Pointer to the pull working data. */
        double                           potentialOffset_;     /**< The offset of the bias potential which changes due to bias updates. */
};

}      // namespace gmx

#endif /* GMX_AWH_H */
