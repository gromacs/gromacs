/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
//
// Created by Eric Irrgang on 8/17/18.
//

#ifndef GROMACS_MD_IMPL_H
#define GROMACS_MD_IMPL_H

#include "md.h"

namespace gmx
{

class MDIntegrator::Impl
{
    public:
        /*!
         * \brief Extract non-owning fields from parameter structure.
         *
         * \param container previously initialized parameter structure.
         */
        explicit Impl(const IntegratorParamsContainer &container);

        ~Impl();

        /*!
         * \brief Provide the run() behavior for MDIntegrator.
         */
        void run();

        //! Handles logging.
        FILE                            *fplog;
        //! Handles communication.
        t_commrec                       *cr;
        //! Coordinates multi-simulations.
        const gmx_multisim_t            *ms;
        //! Handles logging.
        const MDLogger                  &mdlog;
        //! Count of input file options.
        int                              nfile;
        //! Content of input file options.
        const t_filenm                  *fnm;
        //! Handles writing text output.
        const gmx_output_env_t          *oenv;
        //! Contains command-line options to mdrun.
        const MdrunOptions              &mdrunOptions;
        //! Handles virtual sites.
        gmx_vsite_t                     *vsite;
        //! Handles constraints.
        Constraints                     *constr;
        //! Handles enforced rotation.
        gmx_enfrot                      *enforcedRotation;
        //! Handles box deformation.
        BoxDeformation                  *deform;
        //! Handles writing output files.
        IMDOutputProvider               *outputProvider;
        //! Contains user input mdp options.
        t_inputrec                      *inputrec;
        //! Full system topology.
        gmx_mtop_t                      *top_global;
        //! Helper struct for force calculations.
        t_fcdata                        *fcd;
        //! Full simulation state (only non-nullptr on master rank).
        t_state                         *state_global;
        //! History of simulation observables.
        ObservablesHistory              *observablesHistory;
        //! Atom parameters for this domain.
        MDAtoms                         *mdAtoms;
        //! Manages flop accounting.
        t_nrnb                          *nrnb;
        //! Manages wall cycle accounting.
        gmx_wallcycle                   *wcycle;
        //! Parameters for force calculations.
        t_forcerec                      *fr;
        //! Parameters for replica exchange algorihtms.
        const ReplicaExchangeParameters &replExParams;
        //! Parameters for membrane embedding.
        gmx_membed_t                    *membed;
        //! Manages wall time accounting.
        gmx_walltime_accounting         *walltime_accounting;


        std::unique_ptr<md::Context> context_ {nullptr};
};

}      // end namespace gmx

#endif //GROMACS_MD_IMPL_H
