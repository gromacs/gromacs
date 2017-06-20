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
/*! \libinternal \brief Declares the integrator type for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_INTEGRATOR_H
#define GMX_MDLIB_INTEGRATOR_H

#include <cstdio>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class energyhistory_t;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_output_env_t;
struct ObservablesHistory;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
class t_state;

namespace gmx
{

class IMDOutputProvider;
class MDLogger;

/*! \libinternal \brief Integrator parameters container.
 *
 * A container for parameters to integrator_t functions. Allows simpler
 * maintenance of multiple integrator implementations. Currently, the container
 * *does not own* most of its data and must not construct or destruct pointed-to
 * objects, nor should it be allowed to live for more than a line or two.
 * Further encapsulation requires clarification of data ownership and
 * access policies.
 *
 * \ingroup module_mdlib
 */
class IntegratorParams
{
    public:
        /* \libinternal \brief Integrator algorithm implementation parameters.
         *
         * \param[in] fplog               Log file for output
         * \param[in] commrec                  Communication record
         * \param[in] mdlog               Log writer for important output
         * \param[in] nfile               Number of files
         * \param[in] fnm                 Filename structure array
         * \param[in] oenv                Output information
         * \param[in] bVerbose            Verbose output or not
         * \param[in] nstglobalcomm       How often global communication is done
         * \param[in] vsite               Virtual site information
         * \param[in] constr              Constraint information
         * \param[in] stepout             How often we writen to the console
         * \param[in] outputProvider      Additional output provider
         * \param[in] inputrec            Input record with mdp options
         * \param[in] top_global          Molecular topology for the whole system
         * \param[in] fcd                 Force and constraint data
         * \param[in] state_global        The state (x, v, f, box etc.) of the whole system
         * \param[in] observablesHistory  The observables statistics history
         * \param[in] mdatoms             Structure containing atom information
         * \param[in] nrnb                Accounting for floating point operations
         * \param[in] wcycle              Wall cycle timing information
         * \param[in] fr                  Force record with cut-off information and more
         * \param[in] replExParams        Parameters for the replica exchange algorithm
         * \param[in] membed              Membrane embedding data structure
         * \param[in] cpt_period          How often to checkpoint the simulation
         * \param[in] max_hours           Maximum length of the simulation (wall time)
         * \param[in] imdport             Interactive MD port (socket)
         * \param[in] Flags               Flags to control mdrun
         * \param[in] walltime_accounting More timing information
         * \ingroup module_mdlib
         */
        IntegratorParams(FILE *fplog, t_commrec *commrec, const gmx::MDLogger &mdlog,
                         int nfile, const t_filenm* fnm,
                         const gmx_output_env_t *oenv, gmx_bool bVerbose,
                         int nstglobalcomm,
                         gmx_vsite_t *vsite, gmx_constr_t constr,
                         int stepout,
                         gmx::IMDOutputProvider *outputProvider,
                         t_inputrec *inputrec,
                         gmx_mtop_t *top_global, t_fcdata *fcd,
                         t_state *state_global,
                         ObservablesHistory *observablesHistory,
                         t_mdatoms *mdatoms,
                         t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                         t_forcerec *fr,
                         const ReplicaExchangeParameters &replExParams,
                         gmx_membed_t gmx_unused * membed,
                         real cpt_period, real max_hours,
                         int imdport,
                         unsigned long Flags,
                         gmx_walltime_accounting_t walltime_accounting) :
            fpLog_ {fplog},
        commRec_ {commrec},
        mdLog_ {&mdlog},
        nFile_ {nfile},
        fnm_ {fnm},
        oenv_ {oenv},
        verbose_ {bVerbose},
        nstGlobalComm_ {nstglobalcomm},
        vSite_ {vsite},
        constraints_ {constr},
        stepOut_ {stepout},
        outputProvider_ {outputProvider},
        inputRec_ {inputrec},
        topGlobal_ {top_global},
        fcd_ {fcd},
        stateGlobal_ {state_global},
        observablesHistory_ {observablesHistory},
        mdAtoms_ {mdatoms},
        nrnb_ {nrnb},
        wCycle_ {wcycle},
        forceRec_ {fr},
        replExParams_ {&replExParams},
        membed_ {membed},
        cptPeriod_ {cpt_period},
        maxHours_ {max_hours},
        imdPort_ {imdport},
        flags_ {Flags},
        walltimeAccounting_ {walltime_accounting}
        {};

        // The objects don't "own" their data, so they can be copied or assigned
        //just as easily as they can be constructed, but probably shouldn't be.
        /// Default copy constructor.
        IntegratorParams(const IntegratorParams &) = default;

        // In the first implementation, IntegratorParams doesn't own any of its
        // pointers and must not delete any of them.
        /// Default destructor.
        ~IntegratorParams() = default;

        /// Get log file pointer.
        FILE* fpLog() const { return fpLog_; };

        /// communicator record.
        t_commrec* commRec() const { return commRec_; };

        /// MDLogger instance
        const gmx::MDLogger &mdLog() const { return *mdLog_; };

        /// number of files in filename structure array.
        int nFile() const { return nFile_; };

        /// filename structure array.
        const t_filenm* fnm() const { return fnm_; };

        /// Output environment.
        const gmx_output_env_t* oenv() const { return oenv_; };

        /// Current setting for verbose output.
        gmx_bool verbose() const { return verbose_; };

        /// Interval of global communication.
        int nstGlobalComm() const { return nstGlobalComm_; };

        /// Virtual site information.
        gmx_vsite_t* vSite() const { return vSite_; };

        /// Constraint information.
        gmx_constr_t constraints() const { return constraints_; };

        /// Console output interval.
        int stepOut() const { return stepOut_; };

        /// Additional output provider.
        gmx::IMDOutputProvider* outputProvider() const { return outputProvider_; };

        /// Input record with mdp options.
        t_inputrec* inputRec() const { return inputRec_; };

        /// Molecular topology for the whole system.
        gmx_mtop_t* topGlobal() const { return topGlobal_; };

        /// Force and constraint data.
        t_fcdata* fcd() const { return fcd_; };

        /// The state (x, v, f, box etc.) of the whole system.
        t_state* stateGlobal() const { return stateGlobal_; };

        /// The observables statistics history.
        ObservablesHistory* observablesHistory() const { return observablesHistory_; };

        /// Structure containing atom information.
        t_mdatoms* mdAtoms() const { return mdAtoms_; };

        /// Accounting for floating point operations.
        t_nrnb* nrnb() const { return nrnb_; };

        /// Wall cycle timing information.
        gmx_wallcycle_t wCycle() const { return wCycle_; };

        /// Force record with cut-off information and more.
        t_forcerec* forceRec() const { return forceRec_; };

        /// Parameters for the replica exchange algorithm.
        const ReplicaExchangeParameters &replExParams() const { return *replExParams_; };

        /// Membrane embedding data structure.
        gmx_membed_t* membed() const { return membed_; };

        /// How often to checkpoint the simulation.
        real cptPeriod() const { return cptPeriod_; };

        /// Maximum length of the simulation (wall time).
        real maxHours() const { return maxHours_; };

        /// Interactive MD port (socket)
        int imdPort() const { return imdPort_; };

        /// Flags to control mdrun.
        unsigned long flags() const { return flags_; };

        /// More timing information.
        gmx_walltime_accounting_t walltimeAccounting() const { return walltimeAccounting_; };

    private:
        /// \cond
        // These private members don't need to be in doxygen...
        FILE                            *fpLog_;
        t_commrec                       *commRec_;
        const gmx::MDLogger             *mdLog_;
        int                              nFile_;
        const t_filenm                 * fnm_;
        const gmx_output_env_t          *oenv_;
        gmx_bool                         verbose_;
        int                              nstGlobalComm_;
        gmx_vsite_t                     *vSite_;
        // not sure why constraint struct is passed by value or whether it should be.
        gmx_constr_t                     constraints_;
        int                              stepOut_;
        gmx::IMDOutputProvider          *outputProvider_;
        t_inputrec                      *inputRec_;
        gmx_mtop_t                      *topGlobal_;
        t_fcdata                        *fcd_;
        t_state                         *stateGlobal_;
        ObservablesHistory              *observablesHistory_;
        t_mdatoms                       *mdAtoms_;
        t_nrnb                          *nrnb_;
        gmx_wallcycle_t                  wCycle_;
        t_forcerec                      *forceRec_;
        const ReplicaExchangeParameters* replExParams_;
        gmx_membed_t                    *membed_;
        real                             cptPeriod_;
        real                             maxHours_;
        int                              imdPort_;
        unsigned long                    flags_;
        gmx_walltime_accounting_t        walltimeAccounting_;
        /// \endcond
};

/*! \libinternal \brief Integrator algorithm implementation.
 *
 * Provide a typedef with which to declare functions with the signature of an
 * integrator implementation. This is an intermediate solution. The next step
 * of encapsulation would be turning functions with this signature into functors
 * of classes implementing a common integrator_t-like interface.
 * \ingroup module_mdlib
 */
typedef double integrator_t (const IntegratorParams &params);
}      // namespace gmx

#endif // GMX_MDLIB_INTEGRATOR_H
