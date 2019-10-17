/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Defines the high-level constraint code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "constr.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

namespace gmx
{

/* \brief Impl class for Constraints
 *
 * \todo Members like md, idef are valid only for the lifetime of a
 * domain, which would be good to make clearer in the structure of the
 * code. It should not be possible to call apply() if setConstraints()
 * has not been called. For example, this could be achieved if
 * setConstraints returned a valid object with such a method.  That
 * still requires that we manage the lifetime of that object
 * correctly, however. */
class Constraints::Impl
{
public:
    Impl(const gmx_mtop_t&     mtop_p,
         const t_inputrec&     ir_p,
         pull_t*               pull_work,
         FILE*                 log_p,
         const t_mdatoms&      md_p,
         const t_commrec*      cr_p,
         const gmx_multisim_t* ms,
         t_nrnb*               nrnb,
         gmx_wallcycle*        wcycle_p,
         bool                  pbcHandlingRequired,
         int                   numConstraints,
         int                   numSettles);
    ~Impl();
    void setConstraints(const gmx_localtop_t& top, const t_mdatoms& md);
    bool apply(bool               bLog,
               bool               bEner,
               int64_t            step,
               int                delta_step,
               real               step_scaling,
               rvec*              x,
               rvec*              xprime,
               rvec*              min_proj,
               const matrix       box,
               real               lambda,
               real*              dvdlambda,
               rvec*              v,
               tensor*            vir,
               ConstraintVariable econq);
    //! The total number of constraints.
    int ncon_tot = 0;
    //! The number of flexible constraints.
    int nflexcon = 0;
    //! A list of atoms to constraints for each moleculetype.
    std::vector<t_blocka> at2con_mt;
    //! A list of atoms to settles for each moleculetype
    std::vector<std::vector<int>> at2settle_mt;
    //! LINCS data.
    Lincs* lincsd = nullptr; // TODO this should become a unique_ptr
    //! SHAKE data.
    shakedata* shaked = nullptr;
    //! SETTLE data.
    settledata* settled = nullptr;
    //! The maximum number of warnings.
    int maxwarn = 0;
    //! The number of warnings for LINCS.
    int warncount_lincs = 0;
    //! The number of warnings for SETTLE.
    int warncount_settle = 0;
    //! The essential dynamics data.
    gmx_edsam* ed = nullptr;

    //! Thread-local virial contribution.
    tensor* vir_r_m_dr_th = { nullptr };
    //! Did a settle error occur?
    bool* bSettleErrorHasOccurred = nullptr;

    //! Pointer to the global topology - only used for printing warnings.
    const gmx_mtop_t& mtop;
    //! Parameters for the interactions in this domain.
    const t_idef* idef = nullptr;
    //! Data about atoms in this domain.
    const t_mdatoms& md;
    //! Whether we need to do pbc for handling bonds.
    bool pbcHandlingRequired_ = false;

    //! Logging support.
    FILE* log = nullptr;
    //! Communication support.
    const t_commrec* cr = nullptr;
    //! Multi-sim support.
    const gmx_multisim_t* ms = nullptr;
    //! Pulling code object, if any.
    pull_t* pull_work = nullptr;
    /*!\brief Input options.
     *
     * \todo Replace with IMdpOptions */
    const t_inputrec& ir;
    //! Flop counting support.
    t_nrnb* nrnb = nullptr;
    //! Tracks wallcycle usage.
    gmx_wallcycle* wcycle;
};

Constraints::~Constraints() = default;

int Constraints::numFlexibleConstraints() const
{
    return impl_->nflexcon;
}

bool Constraints::havePerturbedConstraints() const
{
    const gmx_ffparams_t& ffparams = impl_->mtop.ffparams;

    for (size_t i = 0; i < ffparams.functype.size(); i++)
    {
        if ((ffparams.functype[i] == F_CONSTR || ffparams.functype[i] == F_CONSTRNC)
            && ffparams.iparams[i].constr.dA != ffparams.iparams[i].constr.dB)
        {
            return true;
        }
    }

    return false;
}

//! Clears constraint quantities for atoms in nonlocal region.
static void clear_constraint_quantity_nonlocal(gmx_domdec_t* dd, rvec* q)
{
    int nonlocal_at_start, nonlocal_at_end, at;

    dd_get_constraint_range(dd, &nonlocal_at_start, &nonlocal_at_end);

    for (at = nonlocal_at_start; at < nonlocal_at_end; at++)
    {
        clear_rvec(q[at]);
    }
}

void too_many_constraint_warnings(int eConstrAlg, int warncount)
{
    gmx_fatal(
            FARGS,
            "Too many %s warnings (%d)\n"
            "If you know what you are doing you can %s"
            "set the environment variable GMX_MAXCONSTRWARN to -1,\n"
            "but normally it is better to fix the problem",
            (eConstrAlg == econtLINCS) ? "LINCS" : "SETTLE", warncount,
            (eConstrAlg == econtLINCS) ? "adjust the lincs warning threshold in your mdp file\nor " : "\n");
}

//! Writes out coordinates.
static void write_constr_pdb(const char*       fn,
                             const char*       title,
                             const gmx_mtop_t& mtop,
                             int               start,
                             int               homenr,
                             const t_commrec*  cr,
                             const rvec        x[],
                             const matrix      box)
{
    char          fname[STRLEN];
    FILE*         out;
    int           dd_ac0 = 0, dd_ac1 = 0, i, ii, resnr;
    gmx_domdec_t* dd;
    const char *  anm, *resnm;

    dd = nullptr;
    if (DOMAINDECOMP(cr))
    {
        dd = cr->dd;
        dd_get_constraint_range(dd, &dd_ac0, &dd_ac1);
        start  = 0;
        homenr = dd_ac1;
    }

    if (PAR(cr))
    {
        sprintf(fname, "%s_n%d.pdb", fn, cr->sim_nodeid);
    }
    else
    {
        sprintf(fname, "%s.pdb", fn);
    }

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, -1, box);
    int molb = 0;
    for (i = start; i < start + homenr; i++)
    {
        if (dd != nullptr)
        {
            if (i >= dd_numHomeAtoms(*dd) && i < dd_ac0)
            {
                continue;
            }
            ii = dd->globalAtomIndices[i];
        }
        else
        {
            ii = i;
        }
        mtopGetAtomAndResidueName(mtop, ii, &molb, &anm, &resnr, &resnm, nullptr);
        gmx_fprintf_pdb_atomline(out, epdbATOM, ii + 1, anm, ' ', resnm, ' ', resnr, ' ',
                                 10 * x[i][XX], 10 * x[i][YY], 10 * x[i][ZZ], 1.0, 0.0, "");
    }
    fprintf(out, "TER\n");

    gmx_fio_fclose(out);
}

//! Writes out domain contents to help diagnose crashes.
static void dump_confs(FILE*             log,
                       int64_t           step,
                       const gmx_mtop_t& mtop,
                       int               start,
                       int               homenr,
                       const t_commrec*  cr,
                       const rvec        x[],
                       rvec              xprime[],
                       const matrix      box)
{
    char buf[STRLEN], buf2[22];

    char* env = getenv("GMX_SUPPRESS_DUMP");
    if (env)
    {
        return;
    }

    sprintf(buf, "step%sb", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "initial coordinates", mtop, start, homenr, cr, x, box);
    sprintf(buf, "step%sc", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "coordinates after constraining", mtop, start, homenr, cr, xprime, box);
    if (log)
    {
        fprintf(log, "Wrote pdb files with previous and current coordinates\n");
    }
    fprintf(stderr, "Wrote pdb files with previous and current coordinates\n");
}

bool Constraints::apply(bool               bLog,
                        bool               bEner,
                        int64_t            step,
                        int                delta_step,
                        real               step_scaling,
                        rvec*              x,
                        rvec*              xprime,
                        rvec*              min_proj,
                        const matrix       box,
                        real               lambda,
                        real*              dvdlambda,
                        rvec*              v,
                        tensor*            vir,
                        ConstraintVariable econq)
{
    return impl_->apply(bLog, bEner, step, delta_step, step_scaling, x, xprime, min_proj, box,
                        lambda, dvdlambda, v, vir, econq);
}

bool Constraints::Impl::apply(bool               bLog,
                              bool               bEner,
                              int64_t            step,
                              int                delta_step,
                              real               step_scaling,
                              rvec*              x,
                              rvec*              xprime,
                              rvec*              min_proj,
                              const matrix       box,
                              real               lambda,
                              real*              dvdlambda,
                              rvec*              v,
                              tensor*            vir,
                              ConstraintVariable econq)
{
    bool   bOK, bDump;
    int    start, homenr;
    tensor vir_r_m_dr;
    real   scaled_delta_t;
    real   invdt, vir_fac = 0, t;
    int    nsettle;
    t_pbc  pbc, *pbc_null;
    char   buf[22];
    int    nth;

    wallcycle_start(wcycle, ewcCONSTR);

    if (econq == ConstraintVariable::ForceDispl && !EI_ENERGY_MINIMIZATION(ir.eI))
    {
        gmx_incons(
                "constrain called for forces displacements while not doing energy minimization, "
                "can not do this while the LINCS and SETTLE constraint connection matrices are "
                "mass weighted");
    }

    bOK   = TRUE;
    bDump = FALSE;

    start  = 0;
    homenr = md.homenr;

    scaled_delta_t = step_scaling * ir.delta_t;

    /* Prepare time step for use in constraint implementations, and
       avoid generating inf when ir.delta_t = 0. */
    if (ir.delta_t == 0)
    {
        invdt = 0.0;
    }
    else
    {
        invdt = 1.0 / scaled_delta_t;
    }

    if (ir.efep != efepNO && EI_DYNAMICS(ir.eI))
    {
        /* Set the constraint lengths for the step at which this configuration
         * is meant to be. The invmasses should not be changed.
         */
        lambda += delta_step * ir.fepvals->delta_lambda;
    }

    if (vir != nullptr)
    {
        clear_mat(vir_r_m_dr);
    }
    const t_ilist* settle = &idef->il[F_SETTLE];
    nsettle               = settle->nr / (1 + NRAL(F_SETTLE));

    if (nsettle > 0)
    {
        nth = gmx_omp_nthreads_get(emntSETTLE);
    }
    else
    {
        nth = 1;
    }

    /* We do not need full pbc when constraints do not cross update groups
     * i.e. when dd->constraint_comm==NULL.
     * Note that PBC for constraints is different from PBC for bondeds.
     * For constraints there is both forward and backward communication.
     */
    if (ir.ePBC != epbcNONE && (cr->dd || pbcHandlingRequired_)
        && !(cr->dd && cr->dd->constraint_comm == nullptr))
    {
        /* With pbc=screw the screw has been changed to a shift
         * by the constraint coordinate communication routine,
         * so that here we can use normal pbc.
         */
        pbc_null = set_pbc_dd(&pbc, ir.ePBC, DOMAINDECOMP(cr) ? cr->dd->nc : nullptr, FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    /* Communicate the coordinates required for the non-local constraints
     * for LINCS and/or SETTLE.
     */
    if (cr->dd)
    {
        dd_move_x_constraints(cr->dd, box, x, xprime, econq == ConstraintVariable::Positions);

        if (v != nullptr)
        {
            /* We need to initialize the non-local components of v.
             * We never actually use these values, but we do increment them,
             * so we should avoid uninitialized variables and overflows.
             */
            clear_constraint_quantity_nonlocal(cr->dd, v);
        }
    }

    if (lincsd != nullptr)
    {
        bOK = constrain_lincs(bLog || bEner, ir, step, lincsd, md, cr, ms, x, xprime, min_proj, box,
                              pbc_null, lambda, dvdlambda, invdt, v, vir != nullptr, vir_r_m_dr,
                              econq, nrnb, maxwarn, &warncount_lincs);
        if (!bOK && maxwarn < INT_MAX)
        {
            if (log != nullptr)
            {
                fprintf(log, "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtLINCS], gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }

    if (shaked != nullptr)
    {
        bOK = constrain_shake(log, shaked, md.invmass, *idef, ir, x, xprime, min_proj, nrnb, lambda,
                              dvdlambda, invdt, v, vir != nullptr, vir_r_m_dr, maxwarn < INT_MAX, econq);

        if (!bOK && maxwarn < INT_MAX)
        {
            if (log != nullptr)
            {
                fprintf(log, "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtSHAKE], gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }

    if (nsettle > 0)
    {
        bool bSettleErrorHasOccurred0 = false;

        switch (econq)
        {
            case ConstraintVariable::Positions:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (int th = 0; th < nth; th++)
                {
                    try
                    {
                        if (th > 0)
                        {
                            clear_mat(vir_r_m_dr_th[th]);
                        }

                        csettle(settled, nth, th, pbc_null, x[0], xprime[0], invdt, v ? v[0] : nullptr,
                                vir != nullptr, th == 0 ? vir_r_m_dr : vir_r_m_dr_th[th],
                                th == 0 ? &bSettleErrorHasOccurred0 : &bSettleErrorHasOccurred[th]);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                if (v != nullptr)
                {
                    inc_nrnb(nrnb, eNR_CONSTR_V, nsettle * 3);
                }
                if (vir != nullptr)
                {
                    inc_nrnb(nrnb, eNR_CONSTR_VIR, nsettle * 3);
                }
                break;
            case ConstraintVariable::Velocities:
            case ConstraintVariable::Derivative:
            case ConstraintVariable::Force:
            case ConstraintVariable::ForceDispl:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (int th = 0; th < nth; th++)
                {
                    try
                    {
                        int calcvir_atom_end;

                        if (vir == nullptr)
                        {
                            calcvir_atom_end = 0;
                        }
                        else
                        {
                            calcvir_atom_end = md.homenr;
                        }

                        if (th > 0)
                        {
                            clear_mat(vir_r_m_dr_th[th]);
                        }

                        int start_th = (nsettle * th) / nth;
                        int end_th   = (nsettle * (th + 1)) / nth;

                        if (start_th >= 0 && end_th - start_th > 0)
                        {
                            settle_proj(settled, econq, end_th - start_th,
                                        settle->iatoms + start_th * (1 + NRAL(F_SETTLE)), pbc_null,
                                        x, xprime, min_proj, calcvir_atom_end,
                                        th == 0 ? vir_r_m_dr : vir_r_m_dr_th[th]);
                        }
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                /* This is an overestimate */
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                break;
            case ConstraintVariable::Deriv_FlexCon:
                /* Nothing to do, since the are no flexible constraints in settles */
                break;
            default: gmx_incons("Unknown constraint quantity for settle");
        }

        if (vir != nullptr)
        {
            /* Reduce the virial contributions over the threads */
            for (int th = 1; th < nth; th++)
            {
                m_add(vir_r_m_dr, vir_r_m_dr_th[th], vir_r_m_dr);
            }
        }

        if (econq == ConstraintVariable::Positions)
        {
            for (int th = 1; th < nth; th++)
            {
                bSettleErrorHasOccurred0 = bSettleErrorHasOccurred0 || bSettleErrorHasOccurred[th];
            }

            if (bSettleErrorHasOccurred0)
            {
                char buf[STRLEN];
                sprintf(buf,
                        "\nstep "
                        "%" PRId64
                        ": One or more water molecules can not be settled.\n"
                        "Check for bad contacts and/or reduce the timestep if appropriate.\n",
                        step);
                if (log)
                {
                    fprintf(log, "%s", buf);
                }
                fprintf(stderr, "%s", buf);
                warncount_settle++;
                if (warncount_settle > maxwarn)
                {
                    too_many_constraint_warnings(-1, warncount_settle);
                }
                bDump = TRUE;

                bOK = FALSE;
            }
        }
    }

    if (vir != nullptr)
    {
        /* The normal uses of constrain() pass step_scaling = 1.0.
         * The call to constrain() for SD1 that passes step_scaling =
         * 0.5 also passes vir = NULL, so cannot reach this
         * assertion. This assertion should remain until someone knows
         * that this path works for their intended purpose, and then
         * they can use scaled_delta_t instead of ir.delta_t
         * below. */
        assert(gmx_within_tol(step_scaling, 1.0, GMX_REAL_EPS));
        switch (econq)
        {
            case ConstraintVariable::Positions: vir_fac = 0.5 / (ir.delta_t * ir.delta_t); break;
            case ConstraintVariable::Velocities: vir_fac = 0.5 / ir.delta_t; break;
            case ConstraintVariable::Force:
            case ConstraintVariable::ForceDispl: vir_fac = 0.5; break;
            default: gmx_incons("Unsupported constraint quantity for virial");
        }

        if (EI_VV(ir.eI))
        {
            vir_fac *= 2; /* only constraining over half the distance here */
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                (*vir)[i][j] = vir_fac * vir_r_m_dr[i][j];
            }
        }
    }

    if (bDump)
    {
        dump_confs(log, step, mtop, start, homenr, cr, x, xprime, box);
    }

    if (econq == ConstraintVariable::Positions)
    {
        if (ir.bPull && pull_have_constraint(pull_work))
        {
            if (EI_DYNAMICS(ir.eI))
            {
                t = ir.init_t + (step + delta_step) * ir.delta_t;
            }
            else
            {
                t = ir.init_t;
            }
            set_pbc(&pbc, ir.ePBC, box);
            pull_constraint(pull_work, &md, &pbc, cr, ir.delta_t, t, x, xprime, v, *vir);
        }
        if (ed && delta_step > 0)
        {
            /* apply the essential dynamics constraints here */
            do_edsam(&ir, step, cr, xprime, v, box, ed);
        }
    }
    wallcycle_stop(wcycle, ewcCONSTR);

    if (v != nullptr && md.cFREEZE)
    {
        /* Set the velocities of frozen dimensions to zero */

        int gmx_unused numThreads = gmx_omp_nthreads_get(emntUpdate);

#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int i = 0; i < md.homenr; i++)
        {
            int freezeGroup = md.cFREEZE[i];

            for (int d = 0; d < DIM; d++)
            {
                if (ir.opts.nFreeze[freezeGroup][d])
                {
                    v[i][d] = 0;
                }
            }
        }
    }

    return bOK;
}

ArrayRef<real> Constraints::rmsdData() const
{
    if (impl_->lincsd)
    {
        return lincs_rmsdData(impl_->lincsd);
    }
    else
    {
        return {};
    }
}

real Constraints::rmsd() const
{
    if (impl_->lincsd)
    {
        return lincs_rmsd(impl_->lincsd);
    }
    else
    {
        return 0;
    }
}

int Constraints::numConstraintsTotal()
{
    return impl_->ncon_tot;
}

FlexibleConstraintTreatment flexibleConstraintTreatment(bool haveDynamicsIntegrator)
{
    if (haveDynamicsIntegrator)
    {
        return FlexibleConstraintTreatment::Include;
    }
    else
    {
        return FlexibleConstraintTreatment::Exclude;
    }
}

/*! \brief Returns a block struct to go from atoms to constraints
 *
 * The block struct will contain constraint indices with lower indices
 * directly matching the order in F_CONSTR and higher indices matching
 * the order in F_CONSTRNC offset by the number of constraints in F_CONSTR.
 *
 * \param[in]  numAtoms  The number of atoms to construct the list for
 * \param[in]  ilists    The interaction lists, size F_NRE
 * \param[in]  iparams   Interaction parameters, can be null when
 * flexibleConstraintTreatment=Include \param[in]  flexibleConstraintTreatment  The flexible
 * constraint treatment, see enum above \returns a block struct with all constraints for each atom
 */
template<typename T>
static t_blocka makeAtomsToConstraintsList(int                         numAtoms,
                                           const T*                    ilists,
                                           const t_iparams*            iparams,
                                           FlexibleConstraintTreatment flexibleConstraintTreatment)
{
    GMX_ASSERT(flexibleConstraintTreatment == FlexibleConstraintTreatment::Include || iparams != nullptr,
               "With flexible constraint detection we need valid iparams");

    std::vector<int> count(numAtoms);

    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const T&  ilist  = ilists[ftype];
        const int stride = 1 + NRAL(ftype);
        for (int i = 0; i < ilist.size(); i += stride)
        {
            if (flexibleConstraintTreatment == FlexibleConstraintTreatment::Include
                || !isConstraintFlexible(iparams, ilist.iatoms[i]))
            {
                for (int j = 1; j < 3; j++)
                {
                    int a = ilist.iatoms[i + j];
                    count[a]++;
                }
            }
        }
    }

    t_blocka at2con;
    at2con.nr           = numAtoms;
    at2con.nalloc_index = at2con.nr + 1;
    snew(at2con.index, at2con.nalloc_index);
    at2con.index[0] = 0;
    for (int a = 0; a < numAtoms; a++)
    {
        at2con.index[a + 1] = at2con.index[a] + count[a];
        count[a]            = 0;
    }
    at2con.nra      = at2con.index[at2con.nr];
    at2con.nalloc_a = at2con.nra;
    snew(at2con.a, at2con.nalloc_a);

    /* The F_CONSTRNC constraints have constraint numbers
     * that continue after the last F_CONSTR constraint.
     */
    int numConstraints = 0;
    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const T&  ilist  = ilists[ftype];
        const int stride = 1 + NRAL(ftype);
        for (int i = 0; i < ilist.size(); i += stride)
        {
            if (flexibleConstraintTreatment == FlexibleConstraintTreatment::Include
                || !isConstraintFlexible(iparams, ilist.iatoms[i]))
            {
                for (int j = 1; j < 3; j++)
                {
                    int a                                  = ilist.iatoms[i + j];
                    at2con.a[at2con.index[a] + count[a]++] = numConstraints;
                }
            }
            numConstraints++;
        }
    }

    return at2con;
}

t_blocka make_at2con(int                         numAtoms,
                     const t_ilist*              ilist,
                     const t_iparams*            iparams,
                     FlexibleConstraintTreatment flexibleConstraintTreatment)
{
    return makeAtomsToConstraintsList(numAtoms, ilist, iparams, flexibleConstraintTreatment);
}

t_blocka make_at2con(const gmx_moltype_t&           moltype,
                     gmx::ArrayRef<const t_iparams> iparams,
                     FlexibleConstraintTreatment    flexibleConstraintTreatment)
{
    return makeAtomsToConstraintsList(moltype.atoms.nr, moltype.ilist.data(), iparams.data(),
                                      flexibleConstraintTreatment);
}

//! Return the number of flexible constraints in the \c ilist and \c iparams.
template<typename T>
static int countFlexibleConstraintsTemplate(const T* ilist, const t_iparams* iparams)
{
    int nflexcon = 0;
    for (int ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        const int numIatomsPerConstraint = 3;
        for (int i = 0; i < ilist[ftype].size(); i += numIatomsPerConstraint)
        {
            const int type = ilist[ftype].iatoms[i];
            if (iparams[type].constr.dA == 0 && iparams[type].constr.dB == 0)
            {
                nflexcon++;
            }
        }
    }

    return nflexcon;
}

int countFlexibleConstraints(const t_ilist* ilist, const t_iparams* iparams)
{
    return countFlexibleConstraintsTemplate(ilist, iparams);
}

//! Returns the index of the settle to which each atom belongs.
static std::vector<int> make_at2settle(int natoms, const InteractionList& ilist)
{
    /* Set all to no settle */
    std::vector<int> at2s(natoms, -1);

    const int stride = 1 + NRAL(F_SETTLE);

    for (int s = 0; s < ilist.size(); s += stride)
    {
        at2s[ilist.iatoms[s + 1]] = s / stride;
        at2s[ilist.iatoms[s + 2]] = s / stride;
        at2s[ilist.iatoms[s + 3]] = s / stride;
    }

    return at2s;
}

void Constraints::Impl::setConstraints(const gmx_localtop_t& top, const t_mdatoms& md)
{
    idef = &top.idef;

    if (ncon_tot > 0)
    {
        /* With DD we might also need to call LINCS on a domain no constraints for
         * communicating coordinates to other nodes that do have constraints.
         */
        if (ir.eConstrAlg == econtLINCS)
        {
            set_lincs(top.idef, md, EI_DYNAMICS(ir.eI), cr, lincsd);
        }
        if (ir.eConstrAlg == econtSHAKE)
        {
            if (cr->dd)
            {
                // We are using the local topology, so there are only
                // F_CONSTR constraints.
                make_shake_sblock_dd(shaked, &idef->il[F_CONSTR], cr->dd);
            }
            else
            {
                make_shake_sblock_serial(shaked, idef, md);
            }
        }
    }

    if (settled)
    {
        settle_set_constraints(settled, &idef->il[F_SETTLE], md);
    }

    /* Make a selection of the local atoms for essential dynamics */
    if (ed && cr->dd)
    {
        dd_make_local_ed_indices(cr->dd, ed);
    }
}

void Constraints::setConstraints(const gmx_localtop_t& top, const t_mdatoms& md)
{
    impl_->setConstraints(top, md);
}

/*! \brief Makes a per-moleculetype container of mappings from atom
 * indices to constraint indices.
 *
 * Note that flexible constraints are only enabled with a dynamical integrator. */
static std::vector<t_blocka> makeAtomToConstraintMappings(const gmx_mtop_t& mtop,
                                                          FlexibleConstraintTreatment flexibleConstraintTreatment)
{
    std::vector<t_blocka> mapping;
    mapping.reserve(mtop.moltype.size());
    for (const gmx_moltype_t& moltype : mtop.moltype)
    {
        mapping.push_back(make_at2con(moltype, mtop.ffparams.iparams, flexibleConstraintTreatment));
    }
    return mapping;
}

Constraints::Constraints(const gmx_mtop_t&     mtop,
                         const t_inputrec&     ir,
                         pull_t*               pull_work,
                         FILE*                 log,
                         const t_mdatoms&      md,
                         const t_commrec*      cr,
                         const gmx_multisim_t* ms,
                         t_nrnb*               nrnb,
                         gmx_wallcycle*        wcycle,
                         bool                  pbcHandlingRequired,
                         int                   numConstraints,
                         int                   numSettles) :
    impl_(new Impl(mtop, ir, pull_work, log, md, cr, ms, nrnb, wcycle, pbcHandlingRequired, numConstraints, numSettles))
{
}

Constraints::Impl::Impl(const gmx_mtop_t&     mtop_p,
                        const t_inputrec&     ir_p,
                        pull_t*               pull_work,
                        FILE*                 log_p,
                        const t_mdatoms&      md_p,
                        const t_commrec*      cr_p,
                        const gmx_multisim_t* ms_p,
                        t_nrnb*               nrnb_p,
                        gmx_wallcycle*        wcycle_p,
                        bool                  pbcHandlingRequired,
                        int                   numConstraints,
                        int                   numSettles) :
    ncon_tot(numConstraints),
    mtop(mtop_p),
    md(md_p),
    pbcHandlingRequired_(pbcHandlingRequired),
    log(log_p),
    cr(cr_p),
    ms(ms_p),
    pull_work(pull_work),
    ir(ir_p),
    nrnb(nrnb_p),
    wcycle(wcycle_p)
{
    if (numConstraints + numSettles > 0 && ir.epc == epcMTTK)
    {
        gmx_fatal(FARGS, "Constraints are not implemented with MTTK pressure control.");
    }

    nflexcon = 0;
    if (numConstraints > 0)
    {
        at2con_mt = makeAtomToConstraintMappings(mtop, flexibleConstraintTreatment(EI_DYNAMICS(ir.eI)));

        for (const gmx_molblock_t& molblock : mtop.molblock)
        {
            int count = countFlexibleConstraintsTemplate(mtop.moltype[molblock.type].ilist.data(),
                                                         mtop.ffparams.iparams.data());
            nflexcon += molblock.nmol * count;
        }

        if (nflexcon > 0)
        {
            if (log)
            {
                fprintf(log, "There are %d flexible constraints\n", nflexcon);
                if (ir.fc_stepsize == 0)
                {
                    fprintf(log,
                            "\n"
                            "WARNING: step size for flexible constraining = 0\n"
                            "         All flexible constraints will be rigid.\n"
                            "         Will try to keep all flexible constraints at their original "
                            "length,\n"
                            "         but the lengths may exhibit some drift.\n\n");
                    nflexcon = 0;
                }
            }
            if (nflexcon > 0)
            {
                please_cite(log, "Hess2002");
            }
        }

        if (ir.eConstrAlg == econtLINCS)
        {
            lincsd = init_lincs(log, mtop, nflexcon, at2con_mt,
                                DOMAINDECOMP(cr) && ddHaveSplitConstraints(*cr->dd), ir.nLincsIter,
                                ir.nProjOrder);
        }

        if (ir.eConstrAlg == econtSHAKE)
        {
            if (DOMAINDECOMP(cr) && ddHaveSplitConstraints(*cr->dd))
            {
                gmx_fatal(FARGS,
                          "SHAKE is not supported with domain decomposition and constraint that "
                          "cross domain boundaries, use LINCS");
            }
            if (nflexcon)
            {
                gmx_fatal(FARGS,
                          "For this system also velocities and/or forces need to be constrained, "
                          "this can not be done with SHAKE, you should select LINCS");
            }
            please_cite(log, "Ryckaert77a");
            if (ir.bShakeSOR)
            {
                please_cite(log, "Barth95a");
            }

            shaked = shake_init();
        }
    }

    if (numSettles > 0)
    {
        please_cite(log, "Miyamoto92a");

        settled = settle_init(mtop);

        /* Make an atom to settle index for use in domain decomposition */
        for (size_t mt = 0; mt < mtop.moltype.size(); mt++)
        {
            at2settle_mt.emplace_back(
                    make_at2settle(mtop.moltype[mt].atoms.nr, mtop.moltype[mt].ilist[F_SETTLE]));
        }

        /* Allocate thread-local work arrays */
        int nthreads = gmx_omp_nthreads_get(emntSETTLE);
        if (nthreads > 1 && vir_r_m_dr_th == nullptr)
        {
            snew(vir_r_m_dr_th, nthreads);
            snew(bSettleErrorHasOccurred, nthreads);
        }
    }

    maxwarn   = 999;
    char* env = getenv("GMX_MAXCONSTRWARN");
    if (env)
    {
        maxwarn = 0;
        sscanf(env, "%8d", &maxwarn);
        if (maxwarn < 0)
        {
            maxwarn = INT_MAX;
        }
        if (log)
        {
            fprintf(log, "Setting the maximum number of constraint warnings to %d\n", maxwarn);
        }
        if (MASTER(cr))
        {
            fprintf(stderr, "Setting the maximum number of constraint warnings to %d\n", maxwarn);
        }
    }
    warncount_lincs  = 0;
    warncount_settle = 0;
}

Constraints::Impl::~Impl()
{
    for (auto blocka : at2con_mt)
    {
        done_blocka(&blocka);
    }
    if (bSettleErrorHasOccurred != nullptr)
    {
        sfree(bSettleErrorHasOccurred);
    }
    if (vir_r_m_dr_th != nullptr)
    {
        sfree(vir_r_m_dr_th);
    }
    if (settled != nullptr)
    {
        settle_free(settled);
    }
    done_lincs(lincsd);
}

void Constraints::saveEdsamPointer(gmx_edsam* ed)
{
    impl_->ed = ed;
}

ArrayRef<const t_blocka> Constraints::atom2constraints_moltype() const
{
    return impl_->at2con_mt;
}

ArrayRef<const std::vector<int>> Constraints::atom2settle_moltype() const
{
    return impl_->at2settle_mt;
}

void do_constrain_first(FILE*                     fplog,
                        gmx::Constraints*         constr,
                        const t_inputrec*         ir,
                        const t_mdatoms*          md,
                        int                       natoms,
                        ArrayRefWithPadding<RVec> x,
                        ArrayRefWithPadding<RVec> v,
                        const matrix              box,
                        real                      lambda)
{
    int     i, m, start, end;
    int64_t step;
    real    dt = ir->delta_t;
    real    dvdl_dum;
    rvec*   savex;

    auto xRvec = as_rvec_array(x.paddedArrayRef().data());
    auto vRvec = as_rvec_array(v.paddedArrayRef().data());

    /* We need to allocate one element extra, since we might use
     * (unaligned) 4-wide SIMD loads to access rvec entries.
     */
    snew(savex, natoms + 1);

    start = 0;
    end   = md->homenr;

    if (debug)
    {
        fprintf(debug, "vcm: start=%d, homenr=%d, end=%d\n", start, md->homenr, end);
    }
    /* Do a first constrain to reset particles... */
    step = ir->init_step;
    if (fplog)
    {
        char buf[STEPSTRSIZE];
        fprintf(fplog, "\nConstraining the starting coordinates (step %s)\n", gmx_step_str(step, buf));
    }
    dvdl_dum = 0;

    /* constrain the current position */
    constr->apply(TRUE, FALSE, step, 0, 1.0, xRvec, xRvec, nullptr, box, lambda, &dvdl_dum, nullptr,
                  nullptr, gmx::ConstraintVariable::Positions);
    if (EI_VV(ir->eI))
    {
        /* constrain the inital velocity, and save it */
        /* also may be useful if we need the ekin from the halfstep for velocity verlet */
        constr->apply(TRUE, FALSE, step, 0, 1.0, xRvec, vRvec, vRvec, box, lambda, &dvdl_dum,
                      nullptr, nullptr, gmx::ConstraintVariable::Velocities);
    }
    /* constrain the inital velocities at t-dt/2 */
    if (EI_STATE_VELOCITY(ir->eI) && ir->eI != eiVV)
    {
        auto subX = x.paddedArrayRef().subArray(start, end);
        auto subV = v.paddedArrayRef().subArray(start, end);
        for (i = start; (i < end); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                /* Reverse the velocity */
                subV[i][m] = -subV[i][m];
                /* Store the position at t-dt in buf */
                savex[i][m] = subX[i][m] + dt * subV[i][m];
            }
        }
        /* Shake the positions at t=-dt with the positions at t=0
         * as reference coordinates.
         */
        if (fplog)
        {
            char buf[STEPSTRSIZE];
            fprintf(fplog, "\nConstraining the coordinates at t0-dt (step %s)\n", gmx_step_str(step, buf));
        }
        dvdl_dum = 0;
        constr->apply(TRUE, FALSE, step, -1, 1.0, xRvec, savex, nullptr, box, lambda, &dvdl_dum,
                      vRvec, nullptr, gmx::ConstraintVariable::Positions);

        for (i = start; i < end; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                /* Re-reverse the velocities */
                subV[i][m] = -subV[i][m];
            }
        }
    }
    sfree(savex);
}

} // namespace gmx
