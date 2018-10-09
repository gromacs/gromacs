/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "readir.h"

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <string>

#include "gromacs/awh/read-params.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxlib/chargegroup.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/keyvaluetreemdpwriter.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/ikeyvaluetreeerror.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#define MAXPTR 254
#define NOGID  255

/* Resource parameters
 * Do not change any of these until you read the instruction
 * in readinp.h. Some cpp's do not take spaces after the backslash
 * (like the c-shell), which will give you a very weird compiler
 * message.
 */

typedef struct t_inputrec_strings
{
    char tcgrps[STRLEN], tau_t[STRLEN], ref_t[STRLEN],
         acc[STRLEN], accgrps[STRLEN], freeze[STRLEN], frdim[STRLEN],
         energy[STRLEN], user1[STRLEN], user2[STRLEN], vcm[STRLEN], x_compressed_groups[STRLEN],
         couple_moltype[STRLEN], orirefitgrp[STRLEN], egptable[STRLEN], egpexcl[STRLEN],
         wall_atomtype[STRLEN], wall_density[STRLEN], deform[STRLEN], QMMM[STRLEN],
         imd_grp[STRLEN];
    char   fep_lambda[efptNR][STRLEN];
    char   lambda_weights[STRLEN];
    char **pull_grp;
    char **rot_grp;
    char   anneal[STRLEN], anneal_npoints[STRLEN],
           anneal_time[STRLEN], anneal_temp[STRLEN];
    char   QMmethod[STRLEN], QMbasis[STRLEN], QMcharge[STRLEN], QMmult[STRLEN],
           bSH[STRLEN], CASorbitals[STRLEN], CASelectrons[STRLEN], SAon[STRLEN],
           SAoff[STRLEN], SAsteps[STRLEN];

} gmx_inputrec_strings;

static gmx_inputrec_strings *is = nullptr;

void init_inputrec_strings()
{
    if (is)
    {
        gmx_incons("Attempted to call init_inputrec_strings before calling done_inputrec_strings. Only one inputrec (i.e. .mdp file) can be parsed at a time.");
    }
    snew(is, 1);
}

void done_inputrec_strings()
{
    sfree(is);
    is = nullptr;
}


enum {
    egrptpALL,         /* All particles have to be a member of a group.     */
    egrptpALL_GENREST, /* A rest group with name is generated for particles *
                        * that are not part of any group.                   */
    egrptpPART,        /* As egrptpALL_GENREST, but no name is generated    *
                        * for the rest group.                               */
    egrptpONE          /* Merge all selected groups into one group,         *
                        * make a rest group for the remaining particles.    */
};

static const char *constraints[eshNR+1]    = {
    "none", "h-bonds", "all-bonds", "h-angles", "all-angles", nullptr
};

static const char *couple_lam[ecouplamNR+1]    = {
    "vdw-q", "vdw", "q", "none", nullptr
};

static void GetSimTemps(int ntemps, t_simtemp *simtemp, double *temperature_lambdas)
{

    int i;

    for (i = 0; i < ntemps; i++)
    {
        /* simple linear scaling -- allows more control */
        if (simtemp->eSimTempScale == esimtempLINEAR)
        {
            simtemp->temperatures[i] = simtemp->simtemp_low + (simtemp->simtemp_high-simtemp->simtemp_low)*temperature_lambdas[i];
        }
        else if (simtemp->eSimTempScale == esimtempGEOMETRIC)  /* should give roughly equal acceptance for constant heat capacity . . . */
        {
            simtemp->temperatures[i] = simtemp->simtemp_low * std::pow(simtemp->simtemp_high/simtemp->simtemp_low, static_cast<real>((1.0*i)/(ntemps-1)));
        }
        else if (simtemp->eSimTempScale == esimtempEXPONENTIAL)
        {
            simtemp->temperatures[i] = simtemp->simtemp_low + (simtemp->simtemp_high-simtemp->simtemp_low)*(std::expm1(temperature_lambdas[i])/std::expm1(1.0));
        }
        else
        {
            char errorstr[128];
            sprintf(errorstr, "eSimTempScale=%d not defined", simtemp->eSimTempScale);
            gmx_fatal(FARGS, "%s", errorstr);
        }
    }
}



static void _low_check(bool b, const char *s, warninp_t wi)
{
    if (b)
    {
        warning_error(wi, s);
    }
}

static void check_nst(const char *desc_nst, int nst,
                      const char *desc_p, int *p,
                      warninp_t wi)
{
    char buf[STRLEN];

    if (*p > 0 && *p % nst != 0)
    {
        /* Round up to the next multiple of nst */
        *p = ((*p)/nst + 1)*nst;
        sprintf(buf, "%s should be a multiple of %s, changing %s to %d\n",
                desc_p, desc_nst, desc_p, *p);
        warning(wi, buf);
    }
}

static bool ir_NVE(const t_inputrec *ir)
{
    return (EI_MD(ir->eI) && ir->etc == etcNO);
}

static int lcd(int n1, int n2)
{
    int d, i;

    d = 1;
    for (i = 2; (i <= n1 && i <= n2); i++)
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }

    return d;
}

static void process_interaction_modifier(const t_inputrec *ir, int *eintmod)
{
    if (*eintmod == eintmodPOTSHIFT_VERLET)
    {
        if (ir->cutoff_scheme == ecutsVERLET)
        {
            *eintmod = eintmodPOTSHIFT;
        }
        else
        {
            *eintmod = eintmodNONE;
        }
    }
}

void check_ir(const char *mdparin, t_inputrec *ir, t_gromppopts *opts,
              warninp_t wi)
/* Check internal consistency.
 * NOTE: index groups are not set here yet, don't check things
 * like temperature coupling group options here, but in triple_check
 */
{
    /* Strange macro: first one fills the err_buf, and then one can check
     * the condition, which will print the message and increase the error
     * counter.
     */
#define CHECK(b) _low_check(b, err_buf, wi)
    char        err_buf[256], warn_buf[STRLEN];
    int         i, j;
    real        dt_pcoupl;
    t_lambda   *fep    = ir->fepvals;
    t_expanded *expand = ir->expandedvals;

    set_warning_line(wi, mdparin, -1);

    if (ir->coulombtype == eelRF_NEC_UNSUPPORTED)
    {
        sprintf(warn_buf, "%s electrostatics is no longer supported",
                eel_names[eelRF_NEC_UNSUPPORTED]);
        warning_error(wi, warn_buf);
    }

    /* BASIC CUT-OFF STUFF */
    if (ir->rcoulomb < 0)
    {
        warning_error(wi, "rcoulomb should be >= 0");
    }
    if (ir->rvdw < 0)
    {
        warning_error(wi, "rvdw should be >= 0");
    }
    if (ir->rlist < 0 &&
        !(ir->cutoff_scheme == ecutsVERLET && ir->verletbuf_tol > 0))
    {
        warning_error(wi, "rlist should be >= 0");
    }
    sprintf(err_buf, "nstlist can not be smaller than 0. (If you were trying to use the heuristic neighbour-list update scheme for efficient buffering for improved energy conservation, please use the Verlet cut-off scheme instead.)");
    CHECK(ir->nstlist < 0);

    process_interaction_modifier(ir, &ir->coulomb_modifier);
    process_interaction_modifier(ir, &ir->vdw_modifier);

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        warning_note(wi,
                     "The group cutoff scheme is deprecated since GROMACS 5.0 and will be removed in a future "
                     "release when all interaction forms are supported for the verlet scheme. The verlet "
                     "scheme already scales better, and it is compatible with GPUs and other accelerators.");

        if (ir->rlist > 0 && ir->rlist < ir->rcoulomb)
        {
            gmx_fatal(FARGS, "rcoulomb must not be greater than rlist (twin-range schemes are not supported)");
        }
        if (ir->rlist > 0 && ir->rlist < ir->rvdw)
        {
            gmx_fatal(FARGS, "rvdw must not be greater than rlist (twin-range schemes are not supported)");
        }

        if (ir->rlist == 0 && ir->ePBC != epbcNONE)
        {
            warning_error(wi, "Can not have an infinite cut-off with PBC");
        }
    }

    if (ir->cutoff_scheme == ecutsVERLET)
    {
        real rc_max;

        /* Normal Verlet type neighbor-list, currently only limited feature support */
        if (inputrec2nboundeddim(ir) < 3)
        {
            warning_error(wi, "With Verlet lists only full pbc or pbc=xy with walls is supported");
        }

        // We don't (yet) have general Verlet kernels for rcoulomb!=rvdw
        if (ir->rcoulomb != ir->rvdw)
        {
            // Since we have PME coulomb + LJ cut-off kernels with rcoulomb>rvdw
            // for PME load balancing, we can support this exception.
            bool bUsesPmeTwinRangeKernel = (EEL_PME_EWALD(ir->coulombtype) &&
                                            ir->vdwtype == evdwCUT &&
                                            ir->rcoulomb > ir->rvdw);
            if (!bUsesPmeTwinRangeKernel)
            {
                warning_error(wi, "With Verlet lists rcoulomb!=rvdw is not supported (except for rcoulomb>rvdw with PME electrostatics)");
            }
        }

        if (ir->vdwtype == evdwSHIFT || ir->vdwtype == evdwSWITCH)
        {
            if (ir->vdw_modifier == eintmodNONE ||
                ir->vdw_modifier == eintmodPOTSHIFT)
            {
                ir->vdw_modifier = (ir->vdwtype == evdwSHIFT ? eintmodFORCESWITCH : eintmodPOTSWITCH);

                sprintf(warn_buf, "Replacing vdwtype=%s by the equivalent combination of vdwtype=%s and vdw_modifier=%s", evdw_names[ir->vdwtype], evdw_names[evdwCUT], eintmod_names[ir->vdw_modifier]);
                warning_note(wi, warn_buf);

                ir->vdwtype = evdwCUT;
            }
            else
            {
                sprintf(warn_buf, "Unsupported combination of vdwtype=%s and vdw_modifier=%s", evdw_names[ir->vdwtype], eintmod_names[ir->vdw_modifier]);
                warning_error(wi, warn_buf);
            }
        }

        if (!(ir->vdwtype == evdwCUT || ir->vdwtype == evdwPME))
        {
            warning_error(wi, "With Verlet lists only cut-off and PME LJ interactions are supported");
        }
        if (!(ir->coulombtype == eelCUT || EEL_RF(ir->coulombtype) ||
              EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD))
        {
            warning_error(wi, "With Verlet lists only cut-off, reaction-field, PME and Ewald electrostatics are supported");
        }
        if (!(ir->coulomb_modifier == eintmodNONE ||
              ir->coulomb_modifier == eintmodPOTSHIFT))
        {
            sprintf(warn_buf, "coulomb_modifier=%s is not supported with the Verlet cut-off scheme", eintmod_names[ir->coulomb_modifier]);
            warning_error(wi, warn_buf);
        }

        if (EEL_USER(ir->coulombtype))
        {
            sprintf(warn_buf, "Coulomb type %s is not supported with the verlet scheme", eel_names[ir->coulombtype]);
            warning_error(wi, warn_buf);
        }

        if (ir->nstlist <= 0)
        {
            warning_error(wi, "With Verlet lists nstlist should be larger than 0");
        }

        if (ir->nstlist < 10)
        {
            warning_note(wi, "With Verlet lists the optimal nstlist is >= 10, with GPUs >= 20. Note that with the Verlet scheme, nstlist has no effect on the accuracy of your simulation.");
        }

        rc_max = std::max(ir->rvdw, ir->rcoulomb);

        if (ir->verletbuf_tol <= 0)
        {
            if (ir->verletbuf_tol == 0)
            {
                warning_error(wi, "Can not have Verlet buffer tolerance of exactly 0");
            }

            if (ir->rlist < rc_max)
            {
                warning_error(wi, "With verlet lists rlist can not be smaller than rvdw or rcoulomb");
            }

            if (ir->rlist == rc_max && ir->nstlist > 1)
            {
                warning_note(wi, "rlist is equal to rvdw and/or rcoulomb: there is no explicit Verlet buffer. The cluster pair list does have a buffering effect, but choosing a larger rlist might be necessary for good energy conservation.");
            }
        }
        else
        {
            if (ir->rlist > rc_max)
            {
                warning_note(wi, "You have set rlist larger than the interaction cut-off, but you also have verlet-buffer-tolerance > 0. Will set rlist using verlet-buffer-tolerance.");
            }

            if (ir->nstlist == 1)
            {
                /* No buffer required */
                ir->rlist = rc_max;
            }
            else
            {
                if (EI_DYNAMICS(ir->eI))
                {
                    if (inputrec2nboundeddim(ir) < 3)
                    {
                        warning_error(wi, "The box volume is required for calculating rlist from the energy drift with verlet-buffer-tolerance > 0. You are using at least one unbounded dimension, so no volume can be computed. Either use a finite box, or set rlist yourself together with verlet-buffer-tolerance = -1.");
                    }
                    /* Set rlist temporarily so we can continue processing */
                    ir->rlist = rc_max;
                }
                else
                {
                    /* Set the buffer to 5% of the cut-off */
                    ir->rlist = (1.0 + verlet_buffer_ratio_nodynamics)*rc_max;
                }
            }
        }
    }

    /* GENERAL INTEGRATOR STUFF */
    if (!EI_MD(ir->eI))
    {
        if (ir->etc != etcNO)
        {
            if (EI_RANDOM(ir->eI))
            {
                sprintf(warn_buf, "Setting tcoupl from '%s' to 'no'. %s handles temperature coupling implicitly. See the documentation for more information on which parameters affect temperature for %s.", etcoupl_names[ir->etc], ei_names[ir->eI], ei_names[ir->eI]);
            }
            else
            {
                sprintf(warn_buf, "Setting tcoupl from '%s' to 'no'. Temperature coupling does not apply to %s.", etcoupl_names[ir->etc], ei_names[ir->eI]);
            }
            warning_note(wi, warn_buf);
        }
        ir->etc = etcNO;
    }
    if (ir->eI == eiVVAK)
    {
        sprintf(warn_buf, "Integrator method %s is implemented primarily for validation purposes; for molecular dynamics, you should probably be using %s or %s", ei_names[eiVVAK], ei_names[eiMD], ei_names[eiVV]);
        warning_note(wi, warn_buf);
    }
    if (!EI_DYNAMICS(ir->eI))
    {
        if (ir->epc != epcNO)
        {
            sprintf(warn_buf, "Setting pcoupl from '%s' to 'no'. Pressure coupling does not apply to %s.", epcoupl_names[ir->epc], ei_names[ir->eI]);
            warning_note(wi, warn_buf);
        }
        ir->epc = epcNO;
    }
    if (EI_DYNAMICS(ir->eI))
    {
        if (ir->nstcalcenergy < 0)
        {
            ir->nstcalcenergy = ir_optimal_nstcalcenergy(ir);
            if (ir->nstenergy != 0 && ir->nstenergy < ir->nstcalcenergy)
            {
                /* nstcalcenergy larger than nstener does not make sense.
                 * We ideally want nstcalcenergy=nstener.
                 */
                if (ir->nstlist > 0)
                {
                    ir->nstcalcenergy = lcd(ir->nstenergy, ir->nstlist);
                }
                else
                {
                    ir->nstcalcenergy = ir->nstenergy;
                }
            }
        }
        else if ( (ir->nstenergy > 0 && ir->nstcalcenergy > ir->nstenergy) ||
                  (ir->efep != efepNO && ir->fepvals->nstdhdl > 0 &&
                   (ir->nstcalcenergy > ir->fepvals->nstdhdl) ) )

        {
            const char *nsten    = "nstenergy";
            const char *nstdh    = "nstdhdl";
            const char *min_name = nsten;
            int         min_nst  = ir->nstenergy;

            /* find the smallest of ( nstenergy, nstdhdl ) */
            if (ir->efep != efepNO && ir->fepvals->nstdhdl > 0 &&
                (ir->nstenergy == 0 || ir->fepvals->nstdhdl < ir->nstenergy))
            {
                min_nst  = ir->fepvals->nstdhdl;
                min_name = nstdh;
            }
            /* If the user sets nstenergy small, we should respect that */
            sprintf(warn_buf,
                    "Setting nstcalcenergy (%d) equal to %s (%d)",
                    ir->nstcalcenergy, min_name, min_nst);
            warning_note(wi, warn_buf);
            ir->nstcalcenergy = min_nst;
        }

        if (ir->epc != epcNO)
        {
            if (ir->nstpcouple < 0)
            {
                ir->nstpcouple = ir_optimal_nstpcouple(ir);
            }
        }

        if (ir->nstcalcenergy > 0)
        {
            if (ir->efep != efepNO)
            {
                /* nstdhdl should be a multiple of nstcalcenergy */
                check_nst("nstcalcenergy", ir->nstcalcenergy,
                          "nstdhdl", &ir->fepvals->nstdhdl, wi);
                /* nstexpanded should be a multiple of nstcalcenergy */
                check_nst("nstcalcenergy", ir->nstcalcenergy,
                          "nstexpanded", &ir->expandedvals->nstexpanded, wi);
            }
            /* for storing exact averages nstenergy should be
             * a multiple of nstcalcenergy
             */
            check_nst("nstcalcenergy", ir->nstcalcenergy,
                      "nstenergy", &ir->nstenergy, wi);
        }
    }

    if (ir->nsteps == 0 && !ir->bContinuation)
    {
        warning_note(wi, "For a correct single-point energy evaluation with nsteps = 0, use continuation = yes to avoid constraining the input coordinates.");
    }

    /* LD STUFF */
    if ((EI_SD(ir->eI) || ir->eI == eiBD) &&
        ir->bContinuation && ir->ld_seed != -1)
    {
        warning_note(wi, "You are doing a continuation with SD or BD, make sure that ld_seed is different from the previous run (using ld_seed=-1 will ensure this)");
    }

    /* TPI STUFF */
    if (EI_TPI(ir->eI))
    {
        sprintf(err_buf, "TPI only works with pbc = %s", epbc_names[epbcXYZ]);
        CHECK(ir->ePBC != epbcXYZ);
        sprintf(err_buf, "TPI only works with ns = %s", ens_names[ensGRID]);
        CHECK(ir->ns_type != ensGRID);
        sprintf(err_buf, "with TPI nstlist should be larger than zero");
        CHECK(ir->nstlist <= 0);
        sprintf(err_buf, "TPI does not work with full electrostatics other than PME");
        CHECK(EEL_FULL(ir->coulombtype) && !EEL_PME(ir->coulombtype));
        sprintf(err_buf, "TPI does not work (yet) with the Verlet cut-off scheme");
        CHECK(ir->cutoff_scheme == ecutsVERLET);
    }

    /* SHAKE / LINCS */
    if ( (opts->nshake > 0) && (opts->bMorse) )
    {
        sprintf(warn_buf,
                "Using morse bond-potentials while constraining bonds is useless");
        warning(wi, warn_buf);
    }

    if ((EI_SD(ir->eI) || ir->eI == eiBD) &&
        ir->bContinuation && ir->ld_seed != -1)
    {
        warning_note(wi, "You are doing a continuation with SD or BD, make sure that ld_seed is different from the previous run (using ld_seed=-1 will ensure this)");
    }
    /* verify simulated tempering options */

    if (ir->bSimTemp)
    {
        bool bAllTempZero = TRUE;
        for (i = 0; i < fep->n_lambda; i++)
        {
            sprintf(err_buf, "Entry %d for %s must be between 0 and 1, instead is %g", i, efpt_names[efptTEMPERATURE], fep->all_lambda[efptTEMPERATURE][i]);
            CHECK((fep->all_lambda[efptTEMPERATURE][i] < 0) || (fep->all_lambda[efptTEMPERATURE][i] > 1));
            if (fep->all_lambda[efptTEMPERATURE][i] > 0)
            {
                bAllTempZero = FALSE;
            }
        }
        sprintf(err_buf, "if simulated tempering is on, temperature-lambdas may not be all zero");
        CHECK(bAllTempZero == TRUE);

        sprintf(err_buf, "Simulated tempering is currently only compatible with md-vv");
        CHECK(ir->eI != eiVV);

        /* check compatability of the temperature coupling with simulated tempering */

        if (ir->etc == etcNOSEHOOVER)
        {
            sprintf(warn_buf, "Nose-Hoover based temperature control such as [%s] my not be entirelyconsistent with simulated tempering", etcoupl_names[ir->etc]);
            warning_note(wi, warn_buf);
        }

        /* check that the temperatures make sense */

        sprintf(err_buf, "Higher simulated tempering temperature (%g) must be >= than the simulated tempering lower temperature (%g)", ir->simtempvals->simtemp_high, ir->simtempvals->simtemp_low);
        CHECK(ir->simtempvals->simtemp_high <= ir->simtempvals->simtemp_low);

        sprintf(err_buf, "Higher simulated tempering temperature (%g) must be >= zero", ir->simtempvals->simtemp_high);
        CHECK(ir->simtempvals->simtemp_high <= 0);

        sprintf(err_buf, "Lower simulated tempering temperature (%g) must be >= zero", ir->simtempvals->simtemp_low);
        CHECK(ir->simtempvals->simtemp_low <= 0);
    }

    /* verify free energy options */

    if (ir->efep != efepNO)
    {
        fep = ir->fepvals;
        sprintf(err_buf, "The soft-core power is %d and can only be 1 or 2",
                fep->sc_power);
        CHECK(fep->sc_alpha != 0 && fep->sc_power != 1 && fep->sc_power != 2);

        sprintf(err_buf, "The soft-core sc-r-power is %d and can only be 6 or 48",
                static_cast<int>(fep->sc_r_power));
        CHECK(fep->sc_alpha != 0 && fep->sc_r_power != 6.0 && fep->sc_r_power != 48.0);

        sprintf(err_buf, "Can't use positive delta-lambda (%g) if initial state/lambda does not start at zero", fep->delta_lambda);
        CHECK(fep->delta_lambda > 0 && ((fep->init_fep_state > 0) ||  (fep->init_lambda > 0)));

        sprintf(err_buf, "Can't use positive delta-lambda (%g) with expanded ensemble simulations", fep->delta_lambda);
        CHECK(fep->delta_lambda > 0 && (ir->efep == efepEXPANDED));

        sprintf(err_buf, "Can only use expanded ensemble with md-vv (for now)");
        CHECK(!(EI_VV(ir->eI)) && (ir->efep == efepEXPANDED));

        sprintf(err_buf, "Free-energy not implemented for Ewald");
        CHECK(ir->coulombtype == eelEWALD);

        /* check validty of lambda inputs */
        if (fep->n_lambda == 0)
        {
            /* Clear output in case of no states:*/
            sprintf(err_buf, "init-lambda-state set to %d: no lambda states are defined.", fep->init_fep_state);
            CHECK((fep->init_fep_state >= 0) && (fep->n_lambda == 0));
        }
        else
        {
            sprintf(err_buf, "initial thermodynamic state %d does not exist, only goes to %d", fep->init_fep_state, fep->n_lambda-1);
            CHECK((fep->init_fep_state >= fep->n_lambda));
        }

        sprintf(err_buf, "Lambda state must be set, either with init-lambda-state or with init-lambda");
        CHECK((fep->init_fep_state < 0) && (fep->init_lambda < 0));

        sprintf(err_buf, "init-lambda=%g while init-lambda-state=%d. Lambda state must be set either with init-lambda-state or with init-lambda, but not both",
                fep->init_lambda, fep->init_fep_state);
        CHECK((fep->init_fep_state >= 0) && (fep->init_lambda >= 0));



        if ((fep->init_lambda >= 0) && (fep->delta_lambda == 0))
        {
            int n_lambda_terms;
            n_lambda_terms = 0;
            for (i = 0; i < efptNR; i++)
            {
                if (fep->separate_dvdl[i])
                {
                    n_lambda_terms++;
                }
            }
            if (n_lambda_terms > 1)
            {
                sprintf(warn_buf, "If lambda vector states (fep-lambdas, coul-lambdas etc.) are set, don't use init-lambda to set lambda state (except for slow growth). Use init-lambda-state instead.");
                warning(wi, warn_buf);
            }

            if (n_lambda_terms < 2 && fep->n_lambda > 0)
            {
                warning_note(wi,
                             "init-lambda is deprecated for setting lambda state (except for slow growth). Use init-lambda-state instead.");
            }
        }

        for (j = 0; j < efptNR; j++)
        {
            for (i = 0; i < fep->n_lambda; i++)
            {
                sprintf(err_buf, "Entry %d for %s must be between 0 and 1, instead is %g", i, efpt_names[j], fep->all_lambda[j][i]);
                CHECK((fep->all_lambda[j][i] < 0) || (fep->all_lambda[j][i] > 1));
            }
        }

        if ((fep->sc_alpha > 0) && (!fep->bScCoul))
        {
            for (i = 0; i < fep->n_lambda; i++)
            {
                sprintf(err_buf, "For state %d, vdw-lambdas (%f) is changing with vdw softcore, while coul-lambdas (%f) is nonzero without coulomb softcore: this will lead to crashes, and is not supported.", i, fep->all_lambda[efptVDW][i],
                        fep->all_lambda[efptCOUL][i]);
                CHECK((fep->sc_alpha > 0) &&
                      (((fep->all_lambda[efptCOUL][i] > 0.0) &&
                        (fep->all_lambda[efptCOUL][i] < 1.0)) &&
                       ((fep->all_lambda[efptVDW][i] > 0.0) &&
                        (fep->all_lambda[efptVDW][i] < 1.0))));
            }
        }

        if ((fep->bScCoul) && (EEL_PME(ir->coulombtype)))
        {
            real sigma, lambda, r_sc;

            sigma  = 0.34;
            /* Maximum estimate for A and B charges equal with lambda power 1 */
            lambda = 0.5;
            r_sc   = std::pow(lambda*fep->sc_alpha*std::pow(sigma/ir->rcoulomb, fep->sc_r_power) + 1.0, 1.0/fep->sc_r_power);
            sprintf(warn_buf, "With PME there is a minor soft core effect present at the cut-off, proportional to (LJsigma/rcoulomb)^%g. This could have a minor effect on energy conservation, but usually other effects dominate. With a common sigma value of %g nm the fraction of the particle-particle potential at the cut-off at lambda=%g is around %.1e, while ewald-rtol is %.1e.",
                    fep->sc_r_power,
                    sigma, lambda, r_sc - 1.0, ir->ewald_rtol);
            warning_note(wi, warn_buf);
        }

        /*  Free Energy Checks -- In an ideal world, slow growth and FEP would
            be treated differently, but that's the next step */

        for (i = 0; i < efptNR; i++)
        {
            for (j = 0; j < fep->n_lambda; j++)
            {
                sprintf(err_buf, "%s[%d] must be between 0 and 1", efpt_names[i], j);
                CHECK((fep->all_lambda[i][j] < 0) || (fep->all_lambda[i][j] > 1));
            }
        }
    }

    if ((ir->bSimTemp) || (ir->efep == efepEXPANDED))
    {
        fep    = ir->fepvals;

        /* checking equilibration of weights inputs for validity */

        sprintf(err_buf, "weight-equil-number-all-lambda (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_n_at_lam, elmceq_names[elmceqNUMATLAM]);
        CHECK((expand->equil_n_at_lam > 0) && (expand->elmceq != elmceqNUMATLAM));

        sprintf(err_buf, "weight-equil-number-samples (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_samples, elmceq_names[elmceqSAMPLES]);
        CHECK((expand->equil_samples > 0) && (expand->elmceq != elmceqSAMPLES));

        sprintf(err_buf, "weight-equil-number-steps (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_steps, elmceq_names[elmceqSTEPS]);
        CHECK((expand->equil_steps > 0) && (expand->elmceq != elmceqSTEPS));

        sprintf(err_buf, "weight-equil-wl-delta (%d) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_samples, elmceq_names[elmceqWLDELTA]);
        CHECK((expand->equil_wl_delta > 0) && (expand->elmceq != elmceqWLDELTA));

        sprintf(err_buf, "weight-equil-count-ratio (%f) is ignored if lmc-weights-equil is not equal to %s",
                expand->equil_ratio, elmceq_names[elmceqRATIO]);
        CHECK((expand->equil_ratio > 0) && (expand->elmceq != elmceqRATIO));

        sprintf(err_buf, "weight-equil-number-all-lambda (%d) must be a positive integer if lmc-weights-equil=%s",
                expand->equil_n_at_lam, elmceq_names[elmceqNUMATLAM]);
        CHECK((expand->equil_n_at_lam <= 0) && (expand->elmceq == elmceqNUMATLAM));

        sprintf(err_buf, "weight-equil-number-samples (%d) must be a positive integer if lmc-weights-equil=%s",
                expand->equil_samples, elmceq_names[elmceqSAMPLES]);
        CHECK((expand->equil_samples <= 0) && (expand->elmceq == elmceqSAMPLES));

        sprintf(err_buf, "weight-equil-number-steps (%d) must be a positive integer if lmc-weights-equil=%s",
                expand->equil_steps, elmceq_names[elmceqSTEPS]);
        CHECK((expand->equil_steps <= 0) && (expand->elmceq == elmceqSTEPS));

        sprintf(err_buf, "weight-equil-wl-delta (%f) must be > 0 if lmc-weights-equil=%s",
                expand->equil_wl_delta, elmceq_names[elmceqWLDELTA]);
        CHECK((expand->equil_wl_delta <= 0) && (expand->elmceq == elmceqWLDELTA));

        sprintf(err_buf, "weight-equil-count-ratio (%f) must be > 0 if lmc-weights-equil=%s",
                expand->equil_ratio, elmceq_names[elmceqRATIO]);
        CHECK((expand->equil_ratio <= 0) && (expand->elmceq == elmceqRATIO));

        sprintf(err_buf, "lmc-weights-equil=%s only possible when lmc-stats = %s or lmc-stats %s",
                elmceq_names[elmceqWLDELTA], elamstats_names[elamstatsWL], elamstats_names[elamstatsWWL]);
        CHECK((expand->elmceq == elmceqWLDELTA) && (!EWL(expand->elamstats)));

        sprintf(err_buf, "lmc-repeats (%d) must be greater than 0", expand->lmc_repeats);
        CHECK((expand->lmc_repeats <= 0));
        sprintf(err_buf, "minimum-var-min (%d) must be greater than 0", expand->minvarmin);
        CHECK((expand->minvarmin <= 0));
        sprintf(err_buf, "weight-c-range (%d) must be greater or equal to 0", expand->c_range);
        CHECK((expand->c_range < 0));
        sprintf(err_buf, "init-lambda-state (%d) must be zero if lmc-forced-nstart (%d)> 0 and lmc-move != 'no'",
                fep->init_fep_state, expand->lmc_forced_nstart);
        CHECK((fep->init_fep_state != 0) && (expand->lmc_forced_nstart > 0) && (expand->elmcmove != elmcmoveNO));
        sprintf(err_buf, "lmc-forced-nstart (%d) must not be negative", expand->lmc_forced_nstart);
        CHECK((expand->lmc_forced_nstart < 0));
        sprintf(err_buf, "init-lambda-state (%d) must be in the interval [0,number of lambdas)", fep->init_fep_state);
        CHECK((fep->init_fep_state < 0) || (fep->init_fep_state >= fep->n_lambda));

        sprintf(err_buf, "init-wl-delta (%f) must be greater than or equal to 0", expand->init_wl_delta);
        CHECK((expand->init_wl_delta < 0));
        sprintf(err_buf, "wl-ratio (%f) must be between 0 and 1", expand->wl_ratio);
        CHECK((expand->wl_ratio <= 0) || (expand->wl_ratio >= 1));
        sprintf(err_buf, "wl-scale (%f) must be between 0 and 1", expand->wl_scale);
        CHECK((expand->wl_scale <= 0) || (expand->wl_scale >= 1));

        /* if there is no temperature control, we need to specify an MC temperature */
        if (!integratorHasReferenceTemperature(ir) && (expand->elmcmove != elmcmoveNO) && (expand->mc_temp <= 0.0))
        {
            sprintf(err_buf, "If there is no temperature control, and lmc-mcmove!='no', mc_temp must be set to a positive number");
            warning_error(wi, err_buf);
        }
        if (expand->nstTij > 0)
        {
            sprintf(err_buf, "nstlog must be non-zero");
            CHECK(ir->nstlog == 0);
            sprintf(err_buf, "nst-transition-matrix (%d) must be an integer multiple of nstlog (%d)",
                    expand->nstTij, ir->nstlog);
            CHECK((expand->nstTij % ir->nstlog) != 0);
        }
    }

    /* PBC/WALLS */
    sprintf(err_buf, "walls only work with pbc=%s", epbc_names[epbcXY]);
    CHECK(ir->nwall && ir->ePBC != epbcXY);

    /* VACUUM STUFF */
    if (ir->ePBC != epbcXYZ && ir->nwall != 2)
    {
        if (ir->ePBC == epbcNONE)
        {
            if (ir->epc != epcNO)
            {
                warning(wi, "Turning off pressure coupling for vacuum system");
                ir->epc = epcNO;
            }
        }
        else
        {
            sprintf(err_buf, "Can not have pressure coupling with pbc=%s",
                    epbc_names[ir->ePBC]);
            CHECK(ir->epc != epcNO);
        }
        sprintf(err_buf, "Can not have Ewald with pbc=%s", epbc_names[ir->ePBC]);
        CHECK(EEL_FULL(ir->coulombtype));

        sprintf(err_buf, "Can not have dispersion correction with pbc=%s",
                epbc_names[ir->ePBC]);
        CHECK(ir->eDispCorr != edispcNO);
    }

    if (ir->rlist == 0.0)
    {
        sprintf(err_buf, "can only have neighborlist cut-off zero (=infinite)\n"
                "with coulombtype = %s or coulombtype = %s\n"
                "without periodic boundary conditions (pbc = %s) and\n"
                "rcoulomb and rvdw set to zero",
                eel_names[eelCUT], eel_names[eelUSER], epbc_names[epbcNONE]);
        CHECK(((ir->coulombtype != eelCUT) && (ir->coulombtype != eelUSER)) ||
              (ir->ePBC     != epbcNONE) ||
              (ir->rcoulomb != 0.0)      || (ir->rvdw != 0.0));

        if (ir->nstlist > 0)
        {
            warning_note(wi, "Simulating without cut-offs can be (slightly) faster with nstlist=0, nstype=simple and only one MPI rank");
        }
    }

    /* COMM STUFF */
    if (ir->nstcomm == 0)
    {
        ir->comm_mode = ecmNO;
    }
    if (ir->comm_mode != ecmNO)
    {
        if (ir->nstcomm < 0)
        {
            warning(wi, "If you want to remove the rotation around the center of mass, you should set comm_mode = Angular instead of setting nstcomm < 0. nstcomm is modified to its absolute value");
            ir->nstcomm = abs(ir->nstcomm);
        }

        if (ir->nstcalcenergy > 0 && ir->nstcomm < ir->nstcalcenergy)
        {
            warning_note(wi, "nstcomm < nstcalcenergy defeats the purpose of nstcalcenergy, setting nstcomm to nstcalcenergy");
            ir->nstcomm = ir->nstcalcenergy;
        }

        if (ir->comm_mode == ecmANGULAR)
        {
            sprintf(err_buf, "Can not remove the rotation around the center of mass with periodic molecules");
            CHECK(ir->bPeriodicMols);
            if (ir->ePBC != epbcNONE)
            {
                warning(wi, "Removing the rotation around the center of mass in a periodic system, this can lead to artifacts. Only use this on a single (cluster of) molecules. This cluster should not cross periodic boundaries.");
            }
        }
    }

    if (EI_STATE_VELOCITY(ir->eI) && !EI_SD(ir->eI) && ir->ePBC == epbcNONE && ir->comm_mode != ecmANGULAR)
    {
        sprintf(warn_buf, "Tumbling and flying ice-cubes: We are not removing rotation around center of mass in a non-periodic system. You should probably set comm_mode = ANGULAR or use integrator = %s.", ei_names[eiSD1]);
        warning_note(wi, warn_buf);
    }

    /* TEMPERATURE COUPLING */
    if (ir->etc == etcYES)
    {
        ir->etc = etcBERENDSEN;
        warning_note(wi, "Old option for temperature coupling given: "
                     "changing \"yes\" to \"Berendsen\"\n");
    }

    if ((ir->etc == etcNOSEHOOVER) || (ir->epc == epcMTTK))
    {
        if (ir->opts.nhchainlength < 1)
        {
            sprintf(warn_buf, "number of Nose-Hoover chains (currently %d) cannot be less than 1,reset to 1\n", ir->opts.nhchainlength);
            ir->opts.nhchainlength = 1;
            warning(wi, warn_buf);
        }

        if (ir->etc == etcNOSEHOOVER && !EI_VV(ir->eI) && ir->opts.nhchainlength > 1)
        {
            warning_note(wi, "leapfrog does not yet support Nose-Hoover chains, nhchainlength reset to 1");
            ir->opts.nhchainlength = 1;
        }
    }
    else
    {
        ir->opts.nhchainlength = 0;
    }

    if (ir->eI == eiVVAK)
    {
        sprintf(err_buf, "%s implemented primarily for validation, and requires nsttcouple = 1 and nstpcouple = 1.",
                ei_names[eiVVAK]);
        CHECK((ir->nsttcouple != 1) || (ir->nstpcouple != 1));
    }

    if (ETC_ANDERSEN(ir->etc))
    {
        sprintf(err_buf, "%s temperature control not supported for integrator %s.", etcoupl_names[ir->etc], ei_names[ir->eI]);
        CHECK(!(EI_VV(ir->eI)));

        if (ir->nstcomm > 0 && (ir->etc == etcANDERSEN))
        {
            sprintf(warn_buf, "Center of mass removal not necessary for %s.  All velocities of coupled groups are rerandomized periodically, so flying ice cube errors will not occur.", etcoupl_names[ir->etc]);
            warning_note(wi, warn_buf);
        }

        sprintf(err_buf, "nstcomm must be 1, not %d for %s, as velocities of atoms in coupled groups are randomized every time step", ir->nstcomm, etcoupl_names[ir->etc]);
        CHECK(ir->nstcomm > 1 && (ir->etc == etcANDERSEN));
    }

    if (ir->etc == etcBERENDSEN)
    {
        sprintf(warn_buf, "The %s thermostat does not generate the correct kinetic energy distribution. You might want to consider using the %s thermostat.",
                ETCOUPLTYPE(ir->etc), ETCOUPLTYPE(etcVRESCALE));
        warning_note(wi, warn_buf);
    }

    if ((ir->etc == etcNOSEHOOVER || ETC_ANDERSEN(ir->etc))
        && ir->epc == epcBERENDSEN)
    {
        sprintf(warn_buf, "Using Berendsen pressure coupling invalidates the "
                "true ensemble for the thermostat");
        warning(wi, warn_buf);
    }

    /* PRESSURE COUPLING */
    if (ir->epc == epcISOTROPIC)
    {
        ir->epc = epcBERENDSEN;
        warning_note(wi, "Old option for pressure coupling given: "
                     "changing \"Isotropic\" to \"Berendsen\"\n");
    }

    if (ir->epc != epcNO)
    {
        dt_pcoupl = ir->nstpcouple*ir->delta_t;

        sprintf(err_buf, "tau-p must be > 0 instead of %g\n", ir->tau_p);
        CHECK(ir->tau_p <= 0);

        if (ir->tau_p/dt_pcoupl < pcouple_min_integration_steps(ir->epc) - 10*GMX_REAL_EPS)
        {
            sprintf(warn_buf, "For proper integration of the %s barostat, tau-p (%g) should be at least %d times larger than nstpcouple*dt (%g)",
                    EPCOUPLTYPE(ir->epc), ir->tau_p, pcouple_min_integration_steps(ir->epc), dt_pcoupl);
            warning(wi, warn_buf);
        }

        sprintf(err_buf, "compressibility must be > 0 when using pressure"
                " coupling %s\n", EPCOUPLTYPE(ir->epc));
        CHECK(ir->compress[XX][XX] < 0 || ir->compress[YY][YY] < 0 ||
              ir->compress[ZZ][ZZ] < 0 ||
              (trace(ir->compress) == 0 && ir->compress[YY][XX] <= 0 &&
               ir->compress[ZZ][XX] <= 0 && ir->compress[ZZ][YY] <= 0));

        if (epcPARRINELLORAHMAN == ir->epc && opts->bGenVel)
        {
            sprintf(warn_buf,
                    "You are generating velocities so I am assuming you "
                    "are equilibrating a system. You are using "
                    "%s pressure coupling, but this can be "
                    "unstable for equilibration. If your system crashes, try "
                    "equilibrating first with Berendsen pressure coupling. If "
                    "you are not equilibrating the system, you can probably "
                    "ignore this warning.",
                    epcoupl_names[ir->epc]);
            warning(wi, warn_buf);
        }
    }

    if (EI_VV(ir->eI))
    {
        if (ir->epc > epcNO)
        {
            if ((ir->epc != epcBERENDSEN) && (ir->epc != epcMTTK))
            {
                warning_error(wi, "for md-vv and md-vv-avek, can only use Berendsen and Martyna-Tuckerman-Tobias-Klein (MTTK) equations for pressure control; MTTK is equivalent to Parrinello-Rahman.");
            }
        }
    }
    else
    {
        if (ir->epc == epcMTTK)
        {
            warning_error(wi, "MTTK pressure coupling requires a Velocity-verlet integrator");
        }
    }

    /* ELECTROSTATICS */
    /* More checks are in triple check (grompp.c) */

    if (ir->coulombtype == eelSWITCH)
    {
        sprintf(warn_buf, "coulombtype = %s is only for testing purposes and can lead to serious "
                "artifacts, advice: use coulombtype = %s",
                eel_names[ir->coulombtype],
                eel_names[eelRF_ZERO]);
        warning(wi, warn_buf);
    }

    if (EEL_RF(ir->coulombtype) && ir->epsilon_rf == 1 && ir->epsilon_r != 1)
    {
        sprintf(warn_buf, "epsilon-r = %g and epsilon-rf = 1 with reaction field, proceeding assuming old format and exchanging epsilon-r and epsilon-rf", ir->epsilon_r);
        warning(wi, warn_buf);
        ir->epsilon_rf = ir->epsilon_r;
        ir->epsilon_r  = 1.0;
    }

    if (ir->epsilon_r == 0)
    {
        sprintf(err_buf,
                "It is pointless to use long-range electrostatics with infinite relative permittivity."
                "Since you are effectively turning of electrostatics, a plain cutoff will be much faster.");
        CHECK(EEL_FULL(ir->coulombtype));
    }

    if (getenv("GMX_DO_GALACTIC_DYNAMICS") == nullptr)
    {
        sprintf(err_buf, "epsilon-r must be >= 0 instead of %g\n", ir->epsilon_r);
        CHECK(ir->epsilon_r < 0);
    }

    if (EEL_RF(ir->coulombtype))
    {
        /* reaction field (at the cut-off) */

        if (ir->coulombtype == eelRF_ZERO && ir->epsilon_rf != 0)
        {
            sprintf(warn_buf, "With coulombtype = %s, epsilon-rf must be 0, assuming you meant epsilon_rf=0",
                    eel_names[ir->coulombtype]);
            warning(wi, warn_buf);
            ir->epsilon_rf = 0.0;
        }

        sprintf(err_buf, "epsilon-rf must be >= epsilon-r");
        CHECK((ir->epsilon_rf < ir->epsilon_r && ir->epsilon_rf != 0) ||
              (ir->epsilon_r == 0));
        if (ir->epsilon_rf == ir->epsilon_r)
        {
            sprintf(warn_buf, "Using epsilon-rf = epsilon-r with %s does not make sense",
                    eel_names[ir->coulombtype]);
            warning(wi, warn_buf);
        }
    }
    /* Allow rlist>rcoulomb for tabulated long range stuff. This just
     * means the interaction is zero outside rcoulomb, but it helps to
     * provide accurate energy conservation.
     */
    if (ir_coulomb_might_be_zero_at_cutoff(ir))
    {
        if (ir_coulomb_switched(ir))
        {
            sprintf(err_buf,
                    "With coulombtype = %s rcoulomb_switch must be < rcoulomb. Or, better: Use the potential modifier options!",
                    eel_names[ir->coulombtype]);
            CHECK(ir->rcoulomb_switch >= ir->rcoulomb);
        }
    }
    else if (ir->coulombtype == eelCUT || EEL_RF(ir->coulombtype))
    {
        if (ir->cutoff_scheme == ecutsGROUP && ir->coulomb_modifier == eintmodNONE)
        {
            sprintf(err_buf, "With coulombtype = %s, rcoulomb should be >= rlist unless you use a potential modifier",
                    eel_names[ir->coulombtype]);
            CHECK(ir->rlist > ir->rcoulomb);
        }
    }

    if (ir->coulombtype == eelSWITCH || ir->coulombtype == eelSHIFT)
    {
        sprintf(err_buf,
                "Explicit switch/shift coulomb interactions cannot be used in combination with a secondary coulomb-modifier.");
        CHECK( ir->coulomb_modifier != eintmodNONE);
    }
    if (ir->vdwtype == evdwSWITCH || ir->vdwtype == evdwSHIFT)
    {
        sprintf(err_buf,
                "Explicit switch/shift vdw interactions cannot be used in combination with a secondary vdw-modifier.");
        CHECK( ir->vdw_modifier != eintmodNONE);
    }

    if (ir->coulombtype == eelSWITCH || ir->coulombtype == eelSHIFT ||
        ir->vdwtype == evdwSWITCH || ir->vdwtype == evdwSHIFT)
    {
        sprintf(warn_buf,
                "The switch/shift interaction settings are just for compatibility; you will get better "
                "performance from applying potential modifiers to your interactions!\n");
        warning_note(wi, warn_buf);
    }

    if (ir->coulombtype == eelPMESWITCH || ir->coulomb_modifier == eintmodPOTSWITCH)
    {
        if (ir->rcoulomb_switch/ir->rcoulomb < 0.9499)
        {
            real percentage  = 100*(ir->rcoulomb-ir->rcoulomb_switch)/ir->rcoulomb;
            sprintf(warn_buf, "The switching range should be 5%% or less (currently %.2f%% using a switching range of %4f-%4f) for accurate electrostatic energies, energy conservation will be good regardless, since ewald_rtol = %g.",
                    percentage, ir->rcoulomb_switch, ir->rcoulomb, ir->ewald_rtol);
            warning(wi, warn_buf);
        }
    }

    if (ir->vdwtype == evdwSWITCH || ir->vdw_modifier == eintmodPOTSWITCH)
    {
        if (ir->rvdw_switch == 0)
        {
            sprintf(warn_buf, "rvdw-switch is equal 0 even though you are using a switched Lennard-Jones potential.  This suggests it was not set in the mdp, which can lead to large energy errors.  In GROMACS, 0.05 to 0.1 nm is often a reasonable vdw switching range.");
            warning(wi, warn_buf);
        }
    }

    if (EEL_FULL(ir->coulombtype))
    {
        if (ir->coulombtype == eelPMESWITCH || ir->coulombtype == eelPMEUSER ||
            ir->coulombtype == eelPMEUSERSWITCH)
        {
            sprintf(err_buf, "With coulombtype = %s, rcoulomb must be <= rlist",
                    eel_names[ir->coulombtype]);
            CHECK(ir->rcoulomb > ir->rlist);
        }
        else if (ir->cutoff_scheme == ecutsGROUP && ir->coulomb_modifier == eintmodNONE)
        {
            if (ir->coulombtype == eelPME || ir->coulombtype == eelP3M_AD)
            {
                sprintf(err_buf,
                        "With coulombtype = %s (without modifier), rcoulomb must be equal to rlist.\n"
                        "For optimal energy conservation,consider using\n"
                        "a potential modifier.", eel_names[ir->coulombtype]);
                CHECK(ir->rcoulomb != ir->rlist);
            }
        }
    }

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        // TODO: Move these checks into the ewald module with the options class
        int orderMin = 3;
        int orderMax = (ir->coulombtype == eelP3M_AD ? 8 : 12);

        if (ir->pme_order < orderMin || ir->pme_order > orderMax)
        {
            sprintf(warn_buf, "With coulombtype = %s, you should have %d <= pme-order <= %d", eel_names[ir->coulombtype], orderMin, orderMax);
            warning_error(wi, warn_buf);
        }
    }

    if (ir->nwall == 2 && EEL_FULL(ir->coulombtype))
    {
        if (ir->ewald_geometry == eewg3D)
        {
            sprintf(warn_buf, "With pbc=%s you should use ewald-geometry=%s",
                    epbc_names[ir->ePBC], eewg_names[eewg3DC]);
            warning(wi, warn_buf);
        }
        /* This check avoids extra pbc coding for exclusion corrections */
        sprintf(err_buf, "wall-ewald-zfac should be >= 2");
        CHECK(ir->wall_ewald_zfac < 2);
    }
    if ((ir->ewald_geometry == eewg3DC) && (ir->ePBC != epbcXY) &&
        EEL_FULL(ir->coulombtype))
    {
        sprintf(warn_buf, "With %s and ewald_geometry = %s you should use pbc = %s",
                eel_names[ir->coulombtype], eewg_names[eewg3DC], epbc_names[epbcXY]);
        warning(wi, warn_buf);
    }
    if ((ir->epsilon_surface != 0) && EEL_FULL(ir->coulombtype))
    {
        if (ir->cutoff_scheme == ecutsVERLET)
        {
            sprintf(warn_buf, "Since molecules/charge groups are broken using the Verlet scheme, you can not use a dipole correction to the %s electrostatics.",
                    eel_names[ir->coulombtype]);
            warning(wi, warn_buf);
        }
        else
        {
            sprintf(warn_buf, "Dipole corrections to %s electrostatics only work if all charge groups that can cross PBC boundaries are dipoles. If this is not the case set epsilon_surface to 0",
                    eel_names[ir->coulombtype]);
            warning_note(wi, warn_buf);
        }
    }

    if (ir_vdw_switched(ir))
    {
        sprintf(err_buf, "With switched vdw forces or potentials, rvdw-switch must be < rvdw");
        CHECK(ir->rvdw_switch >= ir->rvdw);

        if (ir->rvdw_switch < 0.5*ir->rvdw)
        {
            sprintf(warn_buf, "You are applying a switch function to vdw forces or potentials from %g to %g nm, which is more than half the interaction range, whereas switch functions are intended to act only close to the cut-off.",
                    ir->rvdw_switch, ir->rvdw);
            warning_note(wi, warn_buf);
        }
    }
    else if (ir->vdwtype == evdwCUT || ir->vdwtype == evdwPME)
    {
        if (ir->cutoff_scheme == ecutsGROUP && ir->vdw_modifier == eintmodNONE)
        {
            sprintf(err_buf, "With vdwtype = %s, rvdw must be >= rlist unless you use a potential modifier", evdw_names[ir->vdwtype]);
            CHECK(ir->rlist > ir->rvdw);
        }
    }

    if (ir->vdwtype == evdwPME)
    {
        if (!(ir->vdw_modifier == eintmodNONE || ir->vdw_modifier == eintmodPOTSHIFT))
        {
            sprintf(err_buf, "With vdwtype = %s, the only supported modifiers are %s and %s",
                    evdw_names[ir->vdwtype],
                    eintmod_names[eintmodPOTSHIFT],
                    eintmod_names[eintmodNONE]);
            warning_error(wi, err_buf);
        }
    }

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        if (((ir->coulomb_modifier != eintmodNONE && ir->rcoulomb == ir->rlist) ||
             (ir->vdw_modifier != eintmodNONE && ir->rvdw == ir->rlist)))
        {
            warning_note(wi, "With exact cut-offs, rlist should be "
                         "larger than rcoulomb and rvdw, so that there "
                         "is a buffer region for particle motion "
                         "between neighborsearch steps");
        }

        if (ir_coulomb_is_zero_at_cutoff(ir) && ir->rlist <= ir->rcoulomb)
        {
            sprintf(warn_buf, "For energy conservation with switch/shift potentials, rlist should be 0.1 to 0.3 nm larger than rcoulomb.");
            warning_note(wi, warn_buf);
        }
        if (ir_vdw_switched(ir) && (ir->rlist <= ir->rvdw))
        {
            sprintf(warn_buf, "For energy conservation with switch/shift potentials, rlist should be 0.1 to 0.3 nm larger than rvdw.");
            warning_note(wi, warn_buf);
        }
    }

    if (ir->vdwtype == evdwUSER && ir->eDispCorr != edispcNO)
    {
        warning_note(wi, "You have selected user tables with dispersion correction, the dispersion will be corrected to -C6/r^6 beyond rvdw_switch (the tabulated interaction between rvdw_switch and rvdw will not be double counted). Make sure that you really want dispersion correction to -C6/r^6.");
    }

    if (ir->eI == eiLBFGS && (ir->coulombtype == eelCUT || ir->vdwtype == evdwCUT)
        && ir->rvdw != 0)
    {
        warning(wi, "For efficient BFGS minimization, use switch/shift/pme instead of cut-off.");
    }

    if (ir->eI == eiLBFGS && ir->nbfgscorr <= 0)
    {
        warning(wi, "Using L-BFGS with nbfgscorr<=0 just gets you steepest descent.");
    }

    /* ENERGY CONSERVATION */
    if (ir_NVE(ir) && ir->cutoff_scheme == ecutsGROUP)
    {
        if (!ir_vdw_might_be_zero_at_cutoff(ir) && ir->rvdw > 0 && ir->vdw_modifier == eintmodNONE)
        {
            sprintf(warn_buf, "You are using a cut-off for VdW interactions with NVE, for good energy conservation use vdwtype = %s (possibly with DispCorr)",
                    evdw_names[evdwSHIFT]);
            warning_note(wi, warn_buf);
        }
        if (!ir_coulomb_might_be_zero_at_cutoff(ir) && ir->rcoulomb > 0)
        {
            sprintf(warn_buf, "You are using a cut-off for electrostatics with NVE, for good energy conservation use coulombtype = %s or %s",
                    eel_names[eelPMESWITCH], eel_names[eelRF_ZERO]);
            warning_note(wi, warn_buf);
        }
    }

    /* IMPLICIT SOLVENT */
    if (ir->coulombtype == eelGB_NOTUSED)
    {
        sprintf(warn_buf, "Invalid option %s for coulombtype",
                eel_names[ir->coulombtype]);
        warning_error(wi, warn_buf);
    }

    if (ir->bQMMM)
    {
        if (ir->cutoff_scheme != ecutsGROUP)
        {
            warning_error(wi, "QMMM is currently only supported with cutoff-scheme=group");
        }
        if (!EI_DYNAMICS(ir->eI))
        {
            char buf[STRLEN];
            sprintf(buf, "QMMM is only supported with dynamics, not with integrator %s", ei_names[ir->eI]);
            warning_error(wi, buf);
        }
    }

    if (ir->bAdress)
    {
        gmx_fatal(FARGS, "AdResS simulations are no longer supported");
    }
}

/* interpret a number of doubles from a string and put them in an array,
   after allocating space for them.
   str = the input string
   n = the (pre-allocated) number of doubles read
   r = the output array of doubles. */
static void parse_n_real(char *str, int *n, real **r, warninp_t wi)
{
    auto values = gmx::splitString(str);
    *n = values.size();

    snew(*r, *n);
    for (int i = 0; i < *n; i++)
    {
        try
        {
            (*r)[i] = gmx::fromString<real>(values[i]);
        }
        catch (gmx::GromacsException &)
        {
            warning_error(wi, "Invalid value " + values[i] + " in string in mdp file. Expected a real number.");
        }
    }
}


static void do_fep_params(t_inputrec *ir, char fep_lambda[][STRLEN], char weights[STRLEN], warninp_t wi)
{

    int         i, j, max_n_lambda, nweights, nfep[efptNR];
    t_lambda   *fep    = ir->fepvals;
    t_expanded *expand = ir->expandedvals;
    real      **count_fep_lambdas;
    bool        bOneLambda = TRUE;

    snew(count_fep_lambdas, efptNR);

    /* FEP input processing */
    /* first, identify the number of lambda values for each type.
       All that are nonzero must have the same number */

    for (i = 0; i < efptNR; i++)
    {
        parse_n_real(fep_lambda[i], &(nfep[i]), &(count_fep_lambdas[i]), wi);
    }

    /* now, determine the number of components.  All must be either zero, or equal. */

    max_n_lambda = 0;
    for (i = 0; i < efptNR; i++)
    {
        if (nfep[i] > max_n_lambda)
        {
            max_n_lambda = nfep[i];  /* here's a nonzero one.  All of them
                                        must have the same number if its not zero.*/
            break;
        }
    }

    for (i = 0; i < efptNR; i++)
    {
        if (nfep[i] == 0)
        {
            ir->fepvals->separate_dvdl[i] = FALSE;
        }
        else if (nfep[i] == max_n_lambda)
        {
            if (i != efptTEMPERATURE)  /* we treat this differently -- not really a reason to compute the derivative with
                                          respect to the temperature currently */
            {
                ir->fepvals->separate_dvdl[i] = TRUE;
            }
        }
        else
        {
            gmx_fatal(FARGS, "Number of lambdas (%d) for FEP type %s not equal to number of other types (%d)",
                      nfep[i], efpt_names[i], max_n_lambda);
        }
    }
    /* we don't print out dhdl if the temperature is changing, since we can't correctly define dhdl in this case */
    ir->fepvals->separate_dvdl[efptTEMPERATURE] = FALSE;

    /* the number of lambdas is the number we've read in, which is either zero
       or the same for all */
    fep->n_lambda = max_n_lambda;

    /* allocate space for the array of lambda values */
    snew(fep->all_lambda, efptNR);
    /* if init_lambda is defined, we need to set lambda */
    if ((fep->init_lambda > 0) && (fep->n_lambda == 0))
    {
        ir->fepvals->separate_dvdl[efptFEP] = TRUE;
    }
    /* otherwise allocate the space for all of the lambdas, and transfer the data */
    for (i = 0; i < efptNR; i++)
    {
        snew(fep->all_lambda[i], fep->n_lambda);
        if (nfep[i] > 0)  /* if it's zero, then the count_fep_lambda arrays
                             are zero */
        {
            for (j = 0; j < fep->n_lambda; j++)
            {
                fep->all_lambda[i][j] = static_cast<double>(count_fep_lambdas[i][j]);
            }
            sfree(count_fep_lambdas[i]);
        }
    }
    sfree(count_fep_lambdas);

    /* "fep-vals" is either zero or the full number. If zero, we'll need to define fep-lambdas for internal
       bookkeeping -- for now, init_lambda */

    if ((nfep[efptFEP] == 0) && (fep->init_lambda >= 0))
    {
        for (i = 0; i < fep->n_lambda; i++)
        {
            fep->all_lambda[efptFEP][i] = fep->init_lambda;
        }
    }

    /* check to see if only a single component lambda is defined, and soft core is defined.
       In this case, turn on coulomb soft core */

    if (max_n_lambda == 0)
    {
        bOneLambda = TRUE;
    }
    else
    {
        for (i = 0; i < efptNR; i++)
        {
            if ((nfep[i] != 0) && (i != efptFEP))
            {
                bOneLambda = FALSE;
            }
        }
    }
    if ((bOneLambda) && (fep->sc_alpha > 0))
    {
        fep->bScCoul = TRUE;
    }

    /* Fill in the others with the efptFEP if they are not explicitly
       specified (i.e. nfep[i] == 0).  This means if fep is not defined,
       they are all zero. */

    for (i = 0; i < efptNR; i++)
    {
        if ((nfep[i] == 0) && (i != efptFEP))
        {
            for (j = 0; j < fep->n_lambda; j++)
            {
                fep->all_lambda[i][j] = fep->all_lambda[efptFEP][j];
            }
        }
    }


    /* make it easier if sc_r_power = 48 by increasing it to the 4th power, to be in the right scale. */
    if (fep->sc_r_power == 48)
    {
        if (fep->sc_alpha > 0.1)
        {
            gmx_fatal(FARGS, "sc_alpha (%f) for sc_r_power = 48 should usually be between 0.001 and 0.004", fep->sc_alpha);
        }
    }

    /* now read in the weights */
    parse_n_real(weights, &nweights, &(expand->init_lambda_weights), wi);
    if (nweights == 0)
    {
        snew(expand->init_lambda_weights, fep->n_lambda); /* initialize to zero */
    }
    else if (nweights != fep->n_lambda)
    {
        gmx_fatal(FARGS, "Number of weights (%d) is not equal to number of lambda values (%d)",
                  nweights, fep->n_lambda);
    }
    if ((expand->nstexpanded < 0) && (ir->efep != efepNO))
    {
        expand->nstexpanded = fep->nstdhdl;
        /* if you don't specify nstexpanded when doing expanded ensemble free energy calcs, it is set to nstdhdl */
    }
    if ((expand->nstexpanded < 0) && ir->bSimTemp)
    {
        expand->nstexpanded = 2*static_cast<int>(ir->opts.tau_t[0]/ir->delta_t);
        /* if you don't specify nstexpanded when doing expanded ensemble simulated tempering, it is set to
           2*tau_t just to be careful so it's not to frequent  */
    }
}


static void do_simtemp_params(t_inputrec *ir)
{

    snew(ir->simtempvals->temperatures, ir->fepvals->n_lambda);
    GetSimTemps(ir->fepvals->n_lambda, ir->simtempvals, ir->fepvals->all_lambda[efptTEMPERATURE]);
}

static void
convertYesNos(warninp_t /*wi*/, gmx::ArrayRef<const std::string> inputs, const char * /*name*/, gmx_bool *outputs)
{
    int i = 0;
    for (const auto &input : inputs)
    {
        outputs[i] = (gmx_strncasecmp(input.c_str(), "Y", 1) == 0);
        ++i;
    }
}

template <typename T> void
convertInts(warninp_t wi, gmx::ArrayRef<const std::string> inputs, const char *name, T *outputs)
{
    int i = 0;
    for (const auto &input : inputs)
    {
        try
        {
            outputs[i] = gmx::fromStdString<T>(input);
        }
        catch (gmx::GromacsException &)
        {
            auto message = gmx::formatString("Invalid value for mdp option %s. %s should only consist of integers separated by spaces.",
                                             name, name);
            warning_error(wi, message);
        }
        ++i;
    }
}

static void
convertReals(warninp_t wi, gmx::ArrayRef<const std::string> inputs, const char *name, real *outputs)
{
    int i = 0;
    for (const auto &input : inputs)
    {
        try
        {
            outputs[i] = gmx::fromString<real>(input);
        }
        catch (gmx::GromacsException &)
        {
            auto message = gmx::formatString("Invalid value for mdp option %s. %s should only consist of real numbers separated by spaces.",
                                             name, name);
            warning_error(wi, message);
        }
        ++i;
    }
}

static void
convertRvecs(warninp_t wi, gmx::ArrayRef<const std::string> inputs, const char *name, rvec *outputs)
{
    int i = 0, d = 0;
    for (const auto &input : inputs)
    {
        try
        {
            outputs[i][d] = gmx::fromString<real>(input);
        }
        catch (gmx::GromacsException &)
        {
            auto message = gmx::formatString("Invalid value for mdp option %s. %s should only consist of real numbers separated by spaces.",
                                             name, name);
            warning_error(wi, message);
        }
        ++d;
        if (d == DIM)
        {
            d = 0;
            ++i;
        }
    }
}

static void do_wall_params(t_inputrec *ir,
                           char *wall_atomtype, char *wall_density,
                           t_gromppopts *opts,
                           warninp_t wi)
{
    opts->wall_atomtype[0] = nullptr;
    opts->wall_atomtype[1] = nullptr;

    ir->wall_atomtype[0] = -1;
    ir->wall_atomtype[1] = -1;
    ir->wall_density[0]  = 0;
    ir->wall_density[1]  = 0;

    if (ir->nwall > 0)
    {
        auto wallAtomTypes = gmx::splitString(wall_atomtype);
        if (wallAtomTypes.size() != size_t(ir->nwall))
        {
            gmx_fatal(FARGS, "Expected %d elements for wall_atomtype, found %zu",
                      ir->nwall, wallAtomTypes.size());
        }
        for (int i = 0; i < ir->nwall; i++)
        {
            opts->wall_atomtype[i] = gmx_strdup(wallAtomTypes[i].c_str());
        }

        if (ir->wall_type == ewt93 || ir->wall_type == ewt104)
        {
            auto wallDensity = gmx::splitString(wall_density);
            if (wallDensity.size() != size_t(ir->nwall))
            {
                gmx_fatal(FARGS, "Expected %d elements for wall-density, found %zu", ir->nwall, wallDensity.size());
            }
            convertReals(wi, wallDensity, "wall-density", ir->wall_density);
            for (int i = 0; i < ir->nwall; i++)
            {
                if (ir->wall_density[i] <= 0)
                {
                    gmx_fatal(FARGS, "wall-density[%d] = %f\n", i, ir->wall_density[i]);
                }
            }
        }
    }
}

static void add_wall_energrps(gmx_groups_t *groups, int nwall, t_symtab *symtab)
{
    int     i;
    t_grps *grps;
    char    str[STRLEN];

    if (nwall > 0)
    {
        srenew(groups->grpname, groups->ngrpname+nwall);
        grps = &(groups->grps[egcENER]);
        srenew(grps->nm_ind, grps->nr+nwall);
        for (i = 0; i < nwall; i++)
        {
            sprintf(str, "wall%d", i);
            groups->grpname[groups->ngrpname] = put_symtab(symtab, str);
            grps->nm_ind[grps->nr++]          = groups->ngrpname++;
        }
    }
}

static void read_expandedparams(std::vector<t_inpfile> *inp,
                                t_expanded *expand, warninp_t wi)
{
    /* read expanded ensemble parameters */
    printStringNewline(inp, "expanded ensemble variables");
    expand->nstexpanded    = get_eint(inp, "nstexpanded", -1, wi);
    expand->elamstats      = get_eeenum(inp, "lmc-stats", elamstats_names, wi);
    expand->elmcmove       = get_eeenum(inp, "lmc-move", elmcmove_names, wi);
    expand->elmceq         = get_eeenum(inp, "lmc-weights-equil", elmceq_names, wi);
    expand->equil_n_at_lam = get_eint(inp, "weight-equil-number-all-lambda", -1, wi);
    expand->equil_samples  = get_eint(inp, "weight-equil-number-samples", -1, wi);
    expand->equil_steps    = get_eint(inp, "weight-equil-number-steps", -1, wi);
    expand->equil_wl_delta = get_ereal(inp, "weight-equil-wl-delta", -1, wi);
    expand->equil_ratio    = get_ereal(inp, "weight-equil-count-ratio", -1, wi);
    printStringNewline(inp, "Seed for Monte Carlo in lambda space");
    expand->lmc_seed            = get_eint(inp, "lmc-seed", -1, wi);
    expand->mc_temp             = get_ereal(inp, "mc-temperature", -1, wi);
    expand->lmc_repeats         = get_eint(inp, "lmc-repeats", 1, wi);
    expand->gibbsdeltalam       = get_eint(inp, "lmc-gibbsdelta", -1, wi);
    expand->lmc_forced_nstart   = get_eint(inp, "lmc-forced-nstart", 0, wi);
    expand->bSymmetrizedTMatrix = (get_eeenum(inp, "symmetrized-transition-matrix", yesno_names, wi) != 0);
    expand->nstTij              = get_eint(inp, "nst-transition-matrix", -1, wi);
    expand->minvarmin           = get_eint(inp, "mininum-var-min", 100, wi); /*default is reasonable */
    expand->c_range             = get_eint(inp, "weight-c-range", 0, wi);    /* default is just C=0 */
    expand->wl_scale            = get_ereal(inp, "wl-scale", 0.8, wi);
    expand->wl_ratio            = get_ereal(inp, "wl-ratio", 0.8, wi);
    expand->init_wl_delta       = get_ereal(inp, "init-wl-delta", 1.0, wi);
    expand->bWLoneovert         = (get_eeenum(inp, "wl-oneovert", yesno_names, wi) != 0);
}

/*! \brief Return whether an end state with the given coupling-lambda
 * value describes fully-interacting VDW.
 *
 * \param[in]  couple_lambda_value  Enumeration ecouplam value describing the end state
 * \return                          Whether VDW is on (i.e. the user chose vdw or vdw-q in the .mdp file)
 */
static bool couple_lambda_has_vdw_on(int couple_lambda_value)
{
    return (couple_lambda_value == ecouplamVDW ||
            couple_lambda_value == ecouplamVDWQ);
}

namespace
{

class MdpErrorHandler : public gmx::IKeyValueTreeErrorHandler
{
    public:
        explicit MdpErrorHandler(warninp_t wi)
            : wi_(wi), mapping_(nullptr)
        {
        }

        void setBackMapping(const gmx::IKeyValueTreeBackMapping &mapping)
        {
            mapping_ = &mapping;
        }

        bool onError(gmx::UserInputError *ex, const gmx::KeyValueTreePath &context) override
        {
            ex->prependContext(gmx::formatString("Error in mdp option \"%s\":",
                                                 getOptionName(context).c_str()));
            std::string message = gmx::formatExceptionMessageToString(*ex);
            warning_error(wi_, message.c_str());
            return true;
        }

    private:
        std::string getOptionName(const gmx::KeyValueTreePath &context)
        {
            if (mapping_ != nullptr)
            {
                gmx::KeyValueTreePath path = mapping_->originalPath(context);
                GMX_ASSERT(path.size() == 1, "Inconsistent mapping back to mdp options");
                return path[0];
            }
            GMX_ASSERT(context.size() == 1, "Inconsistent context for mdp option parsing");
            return context[0];
        }

        warninp_t                            wi_;
        const gmx::IKeyValueTreeBackMapping *mapping_;
};

} // namespace

void get_ir(const char *mdparin, const char *mdparout,
            gmx::MDModules *mdModules, t_inputrec *ir, t_gromppopts *opts,
            WriteMdpHeader writeMdpHeader, warninp_t wi)
{
    char                  *dumstr[2];
    double                 dumdub[2][6];
    int                    i, j, m;
    char                   warn_buf[STRLEN];
    t_lambda              *fep    = ir->fepvals;
    t_expanded            *expand = ir->expandedvals;

    const char            *no_names[] = { "no", nullptr };

    init_inputrec_strings();
    gmx::TextInputFile     stream(mdparin);
    std::vector<t_inpfile> inp = read_inpfile(&stream, mdparin, wi);

    snew(dumstr[0], STRLEN);
    snew(dumstr[1], STRLEN);

    if (-1 == search_einp(inp, "cutoff-scheme"))
    {
        sprintf(warn_buf,
                "%s did not specify a value for the .mdp option "
                "\"cutoff-scheme\". Probably it was first intended for use "
                "with GROMACS before 4.6. In 4.6, the Verlet scheme was "
                "introduced, but the group scheme was still the default. "
                "The default is now the Verlet scheme, so you will observe "
                "different behaviour.", mdparin);
        warning_note(wi, warn_buf);
    }

    /* ignore the following deprecated commands */
    replace_inp_entry(inp, "title", nullptr);
    replace_inp_entry(inp, "cpp", nullptr);
    replace_inp_entry(inp, "domain-decomposition", nullptr);
    replace_inp_entry(inp, "andersen-seed", nullptr);
    replace_inp_entry(inp, "dihre", nullptr);
    replace_inp_entry(inp, "dihre-fc", nullptr);
    replace_inp_entry(inp, "dihre-tau", nullptr);
    replace_inp_entry(inp, "nstdihreout", nullptr);
    replace_inp_entry(inp, "nstcheckpoint", nullptr);
    replace_inp_entry(inp, "optimize-fft", nullptr);
    replace_inp_entry(inp, "adress_type", nullptr);
    replace_inp_entry(inp, "adress_const_wf", nullptr);
    replace_inp_entry(inp, "adress_ex_width", nullptr);
    replace_inp_entry(inp, "adress_hy_width", nullptr);
    replace_inp_entry(inp, "adress_ex_forcecap", nullptr);
    replace_inp_entry(inp, "adress_interface_correction", nullptr);
    replace_inp_entry(inp, "adress_site", nullptr);
    replace_inp_entry(inp, "adress_reference_coords", nullptr);
    replace_inp_entry(inp, "adress_tf_grp_names", nullptr);
    replace_inp_entry(inp, "adress_cg_grp_names", nullptr);
    replace_inp_entry(inp, "adress_do_hybridpairs", nullptr);
    replace_inp_entry(inp, "rlistlong", nullptr);
    replace_inp_entry(inp, "nstcalclr", nullptr);
    replace_inp_entry(inp, "pull-print-com2", nullptr);
    replace_inp_entry(inp, "gb-algorithm", nullptr);
    replace_inp_entry(inp, "nstgbradii", nullptr);
    replace_inp_entry(inp, "rgbradii", nullptr);
    replace_inp_entry(inp, "gb-epsilon-solvent", nullptr);
    replace_inp_entry(inp, "gb-saltconc", nullptr);
    replace_inp_entry(inp, "gb-obc-alpha", nullptr);
    replace_inp_entry(inp, "gb-obc-beta", nullptr);
    replace_inp_entry(inp, "gb-obc-gamma", nullptr);
    replace_inp_entry(inp, "gb-dielectric-offset", nullptr);
    replace_inp_entry(inp, "sa-algorithm", nullptr);
    replace_inp_entry(inp, "sa-surface-tension", nullptr);

    /* replace the following commands with the clearer new versions*/
    replace_inp_entry(inp, "unconstrained-start", "continuation");
    replace_inp_entry(inp, "foreign-lambda", "fep-lambdas");
    replace_inp_entry(inp, "verlet-buffer-drift", "verlet-buffer-tolerance");
    replace_inp_entry(inp, "nstxtcout", "nstxout-compressed");
    replace_inp_entry(inp, "xtc-grps", "compressed-x-grps");
    replace_inp_entry(inp, "xtc-precision", "compressed-x-precision");
    replace_inp_entry(inp, "pull-print-com1", "pull-print-com");

    printStringNewline(&inp, "VARIOUS PREPROCESSING OPTIONS");
    printStringNoNewline(&inp, "Preprocessor information: use cpp syntax.");
    printStringNoNewline(&inp, "e.g.: -I/home/joe/doe -I/home/mary/roe");
    setStringEntry(&inp, "include", opts->include,  nullptr);
    printStringNoNewline(&inp, "e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)");
    setStringEntry(&inp, "define",  opts->define,   nullptr);

    printStringNewline(&inp, "RUN CONTROL PARAMETERS");
    ir->eI = get_eeenum(&inp, "integrator",         ei_names, wi);
    printStringNoNewline(&inp, "Start time and timestep in ps");
    ir->init_t  = get_ereal(&inp, "tinit", 0.0, wi);
    ir->delta_t = get_ereal(&inp, "dt",    0.001, wi);
    ir->nsteps  = get_eint64(&inp, "nsteps",     0, wi);
    printStringNoNewline(&inp, "For exact run continuation or redoing part of a run");
    ir->init_step = get_eint64(&inp, "init-step",  0, wi);
    printStringNoNewline(&inp, "Part index is updated automatically on checkpointing (keeps files separate)");
    ir->simulation_part = get_eint(&inp, "simulation-part", 1, wi);
    printStringNoNewline(&inp, "mode for center of mass motion removal");
    ir->comm_mode = get_eeenum(&inp, "comm-mode",  ecm_names, wi);
    printStringNoNewline(&inp, "number of steps for center of mass motion removal");
    ir->nstcomm = get_eint(&inp, "nstcomm",    100, wi);
    printStringNoNewline(&inp, "group(s) for center of mass motion removal");
    setStringEntry(&inp, "comm-grps",   is->vcm,            nullptr);

    printStringNewline(&inp, "LANGEVIN DYNAMICS OPTIONS");
    printStringNoNewline(&inp, "Friction coefficient (amu/ps) and random seed");
    ir->bd_fric = get_ereal(&inp, "bd-fric",    0.0, wi);
    ir->ld_seed = get_eint64(&inp, "ld-seed",    -1, wi);

    /* Em stuff */
    printStringNewline(&inp, "ENERGY MINIMIZATION OPTIONS");
    printStringNoNewline(&inp, "Force tolerance and initial step-size");
    ir->em_tol      = get_ereal(&inp, "emtol",     10.0, wi);
    ir->em_stepsize = get_ereal(&inp, "emstep", 0.01, wi);
    printStringNoNewline(&inp, "Max number of iterations in relax-shells");
    ir->niter = get_eint(&inp, "niter",      20, wi);
    printStringNoNewline(&inp, "Step size (ps^2) for minimization of flexible constraints");
    ir->fc_stepsize = get_ereal(&inp, "fcstep", 0, wi);
    printStringNoNewline(&inp, "Frequency of steepest descents steps when doing CG");
    ir->nstcgsteep = get_eint(&inp, "nstcgsteep", 1000, wi);
    ir->nbfgscorr  = get_eint(&inp, "nbfgscorr",  10, wi);

    printStringNewline(&inp, "TEST PARTICLE INSERTION OPTIONS");
    ir->rtpi = get_ereal(&inp, "rtpi",   0.05, wi);

    /* Output options */
    printStringNewline(&inp, "OUTPUT CONTROL OPTIONS");
    printStringNoNewline(&inp, "Output frequency for coords (x), velocities (v) and forces (f)");
    ir->nstxout = get_eint(&inp, "nstxout",    0, wi);
    ir->nstvout = get_eint(&inp, "nstvout",    0, wi);
    ir->nstfout = get_eint(&inp, "nstfout",    0, wi);
    printStringNoNewline(&inp, "Output frequency for energies to log file and energy file");
    ir->nstlog        = get_eint(&inp, "nstlog", 1000, wi);
    ir->nstcalcenergy = get_eint(&inp, "nstcalcenergy", 100, wi);
    ir->nstenergy     = get_eint(&inp, "nstenergy",  1000, wi);
    printStringNoNewline(&inp, "Output frequency and precision for .xtc file");
    ir->nstxout_compressed      = get_eint(&inp, "nstxout-compressed",  0, wi);
    ir->x_compression_precision = get_ereal(&inp, "compressed-x-precision", 1000.0, wi);
    printStringNoNewline(&inp, "This selects the subset of atoms for the compressed");
    printStringNoNewline(&inp, "trajectory file. You can select multiple groups. By");
    printStringNoNewline(&inp, "default, all atoms will be written.");
    setStringEntry(&inp, "compressed-x-grps", is->x_compressed_groups, nullptr);
    printStringNoNewline(&inp, "Selection of energy groups");
    setStringEntry(&inp, "energygrps",  is->energy,         nullptr);

    /* Neighbor searching */
    printStringNewline(&inp, "NEIGHBORSEARCHING PARAMETERS");
    printStringNoNewline(&inp, "cut-off scheme (Verlet: particle based cut-offs, group: using charge groups)");
    ir->cutoff_scheme = get_eeenum(&inp, "cutoff-scheme",    ecutscheme_names, wi);
    printStringNoNewline(&inp, "nblist update frequency");
    ir->nstlist = get_eint(&inp, "nstlist",    10, wi);
    printStringNoNewline(&inp, "ns algorithm (simple or grid)");
    ir->ns_type = get_eeenum(&inp, "ns-type",    ens_names, wi);
    printStringNoNewline(&inp, "Periodic boundary conditions: xyz, no, xy");
    ir->ePBC          = get_eeenum(&inp, "pbc",       epbc_names, wi);
    ir->bPeriodicMols = get_eeenum(&inp, "periodic-molecules", yesno_names, wi) != 0;
    printStringNoNewline(&inp, "Allowed energy error due to the Verlet buffer in kJ/mol/ps per atom,");
    printStringNoNewline(&inp, "a value of -1 means: use rlist");
    ir->verletbuf_tol = get_ereal(&inp, "verlet-buffer-tolerance",    0.005, wi);
    printStringNoNewline(&inp, "nblist cut-off");
    ir->rlist = get_ereal(&inp, "rlist",  1.0, wi);
    printStringNoNewline(&inp, "long-range cut-off for switched potentials");

    /* Electrostatics */
    printStringNewline(&inp, "OPTIONS FOR ELECTROSTATICS AND VDW");
    printStringNoNewline(&inp, "Method for doing electrostatics");
    ir->coulombtype      = get_eeenum(&inp, "coulombtype",    eel_names, wi);
    ir->coulomb_modifier = get_eeenum(&inp, "coulomb-modifier",    eintmod_names, wi);
    printStringNoNewline(&inp, "cut-off lengths");
    ir->rcoulomb_switch = get_ereal(&inp, "rcoulomb-switch",    0.0, wi);
    ir->rcoulomb        = get_ereal(&inp, "rcoulomb",   1.0, wi);
    printStringNoNewline(&inp, "Relative dielectric constant for the medium and the reaction field");
    ir->epsilon_r  = get_ereal(&inp, "epsilon-r",  1.0, wi);
    ir->epsilon_rf = get_ereal(&inp, "epsilon-rf", 0.0, wi);
    printStringNoNewline(&inp, "Method for doing Van der Waals");
    ir->vdwtype      = get_eeenum(&inp, "vdw-type",    evdw_names, wi);
    ir->vdw_modifier = get_eeenum(&inp, "vdw-modifier",    eintmod_names, wi);
    printStringNoNewline(&inp, "cut-off lengths");
    ir->rvdw_switch = get_ereal(&inp, "rvdw-switch",    0.0, wi);
    ir->rvdw        = get_ereal(&inp, "rvdw",   1.0, wi);
    printStringNoNewline(&inp, "Apply long range dispersion corrections for Energy and Pressure");
    ir->eDispCorr = get_eeenum(&inp, "DispCorr",  edispc_names, wi);
    printStringNoNewline(&inp, "Extension of the potential lookup tables beyond the cut-off");
    ir->tabext = get_ereal(&inp, "table-extension", 1.0, wi);
    printStringNoNewline(&inp, "Separate tables between energy group pairs");
    setStringEntry(&inp, "energygrp-table", is->egptable,   nullptr);
    printStringNoNewline(&inp, "Spacing for the PME/PPPM FFT grid");
    ir->fourier_spacing = get_ereal(&inp, "fourierspacing", 0.12, wi);
    printStringNoNewline(&inp, "FFT grid size, when a value is 0 fourierspacing will be used");
    ir->nkx = get_eint(&inp, "fourier-nx",         0, wi);
    ir->nky = get_eint(&inp, "fourier-ny",         0, wi);
    ir->nkz = get_eint(&inp, "fourier-nz",         0, wi);
    printStringNoNewline(&inp, "EWALD/PME/PPPM parameters");
    ir->pme_order              = get_eint(&inp, "pme-order",   4, wi);
    ir->ewald_rtol             = get_ereal(&inp, "ewald-rtol", 0.00001, wi);
    ir->ewald_rtol_lj          = get_ereal(&inp, "ewald-rtol-lj", 0.001, wi);
    ir->ljpme_combination_rule = get_eeenum(&inp, "lj-pme-comb-rule", eljpme_names, wi);
    ir->ewald_geometry         = get_eeenum(&inp, "ewald-geometry", eewg_names, wi);
    ir->epsilon_surface        = get_ereal(&inp, "epsilon-surface", 0.0, wi);

    /* Implicit solvation is no longer supported, but we need grompp
       to be able to refuse old .mdp files that would have built a tpr
       to run it. Thus, only "no" is accepted. */
    ir->implicit_solvent = (get_eeenum(&inp, "implicit-solvent", no_names, wi) != 0);

    /* Coupling stuff */
    printStringNewline(&inp, "OPTIONS FOR WEAK COUPLING ALGORITHMS");
    printStringNoNewline(&inp, "Temperature coupling");
    ir->etc                = get_eeenum(&inp, "tcoupl",        etcoupl_names, wi);
    ir->nsttcouple         = get_eint(&inp, "nsttcouple",  -1, wi);
    ir->opts.nhchainlength = get_eint(&inp, "nh-chain-length", 10, wi);
    ir->bPrintNHChains     = (get_eeenum(&inp, "print-nose-hoover-chain-variables", yesno_names, wi) != 0);
    printStringNoNewline(&inp, "Groups to couple separately");
    setStringEntry(&inp, "tc-grps",     is->tcgrps,         nullptr);
    printStringNoNewline(&inp, "Time constant (ps) and reference temperature (K)");
    setStringEntry(&inp, "tau-t",   is->tau_t,      nullptr);
    setStringEntry(&inp, "ref-t",   is->ref_t,      nullptr);
    printStringNoNewline(&inp, "pressure coupling");
    ir->epc        = get_eeenum(&inp, "pcoupl",        epcoupl_names, wi);
    ir->epct       = get_eeenum(&inp, "pcoupltype",       epcoupltype_names, wi);
    ir->nstpcouple = get_eint(&inp, "nstpcouple",  -1, wi);
    printStringNoNewline(&inp, "Time constant (ps), compressibility (1/bar) and reference P (bar)");
    ir->tau_p = get_ereal(&inp, "tau-p",  1.0, wi);
    setStringEntry(&inp, "compressibility", dumstr[0],  nullptr);
    setStringEntry(&inp, "ref-p",       dumstr[1],      nullptr);
    printStringNoNewline(&inp, "Scaling of reference coordinates, No, All or COM");
    ir->refcoord_scaling = get_eeenum(&inp, "refcoord-scaling", erefscaling_names, wi);

    /* QMMM */
    printStringNewline(&inp, "OPTIONS FOR QMMM calculations");
    ir->bQMMM = (get_eeenum(&inp, "QMMM", yesno_names, wi) != 0);
    printStringNoNewline(&inp, "Groups treated Quantum Mechanically");
    setStringEntry(&inp, "QMMM-grps",  is->QMMM,          nullptr);
    printStringNoNewline(&inp, "QM method");
    setStringEntry(&inp, "QMmethod",     is->QMmethod, nullptr);
    printStringNoNewline(&inp, "QMMM scheme");
    ir->QMMMscheme = get_eeenum(&inp, "QMMMscheme",    eQMMMscheme_names, wi);
    printStringNoNewline(&inp, "QM basisset");
    setStringEntry(&inp, "QMbasis",      is->QMbasis, nullptr);
    printStringNoNewline(&inp, "QM charge");
    setStringEntry(&inp, "QMcharge",    is->QMcharge, nullptr);
    printStringNoNewline(&inp, "QM multiplicity");
    setStringEntry(&inp, "QMmult",      is->QMmult, nullptr);
    printStringNoNewline(&inp, "Surface Hopping");
    setStringEntry(&inp, "SH",          is->bSH, nullptr);
    printStringNoNewline(&inp, "CAS space options");
    setStringEntry(&inp, "CASorbitals",      is->CASorbitals,   nullptr);
    setStringEntry(&inp, "CASelectrons",     is->CASelectrons,  nullptr);
    setStringEntry(&inp, "SAon", is->SAon, nullptr);
    setStringEntry(&inp, "SAoff", is->SAoff, nullptr);
    setStringEntry(&inp, "SAsteps", is->SAsteps, nullptr);
    printStringNoNewline(&inp, "Scale factor for MM charges");
    ir->scalefactor = get_ereal(&inp, "MMChargeScaleFactor", 1.0, wi);

    /* Simulated annealing */
    printStringNewline(&inp, "SIMULATED ANNEALING");
    printStringNoNewline(&inp, "Type of annealing for each temperature group (no/single/periodic)");
    setStringEntry(&inp, "annealing",   is->anneal,      nullptr);
    printStringNoNewline(&inp, "Number of time points to use for specifying annealing in each group");
    setStringEntry(&inp, "annealing-npoints", is->anneal_npoints, nullptr);
    printStringNoNewline(&inp, "List of times at the annealing points for each group");
    setStringEntry(&inp, "annealing-time",       is->anneal_time,       nullptr);
    printStringNoNewline(&inp, "Temp. at each annealing point, for each group.");
    setStringEntry(&inp, "annealing-temp",  is->anneal_temp,  nullptr);

    /* Startup run */
    printStringNewline(&inp, "GENERATE VELOCITIES FOR STARTUP RUN");
    opts->bGenVel = (get_eeenum(&inp, "gen-vel",  yesno_names, wi) != 0);
    opts->tempi   = get_ereal(&inp, "gen-temp",    300.0, wi);
    opts->seed    = get_eint(&inp, "gen-seed",     -1, wi);

    /* Shake stuff */
    printStringNewline(&inp, "OPTIONS FOR BONDS");
    opts->nshake = get_eeenum(&inp, "constraints",   constraints, wi);
    printStringNoNewline(&inp, "Type of constraint algorithm");
    ir->eConstrAlg = get_eeenum(&inp, "constraint-algorithm", econstr_names, wi);
    printStringNoNewline(&inp, "Do not constrain the start configuration");
    ir->bContinuation = (get_eeenum(&inp, "continuation", yesno_names, wi) != 0);
    printStringNoNewline(&inp, "Use successive overrelaxation to reduce the number of shake iterations");
    ir->bShakeSOR = (get_eeenum(&inp, "Shake-SOR", yesno_names, wi) != 0);
    printStringNoNewline(&inp, "Relative tolerance of shake");
    ir->shake_tol = get_ereal(&inp, "shake-tol", 0.0001, wi);
    printStringNoNewline(&inp, "Highest order in the expansion of the constraint coupling matrix");
    ir->nProjOrder = get_eint(&inp, "lincs-order", 4, wi);
    printStringNoNewline(&inp, "Number of iterations in the final step of LINCS. 1 is fine for");
    printStringNoNewline(&inp, "normal simulations, but use 2 to conserve energy in NVE runs.");
    printStringNoNewline(&inp, "For energy minimization with constraints it should be 4 to 8.");
    ir->nLincsIter = get_eint(&inp, "lincs-iter", 1, wi);
    printStringNoNewline(&inp, "Lincs will write a warning to the stderr if in one step a bond");
    printStringNoNewline(&inp, "rotates over more degrees than");
    ir->LincsWarnAngle = get_ereal(&inp, "lincs-warnangle", 30.0, wi);
    printStringNoNewline(&inp, "Convert harmonic bonds to morse potentials");
    opts->bMorse = (get_eeenum(&inp, "morse", yesno_names, wi) != 0);

    /* Energy group exclusions */
    printStringNewline(&inp, "ENERGY GROUP EXCLUSIONS");
    printStringNoNewline(&inp, "Pairs of energy groups for which all non-bonded interactions are excluded");
    setStringEntry(&inp, "energygrp-excl", is->egpexcl,     nullptr);

    /* Walls */
    printStringNewline(&inp, "WALLS");
    printStringNoNewline(&inp, "Number of walls, type, atom types, densities and box-z scale factor for Ewald");
    ir->nwall         = get_eint(&inp, "nwall", 0, wi);
    ir->wall_type     = get_eeenum(&inp, "wall-type",   ewt_names, wi);
    ir->wall_r_linpot = get_ereal(&inp, "wall-r-linpot", -1, wi);
    setStringEntry(&inp, "wall-atomtype", is->wall_atomtype, nullptr);
    setStringEntry(&inp, "wall-density",  is->wall_density,  nullptr);
    ir->wall_ewald_zfac = get_ereal(&inp, "wall-ewald-zfac", 3, wi);

    /* COM pulling */
    printStringNewline(&inp, "COM PULLING");
    ir->bPull = (get_eeenum(&inp, "pull", yesno_names, wi) != 0);
    if (ir->bPull)
    {
        snew(ir->pull, 1);
        is->pull_grp = read_pullparams(&inp, ir->pull, wi);
    }

    /* AWH biasing
       NOTE: needs COM pulling input */
    printStringNewline(&inp, "AWH biasing");
    ir->bDoAwh = (get_eeenum(&inp, "awh", yesno_names, wi) != 0);
    if (ir->bDoAwh)
    {
        if (ir->bPull)
        {
            ir->awhParams = gmx::readAndCheckAwhParams(&inp, ir, wi);
        }
        else
        {
            gmx_fatal(FARGS, "AWH biasing is only compatible with COM pulling turned on");
        }
    }

    /* Enforced rotation */
    printStringNewline(&inp, "ENFORCED ROTATION");
    printStringNoNewline(&inp, "Enforced rotation: No or Yes");
    ir->bRot = (get_eeenum(&inp, "rotation", yesno_names, wi) != 0);
    if (ir->bRot)
    {
        snew(ir->rot, 1);
        is->rot_grp = read_rotparams(&inp, ir->rot, wi);
    }

    /* Interactive MD */
    ir->bIMD = FALSE;
    printStringNewline(&inp, "Group to display and/or manipulate in interactive MD session");
    setStringEntry(&inp, "IMD-group", is->imd_grp, nullptr);
    if (is->imd_grp[0] != '\0')
    {
        snew(ir->imd, 1);
        ir->bIMD = TRUE;
    }

    /* Refinement */
    printStringNewline(&inp, "NMR refinement stuff");
    printStringNoNewline(&inp, "Distance restraints type: No, Simple or Ensemble");
    ir->eDisre = get_eeenum(&inp, "disre",     edisre_names, wi);
    printStringNoNewline(&inp, "Force weighting of pairs in one distance restraint: Conservative or Equal");
    ir->eDisreWeighting = get_eeenum(&inp, "disre-weighting", edisreweighting_names, wi);
    printStringNoNewline(&inp, "Use sqrt of the time averaged times the instantaneous violation");
    ir->bDisreMixed = (get_eeenum(&inp, "disre-mixed", yesno_names, wi) != 0);
    ir->dr_fc       = get_ereal(&inp, "disre-fc",  1000.0, wi);
    ir->dr_tau      = get_ereal(&inp, "disre-tau", 0.0, wi);
    printStringNoNewline(&inp, "Output frequency for pair distances to energy file");
    ir->nstdisreout = get_eint(&inp, "nstdisreout", 100, wi);
    printStringNoNewline(&inp, "Orientation restraints: No or Yes");
    opts->bOrire = (get_eeenum(&inp, "orire",   yesno_names, wi) != 0);
    printStringNoNewline(&inp, "Orientation restraints force constant and tau for time averaging");
    ir->orires_fc  = get_ereal(&inp, "orire-fc",  0.0, wi);
    ir->orires_tau = get_ereal(&inp, "orire-tau", 0.0, wi);
    setStringEntry(&inp, "orire-fitgrp", is->orirefitgrp,    nullptr);
    printStringNoNewline(&inp, "Output frequency for trace(SD) and S to energy file");
    ir->nstorireout = get_eint(&inp, "nstorireout", 100, wi);

    /* free energy variables */
    printStringNewline(&inp, "Free energy variables");
    ir->efep = get_eeenum(&inp, "free-energy", efep_names, wi);
    setStringEntry(&inp, "couple-moltype",  is->couple_moltype,  nullptr);
    opts->couple_lam0  = get_eeenum(&inp, "couple-lambda0", couple_lam, wi);
    opts->couple_lam1  = get_eeenum(&inp, "couple-lambda1", couple_lam, wi);
    opts->bCoupleIntra = (get_eeenum(&inp, "couple-intramol", yesno_names, wi) != 0);

    fep->init_lambda    = get_ereal(&inp, "init-lambda", -1, wi); /* start with -1 so
                                                                            we can recognize if
                                                                            it was not entered */
    fep->init_fep_state = get_eint(&inp, "init-lambda-state", -1, wi);
    fep->delta_lambda   = get_ereal(&inp, "delta-lambda", 0.0, wi);
    fep->nstdhdl        = get_eint(&inp, "nstdhdl", 50, wi);
    setStringEntry(&inp, "fep-lambdas", is->fep_lambda[efptFEP], nullptr);
    setStringEntry(&inp, "mass-lambdas", is->fep_lambda[efptMASS], nullptr);
    setStringEntry(&inp, "coul-lambdas", is->fep_lambda[efptCOUL], nullptr);
    setStringEntry(&inp, "vdw-lambdas", is->fep_lambda[efptVDW], nullptr);
    setStringEntry(&inp, "bonded-lambdas", is->fep_lambda[efptBONDED], nullptr);
    setStringEntry(&inp, "restraint-lambdas", is->fep_lambda[efptRESTRAINT], nullptr);
    setStringEntry(&inp, "temperature-lambdas", is->fep_lambda[efptTEMPERATURE], nullptr);
    fep->lambda_neighbors = get_eint(&inp, "calc-lambda-neighbors", 1, wi);
    setStringEntry(&inp, "init-lambda-weights", is->lambda_weights, nullptr);
    fep->edHdLPrintEnergy   = get_eeenum(&inp, "dhdl-print-energy", edHdLPrintEnergy_names, wi);
    fep->sc_alpha           = get_ereal(&inp, "sc-alpha", 0.0, wi);
    fep->sc_power           = get_eint(&inp, "sc-power", 1, wi);
    fep->sc_r_power         = get_ereal(&inp, "sc-r-power", 6.0, wi);
    fep->sc_sigma           = get_ereal(&inp, "sc-sigma", 0.3, wi);
    fep->bScCoul            = (get_eeenum(&inp, "sc-coul", yesno_names, wi) != 0);
    fep->dh_hist_size       = get_eint(&inp, "dh_hist_size", 0, wi);
    fep->dh_hist_spacing    = get_ereal(&inp, "dh_hist_spacing", 0.1, wi);
    fep->separate_dhdl_file = get_eeenum(&inp, "separate-dhdl-file", separate_dhdl_file_names, wi);
    fep->dhdl_derivatives   = get_eeenum(&inp, "dhdl-derivatives", dhdl_derivatives_names, wi);
    fep->dh_hist_size       = get_eint(&inp, "dh_hist_size", 0, wi);
    fep->dh_hist_spacing    = get_ereal(&inp, "dh_hist_spacing", 0.1, wi);

    /* Non-equilibrium MD stuff */
    printStringNewline(&inp, "Non-equilibrium MD stuff");
    setStringEntry(&inp, "acc-grps",    is->accgrps,        nullptr);
    setStringEntry(&inp, "accelerate",  is->acc,            nullptr);
    setStringEntry(&inp, "freezegrps",  is->freeze,         nullptr);
    setStringEntry(&inp, "freezedim",   is->frdim,          nullptr);
    ir->cos_accel = get_ereal(&inp, "cos-acceleration", 0, wi);
    setStringEntry(&inp, "deform",      is->deform,         nullptr);

    /* simulated tempering variables */
    printStringNewline(&inp, "simulated tempering variables");
    ir->bSimTemp                   = (get_eeenum(&inp, "simulated-tempering", yesno_names, wi) != 0);
    ir->simtempvals->eSimTempScale = get_eeenum(&inp, "simulated-tempering-scaling", esimtemp_names, wi);
    ir->simtempvals->simtemp_low   = get_ereal(&inp, "sim-temp-low", 300.0, wi);
    ir->simtempvals->simtemp_high  = get_ereal(&inp, "sim-temp-high", 300.0, wi);

    /* expanded ensemble variables */
    if (ir->efep == efepEXPANDED || ir->bSimTemp)
    {
        read_expandedparams(&inp, expand, wi);
    }

    /* Electric fields */
    {
        gmx::KeyValueTreeObject      convertedValues = flatKeyValueTreeFromInpFile(inp);
        gmx::KeyValueTreeTransformer transform;
        transform.rules()->addRule()
            .keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
        mdModules->initMdpTransform(transform.rules());
        for (const auto &path : transform.mappedPaths())
        {
            GMX_ASSERT(path.size() == 1, "Inconsistent mapping back to mdp options");
            mark_einp_set(inp, path[0].c_str());
        }
        MdpErrorHandler              errorHandler(wi);
        auto                         result
                   = transform.transform(convertedValues, &errorHandler);
        ir->params = new gmx::KeyValueTreeObject(result.object());
        mdModules->adjustInputrecBasedOnModules(ir);
        errorHandler.setBackMapping(result.backMapping());
        mdModules->assignOptionsToModules(*ir->params, &errorHandler);
    }

    /* Ion/water position swapping ("computational electrophysiology") */
    printStringNewline(&inp, "Ion/water position swapping for computational electrophysiology setups");
    printStringNoNewline(&inp, "Swap positions along direction: no, X, Y, Z");
    ir->eSwapCoords = get_eeenum(&inp, "swapcoords", eSwapTypes_names, wi);
    if (ir->eSwapCoords != eswapNO)
    {
        char buf[STRLEN];
        int  nIonTypes;


        snew(ir->swap, 1);
        printStringNoNewline(&inp, "Swap attempt frequency");
        ir->swap->nstswap = get_eint(&inp, "swap-frequency", 1, wi);
        printStringNoNewline(&inp, "Number of ion types to be controlled");
        nIonTypes = get_eint(&inp, "iontypes", 1, wi);
        if (nIonTypes < 1)
        {
            warning_error(wi, "You need to provide at least one ion type for position exchanges.");
        }
        ir->swap->ngrp = nIonTypes + eSwapFixedGrpNR;
        snew(ir->swap->grp, ir->swap->ngrp);
        for (i = 0; i < ir->swap->ngrp; i++)
        {
            snew(ir->swap->grp[i].molname, STRLEN);
        }
        printStringNoNewline(&inp, "Two index groups that contain the compartment-partitioning atoms");
        setStringEntry(&inp, "split-group0", ir->swap->grp[eGrpSplit0].molname, nullptr);
        setStringEntry(&inp, "split-group1", ir->swap->grp[eGrpSplit1].molname, nullptr);
        printStringNoNewline(&inp, "Use center of mass of split groups (yes/no), otherwise center of geometry is used");
        ir->swap->massw_split[0] = (get_eeenum(&inp, "massw-split0", yesno_names, wi) != 0);
        ir->swap->massw_split[1] = (get_eeenum(&inp, "massw-split1", yesno_names, wi) != 0);

        printStringNoNewline(&inp, "Name of solvent molecules");
        setStringEntry(&inp, "solvent-group", ir->swap->grp[eGrpSolvent].molname, nullptr);

        printStringNoNewline(&inp, "Split cylinder: radius, upper and lower extension (nm) (this will define the channels)");
        printStringNoNewline(&inp, "Note that the split cylinder settings do not have an influence on the swapping protocol,");
        printStringNoNewline(&inp, "however, if correctly defined, the permeation events are recorded per channel");
        ir->swap->cyl0r = get_ereal(&inp, "cyl0-r", 2.0, wi);
        ir->swap->cyl0u = get_ereal(&inp, "cyl0-up", 1.0, wi);
        ir->swap->cyl0l = get_ereal(&inp, "cyl0-down", 1.0, wi);
        ir->swap->cyl1r = get_ereal(&inp, "cyl1-r", 2.0, wi);
        ir->swap->cyl1u = get_ereal(&inp, "cyl1-up", 1.0, wi);
        ir->swap->cyl1l = get_ereal(&inp, "cyl1-down", 1.0, wi);

        printStringNoNewline(&inp, "Average the number of ions per compartment over these many swap attempt steps");
        ir->swap->nAverage = get_eint(&inp, "coupl-steps", 10, wi);

        printStringNoNewline(&inp, "Names of the ion types that can be exchanged with solvent molecules,");
        printStringNoNewline(&inp, "and the requested number of ions of this type in compartments A and B");
        printStringNoNewline(&inp, "-1 means fix the numbers as found in step 0");
        for (i = 0; i < nIonTypes; i++)
        {
            int ig = eSwapFixedGrpNR + i;

            sprintf(buf, "iontype%d-name", i);
            setStringEntry(&inp, buf, ir->swap->grp[ig].molname, nullptr);
            sprintf(buf, "iontype%d-in-A", i);
            ir->swap->grp[ig].nmolReq[0] = get_eint(&inp, buf, -1, wi);
            sprintf(buf, "iontype%d-in-B", i);
            ir->swap->grp[ig].nmolReq[1] = get_eint(&inp, buf, -1, wi);
        }

        printStringNoNewline(&inp, "By default (i.e. bulk offset = 0.0), ion/water exchanges happen between layers");
        printStringNoNewline(&inp, "at maximum distance (= bulk concentration) to the split group layers. However,");
        printStringNoNewline(&inp, "an offset b (-1.0 < b < +1.0) can be specified to offset the bulk layer from the middle at 0.0");
        printStringNoNewline(&inp, "towards one of the compartment-partitioning layers (at +/- 1.0).");
        ir->swap->bulkOffset[0] = get_ereal(&inp, "bulk-offsetA", 0.0, wi);
        ir->swap->bulkOffset[1] = get_ereal(&inp, "bulk-offsetB", 0.0, wi);
        if (!(ir->swap->bulkOffset[0] > -1.0 && ir->swap->bulkOffset[0] < 1.0)
            || !(ir->swap->bulkOffset[1] > -1.0 && ir->swap->bulkOffset[1] < 1.0) )
        {
            warning_error(wi, "Bulk layer offsets must be > -1.0 and < 1.0 !");
        }

        printStringNoNewline(&inp, "Start to swap ions if threshold difference to requested count is reached");
        ir->swap->threshold = get_ereal(&inp, "threshold", 1.0, wi);
    }

    /* AdResS is no longer supported, but we need grompp to be able to
       refuse to process old .mdp files that used it. */
    ir->bAdress = (get_eeenum(&inp, "adress", no_names, wi) != 0);

    /* User defined thingies */
    printStringNewline(&inp, "User defined thingies");
    setStringEntry(&inp, "user1-grps",  is->user1,          nullptr);
    setStringEntry(&inp, "user2-grps",  is->user2,          nullptr);
    ir->userint1  = get_eint(&inp, "userint1",   0, wi);
    ir->userint2  = get_eint(&inp, "userint2",   0, wi);
    ir->userint3  = get_eint(&inp, "userint3",   0, wi);
    ir->userint4  = get_eint(&inp, "userint4",   0, wi);
    ir->userreal1 = get_ereal(&inp, "userreal1",  0, wi);
    ir->userreal2 = get_ereal(&inp, "userreal2",  0, wi);
    ir->userreal3 = get_ereal(&inp, "userreal3",  0, wi);
    ir->userreal4 = get_ereal(&inp, "userreal4",  0, wi);
#undef CTYPE

    {
        gmx::TextOutputFile stream(mdparout);
        write_inpfile(&stream, mdparout, &inp, FALSE, writeMdpHeader, wi);

        // Transform module data into a flat key-value tree for output.
        gmx::KeyValueTreeBuilder       builder;
        gmx::KeyValueTreeObjectBuilder builderObject = builder.rootObject();
        mdModules->buildMdpOutput(&builderObject);
        {
            gmx::TextWriter writer(&stream);
            writeKeyValueTreeAsMdp(&writer, builder.build());
        }
        stream.close();
    }

    /* Process options if necessary */
    for (m = 0; m < 2; m++)
    {
        for (i = 0; i < 2*DIM; i++)
        {
            dumdub[m][i] = 0.0;
        }
        if (ir->epc)
        {
            switch (ir->epct)
            {
                case epctISOTROPIC:
                    if (sscanf(dumstr[m], "%lf", &(dumdub[m][XX])) != 1)
                    {
                        warning_error(wi, "Pressure coupling incorrect number of values (I need exactly 1)");
                    }
                    dumdub[m][YY] = dumdub[m][ZZ] = dumdub[m][XX];
                    break;
                case epctSEMIISOTROPIC:
                case epctSURFACETENSION:
                    if (sscanf(dumstr[m], "%lf%lf", &(dumdub[m][XX]), &(dumdub[m][ZZ])) != 2)
                    {
                        warning_error(wi, "Pressure coupling incorrect number of values (I need exactly 2)");
                    }
                    dumdub[m][YY] = dumdub[m][XX];
                    break;
                case epctANISOTROPIC:
                    if (sscanf(dumstr[m], "%lf%lf%lf%lf%lf%lf",
                               &(dumdub[m][XX]), &(dumdub[m][YY]), &(dumdub[m][ZZ]),
                               &(dumdub[m][3]), &(dumdub[m][4]), &(dumdub[m][5])) != 6)
                    {
                        warning_error(wi, "Pressure coupling incorrect number of values (I need exactly 6)");
                    }
                    break;
                default:
                    gmx_fatal(FARGS, "Pressure coupling type %s not implemented yet",
                              epcoupltype_names[ir->epct]);
            }
        }
    }
    clear_mat(ir->ref_p);
    clear_mat(ir->compress);
    for (i = 0; i < DIM; i++)
    {
        ir->ref_p[i][i]    = dumdub[1][i];
        ir->compress[i][i] = dumdub[0][i];
    }
    if (ir->epct == epctANISOTROPIC)
    {
        ir->ref_p[XX][YY] = dumdub[1][3];
        ir->ref_p[XX][ZZ] = dumdub[1][4];
        ir->ref_p[YY][ZZ] = dumdub[1][5];
        if (ir->ref_p[XX][YY] != 0 && ir->ref_p[XX][ZZ] != 0 && ir->ref_p[YY][ZZ] != 0)
        {
            warning(wi, "All off-diagonal reference pressures are non-zero. Are you sure you want to apply a threefold shear stress?\n");
        }
        ir->compress[XX][YY] = dumdub[0][3];
        ir->compress[XX][ZZ] = dumdub[0][4];
        ir->compress[YY][ZZ] = dumdub[0][5];
        for (i = 0; i < DIM; i++)
        {
            for (m = 0; m < i; m++)
            {
                ir->ref_p[i][m]    = ir->ref_p[m][i];
                ir->compress[i][m] = ir->compress[m][i];
            }
        }
    }

    if (ir->comm_mode == ecmNO)
    {
        ir->nstcomm = 0;
    }

    opts->couple_moltype = nullptr;
    if (strlen(is->couple_moltype) > 0)
    {
        if (ir->efep != efepNO)
        {
            opts->couple_moltype = gmx_strdup(is->couple_moltype);
            if (opts->couple_lam0 == opts->couple_lam1)
            {
                warning(wi, "The lambda=0 and lambda=1 states for coupling are identical");
            }
            if (ir->eI == eiMD && (opts->couple_lam0 == ecouplamNONE ||
                                   opts->couple_lam1 == ecouplamNONE))
            {
                warning(wi, "For proper sampling of the (nearly) decoupled state, stochastic dynamics should be used");
            }
        }
        else
        {
            warning_note(wi, "Free energy is turned off, so we will not decouple the molecule listed in your input.");
        }
    }
    /* FREE ENERGY AND EXPANDED ENSEMBLE OPTIONS */
    if (ir->efep != efepNO)
    {
        if (fep->delta_lambda > 0)
        {
            ir->efep = efepSLOWGROWTH;
        }
    }

    if (fep->edHdLPrintEnergy == edHdLPrintEnergyYES)
    {
        fep->edHdLPrintEnergy = edHdLPrintEnergyTOTAL;
        warning_note(wi, "Old option for dhdl-print-energy given: "
                     "changing \"yes\" to \"total\"\n");
    }

    if (ir->bSimTemp && (fep->edHdLPrintEnergy == edHdLPrintEnergyNO))
    {
        /* always print out the energy to dhdl if we are doing
           expanded ensemble, since we need the total energy for
           analysis if the temperature is changing. In some
           conditions one may only want the potential energy, so
           we will allow that if the appropriate mdp setting has
           been enabled. Otherwise, total it is:
         */
        fep->edHdLPrintEnergy = edHdLPrintEnergyTOTAL;
    }

    if ((ir->efep != efepNO) || ir->bSimTemp)
    {
        ir->bExpanded = FALSE;
        if ((ir->efep == efepEXPANDED) || ir->bSimTemp)
        {
            ir->bExpanded = TRUE;
        }
        do_fep_params(ir, is->fep_lambda, is->lambda_weights, wi);
        if (ir->bSimTemp) /* done after fep params */
        {
            do_simtemp_params(ir);
        }

        /* Because sc-coul (=FALSE by default) only acts on the lambda state
         * setup and not on the old way of specifying the free-energy setup,
         * we should check for using soft-core when not needed, since that
         * can complicate the sampling significantly.
         * Note that we only check for the automated coupling setup.
         * If the (advanced) user does FEP through manual topology changes,
         * this check will not be triggered.
         */
        if (ir->efep != efepNO && ir->fepvals->n_lambda == 0 &&
            ir->fepvals->sc_alpha != 0 &&
            (couple_lambda_has_vdw_on(opts->couple_lam0) &&
             couple_lambda_has_vdw_on(opts->couple_lam1)))
        {
            warning(wi, "You are using soft-core interactions while the Van der Waals interactions are not decoupled (note that the sc-coul option is only active when using lambda states). Although this will not lead to errors, you will need much more sampling than without soft-core interactions. Consider using sc-alpha=0.");
        }
    }
    else
    {
        ir->fepvals->n_lambda = 0;
    }

    /* WALL PARAMETERS */

    do_wall_params(ir, is->wall_atomtype, is->wall_density, opts, wi);

    /* ORIENTATION RESTRAINT PARAMETERS */

    if (opts->bOrire && gmx::splitString(is->orirefitgrp).size() != 1)
    {
        warning_error(wi, "ERROR: Need one orientation restraint fit group\n");
    }

    /* DEFORMATION PARAMETERS */

    clear_mat(ir->deform);
    for (i = 0; i < 6; i++)
    {
        dumdub[0][i] = 0;
    }

    double gmx_unused canary;
    int               ndeform = sscanf(is->deform, "%lf %lf %lf %lf %lf %lf %lf",
                                       &(dumdub[0][0]), &(dumdub[0][1]), &(dumdub[0][2]),
                                       &(dumdub[0][3]), &(dumdub[0][4]), &(dumdub[0][5]), &canary);

    if (strlen(is->deform) > 0 && ndeform != 6)
    {
        warning_error(wi, gmx::formatString("Cannot parse exactly 6 box deformation velocities from string '%s'", is->deform).c_str());
    }
    for (i = 0; i < 3; i++)
    {
        ir->deform[i][i] = dumdub[0][i];
    }
    ir->deform[YY][XX] = dumdub[0][3];
    ir->deform[ZZ][XX] = dumdub[0][4];
    ir->deform[ZZ][YY] = dumdub[0][5];
    if (ir->epc != epcNO)
    {
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j <= i; j++)
            {
                if (ir->deform[i][j] != 0 && ir->compress[i][j] != 0)
                {
                    warning_error(wi, "A box element has deform set and compressibility > 0");
                }
            }
        }
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < i; j++)
            {
                if (ir->deform[i][j] != 0)
                {
                    for (m = j; m < DIM; m++)
                    {
                        if (ir->compress[m][j] != 0)
                        {
                            sprintf(warn_buf, "An off-diagonal box element has deform set while compressibility > 0 for the same component of another box vector, this might lead to spurious periodicity effects.");
                            warning(wi, warn_buf);
                        }
                    }
                }
            }
        }
    }

    /* Ion/water position swapping checks */
    if (ir->eSwapCoords != eswapNO)
    {
        if (ir->swap->nstswap < 1)
        {
            warning_error(wi, "swap_frequency must be 1 or larger when ion swapping is requested");
        }
        if (ir->swap->nAverage < 1)
        {
            warning_error(wi, "coupl_steps must be 1 or larger.\n");
        }
        if (ir->swap->threshold < 1.0)
        {
            warning_error(wi, "Ion count threshold must be at least 1.\n");
        }
    }

    sfree(dumstr[0]);
    sfree(dumstr[1]);
}

static int search_QMstring(const char *s, int ng, const char *gn[])
{
    /* same as normal search_string, but this one searches QM strings */
    int i;

    for (i = 0; (i < ng); i++)
    {
        if (gmx_strcasecmp(s, gn[i]) == 0)
        {
            return i;
        }
    }

    gmx_fatal(FARGS, "this QM method or basisset (%s) is not implemented\n!", s);
} /* search_QMstring */

/* We would like gn to be const as well, but C doesn't allow this */
/* TODO this is utility functionality (search for the index of a
   string in a collection), so should be refactored and located more
   centrally. */
int search_string(const char *s, int ng, char *gn[])
{
    int i;

    for (i = 0; (i < ng); i++)
    {
        if (gmx_strcasecmp(s, gn[i]) == 0)
        {
            return i;
        }
    }

    gmx_fatal(FARGS,
              "Group %s referenced in the .mdp file was not found in the index file.\n"
              "Group names must match either [moleculetype] names or custom index group\n"
              "names, in which case you must supply an index file to the '-n' option\n"
              "of grompp.",
              s);
}

static bool do_numbering(int natoms, gmx_groups_t *groups,
                         gmx::ArrayRef<std::string> groupsFromMdpFile,
                         t_blocka *block, char *gnames[],
                         int gtype, int restnm,
                         int grptp, bool bVerbose,
                         warninp_t wi)
{
    unsigned short *cbuf;
    t_grps         *grps = &(groups->grps[gtype]);
    int             j, gid, aj, ognr, ntot = 0;
    const char     *title;
    bool            bRest;
    char            warn_buf[STRLEN];

    title = gtypes[gtype];

    snew(cbuf, natoms);
    /* Mark all id's as not set */
    for (int i = 0; (i < natoms); i++)
    {
        cbuf[i] = NOGID;
    }

    snew(grps->nm_ind, groupsFromMdpFile.size()+1); /* +1 for possible rest group */
    for (int i = 0; i != groupsFromMdpFile.size(); ++i)
    {
        /* Lookup the group name in the block structure */
        gid = search_string(groupsFromMdpFile[i].c_str(), block->nr, gnames);
        if ((grptp != egrptpONE) || (i == 0))
        {
            grps->nm_ind[grps->nr++] = gid;
        }

        /* Now go over the atoms in the group */
        for (j = block->index[gid]; (j < block->index[gid+1]); j++)
        {

            aj = block->a[j];

            /* Range checking */
            if ((aj < 0) || (aj >= natoms))
            {
                gmx_fatal(FARGS, "Invalid atom number %d in indexfile", aj);
            }
            /* Lookup up the old group number */
            ognr = cbuf[aj];
            if (ognr != NOGID)
            {
                gmx_fatal(FARGS, "Atom %d in multiple %s groups (%d and %d)",
                          aj+1, title, ognr+1, i+1);
            }
            else
            {
                /* Store the group number in buffer */
                if (grptp == egrptpONE)
                {
                    cbuf[aj] = 0;
                }
                else
                {
                    cbuf[aj] = i;
                }
                ntot++;
            }
        }
    }

    /* Now check whether we have done all atoms */
    bRest = FALSE;
    if (ntot != natoms)
    {
        if (grptp == egrptpALL)
        {
            gmx_fatal(FARGS, "%d atoms are not part of any of the %s groups",
                      natoms-ntot, title);
        }
        else if (grptp == egrptpPART)
        {
            sprintf(warn_buf, "%d atoms are not part of any of the %s groups",
                    natoms-ntot, title);
            warning_note(wi, warn_buf);
        }
        /* Assign all atoms currently unassigned to a rest group */
        for (j = 0; (j < natoms); j++)
        {
            if (cbuf[j] == NOGID)
            {
                cbuf[j] = grps->nr;
                bRest   = TRUE;
            }
        }
        if (grptp != egrptpPART)
        {
            if (bVerbose)
            {
                fprintf(stderr,
                        "Making dummy/rest group for %s containing %d elements\n",
                        title, natoms-ntot);
            }
            /* Add group name "rest" */
            grps->nm_ind[grps->nr] = restnm;

            /* Assign the rest name to all atoms not currently assigned to a group */
            for (j = 0; (j < natoms); j++)
            {
                if (cbuf[j] == NOGID)
                {
                    cbuf[j] = grps->nr;
                }
            }
            grps->nr++;
        }
    }

    if (grps->nr == 1 && (ntot == 0 || ntot == natoms))
    {
        /* All atoms are part of one (or no) group, no index required */
        groups->ngrpnr[gtype] = 0;
        groups->grpnr[gtype]  = nullptr;
    }
    else
    {
        groups->ngrpnr[gtype] = natoms;
        snew(groups->grpnr[gtype], natoms);
        for (j = 0; (j < natoms); j++)
        {
            groups->grpnr[gtype][j] = cbuf[j];
        }
    }

    sfree(cbuf);

    return (bRest && grptp == egrptpPART);
}

static void calc_nrdf(const gmx_mtop_t *mtop, t_inputrec *ir, char **gnames)
{
    t_grpopts              *opts;
    const gmx_groups_t     *groups;
    pull_params_t          *pull;
    int                     natoms, ai, aj, i, j, d, g, imin, jmin;
    int                    *nrdf2, *na_vcm, na_tot;
    double                 *nrdf_tc, *nrdf_vcm, nrdf_uc, *nrdf_vcm_sub;
    ivec                   *dof_vcm;
    gmx_mtop_atomloop_all_t aloop;
    int                     mol, ftype, as;

    /* Calculate nrdf.
     * First calc 3xnr-atoms for each group
     * then subtract half a degree of freedom for each constraint
     *
     * Only atoms and nuclei contribute to the degrees of freedom...
     */

    opts = &ir->opts;

    groups = &mtop->groups;
    natoms = mtop->natoms;

    /* Allocate one more for a possible rest group */
    /* We need to sum degrees of freedom into doubles,
     * since floats give too low nrdf's above 3 million atoms.
     */
    snew(nrdf_tc, groups->grps[egcTC].nr+1);
    snew(nrdf_vcm, groups->grps[egcVCM].nr+1);
    snew(dof_vcm, groups->grps[egcVCM].nr+1);
    snew(na_vcm, groups->grps[egcVCM].nr+1);
    snew(nrdf_vcm_sub, groups->grps[egcVCM].nr+1);

    for (i = 0; i < groups->grps[egcTC].nr; i++)
    {
        nrdf_tc[i] = 0;
    }
    for (i = 0; i < groups->grps[egcVCM].nr+1; i++)
    {
        nrdf_vcm[i]     = 0;
        clear_ivec(dof_vcm[i]);
        na_vcm[i]       = 0;
        nrdf_vcm_sub[i] = 0;
    }

    snew(nrdf2, natoms);
    aloop = gmx_mtop_atomloop_all_init(mtop);
    const t_atom *atom;
    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
        nrdf2[i] = 0;
        if (atom->ptype == eptAtom || atom->ptype == eptNucleus)
        {
            g = getGroupType(groups, egcFREEZE, i);
            for (d = 0; d < DIM; d++)
            {
                if (opts->nFreeze[g][d] == 0)
                {
                    /* Add one DOF for particle i (counted as 2*1) */
                    nrdf2[i]                              += 2;
                    /* VCM group i has dim d as a DOF */
                    dof_vcm[getGroupType(groups, egcVCM, i)][d]  = 1;
                }
            }
            nrdf_tc [getGroupType(groups, egcTC, i)]  += 0.5*nrdf2[i];
            nrdf_vcm[getGroupType(groups, egcVCM, i)] += 0.5*nrdf2[i];
        }
    }

    as = 0;
    for (const gmx_molblock_t &molb : mtop->molblock)
    {
        const gmx_moltype_t &molt = mtop->moltype[molb.type];
        atom = molt.atoms.atom;
        for (mol = 0; mol < molb.nmol; mol++)
        {
            for (ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
            {
                gmx::ArrayRef<const int> ia = molt.ilist[ftype].iatoms;
                for (i = 0; i < molt.ilist[ftype].size(); )
                {
                    /* Subtract degrees of freedom for the constraints,
                     * if the particles still have degrees of freedom left.
                     * If one of the particles is a vsite or a shell, then all
                     * constraint motion will go there, but since they do not
                     * contribute to the constraints the degrees of freedom do not
                     * change.
                     */
                    ai = as + ia[i + 1];
                    aj = as + ia[i + 2];
                    if (((atom[ia[i + 1]].ptype == eptNucleus) ||
                         (atom[ia[i + 1]].ptype == eptAtom)) &&
                        ((atom[ia[i + 2]].ptype == eptNucleus) ||
                         (atom[ia[i + 2]].ptype == eptAtom)))
                    {
                        if (nrdf2[ai] > 0)
                        {
                            jmin = 1;
                        }
                        else
                        {
                            jmin = 2;
                        }
                        if (nrdf2[aj] > 0)
                        {
                            imin = 1;
                        }
                        else
                        {
                            imin = 2;
                        }
                        imin       = std::min(imin, nrdf2[ai]);
                        jmin       = std::min(jmin, nrdf2[aj]);
                        nrdf2[ai] -= imin;
                        nrdf2[aj] -= jmin;
                        nrdf_tc [getGroupType(groups, egcTC, ai)]  -= 0.5*imin;
                        nrdf_tc [getGroupType(groups, egcTC, aj)]  -= 0.5*jmin;
                        nrdf_vcm[getGroupType(groups, egcVCM, ai)] -= 0.5*imin;
                        nrdf_vcm[getGroupType(groups, egcVCM, aj)] -= 0.5*jmin;
                    }
                    i  += interaction_function[ftype].nratoms+1;
                }
            }
            gmx::ArrayRef<const int> ia = molt.ilist[F_SETTLE].iatoms;
            for (i = 0; i < molt.ilist[F_SETTLE].size(); )
            {
                /* Subtract 1 dof from every atom in the SETTLE */
                for (j = 0; j < 3; j++)
                {
                    ai         = as + ia[i + 1 + j];
                    imin       = std::min(2, nrdf2[ai]);
                    nrdf2[ai] -= imin;
                    nrdf_tc [getGroupType(groups, egcTC, ai)]  -= 0.5*imin;
                    nrdf_vcm[getGroupType(groups, egcVCM, ai)] -= 0.5*imin;
                }
                i  += 4;
            }
            as += molt.atoms.nr;
        }
    }

    if (ir->bPull)
    {
        /* Correct nrdf for the COM constraints.
         * We correct using the TC and VCM group of the first atom
         * in the reference and pull group. If atoms in one pull group
         * belong to different TC or VCM groups it is anyhow difficult
         * to determine the optimal nrdf assignment.
         */
        pull = ir->pull;

        for (i = 0; i < pull->ncoord; i++)
        {
            if (pull->coord[i].eType != epullCONSTRAINT)
            {
                continue;
            }

            imin = 1;

            for (j = 0; j < 2; j++)
            {
                const t_pull_group *pgrp;

                pgrp = &pull->group[pull->coord[i].group[j]];

                if (pgrp->nat > 0)
                {
                    /* Subtract 1/2 dof from each group */
                    ai = pgrp->ind[0];
                    nrdf_tc [getGroupType(groups, egcTC, ai)]  -= 0.5*imin;
                    nrdf_vcm[getGroupType(groups, egcVCM, ai)] -= 0.5*imin;
                    if (nrdf_tc[getGroupType(groups, egcTC, ai)] < 0)
                    {
                        gmx_fatal(FARGS, "Center of mass pulling constraints caused the number of degrees of freedom for temperature coupling group %s to be negative", gnames[groups->grps[egcTC].nm_ind[getGroupType(groups, egcTC, ai)]]);
                    }
                }
                else
                {
                    /* We need to subtract the whole DOF from group j=1 */
                    imin += 1;
                }
            }
        }
    }

    if (ir->nstcomm != 0)
    {
        int ndim_rm_vcm;

        /* We remove COM motion up to dim ndof_com() */
        ndim_rm_vcm = ndof_com(ir);

        /* Subtract ndim_rm_vcm (or less with frozen dimensions) from
         * the number of degrees of freedom in each vcm group when COM
         * translation is removed and 6 when rotation is removed as well.
         */
        for (j = 0; j < groups->grps[egcVCM].nr+1; j++)
        {
            switch (ir->comm_mode)
            {
                case ecmLINEAR:
                case ecmLINEAR_ACCELERATION_CORRECTION:
                    nrdf_vcm_sub[j] = 0;
                    for (d = 0; d < ndim_rm_vcm; d++)
                    {
                        if (dof_vcm[j][d])
                        {
                            nrdf_vcm_sub[j]++;
                        }
                    }
                    break;
                case ecmANGULAR:
                    nrdf_vcm_sub[j] = 6;
                    break;
                default:
                    gmx_incons("Checking comm_mode");
            }
        }

        for (i = 0; i < groups->grps[egcTC].nr; i++)
        {
            /* Count the number of atoms of TC group i for every VCM group */
            for (j = 0; j < groups->grps[egcVCM].nr+1; j++)
            {
                na_vcm[j] = 0;
            }
            na_tot = 0;
            for (ai = 0; ai < natoms; ai++)
            {
                if (getGroupType(groups, egcTC, ai) == i)
                {
                    na_vcm[getGroupType(groups, egcVCM, ai)]++;
                    na_tot++;
                }
            }
            /* Correct for VCM removal according to the fraction of each VCM
             * group present in this TC group.
             */
            nrdf_uc    = nrdf_tc[i];
            nrdf_tc[i] = 0;
            for (j = 0; j < groups->grps[egcVCM].nr+1; j++)
            {
                if (nrdf_vcm[j] > nrdf_vcm_sub[j])
                {
                    nrdf_tc[i] += nrdf_uc*(static_cast<double>(na_vcm[j])/static_cast<double>(na_tot))*
                        (nrdf_vcm[j] - nrdf_vcm_sub[j])/nrdf_vcm[j];
                }
            }
        }
    }
    for (i = 0; (i < groups->grps[egcTC].nr); i++)
    {
        opts->nrdf[i] = nrdf_tc[i];
        if (opts->nrdf[i] < 0)
        {
            opts->nrdf[i] = 0;
        }
        fprintf(stderr,
                "Number of degrees of freedom in T-Coupling group %s is %.2f\n",
                gnames[groups->grps[egcTC].nm_ind[i]], opts->nrdf[i]);
    }

    sfree(nrdf2);
    sfree(nrdf_tc);
    sfree(nrdf_vcm);
    sfree(dof_vcm);
    sfree(na_vcm);
    sfree(nrdf_vcm_sub);
}

static bool do_egp_flag(t_inputrec *ir, gmx_groups_t *groups,
                        const char *option, const char *val, int flag)
{
    /* The maximum number of energy group pairs would be MAXPTR*(MAXPTR+1)/2.
     * But since this is much larger than STRLEN, such a line can not be parsed.
     * The real maximum is the number of names that fit in a string: STRLEN/2.
     */
#define EGP_MAX (STRLEN/2)
    int      j, k, nr;
    char  ***gnames;
    bool     bSet;

    gnames = groups->grpname;

    auto names = gmx::splitString(val);
    if (names.size() % 2 != 0)
    {
        gmx_fatal(FARGS, "The number of groups for %s is odd", option);
    }
    nr   = groups->grps[egcENER].nr;
    bSet = FALSE;
    for (size_t i = 0; i < names.size() / 2; i++)
    {
        j = 0;
        while ((j < nr) &&
               gmx_strcasecmp(names[2*i].c_str(), *(gnames[groups->grps[egcENER].nm_ind[j]])))
        {
            j++;
        }
        if (j == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n",
                      names[2*i].c_str(), option);
        }
        k = 0;
        while ((k < nr) &&
               gmx_strcasecmp(names[2*i+1].c_str(), *(gnames[groups->grps[egcENER].nm_ind[k]])))
        {
            k++;
        }
        if (k == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n",
                      names[2*i+1].c_str(), option);
        }
        if ((j < nr) && (k < nr))
        {
            ir->opts.egp_flags[nr*j+k] |= flag;
            ir->opts.egp_flags[nr*k+j] |= flag;
            bSet = TRUE;
        }
    }

    return bSet;
}


static void make_swap_groups(
        t_swapcoords  *swap,
        t_blocka      *grps,
        char         **gnames)
{
    int          ig = -1, i = 0, gind;
    t_swapGroup *swapg;


    /* Just a quick check here, more thorough checks are in mdrun */
    if (strcmp(swap->grp[eGrpSplit0].molname, swap->grp[eGrpSplit1].molname) == 0)
    {
        gmx_fatal(FARGS, "The split groups can not both be '%s'.", swap->grp[eGrpSplit0].molname);
    }

    /* Get the index atoms of the split0, split1, solvent, and swap groups */
    for (ig = 0; ig < swap->ngrp; ig++)
    {
        swapg      = &swap->grp[ig];
        gind       = search_string(swap->grp[ig].molname, grps->nr, gnames);
        swapg->nat = grps->index[gind+1] - grps->index[gind];

        if (swapg->nat > 0)
        {
            fprintf(stderr, "%s group '%s' contains %d atoms.\n",
                    ig < 3 ? eSwapFixedGrp_names[ig] : "Swap",
                    swap->grp[ig].molname, swapg->nat);
            snew(swapg->ind, swapg->nat);
            for (i = 0; i < swapg->nat; i++)
            {
                swapg->ind[i] = grps->a[grps->index[gind]+i];
            }
        }
        else
        {
            gmx_fatal(FARGS, "Swap group %s does not contain any atoms.", swap->grp[ig].molname);
        }
    }
}


static void make_IMD_group(t_IMD *IMDgroup, char *IMDgname, t_blocka *grps, char **gnames)
{
    int      ig, i;


    ig            = search_string(IMDgname, grps->nr, gnames);
    IMDgroup->nat = grps->index[ig+1] - grps->index[ig];

    if (IMDgroup->nat > 0)
    {
        fprintf(stderr, "Group '%s' with %d atoms can be activated for interactive molecular dynamics (IMD).\n",
                IMDgname, IMDgroup->nat);
        snew(IMDgroup->ind, IMDgroup->nat);
        for (i = 0; i < IMDgroup->nat; i++)
        {
            IMDgroup->ind[i] = grps->a[grps->index[ig]+i];
        }
    }
}

void do_index(const char* mdparin, const char *ndx,
              gmx_mtop_t *mtop,
              bool bVerbose,
              t_inputrec *ir,
              warninp_t wi)
{
    t_blocka     *grps;
    gmx_groups_t *groups;
    int           natoms;
    t_symtab     *symtab;
    t_atoms       atoms_all;
    char          warnbuf[STRLEN], **gnames;
    int           nr;
    real          tau_min;
    int           nstcmin;
    int           i, j, k, restnm;
    bool          bExcl, bTable, bAnneal, bRest;
    char          warn_buf[STRLEN];

    if (bVerbose)
    {
        fprintf(stderr, "processing index file...\n");
    }
    if (ndx == nullptr)
    {
        snew(grps, 1);
        snew(grps->index, 1);
        snew(gnames, 1);
        atoms_all = gmx_mtop_global_atoms(mtop);
        analyse(&atoms_all, grps, &gnames, FALSE, TRUE);
        done_atom(&atoms_all);
    }
    else
    {
        grps = init_index(ndx, &gnames);
    }

    groups = &mtop->groups;
    natoms = mtop->natoms;
    symtab = &mtop->symtab;

    snew(groups->grpname, grps->nr+1);

    for (i = 0; (i < grps->nr); i++)
    {
        groups->grpname[i] = put_symtab(symtab, gnames[i]);
    }
    groups->grpname[i] = put_symtab(symtab, "rest");
    restnm             = i;
    srenew(gnames, grps->nr+1);
    gnames[restnm]   = *(groups->grpname[i]);
    groups->ngrpname = grps->nr+1;

    set_warning_line(wi, mdparin, -1);

    auto temperatureCouplingTauValues       = gmx::splitString(is->tau_t);
    auto temperatureCouplingReferenceValues = gmx::splitString(is->ref_t);
    auto temperatureCouplingGroupNames      = gmx::splitString(is->tcgrps);
    if (temperatureCouplingTauValues.size() != temperatureCouplingGroupNames.size() ||
        temperatureCouplingReferenceValues.size() != temperatureCouplingGroupNames.size())
    {
        gmx_fatal(FARGS, "Invalid T coupling input: %zu groups, %zu ref-t values and "
                  "%zu tau-t values",
                  temperatureCouplingGroupNames.size(),
                  temperatureCouplingReferenceValues.size(),
                  temperatureCouplingTauValues.size());
    }

    const bool useReferenceTemperature = integratorHasReferenceTemperature(ir);
    do_numbering(natoms, groups, temperatureCouplingGroupNames, grps, gnames, egcTC,
                 restnm, useReferenceTemperature ? egrptpALL : egrptpALL_GENREST, bVerbose, wi);
    nr            = groups->grps[egcTC].nr;
    ir->opts.ngtc = nr;
    snew(ir->opts.nrdf, nr);
    snew(ir->opts.tau_t, nr);
    snew(ir->opts.ref_t, nr);
    if (ir->eI == eiBD && ir->bd_fric == 0)
    {
        fprintf(stderr, "bd-fric=0, so tau-t will be used as the inverse friction constant(s)\n");
    }

    if (useReferenceTemperature)
    {
        if (size_t(nr) != temperatureCouplingReferenceValues.size())
        {
            gmx_fatal(FARGS, "Not enough ref-t and tau-t values!");
        }

        tau_min = 1e20;
        convertReals(wi, temperatureCouplingTauValues, "tau-t", ir->opts.tau_t);
        for (i = 0; (i < nr); i++)
        {
            if ((ir->eI == eiBD) && ir->opts.tau_t[i] <= 0)
            {
                sprintf(warn_buf, "With integrator %s tau-t should be larger than 0", ei_names[ir->eI]);
                warning_error(wi, warn_buf);
            }

            if (ir->etc != etcVRESCALE && ir->opts.tau_t[i] == 0)
            {
                warning_note(wi, "tau-t = -1 is the value to signal that a group should not have temperature coupling. Treating your use of tau-t = 0 as if you used -1.");
            }

            if (ir->opts.tau_t[i] >= 0)
            {
                tau_min = std::min(tau_min, ir->opts.tau_t[i]);
            }
        }
        if (ir->etc != etcNO && ir->nsttcouple == -1)
        {
            ir->nsttcouple = ir_optimal_nsttcouple(ir);
        }

        if (EI_VV(ir->eI))
        {
            if ((ir->etc == etcNOSEHOOVER) && (ir->epc == epcBERENDSEN))
            {
                gmx_fatal(FARGS, "Cannot do Nose-Hoover temperature with Berendsen pressure control with md-vv; use either vrescale temperature with berendsen pressure or Nose-Hoover temperature with MTTK pressure");
            }
            if (ir->epc == epcMTTK)
            {
                if (ir->etc != etcNOSEHOOVER)
                {
                    gmx_fatal(FARGS, "Cannot do MTTK pressure coupling without Nose-Hoover temperature control");
                }
                else
                {
                    if (ir->nstpcouple != ir->nsttcouple)
                    {
                        int mincouple = std::min(ir->nstpcouple, ir->nsttcouple);
                        ir->nstpcouple = ir->nsttcouple = mincouple;
                        sprintf(warn_buf, "for current Trotter decomposition methods with vv, nsttcouple and nstpcouple must be equal.  Both have been reset to min(nsttcouple,nstpcouple) = %d", mincouple);
                        warning_note(wi, warn_buf);
                    }
                }
            }
        }
        /* velocity verlet with averaged kinetic energy KE = 0.5*(v(t+1/2) - v(t-1/2)) is implemented
           primarily for testing purposes, and does not work with temperature coupling other than 1 */

        if (ETC_ANDERSEN(ir->etc))
        {
            if (ir->nsttcouple != 1)
            {
                ir->nsttcouple = 1;
                sprintf(warn_buf, "Andersen temperature control methods assume nsttcouple = 1; there is no need for larger nsttcouple > 1, since no global parameters are computed. nsttcouple has been reset to 1");
                warning_note(wi, warn_buf);
            }
        }
        nstcmin = tcouple_min_integration_steps(ir->etc);
        if (nstcmin > 1)
        {
            if (tau_min/(ir->delta_t*ir->nsttcouple) < nstcmin - 10*GMX_REAL_EPS)
            {
                sprintf(warn_buf, "For proper integration of the %s thermostat, tau-t (%g) should be at least %d times larger than nsttcouple*dt (%g)",
                        ETCOUPLTYPE(ir->etc),
                        tau_min, nstcmin,
                        ir->nsttcouple*ir->delta_t);
                warning(wi, warn_buf);
            }
        }
        convertReals(wi, temperatureCouplingReferenceValues, "ref-t", ir->opts.ref_t);
        for (i = 0; (i < nr); i++)
        {
            if (ir->opts.ref_t[i] < 0)
            {
                gmx_fatal(FARGS, "ref-t for group %d negative", i);
            }
        }
        /* set the lambda mc temperature to the md integrator temperature (which should be defined
           if we are in this conditional) if mc_temp is negative */
        if (ir->expandedvals->mc_temp < 0)
        {
            ir->expandedvals->mc_temp = ir->opts.ref_t[0]; /*for now, set to the first reft */
        }
    }

    /* Simulated annealing for each group. There are nr groups */
    auto simulatedAnnealingGroupNames = gmx::splitString(is->anneal);
    if (simulatedAnnealingGroupNames.size() == 1 &&
        gmx_strncasecmp(simulatedAnnealingGroupNames[0].c_str(), "N", 1) == 0)
    {
        simulatedAnnealingGroupNames.resize(0);
    }
    if (!simulatedAnnealingGroupNames.empty() &&
        simulatedAnnealingGroupNames.size() != size_t(nr))
    {
        gmx_fatal(FARGS, "Not enough annealing values: %zu (for %d groups)\n",
                  simulatedAnnealingGroupNames.size(), nr);
    }
    else
    {
        snew(ir->opts.annealing, nr);
        snew(ir->opts.anneal_npoints, nr);
        snew(ir->opts.anneal_time, nr);
        snew(ir->opts.anneal_temp, nr);
        for (i = 0; i < nr; i++)
        {
            ir->opts.annealing[i]      = eannNO;
            ir->opts.anneal_npoints[i] = 0;
            ir->opts.anneal_time[i]    = nullptr;
            ir->opts.anneal_temp[i]    = nullptr;
        }
        if (!simulatedAnnealingGroupNames.empty())
        {
            bAnneal = FALSE;
            for (i = 0; i < nr; i++)
            {
                if (gmx_strncasecmp(simulatedAnnealingGroupNames[i].c_str(), "N", 1) == 0)
                {
                    ir->opts.annealing[i] = eannNO;
                }
                else if (gmx_strncasecmp(simulatedAnnealingGroupNames[i].c_str(), "S", 1) == 0)
                {
                    ir->opts.annealing[i] = eannSINGLE;
                    bAnneal               = TRUE;
                }
                else if (gmx_strncasecmp(simulatedAnnealingGroupNames[i].c_str(), "P", 1) == 0)
                {
                    ir->opts.annealing[i] = eannPERIODIC;
                    bAnneal               = TRUE;
                }
            }
            if (bAnneal)
            {
                /* Read the other fields too */
                auto simulatedAnnealingPoints = gmx::splitString(is->anneal_npoints);
                if (simulatedAnnealingPoints.size() != simulatedAnnealingGroupNames.size())
                {
                    gmx_fatal(FARGS, "Found %zu annealing-npoints values for %zu groups\n",
                              simulatedAnnealingPoints.size(), simulatedAnnealingGroupNames.size());
                }
                convertInts(wi, simulatedAnnealingPoints, "annealing points", ir->opts.anneal_npoints);
                for (k = 0, i = 0; i < nr; i++)
                {
                    if (ir->opts.anneal_npoints[i] == 1)
                    {
                        gmx_fatal(FARGS, "Please specify at least a start and an end point for annealing\n");
                    }
                    snew(ir->opts.anneal_time[i], ir->opts.anneal_npoints[i]);
                    snew(ir->opts.anneal_temp[i], ir->opts.anneal_npoints[i]);
                    k += ir->opts.anneal_npoints[i];
                }

                auto simulatedAnnealingTimes = gmx::splitString(is->anneal_time);
                if (simulatedAnnealingTimes.size() != size_t(k))
                {
                    gmx_fatal(FARGS, "Found %zu annealing-time values, wanted %d\n",
                              simulatedAnnealingTimes.size(), k);
                }
                auto simulatedAnnealingTemperatures = gmx::splitString(is->anneal_temp);
                if (simulatedAnnealingTemperatures.size() != size_t(k))
                {
                    gmx_fatal(FARGS, "Found %zu annealing-temp values, wanted %d\n",
                              simulatedAnnealingTemperatures.size(), k);
                }

                convertReals(wi, simulatedAnnealingTimes, "anneal-time", ir->opts.anneal_time[i]);
                convertReals(wi, simulatedAnnealingTemperatures, "anneal-temp", ir->opts.anneal_temp[i]);
                for (i = 0, k = 0; i < nr; i++)
                {
                    for (j = 0; j < ir->opts.anneal_npoints[i]; j++)
                    {
                        if (j == 0)
                        {
                            if (ir->opts.anneal_time[i][0] > (ir->init_t+GMX_REAL_EPS))
                            {
                                gmx_fatal(FARGS, "First time point for annealing > init_t.\n");
                            }
                        }
                        else
                        {
                            /* j>0 */
                            if (ir->opts.anneal_time[i][j] < ir->opts.anneal_time[i][j-1])
                            {
                                gmx_fatal(FARGS, "Annealing timepoints out of order: t=%f comes after t=%f\n",
                                          ir->opts.anneal_time[i][j], ir->opts.anneal_time[i][j-1]);
                            }
                        }
                        if (ir->opts.anneal_temp[i][j] < 0)
                        {
                            gmx_fatal(FARGS, "Found negative temperature in annealing: %f\n", ir->opts.anneal_temp[i][j]);
                        }
                        k++;
                    }
                }
                /* Print out some summary information, to make sure we got it right */
                for (i = 0, k = 0; i < nr; i++)
                {
                    if (ir->opts.annealing[i] != eannNO)
                    {
                        j = groups->grps[egcTC].nm_ind[i];
                        fprintf(stderr, "Simulated annealing for group %s: %s, %d timepoints\n",
                                *(groups->grpname[j]), eann_names[ir->opts.annealing[i]],
                                ir->opts.anneal_npoints[i]);
                        fprintf(stderr, "Time (ps)   Temperature (K)\n");
                        /* All terms except the last one */
                        for (j = 0; j < (ir->opts.anneal_npoints[i]-1); j++)
                        {
                            fprintf(stderr, "%9.1f      %5.1f\n", ir->opts.anneal_time[i][j], ir->opts.anneal_temp[i][j]);
                        }

                        /* Finally the last one */
                        j = ir->opts.anneal_npoints[i]-1;
                        if (ir->opts.annealing[i] == eannSINGLE)
                        {
                            fprintf(stderr, "%9.1f-     %5.1f\n", ir->opts.anneal_time[i][j], ir->opts.anneal_temp[i][j]);
                        }
                        else
                        {
                            fprintf(stderr, "%9.1f      %5.1f\n", ir->opts.anneal_time[i][j], ir->opts.anneal_temp[i][j]);
                            if (std::fabs(ir->opts.anneal_temp[i][j]-ir->opts.anneal_temp[i][0]) > GMX_REAL_EPS)
                            {
                                warning_note(wi, "There is a temperature jump when your annealing loops back.\n");
                            }
                        }
                    }
                }
            }
        }
    }

    if (ir->bPull)
    {
        make_pull_groups(ir->pull, is->pull_grp, grps, gnames);

        make_pull_coords(ir->pull);
    }

    if (ir->bRot)
    {
        make_rotation_groups(ir->rot, is->rot_grp, grps, gnames);
    }

    if (ir->eSwapCoords != eswapNO)
    {
        make_swap_groups(ir->swap, grps, gnames);
    }

    /* Make indices for IMD session */
    if (ir->bIMD)
    {
        make_IMD_group(ir->imd, is->imd_grp, grps, gnames);
    }

    auto accelerations          = gmx::splitString(is->acc);
    auto accelerationGroupNames = gmx::splitString(is->accgrps);
    if (accelerationGroupNames.size() * DIM != accelerations.size())
    {
        gmx_fatal(FARGS, "Invalid Acceleration input: %zu groups and %zu acc. values",
                  accelerationGroupNames.size(), accelerations.size());
    }
    do_numbering(natoms, groups, accelerationGroupNames, grps, gnames, egcACC,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nr = groups->grps[egcACC].nr;
    snew(ir->opts.acc, nr);
    ir->opts.ngacc = nr;

    convertRvecs(wi, accelerations, "anneal-time", ir->opts.acc);

    auto freezeDims       = gmx::splitString(is->frdim);
    auto freezeGroupNames = gmx::splitString(is->freeze);
    if (freezeDims.size() != DIM * freezeGroupNames.size())
    {
        gmx_fatal(FARGS, "Invalid Freezing input: %zu groups and %zu freeze values",
                  freezeGroupNames.size(), freezeDims.size());
    }
    do_numbering(natoms, groups, freezeGroupNames, grps, gnames, egcFREEZE,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nr             = groups->grps[egcFREEZE].nr;
    ir->opts.ngfrz = nr;
    snew(ir->opts.nFreeze, nr);
    for (i = k = 0; (size_t(i) < freezeGroupNames.size()); i++)
    {
        for (j = 0; (j < DIM); j++, k++)
        {
            ir->opts.nFreeze[i][j] = static_cast<int>(gmx_strncasecmp(freezeDims[k].c_str(), "Y", 1) == 0);
            if (!ir->opts.nFreeze[i][j])
            {
                if (gmx_strncasecmp(freezeDims[k].c_str(), "N", 1) != 0)
                {
                    sprintf(warnbuf, "Please use Y(ES) or N(O) for freezedim only "
                            "(not %s)", freezeDims[k].c_str());
                    warning(wi, warn_buf);
                }
            }
        }
    }
    for (; (i < nr); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            ir->opts.nFreeze[i][j] = 0;
        }
    }

    auto energyGroupNames = gmx::splitString(is->energy);
    do_numbering(natoms, groups, energyGroupNames, grps, gnames, egcENER,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    add_wall_energrps(groups, ir->nwall, symtab);
    ir->opts.ngener = groups->grps[egcENER].nr;
    auto vcmGroupNames = gmx::splitString(is->vcm);
    bRest           =
        do_numbering(natoms, groups, vcmGroupNames, grps, gnames, egcVCM,
                     restnm, vcmGroupNames.empty() ? egrptpALL_GENREST : egrptpPART, bVerbose, wi);
    if (bRest)
    {
        warning(wi, "Some atoms are not part of any center of mass motion removal group.\n"
                "This may lead to artifacts.\n"
                "In most cases one should use one group for the whole system.");
    }

    /* Now we have filled the freeze struct, so we can calculate NRDF */
    calc_nrdf(mtop, ir, gnames);

    auto user1GroupNames = gmx::splitString(is->user1);
    do_numbering(natoms, groups, user1GroupNames, grps, gnames, egcUser1,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    auto user2GroupNames = gmx::splitString(is->user2);
    do_numbering(natoms, groups, user2GroupNames, grps, gnames, egcUser2,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    auto compressedXGroupNames = gmx::splitString(is->x_compressed_groups);
    do_numbering(natoms, groups, compressedXGroupNames, grps, gnames, egcCompressedX,
                 restnm, egrptpONE, bVerbose, wi);
    auto orirefFitGroupNames = gmx::splitString(is->orirefitgrp);
    do_numbering(natoms, groups, orirefFitGroupNames, grps, gnames, egcORFIT,
                 restnm, egrptpALL_GENREST, bVerbose, wi);

    /* QMMM input processing */
    auto qmGroupNames = gmx::splitString(is->QMMM);
    auto qmMethods    = gmx::splitString(is->QMmethod);
    auto qmBasisSets  = gmx::splitString(is->QMbasis);
    if (ir->eI != eiMimic)
    {
        if (qmMethods.size() != qmGroupNames.size() ||
            qmBasisSets.size() != qmGroupNames.size())
        {
            gmx_fatal(FARGS, "Invalid QMMM input: %zu groups %zu basissets"
                      " and %zu methods\n", qmGroupNames.size(),
                      qmBasisSets.size(), qmMethods.size());
        }
        /* group rest, if any, is always MM! */
        do_numbering(natoms, groups, qmGroupNames, grps, gnames, egcQMMM,
                     restnm, egrptpALL_GENREST, bVerbose, wi);
        nr            = qmGroupNames.size(); /*atoms->grps[egcQMMM].nr;*/
        ir->opts.ngQM = qmGroupNames.size();
        snew(ir->opts.QMmethod, nr);
        snew(ir->opts.QMbasis, nr);
        for (i = 0; i < nr; i++)
        {
            /* input consists of strings: RHF CASSCF PM3 .. These need to be
             * converted to the corresponding enum in names.c
             */
            ir->opts.QMmethod[i] = search_QMstring(qmMethods[i].c_str(),
                                                   eQMmethodNR,
                                                   eQMmethod_names);
            ir->opts.QMbasis[i] = search_QMstring(qmBasisSets[i].c_str(),
                                                  eQMbasisNR,
                                                  eQMbasis_names);

        }
        auto qmMultiplicities = gmx::splitString(is->QMmult);
        auto qmCharges        = gmx::splitString(is->QMcharge);
        auto qmbSH            = gmx::splitString(is->bSH);
        snew(ir->opts.QMmult, nr);
        snew(ir->opts.QMcharge, nr);
        snew(ir->opts.bSH, nr);
        convertInts(wi, qmMultiplicities, "QMmult", ir->opts.QMmult);
        convertInts(wi, qmCharges, "QMcharge", ir->opts.QMcharge);
        convertYesNos(wi, qmbSH, "bSH", ir->opts.bSH);

        auto CASelectrons = gmx::splitString(is->CASelectrons);
        auto CASorbitals  = gmx::splitString(is->CASorbitals);
        snew(ir->opts.CASelectrons, nr);
        snew(ir->opts.CASorbitals, nr);
        convertInts(wi, CASelectrons, "CASelectrons", ir->opts.CASelectrons);
        convertInts(wi, CASorbitals, "CASOrbitals", ir->opts.CASorbitals);

        auto SAon    = gmx::splitString(is->SAon);
        auto SAoff   = gmx::splitString(is->SAoff);
        auto SAsteps = gmx::splitString(is->SAsteps);
        snew(ir->opts.SAon, nr);
        snew(ir->opts.SAoff, nr);
        snew(ir->opts.SAsteps, nr);
        convertInts(wi, SAon, "SAon", ir->opts.SAon);
        convertInts(wi, SAoff, "SAoff", ir->opts.SAoff);
        convertInts(wi, SAsteps, "SAsteps", ir->opts.SAsteps);
    }
    else
    {
        /* MiMiC */
        if (qmGroupNames.size() > 1)
        {
            gmx_fatal(FARGS, "Currently, having more than one QM group in MiMiC is not supported");
        }
        /* group rest, if any, is always MM! */
        do_numbering(natoms, groups, qmGroupNames, grps, gnames, egcQMMM,
                     restnm, egrptpALL_GENREST, bVerbose, wi);

        ir->opts.ngQM = qmGroupNames.size();
    }

    /* end of QMMM input */

    if (bVerbose)
    {
        for (i = 0; (i < egcNR); i++)
        {
            fprintf(stderr, "%-16s has %d element(s):", gtypes[i], groups->grps[i].nr);
            for (j = 0; (j < groups->grps[i].nr); j++)
            {
                fprintf(stderr, " %s", *(groups->grpname[groups->grps[i].nm_ind[j]]));
            }
            fprintf(stderr, "\n");
        }
    }

    nr = groups->grps[egcENER].nr;
    snew(ir->opts.egp_flags, nr*nr);

    bExcl = do_egp_flag(ir, groups, "energygrp-excl", is->egpexcl, EGP_EXCL);
    if (bExcl && ir->cutoff_scheme == ecutsVERLET)
    {
        warning_error(wi, "Energy group exclusions are not (yet) implemented for the Verlet scheme");
    }
    if (bExcl && EEL_FULL(ir->coulombtype))
    {
        warning(wi, "Can not exclude the lattice Coulomb energy between energy groups");
    }

    bTable = do_egp_flag(ir, groups, "energygrp-table", is->egptable, EGP_TABLE);
    if (bTable && !(ir->vdwtype == evdwUSER) &&
        !(ir->coulombtype == eelUSER) && !(ir->coulombtype == eelPMEUSER) &&
        !(ir->coulombtype == eelPMEUSERSWITCH))
    {
        gmx_fatal(FARGS, "Can only have energy group pair tables in combination with user tables for VdW and/or Coulomb");
    }

    for (i = 0; (i < grps->nr); i++)
    {
        sfree(gnames[i]);
    }
    sfree(gnames);
    done_blocka(grps);
    sfree(grps);

}



static void check_disre(const gmx_mtop_t *mtop)
{
    if (gmx_mtop_ftype_count(mtop, F_DISRES) > 0)
    {
        const gmx_ffparams_t &ffparams  = mtop->ffparams;
        int                   ndouble   = 0;
        int                   old_label = -1;
        for (int i = 0; i < ffparams.numTypes(); i++)
        {
            int ftype = ffparams.functype[i];
            if (ftype == F_DISRES)
            {
                int label = ffparams.iparams[i].disres.label;
                if (label == old_label)
                {
                    fprintf(stderr, "Distance restraint index %d occurs twice\n", label);
                    ndouble++;
                }
                old_label = label;
            }
        }
        if (ndouble > 0)
        {
            gmx_fatal(FARGS, "Found %d double distance restraint indices,\n"
                      "probably the parameters for multiple pairs in one restraint "
                      "are not identical\n", ndouble);
        }
    }
}

static bool absolute_reference(t_inputrec *ir, gmx_mtop_t *sys,
                               bool posres_only,
                               ivec AbsRef)
{
    int                    d, g, i;
    gmx_mtop_ilistloop_t   iloop;
    int                    nmol;
    t_iparams             *pr;

    clear_ivec(AbsRef);

    if (!posres_only)
    {
        /* Check the COM */
        for (d = 0; d < DIM; d++)
        {
            AbsRef[d] = (d < ndof_com(ir) ? 0 : 1);
        }
        /* Check for freeze groups */
        for (g = 0; g < ir->opts.ngfrz; g++)
        {
            for (d = 0; d < DIM; d++)
            {
                if (ir->opts.nFreeze[g][d] != 0)
                {
                    AbsRef[d] = 1;
                }
            }
        }
    }

    /* Check for position restraints */
    iloop = gmx_mtop_ilistloop_init(sys);
    while (const InteractionLists *ilist = gmx_mtop_ilistloop_next(iloop, &nmol))
    {
        if (nmol > 0 &&
            (AbsRef[XX] == 0 || AbsRef[YY] == 0 || AbsRef[ZZ] == 0))
        {
            for (i = 0; i < (*ilist)[F_POSRES].size(); i += 2)
            {
                pr = &sys->ffparams.iparams[(*ilist)[F_POSRES].iatoms[i]];
                for (d = 0; d < DIM; d++)
                {
                    if (pr->posres.fcA[d] != 0)
                    {
                        AbsRef[d] = 1;
                    }
                }
            }
            for (i = 0; i < (*ilist)[F_FBPOSRES].size(); i += 2)
            {
                /* Check for flat-bottom posres */
                pr = &sys->ffparams.iparams[(*ilist)[F_FBPOSRES].iatoms[i]];
                if (pr->fbposres.k != 0)
                {
                    switch (pr->fbposres.geom)
                    {
                        case efbposresSPHERE:
                            AbsRef[XX] = AbsRef[YY] = AbsRef[ZZ] = 1;
                            break;
                        case efbposresCYLINDERX:
                            AbsRef[YY] = AbsRef[ZZ] = 1;
                            break;
                        case efbposresCYLINDERY:
                            AbsRef[XX] = AbsRef[ZZ] = 1;
                            break;
                        case efbposresCYLINDER:
                        /* efbposres is a synonym for efbposresCYLINDERZ for backwards compatibility */
                        case efbposresCYLINDERZ:
                            AbsRef[XX] = AbsRef[YY] = 1;
                            break;
                        case efbposresX: /* d=XX */
                        case efbposresY: /* d=YY */
                        case efbposresZ: /* d=ZZ */
                            d         = pr->fbposres.geom - efbposresX;
                            AbsRef[d] = 1;
                            break;
                        default:
                            gmx_fatal(FARGS, " Invalid geometry for flat-bottom position restraint.\n"
                                      "Expected nr between 1 and %d. Found %d\n", efbposresNR-1,
                                      pr->fbposres.geom);
                    }
                }
            }
        }
    }

    return (AbsRef[XX] != 0 && AbsRef[YY] != 0 && AbsRef[ZZ] != 0);
}

static void
check_combination_rule_differences(const gmx_mtop_t *mtop, int state,
                                   bool *bC6ParametersWorkWithGeometricRules,
                                   bool *bC6ParametersWorkWithLBRules,
                                   bool *bLBRulesPossible)
{
    int           ntypes, tpi, tpj;
    int          *typecount;
    real          tol;
    double        c6i, c6j, c12i, c12j;
    double        c6, c6_geometric, c6_LB;
    double        sigmai, sigmaj, epsi, epsj;
    bool          bCanDoLBRules, bCanDoGeometricRules;
    const char   *ptr;

    /* A tolerance of 1e-5 seems reasonable for (possibly hand-typed)
     * force-field floating point parameters.
     */
    tol = 1e-5;
    ptr = getenv("GMX_LJCOMB_TOL");
    if (ptr != nullptr)
    {
        double            dbl;
        double gmx_unused canary;

        if (sscanf(ptr, "%lf%lf", &dbl, &canary) != 1)
        {
            gmx_fatal(FARGS, "Could not parse a single floating-point number from GMX_LJCOMB_TOL (%s)", ptr);
        }
        tol = dbl;
    }

    *bC6ParametersWorkWithLBRules         = TRUE;
    *bC6ParametersWorkWithGeometricRules  = TRUE;
    bCanDoLBRules                         = TRUE;
    ntypes                                = mtop->ffparams.atnr;
    snew(typecount, ntypes);
    gmx_mtop_count_atomtypes(mtop, state, typecount);
    *bLBRulesPossible       = TRUE;
    for (tpi = 0; tpi < ntypes; ++tpi)
    {
        c6i  = mtop->ffparams.iparams[(ntypes + 1) * tpi].lj.c6;
        c12i = mtop->ffparams.iparams[(ntypes + 1) * tpi].lj.c12;
        for (tpj = tpi; tpj < ntypes; ++tpj)
        {
            c6j          = mtop->ffparams.iparams[(ntypes + 1) * tpj].lj.c6;
            c12j         = mtop->ffparams.iparams[(ntypes + 1) * tpj].lj.c12;
            c6           = mtop->ffparams.iparams[ntypes * tpi + tpj].lj.c6;
            c6_geometric = std::sqrt(c6i * c6j);
            if (!gmx_numzero(c6_geometric))
            {
                if (!gmx_numzero(c12i) && !gmx_numzero(c12j))
                {
                    sigmai   = gmx::sixthroot(c12i / c6i);
                    sigmaj   = gmx::sixthroot(c12j / c6j);
                    epsi     = c6i * c6i /(4.0 * c12i);
                    epsj     = c6j * c6j /(4.0 * c12j);
                    c6_LB    = 4.0 * std::sqrt(epsi * epsj) * gmx::power6(0.5 * (sigmai + sigmaj));
                }
                else
                {
                    *bLBRulesPossible = FALSE;
                    c6_LB             = c6_geometric;
                }
                bCanDoLBRules = gmx_within_tol(c6_LB, c6, tol);
            }

            if (!bCanDoLBRules)
            {
                *bC6ParametersWorkWithLBRules = FALSE;
            }

            bCanDoGeometricRules = gmx_within_tol(c6_geometric, c6, tol);

            if (!bCanDoGeometricRules)
            {
                *bC6ParametersWorkWithGeometricRules = FALSE;
            }
        }
    }
    sfree(typecount);
}

static void
check_combination_rules(const t_inputrec *ir, const gmx_mtop_t *mtop,
                        warninp_t wi)
{
    bool bLBRulesPossible, bC6ParametersWorkWithGeometricRules, bC6ParametersWorkWithLBRules;

    check_combination_rule_differences(mtop, 0,
                                       &bC6ParametersWorkWithGeometricRules,
                                       &bC6ParametersWorkWithLBRules,
                                       &bLBRulesPossible);
    if (ir->ljpme_combination_rule == eljpmeLB)
    {
        if (!bC6ParametersWorkWithLBRules || !bLBRulesPossible)
        {
            warning(wi, "You are using arithmetic-geometric combination rules "
                    "in LJ-PME, but your non-bonded C6 parameters do not "
                    "follow these rules.");
        }
    }
    else
    {
        if (!bC6ParametersWorkWithGeometricRules)
        {
            if (ir->eDispCorr != edispcNO)
            {
                warning_note(wi, "You are using geometric combination rules in "
                             "LJ-PME, but your non-bonded C6 parameters do "
                             "not follow these rules. "
                             "This will introduce very small errors in the forces and energies in "
                             "your simulations. Dispersion correction will correct total energy "
                             "and/or pressure for isotropic systems, but not forces or surface tensions.");
            }
            else
            {
                warning_note(wi, "You are using geometric combination rules in "
                             "LJ-PME, but your non-bonded C6 parameters do "
                             "not follow these rules. "
                             "This will introduce very small errors in the forces and energies in "
                             "your simulations. If your system is homogeneous, consider using dispersion correction "
                             "for the total energy and pressure.");
            }
        }
    }
}

void triple_check(const char *mdparin, t_inputrec *ir, gmx_mtop_t *sys,
                  warninp_t wi)
{
    char                      err_buf[STRLEN];
    int                       i, m, c, nmol;
    bool                      bCharge, bAcc;
    real                     *mgrp, mt;
    rvec                      acc;
    gmx_mtop_atomloop_block_t aloopb;
    gmx_mtop_atomloop_all_t   aloop;
    ivec                      AbsRef;
    char                      warn_buf[STRLEN];

    set_warning_line(wi, mdparin, -1);

    if (ir->cutoff_scheme == ecutsVERLET &&
        ir->verletbuf_tol > 0 &&
        ir->nstlist > 1 &&
        ((EI_MD(ir->eI) || EI_SD(ir->eI)) &&
         (ir->etc == etcVRESCALE || ir->etc == etcBERENDSEN)))
    {
        /* Check if a too small Verlet buffer might potentially
         * cause more drift than the thermostat can couple off.
         */
        /* Temperature error fraction for warning and suggestion */
        const real T_error_warn    = 0.002;
        const real T_error_suggest = 0.001;
        /* For safety: 2 DOF per atom (typical with constraints) */
        const real nrdf_at         = 2;
        real       T, tau, max_T_error;
        int        i;

        T   = 0;
        tau = 0;
        for (i = 0; i < ir->opts.ngtc; i++)
        {
            T   = std::max(T, ir->opts.ref_t[i]);
            tau = std::max(tau, ir->opts.tau_t[i]);
        }
        if (T > 0)
        {
            /* This is a worst case estimate of the temperature error,
             * assuming perfect buffer estimation and no cancelation
             * of errors. The factor 0.5 is because energy distributes
             * equally over Ekin and Epot.
             */
            max_T_error = 0.5*tau*ir->verletbuf_tol/(nrdf_at*BOLTZ*T);
            if (max_T_error > T_error_warn)
            {
                sprintf(warn_buf, "With a verlet-buffer-tolerance of %g kJ/mol/ps, a reference temperature of %g and a tau_t of %g, your temperature might be off by up to %.1f%%. To ensure the error is below %.1f%%, decrease verlet-buffer-tolerance to %.0e or decrease tau_t.",
                        ir->verletbuf_tol, T, tau,
                        100*max_T_error,
                        100*T_error_suggest,
                        ir->verletbuf_tol*T_error_suggest/max_T_error);
                warning(wi, warn_buf);
            }
        }
    }

    if (ETC_ANDERSEN(ir->etc))
    {
        int i;

        for (i = 0; i < ir->opts.ngtc; i++)
        {
            sprintf(err_buf, "all tau_t must currently be equal using Andersen temperature control, violated for group %d", i);
            CHECK(ir->opts.tau_t[0] != ir->opts.tau_t[i]);
            sprintf(err_buf, "all tau_t must be positive using Andersen temperature control, tau_t[%d]=%10.6f",
                    i, ir->opts.tau_t[i]);
            CHECK(ir->opts.tau_t[i] < 0);
        }

        if (ir->etc == etcANDERSENMASSIVE && ir->comm_mode != ecmNO)
        {
            for (i = 0; i < ir->opts.ngtc; i++)
            {
                int nsteps = gmx::roundToInt(ir->opts.tau_t[i]/ir->delta_t);
                sprintf(err_buf, "tau_t/delta_t for group %d for temperature control method %s must be a multiple of nstcomm (%d), as velocities of atoms in coupled groups are randomized every time step. The input tau_t (%8.3f) leads to %d steps per randomization", i, etcoupl_names[ir->etc], ir->nstcomm, ir->opts.tau_t[i], nsteps);
                CHECK(nsteps % ir->nstcomm != 0);
            }
        }
    }

    if (EI_DYNAMICS(ir->eI) && !EI_SD(ir->eI) && ir->eI != eiBD &&
        ir->comm_mode == ecmNO &&
        !(absolute_reference(ir, sys, FALSE, AbsRef) || ir->nsteps <= 10) &&
        !ETC_ANDERSEN(ir->etc))
    {
        warning(wi, "You are not using center of mass motion removal (mdp option comm-mode), numerical rounding errors can lead to build up of kinetic energy of the center of mass");
    }

    /* Check for pressure coupling with absolute position restraints */
    if (ir->epc != epcNO && ir->refcoord_scaling == erscNO)
    {
        absolute_reference(ir, sys, TRUE, AbsRef);
        {
            for (m = 0; m < DIM; m++)
            {
                if (AbsRef[m] && norm2(ir->compress[m]) > 0)
                {
                    warning(wi, "You are using pressure coupling with absolute position restraints, this will give artifacts. Use the refcoord_scaling option.");
                    break;
                }
            }
        }
    }

    bCharge = FALSE;
    aloopb  = gmx_mtop_atomloop_block_init(sys);
    const t_atom *atom;
    while (gmx_mtop_atomloop_block_next(aloopb, &atom, &nmol))
    {
        if (atom->q != 0 || atom->qB != 0)
        {
            bCharge = TRUE;
        }
    }

    if (!bCharge)
    {
        if (EEL_FULL(ir->coulombtype))
        {
            sprintf(err_buf,
                    "You are using full electrostatics treatment %s for a system without charges.\n"
                    "This costs a lot of performance for just processing zeros, consider using %s instead.\n",
                    EELTYPE(ir->coulombtype), EELTYPE(eelCUT));
            warning(wi, err_buf);
        }
    }
    else
    {
        if (ir->coulombtype == eelCUT && ir->rcoulomb > 0)
        {
            sprintf(err_buf,
                    "You are using a plain Coulomb cut-off, which might produce artifacts.\n"
                    "You might want to consider using %s electrostatics.\n",
                    EELTYPE(eelPME));
            warning_note(wi, err_buf);
        }
    }

    /* Check if combination rules used in LJ-PME are the same as in the force field */
    if (EVDW_PME(ir->vdwtype))
    {
        check_combination_rules(ir, sys, wi);
    }

    /* Generalized reaction field */
    if (ir->opts.ngtc == 0)
    {
        sprintf(err_buf, "No temperature coupling while using coulombtype %s",
                eel_names[eelGRF]);
        CHECK(ir->coulombtype == eelGRF);
    }
    else
    {
        sprintf(err_buf, "When using coulombtype = %s"
                " ref-t for temperature coupling should be > 0",
                eel_names[eelGRF]);
        CHECK((ir->coulombtype == eelGRF) && (ir->opts.ref_t[0] <= 0));
    }

    bAcc = FALSE;
    for (i = 0; (i < sys->groups.grps[egcACC].nr); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            if (fabs(ir->opts.acc[i][m]) > 1e-6)
            {
                bAcc = TRUE;
            }
        }
    }
    if (bAcc)
    {
        clear_rvec(acc);
        snew(mgrp, sys->groups.grps[egcACC].nr);
        aloop = gmx_mtop_atomloop_all_init(sys);
        const t_atom *atom;
        while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
        {
            mgrp[getGroupType(&sys->groups, egcACC, i)] += atom->m;
        }
        mt = 0.0;
        for (i = 0; (i < sys->groups.grps[egcACC].nr); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                acc[m] += ir->opts.acc[i][m]*mgrp[i];
            }
            mt += mgrp[i];
        }
        for (m = 0; (m < DIM); m++)
        {
            if (fabs(acc[m]) > 1e-6)
            {
                const char *dim[DIM] = { "X", "Y", "Z" };
                fprintf(stderr,
                        "Net Acceleration in %s direction, will %s be corrected\n",
                        dim[m], ir->nstcomm != 0 ? "" : "not");
                if (ir->nstcomm != 0 && m < ndof_com(ir))
                {
                    acc[m] /= mt;
                    for (i = 0; (i < sys->groups.grps[egcACC].nr); i++)
                    {
                        ir->opts.acc[i][m] -= acc[m];
                    }
                }
            }
        }
        sfree(mgrp);
    }

    if (ir->efep != efepNO && ir->fepvals->sc_alpha != 0 &&
        !gmx_within_tol(sys->ffparams.reppow, 12.0, 10*GMX_DOUBLE_EPS))
    {
        gmx_fatal(FARGS, "Soft-core interactions are only supported with VdW repulsion power 12");
    }

    if (ir->bPull)
    {
        bool bWarned;

        bWarned = FALSE;
        for (i = 0; i < ir->pull->ncoord && !bWarned; i++)
        {
            if (ir->pull->coord[i].group[0] == 0 ||
                ir->pull->coord[i].group[1] == 0)
            {
                absolute_reference(ir, sys, FALSE, AbsRef);
                for (m = 0; m < DIM; m++)
                {
                    if (ir->pull->coord[i].dim[m] && !AbsRef[m])
                    {
                        warning(wi, "You are using an absolute reference for pulling, but the rest of the system does not have an absolute reference. This will lead to artifacts.");
                        bWarned = TRUE;
                        break;
                    }
                }
            }
        }

        for (i = 0; i < 3; i++)
        {
            for (m = 0; m <= i; m++)
            {
                if ((ir->epc != epcNO && ir->compress[i][m] != 0) ||
                    ir->deform[i][m] != 0)
                {
                    for (c = 0; c < ir->pull->ncoord; c++)
                    {
                        if (ir->pull->coord[c].eGeom == epullgDIRPBC &&
                            ir->pull->coord[c].vec[m] != 0)
                        {
                            gmx_fatal(FARGS, "Can not have dynamic box while using pull geometry '%s' (dim %c)", EPULLGEOM(ir->pull->coord[c].eGeom), 'x'+m);
                        }
                    }
                }
            }
        }
    }

    check_disre(sys);
}

void double_check(t_inputrec *ir, matrix box,
                  bool bHasNormalConstraints,
                  bool bHasAnyConstraints,
                  warninp_t wi)
{
    real        min_size;
    char        warn_buf[STRLEN];
    const char *ptr;

    ptr = check_box(ir->ePBC, box);
    if (ptr)
    {
        warning_error(wi, ptr);
    }

    if (bHasNormalConstraints && ir->eConstrAlg == econtSHAKE)
    {
        if (ir->shake_tol <= 0.0)
        {
            sprintf(warn_buf, "ERROR: shake-tol must be > 0 instead of %g\n",
                    ir->shake_tol);
            warning_error(wi, warn_buf);
        }
    }

    if ( (ir->eConstrAlg == econtLINCS) && bHasNormalConstraints)
    {
        /* If we have Lincs constraints: */
        if (ir->eI == eiMD && ir->etc == etcNO &&
            ir->eConstrAlg == econtLINCS && ir->nLincsIter == 1)
        {
            sprintf(warn_buf, "For energy conservation with LINCS, lincs_iter should be 2 or larger.\n");
            warning_note(wi, warn_buf);
        }

        if ((ir->eI == eiCG || ir->eI == eiLBFGS) && (ir->nProjOrder < 8))
        {
            sprintf(warn_buf, "For accurate %s with LINCS constraints, lincs-order should be 8 or more.", ei_names[ir->eI]);
            warning_note(wi, warn_buf);
        }
        if (ir->epc == epcMTTK)
        {
            warning_error(wi, "MTTK not compatible with lincs -- use shake instead.");
        }
    }

    if (bHasAnyConstraints && ir->epc == epcMTTK)
    {
        warning_error(wi, "Constraints are not implemented with MTTK pressure control.");
    }

    if (ir->LincsWarnAngle > 90.0)
    {
        sprintf(warn_buf, "lincs-warnangle can not be larger than 90 degrees, setting it to 90.\n");
        warning(wi, warn_buf);
        ir->LincsWarnAngle = 90.0;
    }

    if (ir->ePBC != epbcNONE)
    {
        if (ir->nstlist == 0)
        {
            warning(wi, "With nstlist=0 atoms are only put into the box at step 0, therefore drifting atoms might cause the simulation to crash.");
        }
        if (ir->ns_type == ensGRID)
        {
            if (gmx::square(ir->rlist) >= max_cutoff2(ir->ePBC, box))
            {
                sprintf(warn_buf, "ERROR: The cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element. Increase the box size or decrease rlist.\n");
                warning_error(wi, warn_buf);
            }
        }
        else
        {
            min_size = std::min(box[XX][XX], std::min(box[YY][YY], box[ZZ][ZZ]));
            if (2*ir->rlist >= min_size)
            {
                sprintf(warn_buf, "ERROR: One of the box lengths is smaller than twice the cut-off length. Increase the box size or decrease rlist.");
                warning_error(wi, warn_buf);
                if (TRICLINIC(box))
                {
                    fprintf(stderr, "Grid search might allow larger cut-off's than simple search with triclinic boxes.");
                }
            }
        }
    }
}

void check_chargegroup_radii(const gmx_mtop_t *mtop, const t_inputrec *ir,
                             rvec *x,
                             warninp_t wi)
{
    real rvdw1, rvdw2, rcoul1, rcoul2;
    char warn_buf[STRLEN];

    calc_chargegroup_radii(mtop, x, &rvdw1, &rvdw2, &rcoul1, &rcoul2);

    if (rvdw1 > 0)
    {
        printf("Largest charge group radii for Van der Waals: %5.3f, %5.3f nm\n",
               rvdw1, rvdw2);
    }
    if (rcoul1 > 0)
    {
        printf("Largest charge group radii for Coulomb:       %5.3f, %5.3f nm\n",
               rcoul1, rcoul2);
    }

    if (ir->rlist > 0)
    {
        if (rvdw1  + rvdw2  > ir->rlist ||
            rcoul1 + rcoul2 > ir->rlist)
        {
            sprintf(warn_buf,
                    "The sum of the two largest charge group radii (%f) "
                    "is larger than rlist (%f)\n",
                    std::max(rvdw1+rvdw2, rcoul1+rcoul2), ir->rlist);
            warning(wi, warn_buf);
        }
        else
        {
            /* Here we do not use the zero at cut-off macro,
             * since user defined interactions might purposely
             * not be zero at the cut-off.
             */
            if (ir_vdw_is_zero_at_cutoff(ir) &&
                rvdw1 + rvdw2 > ir->rlist - ir->rvdw)
            {
                sprintf(warn_buf, "The sum of the two largest charge group "
                        "radii (%f) is larger than rlist (%f) - rvdw (%f).\n"
                        "With exact cut-offs, better performance can be "
                        "obtained with cutoff-scheme = %s, because it "
                        "does not use charge groups at all.",
                        rvdw1+rvdw2,
                        ir->rlist, ir->rvdw,
                        ecutscheme_names[ecutsVERLET]);
                if (ir_NVE(ir))
                {
                    warning(wi, warn_buf);
                }
                else
                {
                    warning_note(wi, warn_buf);
                }
            }
            if (ir_coulomb_is_zero_at_cutoff(ir) &&
                rcoul1 + rcoul2 > ir->rlist - ir->rcoulomb)
            {
                sprintf(warn_buf, "The sum of the two largest charge group radii (%f) is larger than rlist (%f) - rcoulomb (%f).\n"
                        "With exact cut-offs, better performance can be obtained with cutoff-scheme = %s, because it does not use charge groups at all.",
                        rcoul1+rcoul2,
                        ir->rlist, ir->rcoulomb,
                        ecutscheme_names[ecutsVERLET]);
                if (ir_NVE(ir))
                {
                    warning(wi, warn_buf);
                }
                else
                {
                    warning_note(wi, warn_buf);
                }
            }
        }
    }
}
