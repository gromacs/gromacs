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

#include "readir.h"

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>

#include <cmath>

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
            gmx_fatal(FARGS, errorstr);
        }
    }
}



static void _low_check(gmx_bool b, const char *s, warninp_t wi)
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

static gmx_bool ir_NVE(const t_inputrec *ir)
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

        if (ir->implicit_solvent != eisNO)
        {
            warning_error(wi, "Implicit solvent is not (yet) supported with the with Verlet lists.");
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
        gmx_bool bAllTempZero = TRUE;
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
                (int)fep->sc_r_power);
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
        sprintf(err_buf, "If there is no temperature control, and lmc-mcmove!= 'no',mc_temperature must be set to a positive number");
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

    if (EI_STATE_VELOCITY(ir->eI) && ir->ePBC == epbcNONE && ir->comm_mode != ecmANGULAR)
    {
        warning_note(wi, "Tumbling and or flying ice-cubes: We are not removing rotation around center of mass in a non-periodic system. You should probably set comm_mode = ANGULAR.");
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

    if (ir->epsilon_r != 1 && ir->implicit_solvent == eisGBSA)
    {
        sprintf(warn_buf, "epsilon-r = %g with GB implicit solvent, will use this value for inner dielectric", ir->epsilon_r);
        warning_note(wi, warn_buf);
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
                "It is pointless to use long-range or Generalized Born electrostatics with infinite relative permittivity."
                "Since you are effectively turning of electrostatics, a plain cutoff will be much faster.");
        CHECK(EEL_FULL(ir->coulombtype) || ir->implicit_solvent == eisGBSA);
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

    if (ir->sa_algorithm == esaSTILL)
    {
        sprintf(err_buf, "Still SA algorithm not available yet, use %s or %s instead\n", esa_names[esaAPPROX], esa_names[esaNO]);
        CHECK(ir->sa_algorithm == esaSTILL);
    }

    if (ir->implicit_solvent == eisGBSA)
    {
        sprintf(err_buf, "With GBSA implicit solvent, rgbradii must be equal to rlist.");
        CHECK(ir->rgbradii != ir->rlist);

        if (ir->coulombtype != eelCUT)
        {
            sprintf(err_buf, "With GBSA, coulombtype must be equal to %s\n", eel_names[eelCUT]);
            CHECK(ir->coulombtype != eelCUT);
        }
        if (ir->vdwtype != evdwCUT)
        {
            sprintf(err_buf, "With GBSA, vdw-type must be equal to %s\n", evdw_names[evdwCUT]);
            CHECK(ir->vdwtype != evdwCUT);
        }
        if (ir->nstgbradii < 1)
        {
            sprintf(warn_buf, "Using GBSA with nstgbradii<1, setting nstgbradii=1");
            warning_note(wi, warn_buf);
            ir->nstgbradii = 1;
        }
        if (ir->sa_algorithm == esaNO)
        {
            sprintf(warn_buf, "No SA (non-polar) calculation requested together with GB. Are you sure this is what you want?\n");
            warning_note(wi, warn_buf);
        }
        if (ir->sa_surface_tension < 0 && ir->sa_algorithm != esaNO)
        {
            sprintf(warn_buf, "Value of sa_surface_tension is < 0. Changing it to 2.05016 or 2.25936 kJ/nm^2/mol for Still and HCT/OBC respectively\n");
            warning_note(wi, warn_buf);

            if (ir->gb_algorithm == egbSTILL)
            {
                ir->sa_surface_tension = 0.0049 * CAL2JOULE * 100;
            }
            else
            {
                ir->sa_surface_tension = 0.0054 * CAL2JOULE * 100;
            }
        }
        if (ir->sa_surface_tension == 0 && ir->sa_algorithm != esaNO)
        {
            sprintf(err_buf, "Surface tension set to 0 while SA-calculation requested\n");
            CHECK(ir->sa_surface_tension == 0 && ir->sa_algorithm != esaNO);
        }

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
    char *ptr[MAXPTR];
    char *endptr;
    int   i;
    char  warn_buf[STRLEN];

    *n = str_nelem(str, MAXPTR, ptr);

    snew(*r, *n);
    for (i = 0; i < *n; i++)
    {
        (*r)[i] = strtod(ptr[i], &endptr);
        if (*endptr != 0)
        {
            sprintf(warn_buf, "Invalid value %s in string in mdp file. Expected a real number.", ptr[i]);
            warning_error(wi, warn_buf);
        }
    }
}

static void do_fep_params(t_inputrec *ir, char fep_lambda[][STRLEN], char weights[STRLEN], warninp_t wi)
{

    int         i, j, max_n_lambda, nweights, nfep[efptNR];
    t_lambda   *fep    = ir->fepvals;
    t_expanded *expand = ir->expandedvals;
    real      **count_fep_lambdas;
    gmx_bool    bOneLambda = TRUE;

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
                fep->all_lambda[i][j] = (double)count_fep_lambdas[i][j];
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
        expand->nstexpanded = 2*(int)(ir->opts.tau_t[0]/ir->delta_t);
        /* if you don't specify nstexpanded when doing expanded ensemble simulated tempering, it is set to
           2*tau_t just to be careful so it's not to frequent  */
    }
}


static void do_simtemp_params(t_inputrec *ir)
{

    snew(ir->simtempvals->temperatures, ir->fepvals->n_lambda);
    GetSimTemps(ir->fepvals->n_lambda, ir->simtempvals, ir->fepvals->all_lambda[efptTEMPERATURE]);

    return;
}

static void do_wall_params(t_inputrec *ir,
                           char *wall_atomtype, char *wall_density,
                           t_gromppopts *opts)
{
    int    nstr, i;
    char  *names[MAXPTR];
    double dbl;

    opts->wall_atomtype[0] = nullptr;
    opts->wall_atomtype[1] = nullptr;

    ir->wall_atomtype[0] = -1;
    ir->wall_atomtype[1] = -1;
    ir->wall_density[0]  = 0;
    ir->wall_density[1]  = 0;

    if (ir->nwall > 0)
    {
        nstr = str_nelem(wall_atomtype, MAXPTR, names);
        if (nstr != ir->nwall)
        {
            gmx_fatal(FARGS, "Expected %d elements for wall_atomtype, found %d",
                      ir->nwall, nstr);
        }
        for (i = 0; i < ir->nwall; i++)
        {
            opts->wall_atomtype[i] = gmx_strdup(names[i]);
        }

        if (ir->wall_type == ewt93 || ir->wall_type == ewt104)
        {
            nstr = str_nelem(wall_density, MAXPTR, names);
            if (nstr != ir->nwall)
            {
                gmx_fatal(FARGS, "Expected %d elements for wall-density, found %d", ir->nwall, nstr);
            }
            for (i = 0; i < ir->nwall; i++)
            {
                if (sscanf(names[i], "%lf", &dbl) != 1)
                {
                    gmx_fatal(FARGS, "Could not parse wall-density value from string '%s'", names[i]);
                }
                if (dbl <= 0)
                {
                    gmx_fatal(FARGS, "wall-density[%d] = %f\n", i, dbl);
                }
                ir->wall_density[i] = dbl;
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

static void read_expandedparams(int *ninp_p, t_inpfile **inp_p,
                                t_expanded *expand, warninp_t wi)
{
    int        ninp;
    t_inpfile *inp;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read expanded ensemble parameters */
    CCTYPE ("expanded ensemble variables");
    ITYPE ("nstexpanded", expand->nstexpanded, -1);
    EETYPE("lmc-stats", expand->elamstats, elamstats_names);
    EETYPE("lmc-move", expand->elmcmove, elmcmove_names);
    EETYPE("lmc-weights-equil", expand->elmceq, elmceq_names);
    ITYPE ("weight-equil-number-all-lambda", expand->equil_n_at_lam, -1);
    ITYPE ("weight-equil-number-samples", expand->equil_samples, -1);
    ITYPE ("weight-equil-number-steps", expand->equil_steps, -1);
    RTYPE ("weight-equil-wl-delta", expand->equil_wl_delta, -1);
    RTYPE ("weight-equil-count-ratio", expand->equil_ratio, -1);
    CCTYPE("Seed for Monte Carlo in lambda space");
    ITYPE ("lmc-seed", expand->lmc_seed, -1);
    RTYPE ("mc-temperature", expand->mc_temp, -1);
    ITYPE ("lmc-repeats", expand->lmc_repeats, 1);
    ITYPE ("lmc-gibbsdelta", expand->gibbsdeltalam, -1);
    ITYPE ("lmc-forced-nstart", expand->lmc_forced_nstart, 0);
    EETYPE("symmetrized-transition-matrix", expand->bSymmetrizedTMatrix, yesno_names);
    ITYPE("nst-transition-matrix", expand->nstTij, -1);
    ITYPE ("mininum-var-min", expand->minvarmin, 100); /*default is reasonable */
    ITYPE ("weight-c-range", expand->c_range, 0);      /* default is just C=0 */
    RTYPE ("wl-scale", expand->wl_scale, 0.8);
    RTYPE ("wl-ratio", expand->wl_ratio, 0.8);
    RTYPE ("init-wl-delta", expand->init_wl_delta, 1.0);
    EETYPE("wl-oneovert", expand->bWLoneovert, yesno_names);

    *ninp_p   = ninp;
    *inp_p    = inp;

    return;
}

/*! \brief Return whether an end state with the given coupling-lambda
 * value describes fully-interacting VDW.
 *
 * \param[in]  couple_lambda_value  Enumeration ecouplam value describing the end state
 * \return                          Whether VDW is on (i.e. the user chose vdw or vdw-q in the .mdp file)
 */
static gmx_bool couple_lambda_has_vdw_on(int couple_lambda_value)
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

        virtual bool onError(gmx::UserInputError *ex, const gmx::KeyValueTreePath &context)
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
    char       *dumstr[2];
    double      dumdub[2][6];
    t_inpfile  *inp;
    const char *tmp;
    int         i, j, m, ninp;
    char        warn_buf[STRLEN];
    t_lambda   *fep    = ir->fepvals;
    t_expanded *expand = ir->expandedvals;

    init_inputrec_strings();
    gmx::TextInputFile stream(mdparin);
    inp = read_inpfile(&stream, mdparin, &ninp, wi);

    snew(dumstr[0], STRLEN);
    snew(dumstr[1], STRLEN);

    if (-1 == search_einp(ninp, inp, "cutoff-scheme"))
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
    REM_TYPE("title");
    REM_TYPE("cpp");
    REM_TYPE("domain-decomposition");
    REM_TYPE("andersen-seed");
    REM_TYPE("dihre");
    REM_TYPE("dihre-fc");
    REM_TYPE("dihre-tau");
    REM_TYPE("nstdihreout");
    REM_TYPE("nstcheckpoint");
    REM_TYPE("optimize-fft");
    REM_TYPE("adress_type");
    REM_TYPE("adress_const_wf");
    REM_TYPE("adress_ex_width");
    REM_TYPE("adress_hy_width");
    REM_TYPE("adress_ex_forcecap");
    REM_TYPE("adress_interface_correction");
    REM_TYPE("adress_site");
    REM_TYPE("adress_reference_coords");
    REM_TYPE("adress_tf_grp_names");
    REM_TYPE("adress_cg_grp_names");
    REM_TYPE("adress_do_hybridpairs");
    REM_TYPE("rlistlong");
    REM_TYPE("nstcalclr");
    REM_TYPE("pull-print-com2");

    /* replace the following commands with the clearer new versions*/
    REPL_TYPE("unconstrained-start", "continuation");
    REPL_TYPE("foreign-lambda", "fep-lambdas");
    REPL_TYPE("verlet-buffer-drift", "verlet-buffer-tolerance");
    REPL_TYPE("nstxtcout", "nstxout-compressed");
    REPL_TYPE("xtc-grps", "compressed-x-grps");
    REPL_TYPE("xtc-precision", "compressed-x-precision");
    REPL_TYPE("pull-print-com1", "pull-print-com");

    CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
    CTYPE ("Preprocessor information: use cpp syntax.");
    CTYPE ("e.g.: -I/home/joe/doe -I/home/mary/roe");
    STYPE ("include", opts->include,  nullptr);
    CTYPE ("e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)");
    STYPE ("define",  opts->define,   nullptr);

    CCTYPE ("RUN CONTROL PARAMETERS");
    EETYPE("integrator",  ir->eI,         ei_names);
    CTYPE ("Start time and timestep in ps");
    RTYPE ("tinit",   ir->init_t, 0.0);
    RTYPE ("dt",      ir->delta_t,    0.001);
    STEPTYPE ("nsteps",   ir->nsteps,     0);
    CTYPE ("For exact run continuation or redoing part of a run");
    STEPTYPE ("init-step", ir->init_step,  0);
    CTYPE ("Part index is updated automatically on checkpointing (keeps files separate)");
    ITYPE ("simulation-part", ir->simulation_part, 1);
    CTYPE ("mode for center of mass motion removal");
    EETYPE("comm-mode",   ir->comm_mode,  ecm_names);
    CTYPE ("number of steps for center of mass motion removal");
    ITYPE ("nstcomm", ir->nstcomm,    100);
    CTYPE ("group(s) for center of mass motion removal");
    STYPE ("comm-grps",   is->vcm,            nullptr);

    CCTYPE ("LANGEVIN DYNAMICS OPTIONS");
    CTYPE ("Friction coefficient (amu/ps) and random seed");
    RTYPE ("bd-fric",     ir->bd_fric,    0.0);
    STEPTYPE ("ld-seed",  ir->ld_seed,    -1);

    /* Em stuff */
    CCTYPE ("ENERGY MINIMIZATION OPTIONS");
    CTYPE ("Force tolerance and initial step-size");
    RTYPE ("emtol",       ir->em_tol,     10.0);
    RTYPE ("emstep",      ir->em_stepsize, 0.01);
    CTYPE ("Max number of iterations in relax-shells");
    ITYPE ("niter",       ir->niter,      20);
    CTYPE ("Step size (ps^2) for minimization of flexible constraints");
    RTYPE ("fcstep",      ir->fc_stepsize, 0);
    CTYPE ("Frequency of steepest descents steps when doing CG");
    ITYPE ("nstcgsteep",  ir->nstcgsteep, 1000);
    ITYPE ("nbfgscorr",   ir->nbfgscorr,  10);

    CCTYPE ("TEST PARTICLE INSERTION OPTIONS");
    RTYPE ("rtpi",    ir->rtpi,   0.05);

    /* Output options */
    CCTYPE ("OUTPUT CONTROL OPTIONS");
    CTYPE ("Output frequency for coords (x), velocities (v) and forces (f)");
    ITYPE ("nstxout", ir->nstxout,    0);
    ITYPE ("nstvout", ir->nstvout,    0);
    ITYPE ("nstfout", ir->nstfout,    0);
    CTYPE ("Output frequency for energies to log file and energy file");
    ITYPE ("nstlog",  ir->nstlog, 1000);
    ITYPE ("nstcalcenergy", ir->nstcalcenergy, 100);
    ITYPE ("nstenergy",   ir->nstenergy,  1000);
    CTYPE ("Output frequency and precision for .xtc file");
    ITYPE ("nstxout-compressed", ir->nstxout_compressed,  0);
    RTYPE ("compressed-x-precision", ir->x_compression_precision, 1000.0);
    CTYPE ("This selects the subset of atoms for the compressed");
    CTYPE ("trajectory file. You can select multiple groups. By");
    CTYPE ("default, all atoms will be written.");
    STYPE ("compressed-x-grps", is->x_compressed_groups, nullptr);
    CTYPE ("Selection of energy groups");
    STYPE ("energygrps",  is->energy,         nullptr);

    /* Neighbor searching */
    CCTYPE ("NEIGHBORSEARCHING PARAMETERS");
    CTYPE ("cut-off scheme (Verlet: particle based cut-offs, group: using charge groups)");
    EETYPE("cutoff-scheme",     ir->cutoff_scheme,    ecutscheme_names);
    CTYPE ("nblist update frequency");
    ITYPE ("nstlist", ir->nstlist,    10);
    CTYPE ("ns algorithm (simple or grid)");
    EETYPE("ns-type",     ir->ns_type,    ens_names);
    CTYPE ("Periodic boundary conditions: xyz, no, xy");
    EETYPE("pbc",         ir->ePBC,       epbc_names);
    EETYPE("periodic-molecules", ir->bPeriodicMols, yesno_names);
    CTYPE ("Allowed energy error due to the Verlet buffer in kJ/mol/ps per atom,");
    CTYPE ("a value of -1 means: use rlist");
    RTYPE("verlet-buffer-tolerance", ir->verletbuf_tol,    0.005);
    CTYPE ("nblist cut-off");
    RTYPE ("rlist",   ir->rlist,  1.0);
    CTYPE ("long-range cut-off for switched potentials");

    /* Electrostatics */
    CCTYPE ("OPTIONS FOR ELECTROSTATICS AND VDW");
    CTYPE ("Method for doing electrostatics");
    EETYPE("coulombtype", ir->coulombtype,    eel_names);
    EETYPE("coulomb-modifier",    ir->coulomb_modifier,    eintmod_names);
    CTYPE ("cut-off lengths");
    RTYPE ("rcoulomb-switch", ir->rcoulomb_switch,    0.0);
    RTYPE ("rcoulomb",    ir->rcoulomb,   1.0);
    CTYPE ("Relative dielectric constant for the medium and the reaction field");
    RTYPE ("epsilon-r",   ir->epsilon_r,  1.0);
    RTYPE ("epsilon-rf",  ir->epsilon_rf, 0.0);
    CTYPE ("Method for doing Van der Waals");
    EETYPE("vdw-type",    ir->vdwtype,    evdw_names);
    EETYPE("vdw-modifier",    ir->vdw_modifier,    eintmod_names);
    CTYPE ("cut-off lengths");
    RTYPE ("rvdw-switch", ir->rvdw_switch,    0.0);
    RTYPE ("rvdw",    ir->rvdw,   1.0);
    CTYPE ("Apply long range dispersion corrections for Energy and Pressure");
    EETYPE("DispCorr",    ir->eDispCorr,  edispc_names);
    CTYPE ("Extension of the potential lookup tables beyond the cut-off");
    RTYPE ("table-extension", ir->tabext, 1.0);
    CTYPE ("Separate tables between energy group pairs");
    STYPE ("energygrp-table", is->egptable,   nullptr);
    CTYPE ("Spacing for the PME/PPPM FFT grid");
    RTYPE ("fourierspacing", ir->fourier_spacing, 0.12);
    CTYPE ("FFT grid size, when a value is 0 fourierspacing will be used");
    ITYPE ("fourier-nx",  ir->nkx,         0);
    ITYPE ("fourier-ny",  ir->nky,         0);
    ITYPE ("fourier-nz",  ir->nkz,         0);
    CTYPE ("EWALD/PME/PPPM parameters");
    ITYPE ("pme-order",   ir->pme_order,   4);
    RTYPE ("ewald-rtol",  ir->ewald_rtol, 0.00001);
    RTYPE ("ewald-rtol-lj", ir->ewald_rtol_lj, 0.001);
    EETYPE("lj-pme-comb-rule", ir->ljpme_combination_rule, eljpme_names);
    EETYPE("ewald-geometry", ir->ewald_geometry, eewg_names);
    RTYPE ("epsilon-surface", ir->epsilon_surface, 0.0);

    CCTYPE("IMPLICIT SOLVENT ALGORITHM");
    EETYPE("implicit-solvent", ir->implicit_solvent, eis_names);

    CCTYPE ("GENERALIZED BORN ELECTROSTATICS");
    CTYPE ("Algorithm for calculating Born radii");
    EETYPE("gb-algorithm", ir->gb_algorithm, egb_names);
    CTYPE ("Frequency of calculating the Born radii inside rlist");
    ITYPE ("nstgbradii", ir->nstgbradii, 1);
    CTYPE ("Cutoff for Born radii calculation; the contribution from atoms");
    CTYPE ("between rlist and rgbradii is updated every nstlist steps");
    RTYPE ("rgbradii",  ir->rgbradii, 1.0);
    CTYPE ("Dielectric coefficient of the implicit solvent");
    RTYPE ("gb-epsilon-solvent", ir->gb_epsilon_solvent, 80.0);
    CTYPE ("Salt concentration in M for Generalized Born models");
    RTYPE ("gb-saltconc",  ir->gb_saltconc, 0.0);
    CTYPE ("Scaling factors used in the OBC GB model. Default values are OBC(II)");
    RTYPE ("gb-obc-alpha", ir->gb_obc_alpha, 1.0);
    RTYPE ("gb-obc-beta", ir->gb_obc_beta, 0.8);
    RTYPE ("gb-obc-gamma", ir->gb_obc_gamma, 4.85);
    RTYPE ("gb-dielectric-offset", ir->gb_dielectric_offset, 0.009);
    EETYPE("sa-algorithm", ir->sa_algorithm, esa_names);
    CTYPE ("Surface tension (kJ/mol/nm^2) for the SA (nonpolar surface) part of GBSA");
    CTYPE ("The value -1 will set default value for Still/HCT/OBC GB-models.");
    RTYPE ("sa-surface-tension", ir->sa_surface_tension, -1);

    /* Coupling stuff */
    CCTYPE ("OPTIONS FOR WEAK COUPLING ALGORITHMS");
    CTYPE ("Temperature coupling");
    EETYPE("tcoupl",  ir->etc,        etcoupl_names);
    ITYPE ("nsttcouple", ir->nsttcouple,  -1);
    ITYPE("nh-chain-length",     ir->opts.nhchainlength, 10);
    EETYPE("print-nose-hoover-chain-variables", ir->bPrintNHChains, yesno_names);
    CTYPE ("Groups to couple separately");
    STYPE ("tc-grps",     is->tcgrps,         nullptr);
    CTYPE ("Time constant (ps) and reference temperature (K)");
    STYPE ("tau-t",   is->tau_t,      nullptr);
    STYPE ("ref-t",   is->ref_t,      nullptr);
    CTYPE ("pressure coupling");
    EETYPE("pcoupl",  ir->epc,        epcoupl_names);
    EETYPE("pcoupltype",  ir->epct,       epcoupltype_names);
    ITYPE ("nstpcouple", ir->nstpcouple,  -1);
    CTYPE ("Time constant (ps), compressibility (1/bar) and reference P (bar)");
    RTYPE ("tau-p",   ir->tau_p,  1.0);
    STYPE ("compressibility", dumstr[0],  nullptr);
    STYPE ("ref-p",       dumstr[1],      nullptr);
    CTYPE ("Scaling of reference coordinates, No, All or COM");
    EETYPE ("refcoord-scaling", ir->refcoord_scaling, erefscaling_names);

    /* QMMM */
    CCTYPE ("OPTIONS FOR QMMM calculations");
    EETYPE("QMMM", ir->bQMMM, yesno_names);
    CTYPE ("Groups treated Quantum Mechanically");
    STYPE ("QMMM-grps",  is->QMMM,          nullptr);
    CTYPE ("QM method");
    STYPE("QMmethod",     is->QMmethod, nullptr);
    CTYPE ("QMMM scheme");
    EETYPE("QMMMscheme",  ir->QMMMscheme,    eQMMMscheme_names);
    CTYPE ("QM basisset");
    STYPE("QMbasis",      is->QMbasis, nullptr);
    CTYPE ("QM charge");
    STYPE ("QMcharge",    is->QMcharge, nullptr);
    CTYPE ("QM multiplicity");
    STYPE ("QMmult",      is->QMmult, nullptr);
    CTYPE ("Surface Hopping");
    STYPE ("SH",          is->bSH, nullptr);
    CTYPE ("CAS space options");
    STYPE ("CASorbitals",      is->CASorbitals,   nullptr);
    STYPE ("CASelectrons",     is->CASelectrons,  nullptr);
    STYPE ("SAon", is->SAon, nullptr);
    STYPE ("SAoff", is->SAoff, nullptr);
    STYPE ("SAsteps", is->SAsteps, nullptr);
    CTYPE ("Scale factor for MM charges");
    RTYPE ("MMChargeScaleFactor", ir->scalefactor, 1.0);

    /* Simulated annealing */
    CCTYPE("SIMULATED ANNEALING");
    CTYPE ("Type of annealing for each temperature group (no/single/periodic)");
    STYPE ("annealing",   is->anneal,      nullptr);
    CTYPE ("Number of time points to use for specifying annealing in each group");
    STYPE ("annealing-npoints", is->anneal_npoints, nullptr);
    CTYPE ("List of times at the annealing points for each group");
    STYPE ("annealing-time",       is->anneal_time,       nullptr);
    CTYPE ("Temp. at each annealing point, for each group.");
    STYPE ("annealing-temp",  is->anneal_temp,  nullptr);

    /* Startup run */
    CCTYPE ("GENERATE VELOCITIES FOR STARTUP RUN");
    EETYPE("gen-vel",     opts->bGenVel,  yesno_names);
    RTYPE ("gen-temp",    opts->tempi,    300.0);
    ITYPE ("gen-seed",    opts->seed,     -1);

    /* Shake stuff */
    CCTYPE ("OPTIONS FOR BONDS");
    EETYPE("constraints", opts->nshake,   constraints);
    CTYPE ("Type of constraint algorithm");
    EETYPE("constraint-algorithm",  ir->eConstrAlg, econstr_names);
    CTYPE ("Do not constrain the start configuration");
    EETYPE("continuation", ir->bContinuation, yesno_names);
    CTYPE ("Use successive overrelaxation to reduce the number of shake iterations");
    EETYPE("Shake-SOR", ir->bShakeSOR, yesno_names);
    CTYPE ("Relative tolerance of shake");
    RTYPE ("shake-tol", ir->shake_tol, 0.0001);
    CTYPE ("Highest order in the expansion of the constraint coupling matrix");
    ITYPE ("lincs-order", ir->nProjOrder, 4);
    CTYPE ("Number of iterations in the final step of LINCS. 1 is fine for");
    CTYPE ("normal simulations, but use 2 to conserve energy in NVE runs.");
    CTYPE ("For energy minimization with constraints it should be 4 to 8.");
    ITYPE ("lincs-iter", ir->nLincsIter, 1);
    CTYPE ("Lincs will write a warning to the stderr if in one step a bond");
    CTYPE ("rotates over more degrees than");
    RTYPE ("lincs-warnangle", ir->LincsWarnAngle, 30.0);
    CTYPE ("Convert harmonic bonds to morse potentials");
    EETYPE("morse",       opts->bMorse, yesno_names);

    /* Energy group exclusions */
    CCTYPE ("ENERGY GROUP EXCLUSIONS");
    CTYPE ("Pairs of energy groups for which all non-bonded interactions are excluded");
    STYPE ("energygrp-excl", is->egpexcl,     nullptr);

    /* Walls */
    CCTYPE ("WALLS");
    CTYPE ("Number of walls, type, atom types, densities and box-z scale factor for Ewald");
    ITYPE ("nwall", ir->nwall, 0);
    EETYPE("wall-type",     ir->wall_type,   ewt_names);
    RTYPE ("wall-r-linpot", ir->wall_r_linpot, -1);
    STYPE ("wall-atomtype", is->wall_atomtype, nullptr);
    STYPE ("wall-density",  is->wall_density,  nullptr);
    RTYPE ("wall-ewald-zfac", ir->wall_ewald_zfac, 3);

    /* COM pulling */
    CCTYPE("COM PULLING");
    EETYPE("pull",          ir->bPull, yesno_names);
    if (ir->bPull)
    {
        snew(ir->pull, 1);
        is->pull_grp = read_pullparams(&ninp, &inp, ir->pull, wi);
    }

    /* AWH biasing
       NOTE: needs COM pulling input */
    CCTYPE("AWH biasing");
    EETYPE("awh", ir->bDoAwh, yesno_names);
    if (ir->bDoAwh)
    {
        if (ir->bPull)
        {
            ir->awhParams = gmx::readAndCheckAwhParams(&ninp, &inp, ir, wi);
        }
        else
        {
            gmx_fatal(FARGS, "AWH biasing is only compatible with COM pulling turned on");
        }
    }

    /* Enforced rotation */
    CCTYPE("ENFORCED ROTATION");
    CTYPE("Enforced rotation: No or Yes");
    EETYPE("rotation",       ir->bRot, yesno_names);
    if (ir->bRot)
    {
        snew(ir->rot, 1);
        is->rot_grp = read_rotparams(&ninp, &inp, ir->rot, wi);
    }

    /* Interactive MD */
    ir->bIMD = FALSE;
    CCTYPE("Group to display and/or manipulate in interactive MD session");
    STYPE ("IMD-group", is->imd_grp, nullptr);
    if (is->imd_grp[0] != '\0')
    {
        snew(ir->imd, 1);
        ir->bIMD = TRUE;
    }

    /* Refinement */
    CCTYPE("NMR refinement stuff");
    CTYPE ("Distance restraints type: No, Simple or Ensemble");
    EETYPE("disre",       ir->eDisre,     edisre_names);
    CTYPE ("Force weighting of pairs in one distance restraint: Conservative or Equal");
    EETYPE("disre-weighting", ir->eDisreWeighting, edisreweighting_names);
    CTYPE ("Use sqrt of the time averaged times the instantaneous violation");
    EETYPE("disre-mixed", ir->bDisreMixed, yesno_names);
    RTYPE ("disre-fc",    ir->dr_fc,  1000.0);
    RTYPE ("disre-tau",   ir->dr_tau, 0.0);
    CTYPE ("Output frequency for pair distances to energy file");
    ITYPE ("nstdisreout", ir->nstdisreout, 100);
    CTYPE ("Orientation restraints: No or Yes");
    EETYPE("orire",       opts->bOrire,   yesno_names);
    CTYPE ("Orientation restraints force constant and tau for time averaging");
    RTYPE ("orire-fc",    ir->orires_fc,  0.0);
    RTYPE ("orire-tau",   ir->orires_tau, 0.0);
    STYPE ("orire-fitgrp", is->orirefitgrp,    nullptr);
    CTYPE ("Output frequency for trace(SD) and S to energy file");
    ITYPE ("nstorireout", ir->nstorireout, 100);

    /* free energy variables */
    CCTYPE ("Free energy variables");
    EETYPE("free-energy", ir->efep, efep_names);
    STYPE ("couple-moltype",  is->couple_moltype,  nullptr);
    EETYPE("couple-lambda0", opts->couple_lam0, couple_lam);
    EETYPE("couple-lambda1", opts->couple_lam1, couple_lam);
    EETYPE("couple-intramol", opts->bCoupleIntra, yesno_names);

    RTYPE ("init-lambda", fep->init_lambda, -1); /* start with -1 so
                                                    we can recognize if
                                                    it was not entered */
    ITYPE ("init-lambda-state", fep->init_fep_state, -1);
    RTYPE ("delta-lambda", fep->delta_lambda, 0.0);
    ITYPE ("nstdhdl", fep->nstdhdl, 50);
    STYPE ("fep-lambdas", is->fep_lambda[efptFEP], nullptr);
    STYPE ("mass-lambdas", is->fep_lambda[efptMASS], nullptr);
    STYPE ("coul-lambdas", is->fep_lambda[efptCOUL], nullptr);
    STYPE ("vdw-lambdas", is->fep_lambda[efptVDW], nullptr);
    STYPE ("bonded-lambdas", is->fep_lambda[efptBONDED], nullptr);
    STYPE ("restraint-lambdas", is->fep_lambda[efptRESTRAINT], nullptr);
    STYPE ("temperature-lambdas", is->fep_lambda[efptTEMPERATURE], nullptr);
    ITYPE ("calc-lambda-neighbors", fep->lambda_neighbors, 1);
    STYPE ("init-lambda-weights", is->lambda_weights, nullptr);
    EETYPE("dhdl-print-energy", fep->edHdLPrintEnergy, edHdLPrintEnergy_names);
    RTYPE ("sc-alpha", fep->sc_alpha, 0.0);
    ITYPE ("sc-power", fep->sc_power, 1);
    RTYPE ("sc-r-power", fep->sc_r_power, 6.0);
    RTYPE ("sc-sigma", fep->sc_sigma, 0.3);
    EETYPE("sc-coul", fep->bScCoul, yesno_names);
    ITYPE ("dh_hist_size", fep->dh_hist_size, 0);
    RTYPE ("dh_hist_spacing", fep->dh_hist_spacing, 0.1);
    EETYPE("separate-dhdl-file", fep->separate_dhdl_file,
           separate_dhdl_file_names);
    EETYPE("dhdl-derivatives", fep->dhdl_derivatives, dhdl_derivatives_names);
    ITYPE ("dh_hist_size", fep->dh_hist_size, 0);
    RTYPE ("dh_hist_spacing", fep->dh_hist_spacing, 0.1);

    /* Non-equilibrium MD stuff */
    CCTYPE("Non-equilibrium MD stuff");
    STYPE ("acc-grps",    is->accgrps,        nullptr);
    STYPE ("accelerate",  is->acc,            nullptr);
    STYPE ("freezegrps",  is->freeze,         nullptr);
    STYPE ("freezedim",   is->frdim,          nullptr);
    RTYPE ("cos-acceleration", ir->cos_accel, 0);
    STYPE ("deform",      is->deform,         nullptr);

    /* simulated tempering variables */
    CCTYPE("simulated tempering variables");
    EETYPE("simulated-tempering", ir->bSimTemp, yesno_names);
    EETYPE("simulated-tempering-scaling", ir->simtempvals->eSimTempScale, esimtemp_names);
    RTYPE("sim-temp-low", ir->simtempvals->simtemp_low, 300.0);
    RTYPE("sim-temp-high", ir->simtempvals->simtemp_high, 300.0);

    /* expanded ensemble variables */
    if (ir->efep == efepEXPANDED || ir->bSimTemp)
    {
        read_expandedparams(&ninp, &inp, expand, wi);
    }

    /* Electric fields */
    {
        gmx::KeyValueTreeObject      convertedValues = flatKeyValueTreeFromInpFile(ninp, inp);
        gmx::KeyValueTreeTransformer transform;
        transform.rules()->addRule()
            .keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
        mdModules->initMdpTransform(transform.rules());
        for (const auto &path : transform.mappedPaths())
        {
            GMX_ASSERT(path.size() == 1, "Inconsistent mapping back to mdp options");
            mark_einp_set(ninp, inp, path[0].c_str());
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
    CCTYPE("Ion/water position swapping for computational electrophysiology setups");
    CTYPE("Swap positions along direction: no, X, Y, Z");
    EETYPE("swapcoords", ir->eSwapCoords, eSwapTypes_names);
    if (ir->eSwapCoords != eswapNO)
    {
        char buf[STRLEN];
        int  nIonTypes;


        snew(ir->swap, 1);
        CTYPE("Swap attempt frequency");
        ITYPE("swap-frequency", ir->swap->nstswap, 1);
        CTYPE("Number of ion types to be controlled");
        ITYPE("iontypes", nIonTypes, 1);
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
        CTYPE("Two index groups that contain the compartment-partitioning atoms");
        STYPE("split-group0", ir->swap->grp[eGrpSplit0].molname, nullptr);
        STYPE("split-group1", ir->swap->grp[eGrpSplit1].molname, nullptr);
        CTYPE("Use center of mass of split groups (yes/no), otherwise center of geometry is used");
        EETYPE("massw-split0", ir->swap->massw_split[0], yesno_names);
        EETYPE("massw-split1", ir->swap->massw_split[1], yesno_names);

        CTYPE("Name of solvent molecules");
        STYPE("solvent-group", ir->swap->grp[eGrpSolvent].molname, nullptr);

        CTYPE("Split cylinder: radius, upper and lower extension (nm) (this will define the channels)");
        CTYPE("Note that the split cylinder settings do not have an influence on the swapping protocol,");
        CTYPE("however, if correctly defined, the permeation events are recorded per channel");
        RTYPE("cyl0-r", ir->swap->cyl0r, 2.0);
        RTYPE("cyl0-up", ir->swap->cyl0u, 1.0);
        RTYPE("cyl0-down", ir->swap->cyl0l, 1.0);
        RTYPE("cyl1-r", ir->swap->cyl1r, 2.0);
        RTYPE("cyl1-up", ir->swap->cyl1u, 1.0);
        RTYPE("cyl1-down", ir->swap->cyl1l, 1.0);

        CTYPE("Average the number of ions per compartment over these many swap attempt steps");
        ITYPE("coupl-steps", ir->swap->nAverage, 10);

        CTYPE("Names of the ion types that can be exchanged with solvent molecules,");
        CTYPE("and the requested number of ions of this type in compartments A and B");
        CTYPE("-1 means fix the numbers as found in step 0");
        for (i = 0; i < nIonTypes; i++)
        {
            int ig = eSwapFixedGrpNR + i;

            sprintf(buf, "iontype%d-name", i);
            STYPE(buf, ir->swap->grp[ig].molname, nullptr);
            sprintf(buf, "iontype%d-in-A", i);
            ITYPE(buf, ir->swap->grp[ig].nmolReq[0], -1);
            sprintf(buf, "iontype%d-in-B", i);
            ITYPE(buf, ir->swap->grp[ig].nmolReq[1], -1);
        }

        CTYPE("By default (i.e. bulk offset = 0.0), ion/water exchanges happen between layers");
        CTYPE("at maximum distance (= bulk concentration) to the split group layers. However,");
        CTYPE("an offset b (-1.0 < b < +1.0) can be specified to offset the bulk layer from the middle at 0.0");
        CTYPE("towards one of the compartment-partitioning layers (at +/- 1.0).");
        RTYPE("bulk-offsetA", ir->swap->bulkOffset[0], 0.0);
        RTYPE("bulk-offsetB", ir->swap->bulkOffset[1], 0.0);
        if (!(ir->swap->bulkOffset[0] > -1.0 && ir->swap->bulkOffset[0] < 1.0)
            || !(ir->swap->bulkOffset[1] > -1.0 && ir->swap->bulkOffset[1] < 1.0) )
        {
            warning_error(wi, "Bulk layer offsets must be > -1.0 and < 1.0 !");
        }

        CTYPE("Start to swap ions if threshold difference to requested count is reached");
        RTYPE("threshold", ir->swap->threshold, 1.0);
    }

    /* AdResS is no longer supported, but we need mdrun to be able to refuse to run old AdResS .tpr files */
    EETYPE("adress", ir->bAdress, yesno_names);

    /* User defined thingies */
    CCTYPE ("User defined thingies");
    STYPE ("user1-grps",  is->user1,          nullptr);
    STYPE ("user2-grps",  is->user2,          nullptr);
    ITYPE ("userint1",    ir->userint1,   0);
    ITYPE ("userint2",    ir->userint2,   0);
    ITYPE ("userint3",    ir->userint3,   0);
    ITYPE ("userint4",    ir->userint4,   0);
    RTYPE ("userreal1",   ir->userreal1,  0);
    RTYPE ("userreal2",   ir->userreal2,  0);
    RTYPE ("userreal3",   ir->userreal3,  0);
    RTYPE ("userreal4",   ir->userreal4,  0);
#undef CTYPE

    {
        gmx::TextOutputFile stream(mdparout);
        write_inpfile(&stream, mdparout, ninp, inp, FALSE, writeMdpHeader, wi);

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

    for (i = 0; (i < ninp); i++)
    {
        sfree(inp[i].name);
        sfree(inp[i].value);
    }
    sfree(inp);

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

    do_wall_params(ir, is->wall_atomtype, is->wall_density, opts);

    /* ORIENTATION RESTRAINT PARAMETERS */

    if (opts->bOrire && str_nelem(is->orirefitgrp, MAXPTR, nullptr) != 1)
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

    return -1;

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

    return -1;
}

static gmx_bool do_numbering(int natoms, gmx_groups_t *groups, int ng, char *ptrs[],
                             t_blocka *block, char *gnames[],
                             int gtype, int restnm,
                             int grptp, gmx_bool bVerbose,
                             warninp_t wi)
{
    unsigned short *cbuf;
    t_grps         *grps = &(groups->grps[gtype]);
    int             i, j, gid, aj, ognr, ntot = 0;
    const char     *title;
    gmx_bool        bRest;
    char            warn_buf[STRLEN];

    if (debug)
    {
        fprintf(debug, "Starting numbering %d groups of type %d\n", ng, gtype);
    }

    title = gtypes[gtype];

    snew(cbuf, natoms);
    /* Mark all id's as not set */
    for (i = 0; (i < natoms); i++)
    {
        cbuf[i] = NOGID;
    }

    snew(grps->nm_ind, ng+1); /* +1 for possible rest group */
    for (i = 0; (i < ng); i++)
    {
        /* Lookup the group name in the block structure */
        gid = search_string(ptrs[i], block->nr, gnames);
        if ((grptp != egrptpONE) || (i == 0))
        {
            grps->nm_ind[grps->nr++] = gid;
        }
        if (debug)
        {
            fprintf(debug, "Found gid %d for group %s\n", gid, ptrs[i]);
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

static void calc_nrdf(gmx_mtop_t *mtop, t_inputrec *ir, char **gnames)
{
    t_grpopts              *opts;
    gmx_groups_t           *groups;
    pull_params_t          *pull;
    int                     natoms, ai, aj, i, j, d, g, imin, jmin;
    t_iatom                *ia;
    int                    *nrdf2, *na_vcm, na_tot;
    double                 *nrdf_tc, *nrdf_vcm, nrdf_uc, *nrdf_vcm_sub;
    ivec                   *dof_vcm;
    gmx_mtop_atomloop_all_t aloop;
    int                     mb, mol, ftype, as;
    gmx_molblock_t         *molb;
    gmx_moltype_t          *molt;

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
            g = ggrpnr(groups, egcFREEZE, i);
            for (d = 0; d < DIM; d++)
            {
                if (opts->nFreeze[g][d] == 0)
                {
                    /* Add one DOF for particle i (counted as 2*1) */
                    nrdf2[i]                              += 2;
                    /* VCM group i has dim d as a DOF */
                    dof_vcm[ggrpnr(groups, egcVCM, i)][d]  = 1;
                }
            }
            nrdf_tc [ggrpnr(groups, egcTC, i)]  += 0.5*nrdf2[i];
            nrdf_vcm[ggrpnr(groups, egcVCM, i)] += 0.5*nrdf2[i];
        }
    }

    as = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        atom = molt->atoms.atom;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            for (ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
            {
                ia = molt->ilist[ftype].iatoms;
                for (i = 0; i < molt->ilist[ftype].nr; )
                {
                    /* Subtract degrees of freedom for the constraints,
                     * if the particles still have degrees of freedom left.
                     * If one of the particles is a vsite or a shell, then all
                     * constraint motion will go there, but since they do not
                     * contribute to the constraints the degrees of freedom do not
                     * change.
                     */
                    ai = as + ia[1];
                    aj = as + ia[2];
                    if (((atom[ia[1]].ptype == eptNucleus) ||
                         (atom[ia[1]].ptype == eptAtom)) &&
                        ((atom[ia[2]].ptype == eptNucleus) ||
                         (atom[ia[2]].ptype == eptAtom)))
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
                        nrdf_tc [ggrpnr(groups, egcTC, ai)]  -= 0.5*imin;
                        nrdf_tc [ggrpnr(groups, egcTC, aj)]  -= 0.5*jmin;
                        nrdf_vcm[ggrpnr(groups, egcVCM, ai)] -= 0.5*imin;
                        nrdf_vcm[ggrpnr(groups, egcVCM, aj)] -= 0.5*jmin;
                    }
                    ia += interaction_function[ftype].nratoms+1;
                    i  += interaction_function[ftype].nratoms+1;
                }
            }
            ia = molt->ilist[F_SETTLE].iatoms;
            for (i = 0; i < molt->ilist[F_SETTLE].nr; )
            {
                /* Subtract 1 dof from every atom in the SETTLE */
                for (j = 0; j < 3; j++)
                {
                    ai         = as + ia[1+j];
                    imin       = std::min(2, nrdf2[ai]);
                    nrdf2[ai] -= imin;
                    nrdf_tc [ggrpnr(groups, egcTC, ai)]  -= 0.5*imin;
                    nrdf_vcm[ggrpnr(groups, egcVCM, ai)] -= 0.5*imin;
                }
                ia += 4;
                i  += 4;
            }
            as += molt->atoms.nr;
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
                    nrdf_tc [ggrpnr(groups, egcTC, ai)]  -= 0.5*imin;
                    nrdf_vcm[ggrpnr(groups, egcVCM, ai)] -= 0.5*imin;
                    if (nrdf_tc[ggrpnr(groups, egcTC, ai)] < 0)
                    {
                        gmx_fatal(FARGS, "Center of mass pulling constraints caused the number of degrees of freedom for temperature coupling group %s to be negative", gnames[groups->grps[egcTC].nm_ind[ggrpnr(groups, egcTC, ai)]]);
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
                if (ggrpnr(groups, egcTC, ai) == i)
                {
                    na_vcm[ggrpnr(groups, egcVCM, ai)]++;
                    na_tot++;
                }
            }
            /* Correct for VCM removal according to the fraction of each VCM
             * group present in this TC group.
             */
            nrdf_uc = nrdf_tc[i];
            if (debug)
            {
                fprintf(debug, "T-group[%d] nrdf_uc = %g\n", i, nrdf_uc);
            }
            nrdf_tc[i] = 0;
            for (j = 0; j < groups->grps[egcVCM].nr+1; j++)
            {
                if (nrdf_vcm[j] > nrdf_vcm_sub[j])
                {
                    nrdf_tc[i] += nrdf_uc*((double)na_vcm[j]/(double)na_tot)*
                        (nrdf_vcm[j] - nrdf_vcm_sub[j])/nrdf_vcm[j];
                }
                if (debug)
                {
                    fprintf(debug, "  nrdf_vcm[%d] = %g, nrdf = %g\n",
                            j, nrdf_vcm[j], nrdf_tc[i]);
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

static gmx_bool do_egp_flag(t_inputrec *ir, gmx_groups_t *groups,
                            const char *option, const char *val, int flag)
{
    /* The maximum number of energy group pairs would be MAXPTR*(MAXPTR+1)/2.
     * But since this is much larger than STRLEN, such a line can not be parsed.
     * The real maximum is the number of names that fit in a string: STRLEN/2.
     */
#define EGP_MAX (STRLEN/2)
    int      nelem, i, j, k, nr;
    char    *names[EGP_MAX];
    char  ***gnames;
    gmx_bool bSet;

    gnames = groups->grpname;

    nelem = str_nelem(val, EGP_MAX, names);
    if (nelem % 2 != 0)
    {
        gmx_fatal(FARGS, "The number of groups for %s is odd", option);
    }
    nr   = groups->grps[egcENER].nr;
    bSet = FALSE;
    for (i = 0; i < nelem/2; i++)
    {
        j = 0;
        while ((j < nr) &&
               gmx_strcasecmp(names[2*i], *(gnames[groups->grps[egcENER].nm_ind[j]])))
        {
            j++;
        }
        if (j == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n",
                      names[2*i], option);
        }
        k = 0;
        while ((k < nr) &&
               gmx_strcasecmp(names[2*i+1], *(gnames[groups->grps[egcENER].nm_ind[k]])))
        {
            k++;
        }
        if (k == nr)
        {
            gmx_fatal(FARGS, "%s in %s is not an energy group\n",
                      names[2*i+1], option);
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
              gmx_bool bVerbose,
              t_inputrec *ir,
              warninp_t wi)
{
    t_blocka     *grps;
    gmx_groups_t *groups;
    int           natoms;
    t_symtab     *symtab;
    t_atoms       atoms_all;
    char          warnbuf[STRLEN], **gnames;
    int           nr, ntcg, ntau_t, nref_t, nacc, nofg, nSA, nSA_points, nSA_time, nSA_temp;
    real          tau_min;
    int           nstcmin;
    int           nacg, nfreeze, nfrdim, nenergy, nvcm, nuser;
    char         *ptr1[MAXPTR], *ptr2[MAXPTR], *ptr3[MAXPTR];
    int           i, j, k, restnm;
    gmx_bool      bExcl, bTable, bSetTCpar, bAnneal, bRest;
    int           nQMmethod, nQMbasis, nQMg;
    char          warn_buf[STRLEN];
    char*         endptr;

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

    ntau_t = str_nelem(is->tau_t, MAXPTR, ptr1);
    nref_t = str_nelem(is->ref_t, MAXPTR, ptr2);
    ntcg   = str_nelem(is->tcgrps, MAXPTR, ptr3);
    if ((ntau_t != ntcg) || (nref_t != ntcg))
    {
        gmx_fatal(FARGS, "Invalid T coupling input: %d groups, %d ref-t values and "
                  "%d tau-t values", ntcg, nref_t, ntau_t);
    }

    bSetTCpar = (ir->etc || EI_SD(ir->eI) || ir->eI == eiBD || EI_TPI(ir->eI));
    do_numbering(natoms, groups, ntcg, ptr3, grps, gnames, egcTC,
                 restnm, bSetTCpar ? egrptpALL : egrptpALL_GENREST, bVerbose, wi);
    nr            = groups->grps[egcTC].nr;
    ir->opts.ngtc = nr;
    snew(ir->opts.nrdf, nr);
    snew(ir->opts.tau_t, nr);
    snew(ir->opts.ref_t, nr);
    if (ir->eI == eiBD && ir->bd_fric == 0)
    {
        fprintf(stderr, "bd-fric=0, so tau-t will be used as the inverse friction constant(s)\n");
    }

    if (bSetTCpar)
    {
        if (nr != nref_t)
        {
            gmx_fatal(FARGS, "Not enough ref-t and tau-t values!");
        }

        tau_min = 1e20;
        for (i = 0; (i < nr); i++)
        {
            ir->opts.tau_t[i] = strtod(ptr1[i], &endptr);
            if (*endptr != 0)
            {
                warning_error(wi, "Invalid value for mdp option tau-t. tau-t should only consist of real numbers separated by spaces.");
            }
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
        for (i = 0; (i < nr); i++)
        {
            ir->opts.ref_t[i] = strtod(ptr2[i], &endptr);
            if (*endptr != 0)
            {
                warning_error(wi, "Invalid value for mdp option ref-t. ref-t should only consist of real numbers separated by spaces.");
            }
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
    nSA = str_nelem(is->anneal, MAXPTR, ptr1);
    if (nSA == 1 && (ptr1[0][0] == 'n' || ptr1[0][0] == 'N'))
    {
        nSA = 0;
    }
    if (nSA > 0 && nSA != nr)
    {
        gmx_fatal(FARGS, "Not enough annealing values: %d (for %d groups)\n", nSA, nr);
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
        if (nSA > 0)
        {
            bAnneal = FALSE;
            for (i = 0; i < nr; i++)
            {
                if (ptr1[i][0] == 'n' || ptr1[i][0] == 'N')
                {
                    ir->opts.annealing[i] = eannNO;
                }
                else if (ptr1[i][0] == 's' || ptr1[i][0] == 'S')
                {
                    ir->opts.annealing[i] = eannSINGLE;
                    bAnneal               = TRUE;
                }
                else if (ptr1[i][0] == 'p' || ptr1[i][0] == 'P')
                {
                    ir->opts.annealing[i] = eannPERIODIC;
                    bAnneal               = TRUE;
                }
            }
            if (bAnneal)
            {
                /* Read the other fields too */
                nSA_points = str_nelem(is->anneal_npoints, MAXPTR, ptr1);
                if (nSA_points != nSA)
                {
                    gmx_fatal(FARGS, "Found %d annealing-npoints values for %d groups\n", nSA_points, nSA);
                }
                for (k = 0, i = 0; i < nr; i++)
                {
                    ir->opts.anneal_npoints[i] = strtol(ptr1[i], &endptr, 10);
                    if (*endptr != 0)
                    {
                        warning_error(wi, "Invalid value for mdp option annealing-npoints. annealing should only consist of integers separated by spaces.");
                    }
                    if (ir->opts.anneal_npoints[i] == 1)
                    {
                        gmx_fatal(FARGS, "Please specify at least a start and an end point for annealing\n");
                    }
                    snew(ir->opts.anneal_time[i], ir->opts.anneal_npoints[i]);
                    snew(ir->opts.anneal_temp[i], ir->opts.anneal_npoints[i]);
                    k += ir->opts.anneal_npoints[i];
                }

                nSA_time = str_nelem(is->anneal_time, MAXPTR, ptr1);
                if (nSA_time != k)
                {
                    gmx_fatal(FARGS, "Found %d annealing-time values, wanted %d\n", nSA_time, k);
                }
                nSA_temp = str_nelem(is->anneal_temp, MAXPTR, ptr2);
                if (nSA_temp != k)
                {
                    gmx_fatal(FARGS, "Found %d annealing-temp values, wanted %d\n", nSA_temp, k);
                }

                for (i = 0, k = 0; i < nr; i++)
                {

                    for (j = 0; j < ir->opts.anneal_npoints[i]; j++)
                    {
                        ir->opts.anneal_time[i][j] = strtod(ptr1[k], &endptr);
                        if (*endptr != 0)
                        {
                            warning_error(wi, "Invalid value for mdp option anneal-time. anneal-time should only consist of real numbers separated by spaces.");
                        }
                        ir->opts.anneal_temp[i][j] = strtod(ptr2[k], &endptr);
                        if (*endptr != 0)
                        {
                            warning_error(wi, "Invalid value for anneal-temp. anneal-temp should only consist of real numbers separated by spaces.");
                        }
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
                            if (fabs(ir->opts.anneal_temp[i][j]-ir->opts.anneal_temp[i][0]) > GMX_REAL_EPS)
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

    nacc = str_nelem(is->acc, MAXPTR, ptr1);
    nacg = str_nelem(is->accgrps, MAXPTR, ptr2);
    if (nacg*DIM != nacc)
    {
        gmx_fatal(FARGS, "Invalid Acceleration input: %d groups and %d acc. values",
                  nacg, nacc);
    }
    do_numbering(natoms, groups, nacg, ptr2, grps, gnames, egcACC,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nr = groups->grps[egcACC].nr;
    snew(ir->opts.acc, nr);
    ir->opts.ngacc = nr;

    for (i = k = 0; (i < nacg); i++)
    {
        for (j = 0; (j < DIM); j++, k++)
        {
            ir->opts.acc[i][j] = strtod(ptr1[k], &endptr);
            if (*endptr != 0)
            {
                warning_error(wi, "Invalid value for mdp option accelerate. accelerate should only consist of real numbers separated by spaces.");
            }
        }
    }
    for (; (i < nr); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            ir->opts.acc[i][j] = 0;
        }
    }

    nfrdim  = str_nelem(is->frdim, MAXPTR, ptr1);
    nfreeze = str_nelem(is->freeze, MAXPTR, ptr2);
    if (nfrdim != DIM*nfreeze)
    {
        gmx_fatal(FARGS, "Invalid Freezing input: %d groups and %d freeze values",
                  nfreeze, nfrdim);
    }
    do_numbering(natoms, groups, nfreeze, ptr2, grps, gnames, egcFREEZE,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nr             = groups->grps[egcFREEZE].nr;
    ir->opts.ngfrz = nr;
    snew(ir->opts.nFreeze, nr);
    for (i = k = 0; (i < nfreeze); i++)
    {
        for (j = 0; (j < DIM); j++, k++)
        {
            ir->opts.nFreeze[i][j] = (gmx_strncasecmp(ptr1[k], "Y", 1) == 0);
            if (!ir->opts.nFreeze[i][j])
            {
                if (gmx_strncasecmp(ptr1[k], "N", 1) != 0)
                {
                    sprintf(warnbuf, "Please use Y(ES) or N(O) for freezedim only "
                            "(not %s)", ptr1[k]);
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

    nenergy = str_nelem(is->energy, MAXPTR, ptr1);
    do_numbering(natoms, groups, nenergy, ptr1, grps, gnames, egcENER,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    add_wall_energrps(groups, ir->nwall, symtab);
    ir->opts.ngener = groups->grps[egcENER].nr;
    nvcm            = str_nelem(is->vcm, MAXPTR, ptr1);
    bRest           =
        do_numbering(natoms, groups, nvcm, ptr1, grps, gnames, egcVCM,
                     restnm, nvcm == 0 ? egrptpALL_GENREST : egrptpPART, bVerbose, wi);
    if (bRest)
    {
        warning(wi, "Some atoms are not part of any center of mass motion removal group.\n"
                "This may lead to artifacts.\n"
                "In most cases one should use one group for the whole system.");
    }

    /* Now we have filled the freeze struct, so we can calculate NRDF */
    calc_nrdf(mtop, ir, gnames);

    nuser = str_nelem(is->user1, MAXPTR, ptr1);
    do_numbering(natoms, groups, nuser, ptr1, grps, gnames, egcUser1,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nuser = str_nelem(is->user2, MAXPTR, ptr1);
    do_numbering(natoms, groups, nuser, ptr1, grps, gnames, egcUser2,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nuser = str_nelem(is->x_compressed_groups, MAXPTR, ptr1);
    do_numbering(natoms, groups, nuser, ptr1, grps, gnames, egcCompressedX,
                 restnm, egrptpONE, bVerbose, wi);
    nofg = str_nelem(is->orirefitgrp, MAXPTR, ptr1);
    do_numbering(natoms, groups, nofg, ptr1, grps, gnames, egcORFIT,
                 restnm, egrptpALL_GENREST, bVerbose, wi);

    /* QMMM input processing */
    nQMg          = str_nelem(is->QMMM, MAXPTR, ptr1);
    nQMmethod     = str_nelem(is->QMmethod, MAXPTR, ptr2);
    nQMbasis      = str_nelem(is->QMbasis, MAXPTR, ptr3);
    if ((nQMmethod != nQMg) || (nQMbasis != nQMg))
    {
        gmx_fatal(FARGS, "Invalid QMMM input: %d groups %d basissets"
                  " and %d methods\n", nQMg, nQMbasis, nQMmethod);
    }
    /* group rest, if any, is always MM! */
    do_numbering(natoms, groups, nQMg, ptr1, grps, gnames, egcQMMM,
                 restnm, egrptpALL_GENREST, bVerbose, wi);
    nr            = nQMg; /*atoms->grps[egcQMMM].nr;*/
    ir->opts.ngQM = nQMg;
    snew(ir->opts.QMmethod, nr);
    snew(ir->opts.QMbasis, nr);
    for (i = 0; i < nr; i++)
    {
        /* input consists of strings: RHF CASSCF PM3 .. These need to be
         * converted to the corresponding enum in names.c
         */
        ir->opts.QMmethod[i] = search_QMstring(ptr2[i], eQMmethodNR,
                                               eQMmethod_names);
        ir->opts.QMbasis[i]  = search_QMstring(ptr3[i], eQMbasisNR,
                                               eQMbasis_names);

    }
    str_nelem(is->QMmult, MAXPTR, ptr1);
    str_nelem(is->QMcharge, MAXPTR, ptr2);
    str_nelem(is->bSH, MAXPTR, ptr3);
    snew(ir->opts.QMmult, nr);
    snew(ir->opts.QMcharge, nr);
    snew(ir->opts.bSH, nr);

    for (i = 0; i < nr; i++)
    {
        ir->opts.QMmult[i]   = strtol(ptr1[i], &endptr, 10);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option QMmult. QMmult should only consist of integers separated by spaces.");
        }
        ir->opts.QMcharge[i] = strtol(ptr2[i], &endptr, 10);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option QMcharge. QMcharge should only consist of integers separated by spaces.");
        }
        ir->opts.bSH[i]      = (gmx_strncasecmp(ptr3[i], "Y", 1) == 0);
    }

    str_nelem(is->CASelectrons, MAXPTR, ptr1);
    str_nelem(is->CASorbitals, MAXPTR, ptr2);
    snew(ir->opts.CASelectrons, nr);
    snew(ir->opts.CASorbitals, nr);
    for (i = 0; i < nr; i++)
    {
        ir->opts.CASelectrons[i] = strtol(ptr1[i], &endptr, 10);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option CASelectrons. CASelectrons should only consist of integers separated by spaces.");
        }
        ir->opts.CASorbitals[i]  = strtol(ptr2[i], &endptr, 10);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option CASorbitals. CASorbitals should only consist of integers separated by spaces.");
        }
    }

    str_nelem(is->SAon, MAXPTR, ptr1);
    str_nelem(is->SAoff, MAXPTR, ptr2);
    str_nelem(is->SAsteps, MAXPTR, ptr3);
    snew(ir->opts.SAon, nr);
    snew(ir->opts.SAoff, nr);
    snew(ir->opts.SAsteps, nr);

    for (i = 0; i < nr; i++)
    {
        ir->opts.SAon[i]    = strtod(ptr1[i], &endptr);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option SAon. SAon should only consist of real numbers separated by spaces.");
        }
        ir->opts.SAoff[i]   = strtod(ptr2[i], &endptr);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option SAoff. SAoff should only consist of real numbers separated by spaces.");
        }
        ir->opts.SAsteps[i] = strtol(ptr3[i], &endptr, 10);
        if (*endptr != 0)
        {
            warning_error(wi, "Invalid value for mdp option SAsteps. SAsteps should only consist of integers separated by spaces.");
        }
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



static void check_disre(gmx_mtop_t *mtop)
{
    gmx_ffparams_t *ffparams;
    t_functype     *functype;
    t_iparams      *ip;
    int             i, ndouble, ftype;
    int             label, old_label;

    if (gmx_mtop_ftype_count(mtop, F_DISRES) > 0)
    {
        ffparams  = &mtop->ffparams;
        functype  = ffparams->functype;
        ip        = ffparams->iparams;
        ndouble   = 0;
        old_label = -1;
        for (i = 0; i < ffparams->ntypes; i++)
        {
            ftype = functype[i];
            if (ftype == F_DISRES)
            {
                label = ip[i].disres.label;
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

static gmx_bool absolute_reference(t_inputrec *ir, gmx_mtop_t *sys,
                                   gmx_bool posres_only,
                                   ivec AbsRef)
{
    int                  d, g, i;
    gmx_mtop_ilistloop_t iloop;
    t_ilist             *ilist;
    int                  nmol;
    t_iparams           *pr;

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
    while (gmx_mtop_ilistloop_next(iloop, &ilist, &nmol))
    {
        if (nmol > 0 &&
            (AbsRef[XX] == 0 || AbsRef[YY] == 0 || AbsRef[ZZ] == 0))
        {
            for (i = 0; i < ilist[F_POSRES].nr; i += 2)
            {
                pr = &sys->ffparams.iparams[ilist[F_POSRES].iatoms[i]];
                for (d = 0; d < DIM; d++)
                {
                    if (pr->posres.fcA[d] != 0)
                    {
                        AbsRef[d] = 1;
                    }
                }
            }
            for (i = 0; i < ilist[F_FBPOSRES].nr; i += 2)
            {
                /* Check for flat-bottom posres */
                pr = &sys->ffparams.iparams[ilist[F_FBPOSRES].iatoms[i]];
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
                                   gmx_bool *bC6ParametersWorkWithGeometricRules,
                                   gmx_bool *bC6ParametersWorkWithLBRules,
                                   gmx_bool *bLBRulesPossible)
{
    int           ntypes, tpi, tpj;
    int          *typecount;
    real          tol;
    double        c6i, c6j, c12i, c12j;
    double        c6, c6_geometric, c6_LB;
    double        sigmai, sigmaj, epsi, epsj;
    gmx_bool      bCanDoLBRules, bCanDoGeometricRules;
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

            if (FALSE == bCanDoLBRules)
            {
                *bC6ParametersWorkWithLBRules = FALSE;
            }

            bCanDoGeometricRules = gmx_within_tol(c6_geometric, c6, tol);

            if (FALSE == bCanDoGeometricRules)
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
    gmx_bool bLBRulesPossible, bC6ParametersWorkWithGeometricRules, bC6ParametersWorkWithLBRules;

    check_combination_rule_differences(mtop, 0,
                                       &bC6ParametersWorkWithGeometricRules,
                                       &bC6ParametersWorkWithLBRules,
                                       &bLBRulesPossible);
    if (ir->ljpme_combination_rule == eljpmeLB)
    {
        if (FALSE == bC6ParametersWorkWithLBRules || FALSE == bLBRulesPossible)
        {
            warning(wi, "You are using arithmetic-geometric combination rules "
                    "in LJ-PME, but your non-bonded C6 parameters do not "
                    "follow these rules.");
        }
    }
    else
    {
        if (FALSE == bC6ParametersWorkWithGeometricRules)
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
    gmx_bool                  bCharge, bAcc;
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
                int nsteps = static_cast<int>(ir->opts.tau_t[i]/ir->delta_t + 0.5);
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
        if (ir->coulombtype == eelCUT && ir->rcoulomb > 0 && !ir->implicit_solvent)
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
            mgrp[ggrpnr(&sys->groups, egcACC, i)] += atom->m;
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
        gmx_bool bWarned;

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
                  gmx_bool bHasNormalConstraints,
                  gmx_bool bHasAnyConstraints,
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
