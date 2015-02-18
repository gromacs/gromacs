/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/names.h"

#include "gromacs/legacyheaders/typedefs.h"

/* note: these arrays should correspond to enums in include/types/enums.h */

const char *epbc_names[epbcNR+1] =
{
    "xyz", "no", "xy", "screw", NULL
};

const char *ens_names[ensNR+1] =
{
    "Grid", "Simple", NULL
};

const char *ei_names[eiNR+1] =
{
    "md", "steep", "cg", "bd", "sd2", "nm", "l-bfgs", "tpi", "tpic", "sd", "md-vv", "md-vv-avek", NULL
};

const char *bool_names[BOOL_NR+1] =
{
    "FALSE", "TRUE", NULL
};

const char *yesno_names[BOOL_NR+1] =
{
    "no", "yes", NULL
};

const char *ptype_str[eptNR+1] = {
    "Atom", "Nucleus", "Shell", "Bond", "VSite", NULL
};

const char *ecutscheme_names[ecutsNR+1] = {
    "Verlet", "Group", NULL
};

const char *eel_names[eelNR+1] = {
    "Cut-off", "Reaction-Field", "Generalized-Reaction-Field",
    "PME", "Ewald", "P3M-AD", "Poisson", "Switch", "Shift", "User",
    "Generalized-Born", "Reaction-Field-nec", "Encad-shift",
    "PME-User", "PME-Switch", "PME-User-Switch",
    "Reaction-Field-zero", NULL
};

const char *eewg_names[eewgNR+1] = {
    "3d", "3dc", NULL
};

const char *eljpme_names[eljpmeNR+1] = {
    "Geometric", "Lorentz-Berthelot", NULL
};

const char *evdw_names[evdwNR+1] = {
    "Cut-off", "Switch", "Shift", "User", "Encad-shift",
    "PME", NULL
};

const char *econstr_names[econtNR+1] = {
    "Lincs", "Shake", NULL
};

const char *eintmod_names[eintmodNR+1] = {
    "Potential-shift-Verlet", "Potential-shift", "None", "Potential-switch", "Exact-cutoff", "Force-switch", NULL
};

const char *egrp_nm[egNR+1] = {
    "Coul-SR", "LJ-SR", "Buck-SR", "Coul-LR", "LJ-LR", "Buck-LR",
    "Coul-14", "LJ-14", NULL
};

const char *etcoupl_names[etcNR+1] = {
    "No", "Berendsen", "Nose-Hoover", "yes", "Andersen", "Andersen-massive", "V-rescale", NULL
}; /* yes is alias for berendsen */

const char *epcoupl_names[epcNR+1] = {
    "No", "Berendsen", "Parrinello-Rahman", "Isotropic", "MTTK", NULL
}; /* isotropic is alias for berendsen */

const char *epcoupltype_names[epctNR+1] = {
    "Isotropic", "Semiisotropic", "Anisotropic", "Surface-Tension", NULL
};

const char *erefscaling_names[erscNR+1] = {
    "No", "All", "COM", NULL
};

const char *edisre_names[edrNR+1] = {
    "No", "Simple", "Ensemble", NULL
};

const char *edisreweighting_names[edrwNR+1] = {
    "Conservative", "Equal", NULL
};

const char *enbf_names[eNBF_NR+1] = {
    "", "LJ", "Buckingham", NULL
};

const char *ecomb_names[eCOMB_NR+1] = {
    "", "Geometric", "Arithmetic", "GeomSigEps", NULL
};

const char *gtypes[egcNR+1] = {
    "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
    "User1", "User2", "VCM", "Compressed X", "Or. Res. Fit", "QMMM", NULL
};

const char *esimtemp_names[esimtempNR+1] = {
    "geometric", "exponential", "linear", NULL
};

const char *efep_names[efepNR+1] = {
    "no", "yes", "static", "slow-growth", "expanded", NULL
};

const char *efpt_names[efptNR+1] = {
    "fep-lambdas", "mass-lambdas", "coul-lambdas", "vdw-lambdas", "bonded-lambdas", "restraint-lambdas", "temperature-lambdas", NULL
};

const char *efpt_singular_names[efptNR+1] = {
    "fep-lambda", "mass-lambda", "coul-lambda", "vdw-lambda", "bonded-lambda", "restraint-lambda", "temperature-lambda", NULL
};

const char *edHdLPrintEnergy_names[edHdLPrintEnergyNR+1] = {
    "no", "total", "potential", "yes", NULL
};

const char *elamstats_names[elamstatsNR+1] = {
    "no", "metropolis-transition", "barker-transition", "minvar", "wang-landau", "weighted-wang-landau", NULL
};

const char *elmcmove_names[elmcmoveNR+1] = {
    "no", "metropolis", "barker", "gibbs", "metropolized-gibbs", NULL
};

const char *elmceq_names[elmceqNR+1] = {
    "no", "yes", "wl-delta", "number-all-lambda", "number-steps", "number-samples", "count-ratio", NULL
};

const char *separate_dhdl_file_names[esepdhdlfileNR+1] = {
    "yes", "no", NULL
};

const char *dhdl_derivatives_names[edhdlderivativesNR+1] = {
    "yes", "no", NULL
};

const char *esol_names[esolNR+1] = {
    "No", "SPC", "TIP4p", NULL
};

const char *edispc_names[edispcNR+1] = {
    "No", "EnerPres", "Ener", "AllEnerPres", "AllEner", NULL
};

const char *ecm_names[ecmNR+1] = {
    "Linear", "Angular", "None", NULL
};

const char *eann_names[eannNR+1] = {
    "No", "Single", "Periodic", NULL
};

const char *eis_names[eisNR+1] = {
    "No", "GBSA", NULL
};

const char *egb_names[egbNR+1] = {
    "Still", "HCT", "OBC", NULL
};

const char *esa_names[esaNR+1] = {
    "Ace-approximation", "None", "Still", NULL
};

const char *ewt_names[ewtNR+1] = {
    "9-3", "10-4", "table", "12-6", NULL
};

const char *epull_names[epullNR+1] = {
    "umbrella", "constraint", "constant-force", "flat-bottom", NULL
};

const char *epullg_names[epullgNR+1] = {
    "distance", "direction", "cylinder", "direction-periodic", "direction-relative", NULL
};

const char *erotg_names[erotgNR+1] = {
    "iso", "iso-pf", "pm", "pm-pf", "rm", "rm-pf", "rm2", "rm2-pf", "flex", "flex-t", "flex2", "flex2-t", NULL
};

const char *erotg_fitnames[erotgFitNR+1] = {
    "rmsd", "norm", "potential", NULL
};

const char *eSwapTypes_names[eSwapTypesNR+1] = {
    "no", "X", "Y", "Z", NULL
};

const char *eQMmethod_names[eQMmethodNR+1] = {
    "AM1", "PM3", "RHF",
    "UHF", "DFT", "B3LYP", "MP2", "CASSCF", "B3LYPLAN",
    "DIRECT", NULL
};

const char *eQMbasis_names[eQMbasisNR+1] = {
    "STO3G", "STO-3G", "3-21G",
    "3-21G*", "3-21+G*", "6-21G",
    "6-31G", "6-31G*", "6-31+G*",
    "6-311G", NULL
};

const char *eQMMMscheme_names[eQMMMschemeNR+1] = {
    "normal", "ONIOM", NULL
};

const char *eMultentOpt_names[eMultentOptNR+1] = {
    "multiple_entries", "no", "use_last", NULL
};

const char *eAdresstype_names[eAdressNR+1] = {
    "off", "constant", "xsplit", "sphere", NULL
};

const char *eAdressICtype_names[eAdressICNR+1] = {
    "off", "thermoforce", NULL
};

const char *eAdressSITEtype_names[eAdressSITENR+1] = {
    "com", "cog", "atom", "atomperatom", NULL
};

const char *gmx_nblist_geometry_names[GMX_NBLIST_GEOMETRY_NR+1] = {
    "Particle-Particle", "Water3-Particle", "Water3-Water3", "Water4-Particle", "Water4-Water4", "CG-CG", NULL
};

const char *gmx_nbkernel_elec_names[GMX_NBKERNEL_ELEC_NR+1] =
{
    "None", "Coulomb", "Reaction-Field", "Cubic-Spline-Table", "Generalized-Born", "Ewald", NULL
};

const char *gmx_nbkernel_vdw_names[GMX_NBKERNEL_VDW_NR+1] =
{
    "None", "Lennard-Jones", "Buckingham", "Cubic-Spline-Table", "LJEwald", NULL
};
