/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_names_c = "$Id$";
#include "typedefs.h"
#include "names.h"

/* note: these arrays should correspond to enums in include/types/enums.h */

char *eblock_names[ebNR+1]=
{
  "CGS","MOLS","SBLOCKS",NULL
};

char *epbc_names[epbcNR+1]=
{
  "xyz", "no", NULL
};

char *ens_names[ensNR+1]=
{
  "Grid","Simple", NULL
};

char *ei_names[eiNR+1]=
{
  "md", "steep", "cg", "bd", "sd", "nm", NULL 
};

char *bool_names[BOOL_NR+1]=
{
  "FALSE","TRUE", NULL
};

char *yesno_names[BOOL_NR+1]=
{
  "no","yes", NULL
};

char *ptype_str[eptNR+1] = {
  "Atom", "Nucleus", "Shell", "Bond", "Dummy", NULL
};

char *eel_names[eelNR+1] = {
  "Cut-off", "Reaction-Field", "Generalized-Reaction-Field",
  "PME", "Ewald", "PPPM", "Poisson", "Switch", "Shift", "User", NULL
};

char *eewg_names[eewgNR+1] = {
  "3d", "3dc", NULL
};

char *evdw_names[evdwNR+1] = {
  "Cut-off", "Switch", "Shift", "User", NULL
};

char *eshake_names[estNR+1] = {
  "Lincs", "Shake", NULL
};

char *egrp_nm[egNR+1] = { 
  "Coul-SR","LJ","Buck", "Coul-LR", "LJ-LR", "Coul-14", "LJ-14", NULL
};

char *etcoupl_names[etcNR+1] = {
  "No", "Berendsen", "Nose-Hoover", "yes", NULL
}; /* yes is alias for berendsen */

char *epcoupl_names[epcNR+1] = {
  "No", "Berendsen", "Parrinello-Rahman", "Isotropic", NULL
}; /* isotropic is alias for berendsen */

char *epcoupltype_names[epctNR+1] = {
  "Isotropic", "Semiisotropic", "Anisotropic", "Surface-Tension", NULL
};

char *edisre_names[edrNR+1] = {
  "No", "Simple", "Ensemble", NULL
};

char *edisreweighting_names[edrwNR+1] = {
  "Conservative", "Equal", NULL
};

char *enbf_names[eNBF_NR+1] = {
  "", "LJ", "Buckingham", NULL
};

char *ecomb_names[eCOMB_NR+1] = {
  "", "Arithmetic", "Geometric", "ArithSigEps", NULL
};

char *gtypes[egcNR+1] = {
  "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
  "User1", "User2", "VCM", "XTC", "Or. Res. Fit", NULL
};

char *efep_names[efepNR+1] = {
  "no", "yes", NULL
};

char *esolv_names[esolNR+1] = {
  "General", "MNO Solvent", "Water", "Water-Water", NULL
};

char *edispc_names[edispcNR+1] = {
  "No", "EnerPres", "Ener", NULL
};

char *ecm_names[ecmNR+1] = { 
  "Linear", "Angular", "None", NULL 
};
 
