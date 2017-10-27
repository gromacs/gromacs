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

#ifndef GMX_GMXPREPROCESS_READIR_H
#define GMX_GMXPREPROCESS_READIR_H

#include "gromacs/fileio/readinp.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"

namespace gmx
{
class MDModules;
}

struct gmx_groups_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_params_t;
struct pull_t;
struct t_grpopts;
struct t_inpfile;
struct t_inputrec;
struct t_rot;
struct warninp;
typedef warninp *warninp_t;

enum {
    eshNONE, eshHBONDS, eshALLBONDS, eshHANGLES, eshALLANGLES, eshNR
};

enum {
    ecouplamVDWQ, ecouplamVDW, ecouplamQ, ecouplamNONE, ecouplamNR
};

struct t_gromppopts
{
    int      warnings;
    int      nshake;
    char    *include;
    char    *define;
    gmx_bool bGenVel;
    gmx_bool bGenPairs;
    real     tempi;
    int      seed;
    gmx_bool bOrire;
    gmx_bool bMorse;
    char    *wall_atomtype[2];
    char    *couple_moltype;
    int      couple_lam0;
    int      couple_lam1;
    gmx_bool bCoupleIntra;
};

/*! \brief Initialise object to hold strings parsed from an .mdp file */
void init_inputrec_strings();

/*! \brief Clean up object that holds strings parsed from an .mdp file */
void done_inputrec_strings();

void check_ir(const char *mdparin, t_inputrec *ir, t_gromppopts *opts,
              warninp_t wi);
/* Validate inputrec data.
 * Fatal errors will be added to nerror.
 */
int search_string(const char *s, int ng, char *gn[]);
/* Returns the index of string s in the index groups */

void double_check(t_inputrec *ir, matrix box,
                  gmx_bool bHasNormalConstraints,
                  gmx_bool bHasAnyConstraints,
                  warninp_t wi);
/* Do more checks */

void triple_check(const char *mdparin, t_inputrec *ir, gmx_mtop_t *sys,
                  warninp_t wi);
/* Do even more checks */

void check_chargegroup_radii(const gmx_mtop_t *mtop, const t_inputrec *ir,
                             rvec *x,
                             warninp_t wi);
/* Even more checks, charge group radii vs. cut-off's only. */

void get_ir(const char *mdparin, const char *mdparout,
            gmx::MDModules *mdModules, t_inputrec *ir, t_gromppopts *opts,
            WriteMdpHeader writeMdpHeader, warninp_t wi);
/* Read the input file, and retrieve data for inputrec.
 * More data are read, but the are only evaluated when the next
 * function is called. Also prints the input file back to mdparout.
 */

void do_index(const char* mdparin,
              const char *ndx,
              gmx_mtop_t *mtop,
              gmx_bool    bVerbose,
              t_inputrec *ir,
              warninp_t   wi);
/* Read the index file and assign grp numbers to atoms.
 */

/* Routines In readpull.c */

char **read_pullparams(int *ninp_p, t_inpfile **inp,
                       pull_params_t *pull,
                       warninp_t wi);
/* Reads the pull parameters, returns a list of the pull group names */

void make_pull_groups(pull_params_t *pull,
                      char **pgnames,
                      const t_blocka *grps, char **gnames);
/* Process the pull group parameters after reading the index groups */

void make_pull_coords(pull_params_t *pull);
/* Process the pull coordinates after reading the pull groups */

pull_t *set_pull_init(t_inputrec *ir, const gmx_mtop_t *mtop,
                      rvec *x, matrix box, real lambda,
                      const gmx_output_env_t *oenv);
/* Prints the initial pull group distances in x.
 * If requested, adds the current distance to the initial reference location.
 * Returns the pull_t pull work struct. This should be passed to finish_pull()
 * after all modules have registered their external potentials, if present.
 */

char **read_rotparams(int *ninp_p, t_inpfile **inp, t_rot *rot, warninp_t wi);
/* Reads enforced rotation parameters, returns a list of the rot group names */

void make_rotation_groups(t_rot *rot, char **rotgnames,
                          t_blocka *grps, char **gnames);
/* Process the rotation parameters after reading the index groups */

void set_reference_positions(t_rot *rot, rvec *x, matrix box,
                             const char *fn, gmx_bool bSet, warninp_t wi);

#endif
