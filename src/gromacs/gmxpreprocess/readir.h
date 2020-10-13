/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

#include <string>

#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class MDModules;
struct MdModulesNotifier;
} // namespace gmx

struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_params_t;
struct pull_t;
struct t_blocka;
struct t_grpopts;
struct t_inpfile;
struct t_inputrec;
struct t_pull_group;
struct t_pull_coord;
struct t_rot;
struct warninp;
typedef warninp* warninp_t;

enum
{
    eshNONE,
    eshHBONDS,
    eshALLBONDS,
    eshHANGLES,
    eshALLANGLES,
    eshNR
};

enum
{
    ecouplamVDWQ,
    ecouplamVDW,
    ecouplamQ,
    ecouplamNONE,
    ecouplamNR
};

struct t_gromppopts
{
    int         warnings     = 0;
    int         nshake       = 0;
    char*       include      = nullptr;
    char*       define       = nullptr;
    bool        bGenVel      = false;
    bool        bGenPairs    = false;
    real        tempi        = 0;
    int         seed         = 0;
    int         numMtsLevels = 0;
    std::string mtsLevel2Forces;
    bool        bOrire           = false;
    bool        bMorse           = false;
    char*       wall_atomtype[2] = { nullptr, nullptr };
    char*       couple_moltype   = nullptr;
    int         couple_lam0      = 0;
    int         couple_lam1      = 0;
    bool        bCoupleIntra     = false;
};

/*! \brief Initialise object to hold strings parsed from an .mdp file */
void init_inputrec_strings();

/*! \brief Clean up object that holds strings parsed from an .mdp file */
void done_inputrec_strings();

void check_ir(const char*                   mdparin,
              const gmx::MdModulesNotifier& mdModulesNotifier,
              t_inputrec*                   ir,
              t_gromppopts*                 opts,
              warninp_t                     wi);
/* Validate inputrec data.
 * Fatal errors will be added to nerror.
 */
int search_string(const char* s, int ng, char* gn[]);
/* Returns the index of string s in the index groups */

void double_check(t_inputrec* ir, matrix box, bool bHasNormalConstraints, bool bHasAnyConstraints, warninp_t wi);
/* Do more checks */

void triple_check(const char* mdparin, t_inputrec* ir, gmx_mtop_t* sys, warninp_t wi);
/* Do even more checks */

void get_ir(const char*     mdparin,
            const char*     mdparout,
            gmx::MDModules* mdModules,
            t_inputrec*     ir,
            t_gromppopts*   opts,
            WriteMdpHeader  writeMdpHeader,
            warninp_t       wi);
/* Read the input file, and retrieve data for inputrec.
 * More data are read, but the are only evaluated when the next
 * function is called. Also prints the input file back to mdparout.
 */

void do_index(const char*                   mdparin,
              const char*                   ndx,
              gmx_mtop_t*                   mtop,
              bool                          bVerbose,
              const gmx::MdModulesNotifier& notifier,
              t_inputrec*                   ir,
              warninp_t                     wi);
/* Read the index file and assign grp numbers to atoms.
 */

/* Routines In readpull.c */

std::vector<std::string> read_pullparams(std::vector<t_inpfile>* inp, pull_params_t* pull, warninp_t wi);
/* Reads the pull parameters, returns a list of the pull group names */
void process_pull_groups(gmx::ArrayRef<t_pull_group>      pullGroups,
                         gmx::ArrayRef<const std::string> pullGroupNames,
                         const t_blocka*                  grps,
                         char**                           gnames);
/* Process the pull group parameters after reading the index groups */

void checkPullCoords(gmx::ArrayRef<const t_pull_group> pullGroups,
                     gmx::ArrayRef<const t_pull_coord> pullCoords);
/* Process the pull coordinates after reading the pull groups */

pull_t* set_pull_init(t_inputrec* ir, const gmx_mtop_t* mtop, rvec* x, matrix box, real lambda, warninp_t wi);
/* Prints the initial pull group distances in x.
 * If requested, adds the current distance to the initial reference location.
 * Returns the pull_t pull work struct. This should be passed to finish_pull()
 * after all modules have registered their external potentials, if present.
 */

std::vector<std::string> read_rotparams(std::vector<t_inpfile>* inp, t_rot* rot, warninp_t wi);
/* Reads enforced rotation parameters, returns a list of the rot group names */

void make_rotation_groups(t_rot*                           rot,
                          gmx::ArrayRef<const std::string> rotateGroupNames,
                          t_blocka*                        grps,
                          char**                           gnames);
/* Process the rotation parameters after reading the index groups */

void set_reference_positions(t_rot* rot, rvec* x, matrix box, const char* fn, bool bSet, warninp_t wi);

#endif
