/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_GMXPREPROCESS_READIR_H
#define GMX_GMXPREPROCESS_READIR_H

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class MDModules;
struct MDModulesNotifiers;
} // namespace gmx

struct gmx_mtop_t;
struct gmx_output_eenv_t;
struct IndexGroup;
struct pull_params_t;
struct pull_t;
struct t_grpopts;
struct t_inpfile;
struct t_inputrec;
struct t_pull_group;
struct t_pull_coord;
struct t_rot;
struct warninp;
class WarningHandler;

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
    int                warnings  = 0;
    int                nshake    = 0;
    char*              include   = nullptr;
    char*              define    = nullptr;
    bool               bGenVel   = false;
    real               tempi     = 0;
    bool               bMadeSeed = false;
    int                seed      = 0;
    gmx::GromppMtsOpts mtsOpts;
    bool               bOrire           = false;
    bool               bMorse           = false;
    char*              wall_atomtype[2] = { nullptr, nullptr };
    char*              couple_moltype   = nullptr;
    int                couple_lam0      = 0;
    int                couple_lam1      = 0;
    bool               bCoupleIntra     = false;
    bool               deformInitFlow   = false;
};

/*! \brief Initialise object to hold strings parsed from an .mdp file */
void init_inputrec_strings();

/*! \brief Clean up object that holds strings parsed from an .mdp file */
void done_inputrec_strings();

/*! \brief Performs all validation on \p ir that can be done without index groups and topology
 *
 * Any errors, warnings or notes are added to \p wi
 */
void check_ir(const char*                    mdparin,
              const gmx::MDModulesNotifiers& mdModulesNotifiers,
              t_inputrec*                    ir,
              t_gromppopts*                  opts,
              WarningHandler*                wi);

//! Returns the index of string \p s in \p indexGroups or exit with a verbose fatal error when not found
int getGroupIndex(const std::string& s, gmx::ArrayRef<const IndexGroup> indexGroups);

void double_check(t_inputrec* ir, matrix box, bool bHasNormalConstraints, bool bHasAnyConstraints, WarningHandler* wi);
/* Do more checks */

void triple_check(const char* mdparin, t_inputrec* ir, gmx_mtop_t* sys, WarningHandler* wi);
/* Do even more checks */

void get_ir(const char*     mdparin,
            const char*     mdparout,
            gmx::MDModules* mdModules,
            t_inputrec*     ir,
            t_gromppopts*   opts,
            WriteMdpHeader  writeMdpHeader,
            WarningHandler* wi);
/* Read the input file, and retrieve data for inputrec.
 * More data are read, but the are only evaluated when the next
 * function is called. Also prints the input file back to mdparout.
 */

void do_index(const char*                                 mdparin,
              const std::optional<std::filesystem::path>& ndx,
              gmx_mtop_t*                                 mtop,
              bool                                        bVerbose,
              const gmx::MDModulesNotifiers&              mdModulesNotifiers,
              t_inputrec*                                 ir,
              WarningHandler*                             wi);
/* Read the index file and assign grp numbers to atoms.
 */

/* Routines In readpull.c */

std::vector<std::string> read_pullparams(std::vector<t_inpfile>* inp, pull_params_t* pull, WarningHandler* wi);
/* Reads the pull parameters, returns a list of the pull group names */
void process_pull_groups(gmx::ArrayRef<t_pull_group>      pullGroups,
                         gmx::ArrayRef<const std::string> pullGroupNames,
                         gmx::ArrayRef<const IndexGroup>  indexGroups);
/* Process the pull group parameters after reading the index groups */

void checkPullCoords(gmx::ArrayRef<const t_pull_group> pullGroups,
                     gmx::ArrayRef<const t_pull_coord> pullCoords);
/* Process the pull coordinates after reading the pull groups */

pull_t* set_pull_init(t_inputrec*                    ir,
                      const gmx_mtop_t&              mtop,
                      gmx::ArrayRef<const gmx::RVec> x,
                      matrix                         box,
                      real                           lambda,
                      WarningHandler*                wi);
/* Prints the initial pull group distances in x.
 * If requested, adds the current distance to the initial reference location.
 * Returns the pull_t pull work struct. This should be passed to finish_pull()
 * after all modules have registered their external potentials, if present.
 */

std::vector<std::string> read_rotparams(std::vector<t_inpfile>* inp, t_rot* rot, WarningHandler* wi);
/* Reads enforced rotation parameters, returns a list of the rot group names */

void make_rotation_groups(t_rot*                           rot,
                          gmx::ArrayRef<const std::string> rotateGroupNames,
                          gmx::ArrayRef<const IndexGroup>  indexGroups);
/* Process the rotation parameters after reading the index groups */

void set_reference_positions(t_rot* rot, rvec* x, matrix box, const char* fn, bool bSet, WarningHandler* wi);

#endif
