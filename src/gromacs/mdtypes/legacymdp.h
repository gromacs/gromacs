/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \file
 * \brief Declares data structure for inputrec fields that currently
 * simply map to legacy .mdp fields.
 *
 * \todo Eliminate the need for this handling and data structure.
 *
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_LEGACYMDP_H
#define GMX_MDTYPES_LEGACYMDP_H

#include <cstdio>

#include <memory>
#include <string>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/cstringutil.h"

struct gmx_output_env_t;
struct t_filenm;
struct t_inputrec;

// TODO These will move to InputrecStrings once all mdp options
// are handled with new machinery.
struct t_inputrec_strings
{
    char tcgrps[STRLEN], tau_t[STRLEN], ref_t[STRLEN],
         acc[STRLEN], accgrps[STRLEN], freeze[STRLEN], frdim[STRLEN],
         energy[STRLEN], user1[STRLEN], user2[STRLEN], vcm[STRLEN],
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
           SAoff[STRLEN], SAsteps[STRLEN], bTS[STRLEN], bOPT[STRLEN];

};

struct InputrecStrings
{
    std::string x_compressed_groups;
};

// This struct contains values for mdp options that are handled by
// grompp before it writes out the .tpr file.
// TODO These will move to GromppOptions once all mdp options are handled
// with new machinery.
struct t_gromppopts
{
    int      warnings;
    int      nshake;
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

// See t_gromppopts
// TODO Should it become an IInputRecExtension? Separate from MDModules?
struct GromppOptions
{
    std::string include;
    std::string define;
};

namespace gmx
{

class IKeyValueTreeTransformRules;
class IOptionsContainerWithSections;
class KeyValueTreeObject;
class LegacyMdp;

/*! \brief Creates a transitional module for handling all legacy mdp
 * fields that have not been updated to new-style modules.
 */
std::unique_ptr<LegacyMdp> createLegacyMdpModule(t_inputrec *ir);

/*! \internal
 * \brief Describe legacy mdp fields.
 *
 * TODO
 *
 * The public fields exist to smooth the transition period from
 * LegacyMdp to module handling, and will disappear along with
 * LegacyMdp eventually.
 */
class LegacyMdp : public IInputRecExtension
{
    public:
        LegacyMdp(t_inputrec *ir);
        ~LegacyMdp();

        virtual void initMdpTransform(IKeyValueTreeTransformRules *rules);
        virtual void initMdpOptions(IOptionsContainerWithSections *options);
        virtual void initMdpBackTransform(IKeyValueTreeTransformRules *rules) const;

        virtual void initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                                bool bAppendFiles, const gmx_output_env_t *oenv);
        virtual void finishOutput();
        virtual void initForcerec(t_forcerec *fr);
    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
    public:
        //! Handle to struct for inputrec strings.
        t_inputrec_strings *is_;
        //! Handle to struct for inputrec strings.
        InputrecStrings    *inputrecStrings_;
        //! Handle to struct for grompp options.
        t_gromppopts       *opts_;
        //! Handle to struct for grompp options.
        GromppOptions      *gromppOptions_;
};

} // namespace

#endif
