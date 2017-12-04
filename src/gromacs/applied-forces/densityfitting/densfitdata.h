/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#ifndef GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H
#define GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H

#include "gromacs/utility/real.h"

#include "gromacs/utility/cstringutil.h"
#include "gromacs/mdtypes/md_enums.h"

struct t_blocka;
struct t_inpfile;
struct t_commrec;
struct t_fileio;
struct warninp;
struct gmx_output_env_t;

/* \TODO : very dirty hack until proper conversion to mdpOptions. */
#ifndef t_inputrec_strings
#define t_inputrec_strings
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
    char  *dens_grp;
    char   anneal[STRLEN], anneal_npoints[STRLEN],
           anneal_time[STRLEN], anneal_temp[STRLEN];
    char   QMmethod[STRLEN], QMbasis[STRLEN], QMcharge[STRLEN], QMmult[STRLEN],
           bSH[STRLEN], CASorbitals[STRLEN], CASelectrons[STRLEN], SAon[STRLEN],
           SAoff[STRLEN], SAsteps[STRLEN];
    char densfit_time[STRLEN], densfit_sigma[STRLEN], densfit_k[STRLEN], densfit_temp[STRLEN];

} gmx_inputrec_strings;
#endif

namespace gmx
{

/*! \brief Parameters for density fitting that are known before the simulation.
 *
 * Anything that is known about density fitting and is independent of the course
 * of the simulation belongs in this class.
 * \TODO: grompp should be responsible for all pre-calculation and settings in this classs
 * \TODO: move the reference map into this class
 * \TODO: all members const?
 */
class DensfitData
{
    public:
        DensfitData();
        DensfitData(int npoint, real sigma, real k, real sigma_dist, int nAtoms,
                    int *index, bool bVerbose);
        void makeGroups(char *groupname, t_blocka *grps, char **gnames);
        void broadcast(const t_commrec * cr);
        real currentK(real time) const;
        real currentT(real time) const;
        real currentSigma(real time) const;
        real dist() const;
        int nAtoms() const;
        int *indices() const;
        void printToOut(FILE * fp, const gmx_output_env_t *oenv) const;
        void setNstFitFromEnv(const t_commrec * cr, int nstlist, FILE * fp);
        bool keepAndNumberMaps() const;
        int timePoints() const;
        int nStepsMapOutput() const;
        int nStepsFit() const;
        void print(FILE * fp, int indent) const;
        void do_fio(t_fileio *fio, bool bRead);
        void read_densparams(int *ninp_p, t_inpfile **inp_p, warninp * wi, gmx_inputrec_strings * is);
    private:
        real interpolateLinearly(real time, int npoints, real *time_values,
                                 real *y_values) const;

        int   npoints_;               /* Number of time points for the refinement with
                                         different values for sigma, k and T             */
        real *time_values_;           /* Array (size npoints) of time values (ps)        */
        real *sigma_values_;          /* Same for sigma (nm)                             */
        real *k_values_;              /* Same for k (kJ/mol)                             */
        real *temp_values_;           /* Same for T (K)                                  */
        real  dist_;                  /* Cutoff distance (in sigmas) for density spreading spread each
                                         atom only in the range -dist ... +dist, ignore contributions
                                         further away                                            */
        int      nstfit_;             /* Only recalculate V_fit every nstfit time steps  */
        int      nstout_;             /* Write diagnostic output every nstout time steps */
        int      nstmapout_;          /* Keep simulated maps in these intervals          */
        int      nat_;                /* Number of atoms to be spread on the map         */
        int     *ind_;                /* The global atoms numbers                        */
        gmx_bool bKeepAndNumberMaps_; /* ... or just keep the last one ... */
};

}

#endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H */
