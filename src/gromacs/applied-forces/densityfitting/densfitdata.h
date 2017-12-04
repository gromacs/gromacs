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

#include "densfitmapdata.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"

struct t_blocka;
struct t_inpfile;
struct t_commrec;
struct t_fileio;
struct warninp;
struct gmx_output_env_t;

/* \TODO : very dirty hack until proper conversion to mdpOptions. */
#ifndef t_inputrec_strings
#define t_inputrec_strings
typedef struct t_inputrec_strings {
    char tcgrps[STRLEN], tau_t[STRLEN], ref_t[STRLEN], acc[STRLEN],
         accgrps[STRLEN], freeze[STRLEN], frdim[STRLEN], energy[STRLEN],
         user1[STRLEN], user2[STRLEN], vcm[STRLEN], x_compressed_groups[STRLEN],
         couple_moltype[STRLEN], orirefitgrp[STRLEN], egptable[STRLEN],
         egpexcl[STRLEN], wall_atomtype[STRLEN], wall_density[STRLEN],
         deform[STRLEN], QMMM[STRLEN], imd_grp[STRLEN];
    char   fep_lambda[efptNR][STRLEN];
    char   lambda_weights[STRLEN];
    char **pull_grp;
    char **rot_grp;
    char  *dens_grp;
    char   anneal[STRLEN], anneal_npoints[STRLEN], anneal_time[STRLEN],
           anneal_temp[STRLEN];
    char   QMmethod[STRLEN], QMbasis[STRLEN], QMcharge[STRLEN], QMmult[STRLEN],
           bSH[STRLEN], CASorbitals[STRLEN], CASelectrons[STRLEN], SAon[STRLEN],
           SAoff[STRLEN], SAsteps[STRLEN];
    char densfit_time[STRLEN], densfit_sigma[STRLEN], densfit_k[STRLEN],
         densfit_temp[STRLEN];

} gmx_inputrec_strings;
#endif

namespace gmx
{

class TimeDependentDensfitParameters
{
    public:
        TimeDependentDensfitParameters() = default;
        TimeDependentDensfitParameters(const std::vector<real> &times,
                                       const std::vector<real> &forceConstants,
                                       const std::vector<real> &sigmas);

        TimeDependentDensfitParameters(const char timePointString[STRLEN],
                                       const char forceConstantString[STRLEN],
                                       const char sigmaString[STRLEN], warninp *wi);
        void print(FILE *fp, int indent);
        real currentK(real time) const;
        real currentSigma(real time) const;
        void do_fio(t_fileio *fio, bool bRead);
        void print(FILE *fp, int indent) const;
        const std::vector<real> &timePoints() const;
        const std::vector<real> &forceConstants() const;
        void broadcast(const t_commrec *cr);

    private:
        std::vector<real> timePoints_;     /* npoints of time values */
        std::vector<real> forceConstants_; /* force constants  k (kJ/mol) */
        std::vector<real> sigmas_;         /* Same for sigma (nm) */

        /*! Linearly interpolate function at x0 between given values x,y.
         *
         * Assumes sorted x_values and y_values corresponding to x_values.
         * If x0 > max(x_values) or x0 < min(x_values) return last or first y_value
         * respectively.
         */
        real interpolateLinearly(real x0, const std::vector<real> &x,
                                 const std::vector<real> &y) const;
};

/*! \brief Parameters for density fitting that are known before the simulation.
 *
 * Anything that is known about density fitting and is independent of the course
 * of the simulation belongs in this class.
 * \TODO: grompp should be responsible for all pre-calculation and settings in
 * this classs
 */
class DensfitData
{
    public:
        DensfitData() = default;
        DensfitData(const std::vector<real> &sigma, const std::vector<real> &k,
                    real sigma_dist, int nAtoms, int *index, bool bVerbose);
        void setReferenceMap(const t_mapdata &reference);
        const t_mapdata &referenceMap();
        void makeGroups(char *groupname, t_blocka *grps, char **gnames);
        void broadcast(const t_commrec *cr);
        real dist() const;
        int nAtoms() const;
        int *indices() const;
        void printToOut(FILE *fp, const gmx_output_env_t *oenv) const;
        void setNstFitFromEnv(const t_commrec *cr, int nstlist, FILE *fp);
        bool keepAndNumberMaps() const;
        int nStepsMapOutput() const;
        int nStepsFit() const;
        void print(FILE *fp, int indent) const;
        const TimeDependentDensfitParameters &timeDependent() const;
        void do_fio(t_fileio *fio, bool bRead);
        void read_densparams(int *ninp_p, t_inpfile **inp_p, warninp *wi,
                             gmx_inputrec_strings *is);
        bool areTimeDependent();

    private:
        t_mapdata map_ref_;           /* The reference=experimental map to fit to        */
        TimeDependentDensfitParameters timeDependent_;
        real      dist_;              /* Cutoff distance (in sigmas) for density spreading spread each
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

#endif /* end of include guard: \
          GMX_APPLIEDFORCES_DENSITYFITTTING_DENSFITDATA_H */
