/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_DISPERSIONCORRECTION_H
#define GMX_MDLIB_DISPERSIONCORRECTION_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

struct gmx_mtop_t;
struct interaction_const_t;
struct t_forcerec;
struct t_inputrec;

namespace gmx
{
class MDLogger;
} // namespace gmx

class DispersionCorrection
{
    public:
        DispersionCorrection(const gmx_mtop_t          &mtop,
                             const t_inputrec          &inputrec,
                             bool                       useBuckingham,
                             int                        numAtomTypes,
                             gmx::ArrayRef<const real>  nonbondedForceParameters);

        /* Print dispersion correction information to the log file */
        void print(const gmx::MDLogger &mdlog) const;

        /* Computes and sets energy and virial correction parameters based on the cut-off and LJ-ewald coefficient */
        void setParameters(const interaction_const_t &ic,
                           const char                *tableFileName);

        /* Computes and returns the dispersion correction for the pressure and energy */
        void calculate(const matrix box, real lambda, tensor pres, tensor virial,
                       real *prescorr, real *enercorr, real *dvdlcorr) const;
       
    private:
        int                  eDispCorr_;
        int                  vdwType_;
        int                  eFep_;
        int                  numAtomsForDispersionCorrection_;
        int                  numAtomsForTpi_;
        struct t_forcetable *dispersionCorrectionTable_;

        /* The shift of the shift or user potentials */
        real enershiftsix_;
        real enershifttwelve_;
        /* Integrated differces for energy and virial with cut-off functions */
        real enerdiffsix_;
        real enerdifftwelve_;
        real virdiffsix_;
        real virdifftwelve_;
        /* Constant for long range dispersion correction (average dispersion)
         * for topology A/B ([0]/[1]) */
        std::array<real, 2> avcsix_;
        /* Constant for long range repulsion term. Relative difference of about
         * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
         */
        std::array<real, 2> avctwelve_;
};

#endif
