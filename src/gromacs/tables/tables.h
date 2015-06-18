/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_TABLES_TABLES_H
#define GMX_TABLES_TABLES_H

#include <cstdlib>

#include "gromacs/utility/real.h"

/*! \libinternal \brief
 * Function pointer to fill table.
 *
 * Function pointer used to tell table_spline3_fill_ewald_lr whether it
 * should calculate the grid contribution for electrostatics or LJ.
 * The first parameters is an array of constant used to evaluate
 * the function, the second parameter is the distance.
 */
typedef double (*interaction_potential_function)(const double *, double);

/*! \libinternal\brief
 * Return the long range energy for Ewald Coulomb
 *
 * \param[in] params Parameters to the long range function
 * \param[in] r      Distance
 * \return Long range energy
 */
double v_q_ewald_lr(const double *params, double r);

/*! \libinternal\brief
 * Return the long range energy for Ewald Lennard Jones (r^{-6} only)
 *
 * \param[in] params Parameters to the long range function
 * \param[in] r      Distance
 * \return Long range energy
 */
double v_lj_ewald_lr(const double *params, double r);

namespace gmx
{

//! Enumerated type to differentiate interactions for Ewald summation.
enum EwaldType {
    etCoulomb, etLJ
};

/*! \libinternal\brief Class containing interaction tables.
 *
 * Force/energy interpolation tables, linear in force, quadratic in V
 */
class InteractionTables
{
    public:
        /*! \libinternal \brief Constructor initiating arrays to NULL */
        InteractionTables() : tabScale_(0.0), tabCoulF_(NULL), tabCoulV_(NULL), tabCoulFDV0_(NULL),
                              tabVdwF_(NULL), tabVdwV_(NULL), tabVdwFDV0_(NULL) {};

        /*! \libinternal \brief Destructor feeing memory */
        ~InteractionTables();
        /*! \libinternal\brief Compute and return the table scale
         *
         * \param[in] eeltype       The electrostatics type
         * \param[in] ewaldcoeff_q  The coefficient for Coulomb
         * \param[in] rcoulomb      The Coulomb cut-off (switching distance)
         * \param[in] vdwtype       The Van der Waals type
         * \param[in] ewaldcoeff_lj The coefficient for Van der Waals
         * \param[in] rvdw          The Van der Waals cut-off (switching distance)
         * \return the least scaling factor of the two
         */
        real tabScale(int  eeltype,
                      real ewaldcoeff_q,
                      real rcoulomb,
                      int  vdwtype,
                      real ewaldcoeff_lj,
                      real rvdw);

        /*! \libinternal\brief
         * Fill a table using a spline interpolation.
         *
         * This function interpolates the Ewald mesh potential contribution
         * using a quadratic spline.
         * The force can then be interpolated linearly.
         * \param[in]  ntab       Length of the table
         * \param[in]  dx         Spacing between the points
         * \param[in]  params     Parameters for the function v_lr
         * \param[in]  et         Determines whether Coulomb or LJ tables are filled
         * \param[in]  v_lr       A function pointer that can be used as an
         *                        alternative to the default Ewald long range.
         */
        void fillEwaldSpline3(int                             ntab,
                              double                          dx,
                              double                         *params,
                              EwaldType                       et,
                              interaction_potential_function  v_lr);
        //! Return the Coulomb force table
        real *coulombF() { return tabCoulF_; }
        //! Return the Coulomb potential table
        real *coulombV() { return tabCoulV_; }
        //! Return the Coulomb force + energy table
        real *coulombFDV0() { return tabCoulFDV0_; }
        //! Return the Van der Waals force table
        real *vdwF() { return tabVdwF_; }
        //! Return the Van der Waals potential table
        real *vdwV() { return tabVdwV_; }
        //! Return the Van der Waals force + energy table
        real *vdwFDV0() { return tabVdwFDV0_; }
    private:
        //! Table scaling
        real  tabScale_;
        //! Coulomb force table
        real *tabCoulF_;
        //! Coulomb energy table
        real *tabCoulV_;
        /*! \brief Combined array for Coulomb
         *
         * Coulomb force+energy table, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
         * this is used with single precision x86 SIMD for aligned loads
         */
        real *tabCoulFDV0_;

        //! Vdw force table for LJ-PME
        real *tabVdwF_;
        //! Vdw energy table for LJ-PME
        real *tabVdwV_;
        /*! \brief Combined array for Van der Waals
         *
         * Vdw force+energy table for LJ-PME, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0, this is used with
         * single precision x86 SIMD for aligned loads
         */
        real *tabVdwFDV0_;
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
