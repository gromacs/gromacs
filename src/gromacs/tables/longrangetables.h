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

#include <vector>

/*! \libinternal \brief
 * Function pointer to fill table.
 *
 * Function pointer used to fill the table.
 * The first parameters is an array of constant used to evaluate
 * the function, the second parameter is the distance,
 * the third parameter is the returned energy and the forth parameter the force.
 */
typedef void (*interaction_potential_function)(const double *, double, double *, double *);

/*! \libinternal\brief
 * Compute the long range energy and force for Ewald Coulomb
 *
 * \param[in]  params Parameters to the long range function
 * \param[in]  r      Distance
 * \param[out] e      Energy
 * \param[out] f      Force
 */
void v_q_ewald_lr(const double *params, double r, double *e, double *f);

/*! \libinternal\brief
 * Compute the long range energy and force for Ewald Lennard Jones (r^{-6} only)
 *
 * \param[in]  params Parameters to the long range function
 * \param[in]  r      Distance
 * \param[out] e      Energy
 * \param[out] f      Force
 *
 * \param[in] params Parameters to the long range function
 * \param[in] r      Distance
 * \return Long range energy
 */
void v_lj_ewald_lr(const double *params, double r, double *e, double *f);

namespace gmx
{
class SplineTable;

/*! \libinternal\brief Class containing long range interaction tables.
 *
 * This class generates input for quadratic spline table interpolation.
 * The force table table_f is always generated. The potential and combined
 * tables table_v and table_fdv0 are not generated when NULL is passed.
 * Input can be generated using either an analytical function, passed
 * as function pointer v_could, v_vdw with parameters params_coul and
 * params_vdw or using table inputs.
 * Table input can be passed as a potential alone, using table_v_coul
 * and table_v_vdw or
 * combined with a force table(s) table_f_coul, table_f_vdw.
 * A stride is used (which can be 1):
 * V(i*dx) = table_in_v(i*stride*dx). When no force table input is passed,
 * the stride has to be a multiple of 2.
 * Note that the table input should be smooth, e.g. not containing noise
 * from an (iterative) Boltzmann inversion procedure.
 *
 * The scale (1/spacing) for third order spline interpolation
 * of the Ewald mesh contribution which needs to be subtracted
 * from the non-bonded interactions.
 * Since there is currently only one spacing for Coulomb and LJ,
 * the finest spacing is used if both Ewald types are present.
 *
 * The class implements completely separate tables for the energy
 * and force for both Coulomb and Van der Waals and a combined
 * table which is more efficient to evaluate in the mdrun kernels.
 *
 * The current setup uses the same spacing to provide slightly
 * faster kernels with both Coulomb and LJ Ewald, especially
 * when interleaving both tables (currently not implemented).
 */
class LongRangeInteractionTables
{
    public:
        /*! \brief Constructor
         *
         * \param[in] type_coul       The electrostatics type
         * \param[in] ewaldcoeff_coul The coefficient for Coulomb
         * \param[in] r_coul          The Coulomb cut-off (switching distance)
         * \param[in] params_coul     Parameters for the function v_qq
         * \param[in] v_coul          A function pointer that can be used
         *                            for computing energies (may be NULL)
         * \param[in] table_v_coul    Table containing energies for Coulomb
         *                            (may be NULL)
         * \param[in] table_f_coul    Table containing forces for Coulomb
         *                            (may be NULL)
         * \param[in] type_vdw        The Van der Waals type
         * \param[in] ewaldcoeff_vdw  The coefficient for Van der Waals
         * \param[in] r_vdw           The Van der Waals cut-off (switching distance)
         * \param[in] params_vdw      Parameters for the function v_vdw
         * \param[in] v_vdw           A function pointer that can be used as an
         * \param[in] table_v_vdw     Table containing energies for Van der
         *                            Waals (may be NULL)
         * \param[in] table_f_vdw     Table containing forces for Van der
         *                            Waals (may be NULL)
         * \param[in] table_length    Length of input tables.
         * \param[in] table_dx        Spacing between data points in
         *                            input tables.
         * \param[in] r_tab           Extra length (nm) for tables beyond cut-offs
         */
        LongRangeInteractionTables(int                               type_coul,
                                   double                            ewaldcoeff_coul,
                                   double                            r_coul,
                                   double                           *params_coul,
                                   interaction_potential_function    v_coul,
                                   const double                     *table_v_coul,
                                   const double                     *table_f_coul,
                                   int                               type_vdw,
                                   double                            ewaldcoeff_vdw,
                                   double                            r_vdw,
                                   double                           *params_vdw,
                                   interaction_potential_function    v_vdw,
                                   const double                     *table_v_vdw,
                                   const double                     *table_f_vdw,
                                   unsigned int                      table_length,
                                   double                            table_dx,
                                   double                            r_tab);

        /*! \brief Copy Constructor since Jenkins wants one
         *
         * \param[in] L The instance to be copied
         */
        LongRangeInteractionTables(const LongRangeInteractionTables &L);
        /*! \brief Return table spacing (nm) */
        double tabSpacing() const { return tabSpacing_; }
        /*! \brief Return number of points in the table */
        unsigned int tabSize() const { return tabSize_; }
        //! Return the Coulomb table class
        const SplineTable *coulomb() const { return tableCoulomb_; }
        //! Return the Van der Waals table class
        const SplineTable *vanderWaals() const { return tableVdw_; }
    private:
        /*! \brief Internal function to compute tabSpacing_ and tabLength_
         *
         * \param[in] type_coul       The electrostatics type
         * \param[in] ewaldcoeff_coul The coefficient for Coulomb
         * \param[in] r_coul          The Coulomb cut-off (switching distance)
         * \param[in] type_vdw        The Van der Waals type
         * \param[in] ewaldcoeff_vdw  The coefficient for Van der Waals
         * \param[in] r_vdw           The Van der Waals cut-off (switching distance)
         * \param[in] r_tab           Extra length (nm) for tables beyond cut-off
         */
        void initTabScaleLength(int    type_coul,
                                double ewaldcoeff_coul,
                                double r_coul,
                                int    type_vdw,
                                double ewaldcoeff_vdw,
                                double r_vdw,
                                double r_tab);
        //! Table scaling (nm)
        double              tabSpacing_;
        //! Table size in number of points
        unsigned int        tabSize_;
        /*! \brief Combined array class for Coulomb
         *
         * Coulomb force+energy table, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
         * this is used with single precision x86 SIMD for aligned loads
         */
        SplineTable         *tableCoulomb_;

        /*! \brief Combined array class for Van der Waals
         *
         * Vdw force+energy table for LJ-PME, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0, this is used with
         * single precision x86 SIMD for aligned loads
         */
        SplineTable         *tableVdw_;
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
