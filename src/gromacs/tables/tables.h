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
/*! \libinternal\brief Class containing interaction tables.
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
class InteractionTables
{
    public:
        /*! \brief Constructor
         *
         * \param[in] type_coul       The electrostatics type
         * \param[in] ewaldcoeff_coul The coefficient for Coulomb
         * \param[in] r_coul          The Coulomb cut-off (switching distance)
         * \param[in] params_could    Parameters for the function v_qq
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
         * \param[in] table_in_length Length of input tables.
         * \param[in] table_in_stride Distance between usefule data points in
         *                            input tables.
         */
        InteractionTables(int                               type_coul,
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
                          unsigned int                      table_in_length,
                          unsigned int                      table_in_stride);

        /*! \brief Return table spacing (nm) */
        double tabScale() const { return tabScale_; }
        /*! \brief Return table length (nm) */
        double tabLength() const { return tabLength_; }

        //! Return the Coulomb force table
        const std::vector<double> &coulombF() const { return tabCoulF_; }
        //! Return the Coulomb potential table
        const std::vector<double> &coulombV() const { return tabCoulV_; }
        //! Return the Coulomb force + energy table
        const std::vector<double> &coulombFDV0() const { return tabCoulFDV0_; }
        //! Return the Van der Waals force table
        const std::vector<double> &vdwF() const { return tabVdwF_; }
        //! Return the Van der Waals potential table
        const std::vector<double> &vdwV() const { return tabVdwV_; }
        //! Return the Van der Waals force + energy table
        const std::vector<double> &vdwFDV0() const { return tabVdwFDV0_; }
    private:
        /*! \brief Internal function to compute tabScale_ and tabLength_
         *
         * \param[in] type_coul       The electrostatics type
         * \param[in] ewaldcoeff_coul The coefficient for Coulomb
         * \param[in] r_coul          The Coulomb cut-off (switching distance)
         * \param[in] type_vdw        The Van der Waals type
         * \param[in] ewaldcoeff_vdw  The coefficient for Van der Waals
         * \param[in] r_vdw           The Van der Waals cut-off (switching distance)
         */
        void initTabScaleLength(int    type_coul,
                                double ewaldcoeff_coul,
                                double r_coul,
                                int    type_vdw,
                                double ewaldcoeff_vdw,
                                double r_vdw);
        //! Table scaling (nm)
        double              tabScale_;
        //! Table length (nm)
        double              tabLength_;
        //! Coulomb force table
        std::vector<double> tabCoulF_;
        //! Coulomb energy table
        std::vector<double> tabCoulV_;
        /*! \brief Combined array for Coulomb
         *
         * Coulomb force+energy table, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
         * this is used with single precision x86 SIMD for aligned loads
         */
        std::vector<double> tabCoulFDV0_;

        //! Vdw force table for LJ-PME
        std::vector<double> tabVdwF_;
        //! Vdw energy table for LJ-PME
        std::vector<double> tabVdwV_;
        /*! \brief Combined array for Van der Waals
         *
         * Vdw force+energy table for LJ-PME, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0, this is used with
         * single precision x86 SIMD for aligned loads
         */
        std::vector<double> tabVdwFDV0_;
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
