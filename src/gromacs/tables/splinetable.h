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
#ifndef GMX_TABLES_SPLINETABLE_H
#define GMX_TABLES_SPLINETABLE_H

#include <cstdio>
#include <cstdlib>

#include <vector>

#include "gromacs/utility/classhelpers.h"

/*! \libinternal \brief
 * Function pointer to fill table.
 *
 * Function pointer used to fill the table.
 * The first parameters is an array of constant used to evaluate
 * the function, the second parameter is the distance,
 * the third parameter is the returned energy and the forth parameter the force.
 */
typedef void (*interaction_potential_function)(const double *, double, double *, double *);

namespace gmx
{
/*! \libinternal\brief Class containing one spline interaction table.
 *
 * This class generates input for quadratic spline table interpolation
 * from input tables.
 * If both energy and force tables are input the number of output
 * points is equal to the input.
 * When no force table input is passed the number of points in the
 * table is halved because an extra interpolation is needed.
 * In general it is more accurate to provide both energy and force tables.
 *
 * Functions using the tables should use methods tabSpacing() and tabSize()
 * to get the actual spacing and number of points.
 *
 * Note that the table input should be smooth, e.g. not containing noise
 * from an (iterative) Boltzmann inversion procedure.
 *
 * The class implements completely separate tables for the energy
 * and force for both Coulomb and Van der Waals and a combined
 * table which is more efficient to evaluate in the mdrun kernels.
 */
class QuadraticSplineTable
{
    public:
        /*! \brief Constructor using table input
         *
         * \param[in] ntable     Number of data points in the table
         * \param[in] dx         Space between data points
         * \param[in] table_v    Table containing energies
         * \param[in] table_f    Table containing forces (may be NULL)
         */
        QuadraticSplineTable(unsigned int  ntable,
                             double        dx,
                             const double *table_v,
                             const double *table_f);

        /*! \brief Constructor using analytical function
         *
         * This utility generates a quadratic spline table and fills it.
         * Input is generated using an analytical function, passed
         * as function pointer v_ana with parameters params.
         *
         * The scale (1/spacing) for third order spline interpolation
         * of the Ewald mesh contribution which needs to be subtracted
         * from the non-bonded interactions.
         *
         * \param[in] ntable  Length of input tables.
         * \param[in] dx      Spacing between data points in
         *                    input tables.
         * \param[in] v_ana   A function pointer that can be used
         *                    for computing energies (may be NULL)
         * \param[in] params  Parameters for the function v_analytical
         * \return Pointer to spline table. The calling routine should delete it.
         */
        QuadraticSplineTable(unsigned int                    ntable,
                             double                          dx,
                             interaction_potential_function  v_ana,
                             const double                   *params);

        /*! \brief Return the force at a distance r */
        double force(double r) const;

        /*! \brief Return the energy at a distance r */
        double energy(double r) const;

        /*! \brief Return table spacing (nm) */
        double spacing() const { return tabSpacing_; }

        /*! \brief Return number of points in the table */
        unsigned int size() const { return tabSize_; }

        //! Return the force table
        const std::vector<double> &F() const { return tableF_; }
        //! Return the energy table
        const std::vector<double> &V() const { return tableV_; }
        //! Return the force + energy table
        const std::vector<double> &FDV0() const { return tableFDV0_; }
        //! Dump the contents of the energy and force arrays to a FILE
        void dump(FILE *fp) const;
    private:
        /*! \brief Do the actual filling of the tables
         *
         * \param[in] stride      Distance between data points
         * \param[in] table_in_v  Table containing energies
         * \param[in] table_in_f  Table containing forces (may be NULL)
         */
        void fillTable(unsigned int  stride,
                       const double *table_in_v,
                       const double *table_in_f);
        //! Table scaling (nm)
        double              tabSpacing_;
        //! Length of table (number of points)
        unsigned int        tabSize_;
        //! Force table
        std::vector<double> tableF_;
        //! Energy table
        std::vector<double> tableV_;
        /*! \brief Combined array for force and energy
         *
         * Force+energy table, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
         * this is used with single precision x86 SIMD for aligned loads
         */
        std::vector<double> tableFDV0_;

        GMX_DISALLOW_COPY_AND_ASSIGN(QuadraticSplineTable);
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
