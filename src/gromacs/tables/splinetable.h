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

#include <cstdlib>

#include <vector>

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
class SplineTable
{
    public:
        /*! \brief Constructor
         *
         * \param[in] ntable     Number of data points in the table
         * \param[in] dx         Space between data points
         * \param[in] table_v    Table containing energies
         * \param[in] table_f    Table containing forces (may be NULL)
         */
        SplineTable(unsigned int  ntable,
                    double        dx,
                    const double *table_v,
                    const double *table_f);

        /*! \brief Return the force at a distance r */
        double force(double r);

        /*! \brief Return the energy at a distance r */
        double energy(double r);

        /*! \brief Return table spacing (nm) */
        double tabSpacing() const { return tabSpacing_; }

        /*! \brief Return number of points in the table */
        unsigned int tabSize() const { return tableFDV0_.size()/4; }

        //! Return the force table
        const std::vector<double> &F() const { return tableFDV0_; }
        //! Return the energy table
        const std::vector<double> &V() const { return tableFDV0_; }
        //! Return the force + energy table
        const std::vector<double> &FDV0() const { return tableFDV0_; }
    private:
        //! Table scaling (nm)
        double              tabSpacing_;
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
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
