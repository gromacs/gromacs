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
/*! \libinternal\file
 * \brief Declares class for spline interpolation tables
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 */
#ifndef GMX_TABLES_SPLINETABLE_H
#define GMX_TABLES_SPLINETABLE_H

#include <vector>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

/*! \brief Function pointer to fill table.
 *
 * Function pointer used to fill the table.
 * \param[in]  param A vector of constants used to evaluate the function
 * \param[in]  r     The distance,
 * \param[out] e     The returned energy
 * \param[out] f     The force.
 */
typedef void (*interactionPotentialFunction)(const std::vector<double> &param,
                                             double                     r,
                                             double                    *e,
                                             double                    *f);

namespace gmx
{
/*! \libinternal
 * \brief Class containing quadratic spline interaction table.
 *
 * This class generates input for quadratic spline table interpolation
 * from input tables.
 * If both energy and force tables are input the number of output
 * points is equal to the input. When no force table input is passed
 * the number of points in the table is halved because an extra interpolation
 * is needed. In general it is more accurate to provide both energy and force
 * tables.
 *
 * Functions using the tables should use methods tabSpacing() and tabSize()
 * to get the actual spacing and number of points.
 *
 * Note that the table input should be smooth, e.g. not containing noise
 * e.g. from an (iterative) Boltzmann inversion procedure.
 *
 * The class implements completely separate tables for the energy
 * and force for both Coulomb and Van der Waals and a combined
 * table which is more efficient to evaluate in the mdrun kernels.
 */
class QuadraticSplineTable
{
    public:
        /*! \brief Constructor that fills table using user table(s).
         *
         * Create a spline interpolation table from user provided
         * energies and, optionally, forces. The number of points and
         * the spacing in the final table is determined by whether forces
         * are present. If they are, the number of points is the same as
         * on input, otherwise the number is halved. Similarly spacing
         * will be doubled when no forces are present. This behavior is
         * different from the other constructor below.
         *
         * \param[in] dx     Space between data points in input
         * \param[in] tableV Table containing energies
         * \param[in] tableF Table containing forces (may be empty)
         */
        QuadraticSplineTable(double                     dx,
                             const std::vector<double> &tableV,
                             const std::vector<double> &tableF);

        /*! \brief Constructor that fills table using analytical function
         *
         * This constructor generates a quadratic spline table and fills it.
         * Input is generated using an analytical function, passed
         * as function pointer vAnalytical with parameters params.
         * \param[in] nDataPoints Length of input tables.
         * \param[in] dx          Spacing between data points in
         *                        input tables.
         * \param[in] vAnalytical A function pointer that can be used
         *                        for computing energies
         * \param[in] params      Parameters for the function vAnalytical
         */
        QuadraticSplineTable(std::size_t                   nDataPoints,
                             double                        dx,
                             interactionPotentialFunction  vAnalytical,
                             const std::vector<double>    &params);

        /*! \brief Return the force at a distance r */
        real force(real r) const;

        /*! \brief Return the energy at a distance r */
        real energy(real r) const;

        /*! \brief Return table spacing (nm) */
        real spacing() const { return tabSpacing_; }

        /*! \brief Return number of points in the table */
        std::size_t size() const { return tabSize_; }

        //! Return the force table
        const std::vector<real> &F() const { return tableF_; }

        //! Return the energy table
        const std::vector<real> &V() const { return tableV_; }

        //! Return the force + energy table
        const std::vector<real> &FDV0() const { return tableFDV0_; }

    private:
        /*! \brief Do the actual filling of the tables
         *
         * \param[in] stride    Distance between data points
         * \param[in] tableInV  Table containing energies
         * \param[in] tableInF  Table containing forces (may be NULL)
         */
        void fillTable(std::size_t                stride,
                       const std::vector<double> &tableInV,
                       const std::vector<double> &tableInF);
        //! Table scaling (nm)
        double              tabSpacing_;
        //! Length of table (number of points)
        std::size_t         tabSize_;
        //! Force table
        std::vector<real>   tableF_;
        //! Energy table
        std::vector<real>   tableV_;
        /*! \brief Combined array for force and energy
         *
         * Force+energy table, size of array is 4 times longer
         * than the other tables.
         * Quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
         * this is used with single precision x86 SIMD for aligned loads
         */
        std::vector<real> tableFDV0_;

        GMX_DISALLOW_COPY_AND_ASSIGN(QuadraticSplineTable);
};

}      // namespace gmx

#endif /* GMX_TABLES_TABLES_H */
