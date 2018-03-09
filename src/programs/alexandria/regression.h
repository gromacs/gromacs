/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#ifndef ALEXANDRIA_REGRESSION_H
#define ALEXANDRIA_REGRESSION_H

#include <vector>

#include "gromacs/math/vectypes.h"

class MatrixWrapper
{
    public:
        /*! \brief Constructor
         *
         * \param[in] ncolumn Number of columns in the matrix
         * \param[in] nrow    Number of rows
         */
        MatrixWrapper(int ncolumn, int nrow);

        //! \brief Destructor in charge of freeing memory
        ~MatrixWrapper();

        //!\brief Return number of columns
        int nColumn() const { return ncolumn_; }

        //!\brief Return number of rows
        int nRow() const { return nrow_; }

        /*! \brief Set a value in the matrix
         *
         * \param[in] col    The column
         * \param[in] row    The row
         * \param[in] value  The value
         */
        void set(int col, int row, double value);

        /*! \brief Set a row of values in the matrix
         *
         * \param[in] row    The row number
         * \param[in] value  The values. Array should be ncolumn long.
         */
        void setRow(int row, const double value[]);

        /*! \brief Get a value from the matrix
         *
         * \param[in] col    The column
         * \param[in] row    The row
         * \returns The value
         */
        double get(int col, int row) const;

        /*! \brief Solve a matrix equation A solution = rhs
         *
         * \param[in] rhs  Vector of right hand side values
         * \param[out] solution Pointer to vector of solution
         */
        void solve(std::vector<double> rhs, std::vector<double> *solution);
    private:
        // Number of rows
        int      nrow_;
        // Number of columns
        int      ncolumn_;
        // The internal data structure
        double **a_;
};

void kabsch_rotation(tensor p, tensor q, tensor rotated_p);

#endif
