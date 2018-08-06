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
//
// Created by Eric Irrgang on 10/30/17.
//

#ifndef GROMACS_RESTRAINTFUNCTOR_H
#define GROMACS_RESTRAINTFUNCTOR_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

/*! \libinternal \file
 * \brief Provide a functor class to clarify calculation input, output, and point in time.
 */

struct pull_t;
struct t_mdatoms;
struct t_pbc;
struct t_commrec;


namespace gmx
{

namespace restraint
{

/*!
 * \brief Implement class for applying the restraint calculations.
 *
 * This class implements the stage of MD integration at which the contributions from restraints
 * (such as for structural refinement or free-energy calculation) are calculated.
 *
 * Example:
 *
 *     void MdIntegrator::applyConstraints()
 *     {
 *          auto calculator = RestraintFunctor(pull_work_, mdatoms_, pbc_, commRec_, time_, lambda_, position_);
 *          calculator.calculate(&energy_, &force_, &virial_, &work_);
 *     };
 */
class RestraintFunctor
{
    public:
        /*!
         * \brief Construct the functor with necessary input parameters.
         * \param pull
         * \param atoms
         * \param pbc
         * \param cr
         * \param t
         * \param lambda
         * \param x
         */
        RestraintFunctor(const pull_t &pull, const t_mdatoms &atoms, const t_pbc &pbc, const t_commrec &cr, double t, real lambda, const RVec &x);

        void calculate(real* energy, RVec* force, tensor* virial, real* work);

    private:

};

}      // end namespace gmx::restraint

}      // end namespace gmx

#endif //GROMACS_RESTRAINTFUNCTOR_H
