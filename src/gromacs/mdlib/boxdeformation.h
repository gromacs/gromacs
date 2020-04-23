/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declares interface to box deformation code.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDRUN_BOXDEFORMATION_H
#define GMX_MDRUN_BOXDEFORMATION_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxmpi.h"

struct t_inputrec;
enum class DDRole;
enum class NumRanks;

namespace gmx
{

class BoxDeformation
{
public:
    //! Trivial constructor.
    BoxDeformation(double timeStep, int64_t initialStep, const tensor& deformationTensor, const matrix& referenceBox);

    //! Deform \c x and \c box at this \c step;
    void apply(ArrayRef<RVec> x, matrix box, int64_t step);

private:
    //! The integrator time step.
    double timeStep_;
    //! The initial step number (from the .tpr, which permits checkpointing to work correctly).
    int64_t initialStep_;
    //! Non-zero elements provide a scaling factor for deformation in that box dimension.
    tensor deformationTensor_;
    //! The initial box, ie from the .tpr file.
    matrix referenceBox_;

    GMX_DISALLOW_COPY_AND_ASSIGN(BoxDeformation);
};

/*! \brief Factory function for box deformation module.
 *
 * If the \c inputrec specifies the use of box deformation during the
 * update phase, communicates the \c initialBox from SIMMASTER to
 * other ranks, and constructs and returns an object to manage that
 * update.
 *
 * \throws NotImplementedError if the \c inputrec specifies an
 * unsupported combination.
 */
std::unique_ptr<BoxDeformation> prepareBoxDeformation(const matrix&     initialBox,
                                                      DDRole            ddRole,
                                                      NumRanks          numRanks,
                                                      MPI_Comm          communicator,
                                                      const t_inputrec& inputrec);

} // namespace gmx

#endif
