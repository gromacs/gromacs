/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2016,2019, by the GROMACS development team, led by
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
/*! \libinternal
 * \file
 * \brief
 * Declares routine for computing many correlation functions using OpenMP
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_MANYAUTOCORRELATION_H
#define GMX_MANYAUTOCORRELATION_H

#include <vector>

#include "gromacs/fft/fft.h"
#include "gromacs/utility/real.h"

/*! \brief
 * Perform many autocorrelation calculations.
 *
 * This routine performs many autocorrelation function calculations using FFTs.
 * The GROMACS FFT library wrapper is employed. On return the c vector contain
 * a symmetric function that is useful for further FFT:ing, for instance in order to
 * compute spectra.
 *
 * The vectors c[i] should all have the same length, but this is not checked for.
 *
 * The c arrays will be extend and filled with zero beyond ndata before
 * computing the correlation.
 *
 * The functions uses OpenMP parallellization.
 *
 * \param[inout] c Data array
 * \return fft error code, or zero if everything went fine (see fft/fft.h)
 * \throws gmx::InconsistentInputError if the input is inconsistent.
 */
int many_auto_correl(std::vector<std::vector<real>>* c);

#endif
