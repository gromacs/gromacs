/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <m.ghahremanpour@hotmail.com>
 */

#include "optparam.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

#include "gromacs/random.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

OptParam::OptParam(){};

void OptParam::Init(const char *xvgconv, const char *xvgepot, const gmx_output_env_t *oenv, real seed,
                    real step, int maxiter, int nprint, real temperature, gmx_bool bBound)
{

    xvgconv_     = xvgconv;
    xvgepot_     = xvgepot;
    oenv_        = oenv;
    bBound_      = bBound;
    seed_        = seed;
    step_        = step;
    maxiter_     = maxiter;
    nprint_      = nprint;
    temperature_ = temperature;
}

void OptParam::setSeed(real seed)
{
    seed_ = seed;
}

void OptParam::setMaxiter(int maxiter)
{
    maxiter_ = maxiter;
}

void OptParam::setNprint(int nprint)
{
    nprint_ = nprint;
}

void OptParam::setStep(real step)
{
    step_ = step;
}

void OptParam::setTemperature(real temperature)
{
    temperature_ = temperature;
}
}
