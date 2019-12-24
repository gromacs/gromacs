/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
#ifndef CHARGEMODEL_H
#define CHARGEMODEL_H

#include "gmxpre.h"

#include <string>
#include <map>

#include "gromacs/utility/basedefinitions.h"

namespace alexandria
{

/*! \brief
 * Enumerated type holding the charge models used in PolData
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeModel {
    eqdESP_p    = 0,
    eqdESP_pg   = 1,
    eqdESP_ps   = 2,
    eqdESP_pp   = 3,
    eqdACM_g    = 4,
    eqdACM_pg   = 5,
    eqdACM_ps   = 6,
    eqdYang     = 7,
    eqdBultinck = 8,
    eqdRappe    = 9,
    eqdNR       = 10
};

/*! \brief
 * Enumerated type holding the charge generation algorithms
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeGenerationAlgorithm {
    eqgNONE, 
    eqgACM, 
    eqgESP
};

static bool gmx_unused getEemtypePolarizable(ChargeModel eem)
{
    return (eem == eqdACM_pg ||
            eem == eqdESP_pp || eem == eqdESP_pg || eem == eqdESP_ps);
}

static bool gmx_unused getEemtypeDistributed(ChargeModel eem)
{
    return (eem == eqdACM_g  || eem == eqdACM_pg || eem == eqdACM_ps ||
            eem == eqdESP_pg || eem == eqdESP_ps ||
            eem == eqdYang   || eem == eqdRappe);
}

static bool gmx_unused getEemtypeSlater(ChargeModel eem)
{
    return (eem == eqdACM_ps || eem == eqdESP_ps ||
            eem == eqdYang   || eem == eqdRappe);
}

static bool gmx_unused getEemtypeGaussian(ChargeModel eem)
{
    return (eem == eqdACM_g || eem == eqdACM_pg || eem == eqdESP_pg);
}

/* Return the charge generation algorithm */
static ChargeGenerationAlgorithm gmx_unused chargeGenerationAlgorithm(ChargeModel eem)
{
    if (eem == eqdESP_p  || eem == eqdESP_pg ||
        eem == eqdESP_ps || eem == eqdESP_pp)
    {
        return eqgESP;
    }
    else
    {
        return eqgACM;
    }
}

//! \brief Return the string corresping to eem
const char *getEemtypeName(ChargeModel eem);

//! \brief Return the ChargeModel corresponding to name
ChargeModel name2eemtype(const std::string &name);

} // namespace aleaxndria
#endif
