/*
 * InteractionType.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include "InteractionType.h"

namespace fda {

std::string to_string(InteractionType i)
{
    switch(i) {
        case InteractionType_NONE:
            return "none";
        case InteractionType_BOND:
            return "bond";
        case InteractionType_ANGLE:
            return "angle";
        case InteractionType_DIHEDRAL:
            return "dihedral";
        case InteractionType_POLAR:
            return "polar";
        case InteractionType_COULOMB:
            return "coulomb";
        case InteractionType_LJ:
            return "lj";
        case InteractionType_NB14:
            return "nb14";
        case InteractionType_BONDED:
            return "bonded";
        case InteractionType_NONBONDED:
            return "nonbonded";
        case InteractionType_ALL:
            return "all";
        default:
            return "invalid";
    }
}

InteractionType from_string(std::string const& s)
{
    std::string s_lower_case(s);
    std::transform(s_lower_case.begin(), s_lower_case.end(), s_lower_case.begin(), tolower);
    if (s_lower_case == "none")
        return InteractionType_NONE;
    if (s_lower_case == "bond")
        return InteractionType_BOND;
    else if (s_lower_case == "angle")
        return InteractionType_ANGLE;
    else if (s_lower_case == "dihedral")
        return InteractionType_DIHEDRAL;
    else if (s_lower_case == "polar")
        return InteractionType_POLAR;
    else if (s_lower_case == "coulomb")
        return InteractionType_COULOMB;
    else if (s_lower_case == "lj")
        return InteractionType_LJ;
    else if (s_lower_case == "nb14")
        return InteractionType_NB14;
    else if (s_lower_case == "bonded")
        return InteractionType_BONDED;
    else if (s_lower_case == "nonbonded")
        return InteractionType_NONBONDED;
    else if (s_lower_case == "all")
        return InteractionType_ALL;
    else
        throw std::runtime_error("InteractionType " + s + "unknown");
}

} // namespace fda
