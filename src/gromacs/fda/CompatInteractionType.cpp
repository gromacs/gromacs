/*
 * CompatInteractionType.h
 *
 *  Created on: Jan 10, 2017
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <stdexcept>
#include "CompatInteractionType.h"

namespace fda {

CompatInteractionType to_compat(InteractionType i)
{
    switch(i) {
        case InteractionType_BOND:
            return CompatInteractionType::BOND;
        case InteractionType_ANGLE:
            return CompatInteractionType::ANGLE;
        case InteractionType_DIHEDRAL:
            return CompatInteractionType::DIHEDRAL;
        case InteractionType_POLAR:
            return CompatInteractionType::POLAR;
        case InteractionType_COULOMB:
            return CompatInteractionType::COULOMB;
        case InteractionType_LJ:
            return CompatInteractionType::LJ;
        case InteractionType_NB14:
            return CompatInteractionType::DIHEDRAL;
        default:
        	return CompatInteractionType::NONE;
    }
}

InteractionType from_compat(CompatInteractionType i)
{
    switch(i) {
        case CompatInteractionType::BOND:
            return InteractionType_BOND;
        case CompatInteractionType::ANGLE:
            return InteractionType_ANGLE;
        case CompatInteractionType::DIHEDRAL:
            return InteractionType_DIHEDRAL;
        case CompatInteractionType::POLAR:
            return InteractionType_POLAR;
        case CompatInteractionType::COULOMB:
            return InteractionType_COULOMB;
        case CompatInteractionType::LJ:
            return InteractionType_LJ;
        default:
            throw std::runtime_error("Is not a compat interaction");
    }
}

} // namespace fda
