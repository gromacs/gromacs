/*
 * PureInteractionType.h
 *
 *  Created on: Nov 24, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <stdexcept>
#include "PureInteractionType.h"

namespace fda {

PureInteractionType to_pure(InteractionType i)
{
    switch(i) {
        case InteractionType_BOND:
            return PureInteractionType::BOND;
        case InteractionType_ANGLE:
            return PureInteractionType::ANGLE;
        case InteractionType_DIHEDRAL:
            return PureInteractionType::DIHEDRAL;
        case InteractionType_POLAR:
            return PureInteractionType::POLAR;
        case InteractionType_COULOMB:
            return PureInteractionType::COULOMB;
        case InteractionType_LJ:
            return PureInteractionType::LJ;
        case InteractionType_NB14:
            return PureInteractionType::NB14;
        default:
            throw std::runtime_error("Is not a pure interaction");
    }
}

InteractionType from_pure(PureInteractionType i)
{
    switch(i) {
        case PureInteractionType::BOND:
            return InteractionType_BOND;
        case PureInteractionType::ANGLE:
            return InteractionType_ANGLE;
        case PureInteractionType::DIHEDRAL:
            return InteractionType_DIHEDRAL;
        case PureInteractionType::POLAR:
            return InteractionType_POLAR;
        case PureInteractionType::COULOMB:
            return InteractionType_COULOMB;
        case PureInteractionType::LJ:
            return InteractionType_LJ;
        case PureInteractionType::NB14:
            return InteractionType_NB14;
        default:
            throw std::runtime_error("Is not a pure interaction");
    }
}

} // namespace fda
