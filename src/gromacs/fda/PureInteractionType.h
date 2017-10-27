/*
 * PureInteractionType.h
 *
 *  Created on: Nov 24, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_
#define SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_

#include <type_traits>
#include "InteractionType.h"

namespace fda {

/// Used for indexing detailed forces array!
/// Don't change the numbers!
enum class PureInteractionType : int
{
    BOND        = 0,
    ANGLE       = 1,
    DIHEDRAL    = 2,
    POLAR       = 3,
    COULOMB     = 4,
    LJ          = 5,
    NB14        = 6,
    NUMBER      = 7
};

using T = std::underlying_type<PureInteractionType>::type;

constexpr T to_index(PureInteractionType e)
{
    return static_cast<T>(e);
}

/// Conversion from InteractionType into PureInteractionType
PureInteractionType to_pure(InteractionType i);

/// Conversion from PureInteractionType into InteractionType
InteractionType from_pure(PureInteractionType i);

} // namespace fda

#endif /* SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_ */
