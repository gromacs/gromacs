/*
 * CompatInteractionType.h
 *
 *  Created on: Jan 10, 2017
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_COMPATINTERACTIONTYPE_H_
#define SRC_GROMACS_FDA_COMPATINTERACTIONTYPE_H_

#include <type_traits>
#include "InteractionType.h"

namespace fda {

/// Used for indexing detailed forces array!
/// Don't change the numbers!
enum class CompatInteractionType : int
{
	NONE        = 0,
    BOND        = 1,
    ANGLE       = 2,
    DIHEDRAL    = 3,
    POLAR       = 4,
    LJ          = 5,
    COULOMB     = 6,
    ALL         = 7
};

using T = std::underlying_type<CompatInteractionType>::type;

constexpr T to_index(CompatInteractionType e)
{
    return static_cast<T>(e);
}

/// Conversion from InteractionType into CompatInteractionType
CompatInteractionType to_compat(InteractionType i);

/// Conversion from CompatInteractionType into InteractionType
InteractionType from_compat(CompatInteractionType i);

} // namespace fda

#endif /* SRC_GROMACS_FDA_COMPATINTERACTIONTYPE_H_ */
