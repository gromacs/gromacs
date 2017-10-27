/*
 * InteractionType.h
 *
 *  Created on: Nov 18, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_INTERACTIONTYPE_H_
#define SRC_GROMACS_FDA_INTERACTIONTYPE_H_

#include <string>

namespace fda {

using InteractionType = int;

static const int InteractionType_NONE      =      0;
static const int InteractionType_BOND      = 1 << 0;
static const int InteractionType_ANGLE     = 1 << 1;
static const int InteractionType_DIHEDRAL  = 1 << 2;
static const int InteractionType_POLAR     = 1 << 3;
static const int InteractionType_COULOMB   = 1 << 4;
static const int InteractionType_LJ        = 1 << 5;
static const int InteractionType_NB14      = 1 << 6;
static const int InteractionType_BONDED    = InteractionType_BOND + InteractionType_ANGLE + InteractionType_DIHEDRAL;
static const int InteractionType_NONBONDED = InteractionType_COULOMB + InteractionType_LJ + InteractionType_NB14;
static const int InteractionType_ALL       = InteractionType_BONDED + InteractionType_NONBONDED + InteractionType_POLAR;

/// Convert string into InteractionType
std::string to_string(InteractionType i);

/// Convert InteractionType into string
InteractionType from_string(std::string const& s);

} // namespace fda

#endif /* SRC_GROMACS_FDA_INTERACTIONTYPE_H_ */
