/*
 * ParticleType.h
 *
 *  Created on: Apr 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef PARTICLETYPE_H_
#define PARTICLETYPE_H_

#include "EnumParser.h"
#include <string>
#include <iostream>

namespace fda_analysis {

enum ParticleType { ATOM, RESIDUE };

template <>
inline EnumParser<ParticleType>::EnumParser()
{
    map("atom", ATOM);
    map("residue", RESIDUE);
}

} // namespace fda_analysis

#endif /* PARTICLETYPE_H_ */
