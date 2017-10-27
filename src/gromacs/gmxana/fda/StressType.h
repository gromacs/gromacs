/*
 * StressType.h
 *
 *  Created on: Apr 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef STRESSTYPE_H_
#define STRESSTYPE_H_

#include "EnumParser.h"
#include <string>
#include <iostream>

namespace fda_analysis {

enum StressType { PUNCTUAL, VIRIAL, VIRIAL_VON_MISES };

template <>
inline EnumParser<StressType>::EnumParser()
{
    map("punctual", PUNCTUAL);
    map("virial", VIRIAL);
    map("virial_von_mises", VIRIAL_VON_MISES);
}

} // namespace fda_analysis

#endif /* STRESSTYPE_H_ */
