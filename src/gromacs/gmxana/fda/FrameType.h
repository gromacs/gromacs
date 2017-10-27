/*
 * FrameType.h
 *
 *  Created on: Apr 2, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef FRAMETYPE_H_
#define FRAMETYPE_H_

#include "EnumParser.h"
#include <string>
#include <iostream>

namespace fda_analysis {

enum FrameType { SINGLE, ALL, AVERAGE, SKIP };

template <>
inline EnumParser<FrameType>::EnumParser()
{
    map("single", SINGLE);
    map("all", ALL);
    map("average", AVERAGE);
    map("skip", SKIP);
}

} // namespace fda_analysis

#endif /* FRAMETYPE_H_ */
