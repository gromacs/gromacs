/*
 * ResultFormat.h
 *
 *  Created on: Apr 14, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef RESULTFORMAT_H_
#define RESULTFORMAT_H_

#include "EnumParser.h"
#include <string>
#include <iostream>

namespace fda_analysis {

enum ResultFormat { UNKNOWN, PDB, DIMACS };

template <>
inline EnumParser<ResultFormat>::EnumParser()
{
    map("unknown", UNKNOWN);
    map("pdb", PDB);
    map("dimacs", DIMACS);
}

} // namespace fda_analysis

#endif /* RESULTFORMAT_H_ */
