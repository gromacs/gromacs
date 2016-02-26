/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_SQLITE3_H
#define MOLPROP_SQLITE3_H

#include "molprop.h"

void ReadSqlite3(const char                       *sqlite_file,
                 std::vector<alexandria::MolProp> &mp);

#endif
