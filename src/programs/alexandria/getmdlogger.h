#ifndef ALEXANDRIA_GETMDLOGGER_H
#define ALEXANDRIA_GETMDLOGGER_H
#include <cstdio>

#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/logger.h"

gmx::MDLogger getMdLogger(const t_commrec *cr,
                          FILE            *fplog);
                          
#endif
