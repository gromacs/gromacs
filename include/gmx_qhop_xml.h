#ifndef _GMX_QHOP_XML_H
#define _GMX_QHOP_XML_H

#include "types/gmx_qhop_types.h"
/*#include "gmx_qhop_parm.h"
#include "gmx_qhop_db.h"*/


extern qhop_db_t qhops_read(char *fn);

extern void qhops_write(char *fn,qhop_db_t);

#endif
