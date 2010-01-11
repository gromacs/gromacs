#ifndef _GMX_QHOP_XML_H
#define _GMX_QHOP_XML_H

#include "types/gmx_qhop_types.h"
/*#include "gmx_qhop_parm.h"
#include "gmx_qhop_db.h"*/


extern void qhops_read(char *fn, qhop_db_t qdb);

extern void qhops_write(char *fn,int nqhop,qhop_t qht[]);

#endif
