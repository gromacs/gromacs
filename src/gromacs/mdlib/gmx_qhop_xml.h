#ifndef _GMX_QHOP_XML_H
#define _GMX_QHOP_XML_H

#include "gmx_qhop_parm.h"

extern gmx_qhop_t *gmx_qhops_read(char *fn,int *nqhop);

extern void gmx_qhops_write(char *fn,int nqhop,gmx_qhop_t qht);

#endif
