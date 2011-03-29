#ifndef _TITRATION_DB_XML_H
#define _TITRATION_DB_XML_H

#include "titrationrec.h"

extern qhop_db_t qhops_read(char *fn);

extern void qhops_write(char *fn,qhop_db_t);

#endif
