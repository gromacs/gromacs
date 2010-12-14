#include <string.h>
#include "types/qhoprec.h"
#include "qhop.h"

extern t_qhoprec *mk_qhoprec(void)
{
  t_qhoprec *qr;

  snew(qr,1);
  memset((void *)qr, 0, sizeof(t_qhoprec)); /* Just in case */
  return (qr);
}  /* mk_qhoprec */
