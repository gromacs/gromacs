#include <math.h>
#include "typedefs.h"
#include "gstat.h"
#include "vec.h"

void init_lsq(t_lsq *lsq)
{
  lsq->yy=lsq->yx=lsq->xx=lsq->sx=lsq->sy=0.0;
  lsq->np=0;
}

int npoints_lsq(t_lsq *lsq)
{
  return lsq->np;
}

void done_lsq(t_lsq *lsq)
{
  init_lsq(lsq);
}

void add_lsq(t_lsq *lsq,real x,real y)
{
  lsq->yy+=y*y;
  lsq->yx+=y*x;
  lsq->xx+=x*x;
  lsq->sx+=x;
  lsq->sy+=y;
  lsq->np++;
}

void get_lsq_ab(t_lsq *lsq,real *a,real *b)
{
  double yx,xx,sx,sy;
  
  yx=lsq->yx/lsq->np;
  xx=lsq->xx/lsq->np;
  sx=lsq->sx/lsq->np;
  sy=lsq->sy/lsq->np;
  
  (*a)=(yx-sx*sy)/(xx-sx*sx);
  (*b)=(sy)-(*a)*(sx);
}

real aver_lsq(t_lsq *lsq)
{
  if (lsq->np == 0)
    fatal_error(0,"No points in distribution\n");
  
  return (lsq->sy/lsq->np);
}

real sigma_lsq(t_lsq *lsq)
{
  if (lsq->np == 0)
    fatal_error(0,"No points in distribution\n");
    
  return sqrt(lsq->yy/lsq->np - sqr(lsq->sy/lsq->np));
}

real error_lsq(t_lsq *lsq)
{
  return sigma_lsq(lsq)/sqrt(lsq->np);
}

