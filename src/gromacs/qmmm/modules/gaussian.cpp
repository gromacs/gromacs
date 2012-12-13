#include "gaussian.h"

class QmmmGaussian
{
  QmmmGaussian(const t_commrec *cr,
               const matrix box,
               const gmx_mtop_t *mtop,
               const t_inputrec *ir,
               const t_forcerec *fr)
  {
  }

  ~QmmmGaussian()
  {
  }

  void update(const t_commrec *cr,
              const t_forcerec *fr,
              const rvec x[],
              const t_mdatoms *md,
              const matrix box,
              const gmx_localtop_t *top)
  {
  }

  real calculate(const t_commrec *cr,
                 const rvec x[],
                 rvec f[],
                 const t_forcerec *fr,
                 const t_mdatoms *md)
  {
    return 0;
  }
};
