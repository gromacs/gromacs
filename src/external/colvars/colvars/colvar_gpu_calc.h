#ifndef COLVAR_GPU_CALC_H
#define COLVAR_GPU_CALC_H

#include "colvar_gpu_support.h"
#include "colvarmodule.h"
#include "colvar.h"

namespace colvars_gpu {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
class colvarcomp_gpu_graph {
public:
  colvarcomp_gpu_graph(colvar::cvc* cvc_in, colvarmodule* cvmodule_in): cvc(cvc_in), cvmodule(cvmodule_in) {}
  virtual ~colvarcomp_gpu_graph() {}
  int reset_graphs();
  int calc_value_gpu(bool requires_lattice);
  int calc_gradients_gpu(bool requires_lattice);
  int calc_force_invgrads_gpu(bool requires_lattice);
  int calc_Jacobian_derivative_gpu(bool requires_lattice);
protected:
  gpu_graph_t graph_total_force;
  gpu_graph_t graph_calc_value;
  gpu_graph_t graph_calc_gradients;
  gpu_graph_t graph_calc_Jacobian_derivative;
  colvar::cvc* cvc;
  colvarmodule* cvmodule;
  /// \brief Calculate the variable on GPU
  virtual int add_calc_value_node() { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Calculate the atomic gradients, to be reused later in
  /// order to apply forces on GPU
  virtual int add_calc_gradients_node() { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Calculate the total force from the system using the
  /// inverse atomic gradients on GPU
  virtual int add_calc_force_invgrads_node() { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Calculate the divergence of the inverse atomic gradients on GPU
  virtual int add_calc_Jacobian_derivative_node() { return COLVARS_NOT_IMPLEMENTED; }
};
#endif
}

#endif // COLVAR_GPU_CALC_H
