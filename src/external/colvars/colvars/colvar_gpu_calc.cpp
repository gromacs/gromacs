#include "colvar_gpu_calc.h"
#include "colvarcomp.h"

using namespace colvars_gpu;

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int colvarcomp_gpu_graph::reset_graphs() {
  int error_code = COLVARS_OK;
  error_code |= graph_total_force.reset();
  error_code |= graph_calc_value.reset();
  error_code |= graph_calc_gradients.reset();
  error_code |= graph_calc_Jacobian_derivative.reset();
  return error_code;
}

int colvarcomp_gpu_graph::calc_value_gpu(bool requires_lattice) {
  int error_code = COLVARS_OK;
  if (!graph_calc_value.graph_exec_initialized) {
    error_code |= checkGPUError(cudaGraphCreate(&graph_calc_value.graph, 0));
    error_code |= add_calc_value_node();
    error_code |= graph_calc_value.init_graph_exec(cvmodule, cvc->get_stream());
    if (cvmodule->debug()) {
      const std::string filename = cvmodule->output_prefix() + "_" + cvc->name + "_calc_value.dot";
      error_code |= graph_calc_value.dump_graph(filename);
    }
  }
  if (graph_calc_value.graph_exec_initialized) {
    if (requires_lattice) {
      error_code |= checkGPUError(cudaStreamWaitEvent(
        cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::update_lattice)));
    }
    auto children = cvc->get_children();
    for (auto it = children.begin(); it != children.end(); ++it) {
      switch ((*it)->get_object_type()) {
        case colvardeps::object_t::colvarcomp: {
          // This branch is for future use, since for the time being Colvars doesn't
          // support nested CVCs as children.
          colvar::cvc* child = dynamic_cast<colvar::cvc*>(*it);
          error_code |= checkGPUError(cudaStreamWaitEvent(cvc->get_stream(), child->get_event(colvar::cvc::event_type::calc_value)));
          break;
        }
        case colvardeps::object_t::atom_group: {
          cvm::atom_group* child = dynamic_cast<cvm::atom_group*>(*it);
          // Wait for the atom group GPU calculations
          error_code |= checkGPUError(cudaStreamWaitEvent(
            cvc->get_stream(), child->get_gpu_atom_group()->get_event(
              colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
          break;
        }
        default: {
          return cvmodule->error("Unsupported colvardeps type in calc_value_gpu", COLVARS_BUG_ERROR);
        }
      }
    }
    // Launch the graph
    error_code |= checkGPUError(cudaGraphLaunch(graph_calc_value.graph_exec, cvc->get_stream()));
    // Record the event
    error_code |= checkGPUError(cudaEventRecord(cvc->get_event(colvar::cvc::event_type::calc_value), cvc->get_stream()));
  }
  return error_code;
}

int colvarcomp_gpu_graph::calc_gradients_gpu(bool requires_lattice) {
  int error_code = COLVARS_OK;
  if (!graph_calc_gradients.graph_exec_initialized) {
    error_code |= checkGPUError(cudaGraphCreate(&graph_calc_gradients.graph, 0));
    error_code |= add_calc_gradients_node();
    error_code |= graph_calc_gradients.init_graph_exec(cvmodule, cvc->get_stream());
    if (cvmodule->debug()) {
      const std::string filename = cvmodule->output_prefix() + "_" + cvc->name + "_calc_gradients.dot";
      error_code |= graph_calc_gradients.dump_graph(filename);
    }
  }
  if (graph_calc_gradients.graph_exec_initialized) {
    // We can assume that the calc_value_gpu() has been called, so we can directly
    // push the graph into the stream.
    error_code |= checkGPUError(cudaGraphLaunch(graph_calc_gradients.graph_exec, cvc->get_stream()));
    // Record the event
    error_code |= checkGPUError(cudaEventRecord(cvc->get_event(colvar::cvc::event_type::calc_gradients), cvc->get_stream()));
    // NOTE: Because gradients are stored in atom groups if f_cvc_explicit_gradient,
    // graph_calc_gradients may changed the d_atoms_grad, so I have to make the
    // atom groups streams wait for calc_gradients before any successive apply_force.
    for (auto agi = cvc->atom_groups.begin(); agi != cvc->atom_groups.end(); agi++) {
      error_code |= checkGPUError(cudaStreamWaitEvent(
        (*agi)->get_stream(), cvc->get_event(colvar::cvc::event_type::calc_gradients)));
    }
  }
  return error_code;
}

int colvarcomp_gpu_graph::calc_Jacobian_derivative_gpu(bool requires_lattice) {
  int error_code = COLVARS_OK;
  if (!graph_calc_Jacobian_derivative.graph_exec_initialized) {
    error_code |= checkGPUError(cudaGraphCreate(&graph_calc_Jacobian_derivative.graph, 0));
    error_code |= add_calc_Jacobian_derivative_node();
    error_code |= graph_calc_Jacobian_derivative.init_graph_exec(cvmodule, cvc->get_stream());
    if (cvmodule->debug()) {
      const std::string filename = cvmodule->output_prefix() + "_" + cvc->name + "_calc_Jacobian_derivative.dot";
      error_code |= graph_calc_Jacobian_derivative.dump_graph(filename);
    }
  }
  if (graph_calc_Jacobian_derivative.graph_exec_initialized) {
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::copy_atoms)));
    error_code |= checkGPUError(cudaGraphLaunch(graph_calc_Jacobian_derivative.graph_exec, cvc->get_stream()));
    error_code |= checkGPUError(cudaEventRecord(cvc->get_event(colvar::cvc::event_type::calc_Jacobian_derivative), cvc->get_stream()));
  }
  return error_code;
}

int colvarcomp_gpu_graph::calc_force_invgrads_gpu(bool requires_lattice) {
  int error_code = COLVARS_OK;
  if (!graph_total_force.graph_exec_initialized) {
    error_code |= checkGPUError(cudaGraphCreate(&graph_total_force.graph, 0));
    error_code |= add_calc_force_invgrads_node();
    error_code |= graph_total_force.init_graph_exec(cvmodule, cvc->get_stream());
    if (cvmodule->debug()) {
      const std::string filename = cvmodule->output_prefix() + "_" + cvc->name + "_calc_force_invgrads.dot";
      error_code |= graph_total_force.dump_graph(filename);
    }
  }
  if (graph_total_force.graph_exec_initialized) {
    // In case of "cvmodule->proxy->total_forces_valid() && (!is_enabled(f_cv_total_force_current_step))",
    // calc_cvc_total_force could be called at the first step, so we need
    // to wait for the default event. Even if it is not the case, it is always OK
    // to wait. If total forces have to be calculated after gradients and atom
    // positions, the CUDA graphs/kernels execution order is implicitly gauranteed
    // by the calling order of calc_cvc_values, calc_cvc_gradients, calc_cvc_Jacobians
    // and calc_cvc_total_force, as long as these operations are using the same
    // m_stream.
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::copy_atoms)));
    // Push the graph to stream
    error_code |= checkGPUError(cudaGraphLaunch(graph_total_force.graph_exec, cvc->get_stream()));
    // Record the event
    error_code |= checkGPUError(cudaEventRecord(cvc->get_event(colvar::cvc::event_type::calc_force_invgrads), cvc->get_stream()));
  }
  return error_code;
}
#endif
