#include "colvar_gpu_support.h"
#include "colvarmodule.h"

#include <cstring>

#if defined(__has_include)
#if __cplusplus >= 202302L
# if __has_include(<stacktrace>)
#  define use_cpp_stacktrace
#  include <stacktrace>
# endif
#endif
#endif

namespace colvars_gpu {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int gpuAssert(cudaError_t code, const char *file, int line)
{
  if (code != cudaSuccess) {
    std::string error =
      std::string("GPUassert: ") +
      cudaGetErrorString(code) + " " + file + ", line " + cvm::to_str(line);
#ifdef use_cpp_stacktrace
#ifdef __cpp_lib_stacktrace
      error += "\nBacktrace:\n" + std::to_string(std::stacktrace::current()) + "\n";
#endif
#endif
    return cvm::error_static(error, COLVARS_ERROR);
  }
  return COLVARS_OK;
}

int add_clear_array_node_impl(
  void* dst, const size_t num_elements, const size_t sizeofT,
  cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
    /**< Size of each element in bytes. Must be 1, 2, or 4. */
  const size_t elementSize =
    (sizeofT % 4 == 0) ? 4 :
    ((sizeofT % 2 == 0) ? 2 : 1);
  const size_t width = num_elements * (sizeofT / elementSize);
  cudaMemsetParams memsetParams = {0};
  memsetParams.dst         = (void*)dst;
  memsetParams.value       = 0;
  memsetParams.elementSize = elementSize;
  memsetParams.width       = width;
  memsetParams.height      = 1;
  if (cvm::debug()) {
    cvm::log_static(
      "Add a memset clear node: ptr = " + cvm::to_str(dst) +
      " width = " + cvm::to_str(width) + " elementSize = " +
      cvm::to_str(elementSize) + "\n");
  }
  return checkGPUError(cudaGraphAddMemsetNode(
    &node_out, graph, dependencies.data(),
    dependencies.size(), &memsetParams));
}

int add_copy_node_impl(
  const void* src, void* dst, const size_t num_elements, const size_t sizeofT,
  cudaMemcpyKind kind, cudaGraphNode_t& node_out, cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  if (cvm::debug()) {
    cvm::log_static(
      "Add a memcpy node: src = " + cvm::to_str(src) +
      " dst = " + cvm::to_str(dst) + " num_elements = " +
      cvm::to_str(num_elements) + "\n");
  }
  return checkGPUError(cudaGraphAddMemcpyNode1D(
    &node_out, graph, dependencies.data(), dependencies.size(),
    dst, src, sizeofT * num_elements, kind));
}

int prepare_dependencies(
  const std::vector<std::pair<std::string, bool>>& node_names,
  std::vector<cudaGraphNode_t>& dependencies,
  const std::unordered_map<std::string, cudaGraphNode_t>& map,
  const std::string& caller_operation_name) {
  int error_code = COLVARS_OK;
  for (auto it = node_names.begin(); it != node_names.end(); ++it) {
    const std::string node_name = it->first;
    const bool allow_not_found = it->second;
    if (auto search = map.find(node_name); search != map.end()) {
      dependencies.push_back(search->second);
      if (cvm::debug()) {
        cvm::log_static("Operation " + caller_operation_name +
                " depends on node\" " + node_name + "\"\n");
      }
    } else {
      if (!allow_not_found) {
        error_code |= cvm::error_static("BUG: cannot find node " + node_name + "\n");
      } else {
        if (cvm::debug()) {
          cvm::log_static("Operation " + caller_operation_name +
                  " cannot depend on node\" " + node_name + "\"\n");
        }
      }
    }
  }
  return error_code;
}

gpu_graph_t::~gpu_graph_t() {
  (void)reset();
}

int gpu_graph_t::reset() {
  int error_code = COLVARS_OK;
  if (graph) {
    error_code |= checkGPUError(cudaGraphDestroy(graph));
    graph = nullptr;
  }
  if (graph_exec) {
    error_code |= checkGPUError(cudaGraphExecDestroy(graph_exec));
    graph_exec = nullptr;
  }
  nodes.clear();
  graph_exec_initialized = false;
  return error_code;
}

int gpu_graph_t::init_graph_exec(colvarmodule* cvmodule, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (graph_exec_initialized) return error_code;
  if (!nodes.empty()) {
    cudaGraphInstantiateParams params{0};
    params.flags = cudaGraphInstantiateFlagUpload;
    params.uploadStream = stream;
    error_code |= checkGPUError(cudaGraphInstantiateWithParams(&graph_exec, graph, &params));
    if (params.result_out == cudaGraphInstantiateSuccess) {
      graph_exec_initialized = true;
    } else {
      // Search for the error node
      for (const auto& node: nodes) {
        if (params.errNode_out == node.second) {
          cvmodule->log("Error on initialize node \"" + node.first + "\".\n");
        }
      }
      error_code |= cvmodule->error("Failed to instantiate CUDA graph.", COLVARS_BUG_ERROR);
    }
  }
  return error_code;
}

int gpu_graph_t::dump_graph(const std::string& filename) {
  cudaGraphDebugDotFlags dotFlags = cudaGraphDebugDotFlagsVerbose;
  return checkGPUError(cudaGraphDebugDotPrint(graph, filename.c_str(), dotFlags));
}
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}
