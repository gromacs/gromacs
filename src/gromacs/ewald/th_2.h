template <typename T>
struct gpu_vector {
  T *host_ptr();
  T *dev_ptr();
  void cpy_to_dev();
  void cpy_to_host();
};

struct thread_vectors {
};
