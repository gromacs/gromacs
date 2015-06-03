#ifndef GMX_EWALD_TH_V2_H
#define GMX_EWALD_TH_V2_H

#include <vector>

/**
 * Thread vectors
 */

enum vector_loc {
  LOC_CPU,
  LOC_PINNED,
  LOC_GPU,
};

template <typename T>
struct simple_vector {
  void *p;
  size_t size;
  void loc_alloc(vector_loc loc) {
    switch (loc) {
    case LOC_CPU: p = new T(size); break;
    case LOC_PINNED: p = cudaHostMalloc(sizeof(T) * size, cudaMallocPinned); break; // FIX FIX FIX
    case LOC_GPU: p = cudaMalloc(sizeof(T) * size); break; // FIX FIX
    }
  }
  void loc_free(vector_loc loc) {
    switch (loc) {
    case LOC_CPU: delete[] (T *) p; break;
    case LOC_PINNED: cudaFree(p); break;
    case LOC_GPU: p = cudaFree(p); break;
    }
    p = NULL;
  }
};

struct thread_vectors;

struct local_vectors {
  thread_vectors &tv;
  int thread;
  local_vectors(thread_vectors &tv, int thread) : tv(tv), thread(thread) { }

  void *&vptr(bool device, int id);

  template <typename V>
  V &create_or_recreate(void *&p, int size, vector_loc loc) {
    if (p == NULL || (int) ((V *) p)->size < size) {
      if (p == NULL) {
	if (size > 0) {
	  p = (void *) new V();
          ((V *) p)->size = size;
          ((V *) p)->loc_alloc(loc);
	}
      } else {
	if (size > 0) {
          ((V *) p)->loc_free(loc);
	  delete (V *) p;
	  p = (void *) new V();
          ((V *) p)->size = size;
          ((V *) p)->loc_alloc(loc);
	} else {
          ((V *) p)->loc_free(loc);
	  delete (V *) p;
	  p = NULL;
	}
      }
    }
    return *((V *) p);
  }

  template <typename T>
  T *host(int id, int size) {
    size = size * 2 + 16; // over-allocate
    return create_or_recreate< simple_vector<T> >(vptr(false, id), size,
						  LOC_CPU);
  };

  template <typename T>
  T *pinned(int id, int size) {
    size = size * 2 + 16; // over-allocate
    return create_or_recreate< simple_vector<T> >(vptr(false, id), size,
						  LOC_PINNED);
  }

  template <typename T>
  T *gpu(int id, int size) {
    size = size * 2 + 16; // over-allocate
    return create_or_recreate< simple_vector<T> >(vptr(false, id), size,
						  LOC_GPU);
  }

  template <typename T>
  T *cpy_to_gpu(int id, int size) {
    T *cpu = pinned(id, size);
    T *gpu = gpu(id, size);
    cudaMemcpy(gpu, cpu, sizeof(T) * size); // TODO stream
    return gpu;
  }

  template <typename T>
  T *cpy_to_host(int id, int size) {
    T *gpu = gpu(id, size);
    T *cpu = pinned(id, size);
    cudaMemcpy(cpu, gpu, sizeof(T) * size); // TODO stream
    return cpu;
  }
};

struct thread_vectors {
  int nthread, nid;
  std::vector<void *> v;

  thread_vectors(int nthread, int nid) {
    ensure(nthread, nid);
  }

  void ensure(int nthread, int nid) {
    v.resize(nthread * nid * 2);
    this->nthread = nthread;
    this->nid = nid;
  }

  local_vectors local(int thread) {
    return local_vectors(*this, thread);
  }
};

inline void *&local_vectors::vptr(bool device, int id) {
  int idx = (thread + (device ? tv.nthread : 0)) * tv.nid + id;
  return tv.v[idx];
}

#endif // GMX_EWALD_TH_V2_H
