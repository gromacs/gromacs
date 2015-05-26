#ifndef GMX_EWALD_TH_V_H
#define GMX_EWALD_TH_V_H

#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

/**
 * Thread vectors
 */

struct thread_vectors;

struct local_vectors {
  thread_vectors &tv;
  int thread;
  local_vectors(thread_vectors &tv, int thread) : tv(tv), thread(thread) { }

  void *&vptr(bool device, int id);

  template <typename V>
  V &create_or_resize(void *&p, int size) {
    if (p == NULL || (int) ((V *) p)->size() < size) {
      if (p == NULL) {
	if (size > 0) {
	  p = (void *) new V(size);
	}
      } else {
	if (size > 0) {
	  // unused warning, fixed in thrust 1.7.1
	  // ( https://github.com/thrust/thrust/releases/tag/1.7.1 )
	  //((V *) p)->resize(size);
	  delete (V *) p;
	  p = (void *) new V(size);
	} else {
	  delete (V *) p;
	  p = NULL;
	}
      }
    }
    return *((V *) p);
  }

  template <typename T>
  thrust::host_vector<T> &host(int id, int size) {
    size = size * 2 + 16; // over-allocate
    return create_or_resize< thrust::host_vector<T> >(vptr(false, id), size);
  }

  template <typename T>
  thrust::device_vector<T> &device(int id, int size) {
    size = size * 2 + 16; // over-allocate
    return create_or_resize< thrust::device_vector<T> >(vptr(true, id), size);
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


#endif // GMX_EWALD_TH_V_H
