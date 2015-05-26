#include <cuda.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "thread_mpi/mutex.h"

const bool check_verbose = false;
static tMPI::mutex print_mutex;

template <typename T>
void check(const char *name, T *data, T *expected, int size, gmx_bool bDevice)
{
  print_mutex.lock();
  bool bDiff = false;
  for (int i = 0; i < size; ++i) {
    T cpu_v = expected[i];
    T gpu_v;
    if (bDevice) {
      cudaMemcpy(&gpu_v, &data[i], sizeof(T), cudaMemcpyDeviceToHost);
    } else {
      gpu_v = data[i];
    }
    T diff = gpu_v - cpu_v;
    if (check_verbose) {
      fprintf(stderr, " %d:%f(%f)", i, (double) cpu_v, (double) diff);
    }
    if (diff != 0) {
      if (!bDiff) {
	fprintf(stderr, "%s\n", name);
	bDiff = true;
      }
      T absdiff = diff > 0 ? diff : -diff;
      T abscpu_v = cpu_v > 0 ? cpu_v : -cpu_v;
      T reldiff = absdiff / (abscpu_v > 1e-11 ? abscpu_v : 1e-11);
      if (reldiff > .000001) {
	fprintf(stderr, "%.0fppm", (double) (reldiff * 1e6));
	if (reldiff > .0001) {
	  fprintf(stderr, " value %f vs %f ", (double) cpu_v, (double) gpu_v);
	}
      } else {
	fprintf(stderr, "~");
      }
    }
  }
  if (bDiff) {
    fprintf(stderr, "\n");
  }
  print_mutex.unlock();
}

void check_int(const char *name, int *data, int *expected, int size, gmx_bool bDevice)
{
  check(name, data, expected, size, bDevice);
}

void check_real(const char *name, real *data, real *expected, int size, gmx_bool bDevice)
{
  check(name, data, expected, size, bDevice);
}

void print_lock() {
  print_mutex.lock();
}

void print_unlock() {
  print_mutex.lock();
}
