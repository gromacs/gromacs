#ifndef GMX_EWALD_CHECK_H
#define GMX_EWALD_CHECK_H

struct gpu_flags;

bool run_on_cpu(const gpu_flags &flags);
bool run_on_gpu(const gpu_flags &flags);
bool check_vs_cpu(const gpu_flags &flags);
bool check_vs_cpu_verbose(const gpu_flags &flags);


void check_int(const char *name, int *data, int *expected, int size, gmx_bool bDevice);
void check_real(const char *name, real *data, real *expected, int size, gmx_bool bDevice);
void print_lock();
void print_unlock();

#endif // GMX_EWALD_CHECK_H
