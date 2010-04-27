#ifndef _GMX_GPU_UTILS_H_
#define _GMX_GPU_UTILS_H_

#ifndef __cplusplus
extern "C" {
#endif

int do_quick_memtest(int /*dev_id*/);

int do_full_memtest(int /*dev_id*/);

int do_timed_memtest(int /*dev_id*/, int /*time_limit*/);

int is_supported_cuda_gpu(int /*dev_id*/, char* /*gpu_name*/);

#ifndef __cplusplus
}  /* extern "C" */
#endif

#endif // _GMX_GPU_UTILS_H_

