#ifndef _GMX_GPU_UTILS_H_
#define _GMX_GPU_UTILS_H_

#ifndef __cplusplus
extern "C" {
#endif

/** Runs a quick memory test and returns 0 in case if no error is detected. */
int do_quick_memtest(int /*dev_id*/);

/** Runs a full memory test and returns 0 in case if no error is detected. */
int do_full_memtest(int /*dev_id*/);

/** Runs a time constrained memory test and returns 0 in case if no error is detected. */
int do_timed_memtest(int /*dev_id*/, int /*time_limit*/);

/** Checks whether the GPU with the given device id is supported or not. */
int is_supported_cuda_gpu(int /*dev_id*/, char* /*gpu_name*/);

#ifndef __cplusplus
}  /* extern "C" */
#endif

#endif // _GMX_GPU_UTILS_H_

