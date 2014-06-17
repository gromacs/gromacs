/* This source file has the sole purpose to force C++ linking of the gpu_utils
 * static archive, otherwise the exception handling code generated inside
 * memtestG80 will cause undefined symbols at linking. */

/* Avoid warnings about empty object files */
int gpu_utils_dummy;
