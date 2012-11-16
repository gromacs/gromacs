#ifdef USE_VISIBILITY  /* off by default */
#define GMX_LIBGMX_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBMD_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBGMXANA_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBGMXPREPROCESS_EXPORT __attribute__((__visibility__("default")))
#else
#define GMX_LIBGMX_EXPORT
#define GMX_LIBMD_EXPORT
#define GMX_LIBGMXANA_EXPORT
#define GMX_LIBGMXPREPROCESS_EXPORT
#endif
