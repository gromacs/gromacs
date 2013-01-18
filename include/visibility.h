#ifdef USE_VISIBILITY  /* off by default */
#if defined _WIN32 || defined __CYGWIN__ || defined WINDOWS
#ifdef gmx_EXPORTS
#define GMX_LIBGMX_EXPORT __declspec(dllexport)
#else
#define GMX_LIBGMX_EXPORT __declspec(dllimport)
#endif
#ifdef md_EXPORTS
#define GMX_LIBMD_EXPORT __declspec(dllexport)
#else
#define GMX_LIBMD_EXPORT __declspec(dllimport)
#endif
#ifdef gmxana_EXPORTS
#define GMX_LIBGMXANA_EXPORT __declspec(dllexport)
#else
#define GMX_LIBGMXANA_EXPORT __declspec(dllimport)
#endif
#ifdef gmxpreprocess_EXPORTS
#define GMX_LIBGMXPREPROCESS_EXPORT __declspec(dllexport)
#else
#define GMX_LIBGMXPREPROCESS_EXPORT __declspec(dllimport)
#endif
#else /* Unix */
#define GMX_LIBGMX_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBMD_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBGMXANA_EXPORT __attribute__((__visibility__("default")))
#define GMX_LIBGMXPREPROCESS_EXPORT __attribute__((__visibility__("default")))
#endif
#else /* no USE_VISIBILITY */
#define GMX_LIBGMX_EXPORT
#define GMX_LIBMD_EXPORT
#define GMX_LIBGMXANA_EXPORT
#define GMX_LIBGMXPREPROCESS_EXPORT
#endif /* USE_VISIBILITY */
