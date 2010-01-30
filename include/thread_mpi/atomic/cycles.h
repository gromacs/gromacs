
/* check for cycle counters on supported platforms */
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__)  || defined(__PGIC__)) && (defined(__i386__) || defined(__x86_64__)))
#define TMPI_CYCLE_COUNT
/* x86 or x86-64 with GCC inline assembly */
typedef unsigned long long tmpi_cycles_t;

static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* x86 with GCC inline assembly - pentium TSC register */
    tmpi_cycles_t   cycle;
    unsigned       low,high;

    __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));

    cycle = ((unsigned long long)low) | (((unsigned long long)high)<<32);

    return cycle;
}
#elif (defined(__INTEL_COMPILER) && defined(__ia64__))
#define TMPI_CYCLE_COUNT
typedef unsigned long tmpi_cycles_t;
static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* Intel compiler on ia64 */
    return __getReg(_IA64_REG_AR_ITC);
}
#elif defined(__GNUC__) && defined(__ia64__)
#define TMPI_CYCLE_COUNT
typedef unsigned long tmpi_cycles_t;
static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* ia64 with GCC inline assembly */
    tmpi_cycles_t ret;
    __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(ret));
    return ret;
}
#elif defined(_MSC_VER)
#define TMPI_CYCLE_COUNT
typedef __int64 tmpi_cycles_t;
static __inline tmpi_cycles_t tmpi_cycles_read(void)
{
    return __rdtsc();
}
#endif

