#if HAVE_MM_MALLOC_H
#    include <mm_malloc.h>
#elif HAVE_MALLOC_H
#    include <malloc.h>
#elif HAVE_XMMINTRIN_H
#    include <xmmintrin.h>
#endif

int
main()
{
    void *p = _mm_malloc(8*sizeof(float),64);
    _mm_free(p);
}
