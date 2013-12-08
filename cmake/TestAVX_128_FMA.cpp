#include <immintrin.h>
#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h>
#endif
#ifdef HAVE_INTRIN_H
#include <intrin.h>
#endif

int main()
{
    __m128 x  = _mm_set1_ps(0.5);
    x = _mm_macc_ps(x,x,x);
    return 0;
}
