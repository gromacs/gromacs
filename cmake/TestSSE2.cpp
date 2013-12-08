#include <xmmintrin.h>

int main()
{
    __m128 x  = _mm_set1_ps(0.5);
    x = _mm_rsqrt_ps(x);
    return 0;
}
