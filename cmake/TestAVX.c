#include <immintrin.h>

int main()
{
    __m256 x  = _mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
    x = _mm256_rsqrt_ps(x);
    return 0;
}
