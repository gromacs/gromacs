#include <immintrin.h>

int main()
{
    __m256 x  = _mm256_set1_ps(0.5);
    x = _mm256_fmadd_ps(x,x,x);
    return 0;
}
