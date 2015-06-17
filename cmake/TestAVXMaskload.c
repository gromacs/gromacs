#include<immintrin.h>
int main()
{
    __m256d a;
    __m256i mask;
    double  d[4]={1,2,3,4};

    a = _mm256_setzero_pd();
    mask = _mm256_castpd_si256(a);

#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
    a = _mm256_maskload_pd(d,_mm256_castsi256_pd(mask));
#else
    a = _mm256_maskload_pd(d,mask);
#endif
    return 0;
}

