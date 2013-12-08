#include <smmintrin.h>

int main()
{
    __m128 x  = _mm_set1_ps(0.5);
    x = _mm_dp_ps(x,x,0x77);
    return 0;
}
