#include <cctype>
#include <cmath>

#include "../gmx_lapack.h"


float F77_FUNC(slange, SLANGE)(const char* norm, int* m, int* n, float* a, int* lda, float* work)
{
    const char ch = std::toupper(*norm);
    float      dtemp, sum, max, val, scale;
    int        i, j;

    switch (ch)
    {
        case 'M':
            max = 0.0;
            for (j = 0; j < *n; j++)
                for (i = 0; i < *m; i++)
                {
                    dtemp = std::abs(a[j * (*lda) + i]);
                    if (dtemp > max)
                        max = dtemp;
                }
            val = max;
            break;

        case 'O':
        case '1':
            max = 0.0;
            for (j = 0; j < *n; j++)
            {
                sum = 0.0;
                for (i = 0; i < *m; i++)
                    sum += std::abs(a[j * (*lda) + i]);
                if (sum > max)
                    max = sum;
            }
            val = max;
            break;

        case 'I':
            for (i = 0; i < *m; i++)
                work[i] = 0.0;
            for (j = 0; j < *n; j++)
                for (i = 0; i < *m; i++)
                    work[i] += std::abs(a[j * (*lda) + i]);
            max = 0;
            for (i = 0; i < *m; i++)
                if (work[i] > max)
                    max = work[i];
            val = max;
            break;

        case 'F':
        case 'E':
            scale = 0.0;
            sum   = 1.0;
            i     = 1;
            for (j = 0; j < *n; j++)
                F77_FUNC(slassq, SLASSQ)(m, &(a[j * (*lda) + 0]), &i, &scale, &sum);
            val = scale * std::sqrt(sum);
            break;

        default: val = 0.0; break;
    }
    return val;
}
