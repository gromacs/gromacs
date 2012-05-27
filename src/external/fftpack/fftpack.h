#include "types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

#define Treal real

    void cffti1(int n, Treal wa[], int ifac[]);
    void cfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[], int isign);
    void rffti1(int n, Treal wa[], int ifac[]);
    void rfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[]);
    void rfftb1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[]);

#ifdef __cplusplus
}
#endif
