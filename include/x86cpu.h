#define VENDOR_AMD   0x68747541
#define VENDOR_INTEL 0x756e6547
#define FLAGS_SUPPORT_SSE 0x02000000
#define FLAGS_SUPPORT_EXT_3DNOW 0xc0000000


/* Assembly routines in gmxcpuid.s */
void gmxcpuid(int,unsigned long *,unsigned long *,unsigned long *,unsigned long *);
void checksse();
void check3dnow();

/* Special optimization levels returned by check_x86cpu */
#define X86_NORMAL     0
#define X86_SSE        1
#define X86_3DNOW      2

int check_x86cpu(FILE *log);

