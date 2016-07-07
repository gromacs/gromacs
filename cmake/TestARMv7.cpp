int main()
{
#ifdef __ARM_ARCH_7A__
    return 0;
#else
#error This compiler is not targetting 32-bit ARMv7
#endif
}
