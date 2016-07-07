int main()
{
#ifdef __aarch64__
    return 0;
#else
#error This compiler is not targetting 32-bit ARMv7
#endif
}
