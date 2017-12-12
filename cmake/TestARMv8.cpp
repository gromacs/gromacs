int main()
{
#ifdef __aarch64__
    return 0;
#else
#error This compiler is not targetting 64-bit ARMv8
#endif
}
