int main()
{
#ifdef __aarch64__
    return 0;
#else
#error This compiler is not targetting AArch64
#endif
}
