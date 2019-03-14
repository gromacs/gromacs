int main()
{
#if defined (__i386__) || defined (__x86_64__) || defined (_M_IX86) || defined (_M_X64)
    return 0;
#else
#error This is not x86
#endif
}
