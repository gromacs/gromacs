int main()
{
#if defined (__FUJITSU) && ( defined(__sparc) || defined(__sparcv9) ) && ( defined(__LP64__) || defined(__arch64) )
    return 0;
#else
#error This compiler is not targetting Fujitsu Sparc64
#endif
}
