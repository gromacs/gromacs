int main()
{
#ifdef __MIC__
    return 0;
#else
#error This compiler is not targetting MIC
#endif
}
