int main()
{
#ifdef __bgq__
    return 0;
#else
#error This compiler is not targetting BlueGene/Q
#endif
}
