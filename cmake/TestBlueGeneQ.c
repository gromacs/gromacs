int main()
{
    /* In theory, gcc supports BlueGeneQ, but inpractice its
       performance is miles behind IBM's compiler, so for now we only
       worry about the predefined macro of XLC */
#ifdef __bgq__
    return 0;
#else
#error This compiler is not targetting BlueGene/Q
#endif
}
