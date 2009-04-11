
int
main()
{
#if defined __QK_USER__
  return 0;
#else
#  error not catamount
#endif
}
