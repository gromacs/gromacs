extern void test_tng(void);
extern void test_zlib(void);

int main(int argc, char* argv[])
{
    test_tng();
    test_zlib();
    return 0;
}
