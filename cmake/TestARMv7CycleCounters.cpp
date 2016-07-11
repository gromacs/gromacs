int main()
{
#if defined(__ARM_ARCH_7A__) && defined(__GNUC__)
    unsigned int cycles_lo, cycles_hi;
    asm volatile("mrrc p15, 1, %0, %1, c14" : "=r" (cycles_lo), "=r" (cycles_hi));

    // Return 0 (success) if low or high 32 bits contained anything non-trivial
    return !(cycles_lo > 0 || cycles_hi > 0);
#else
#error This architecture/compiler does not support ARMv7 32-bit cycle counters
#endif
}
