/*
 * $Id$
 *
 */

/* IA64 provides 4 different CPUID registers,
 * CPUID[x], where is from 0 to at least 4.
 * The argument to this function is the index of the desired
 * register, whose contents will be returned as a 64-bit integer.
 *
 * Normally we are interested in the CPU generation numbers in
 * CPUID[3]. The contents of this register is:
 *
 * Bits 63-40: Reserved
 * Bits 39-32: Architecture revision.
 * Bits 31-24: Processor family number.
 * Bits 23-16: Processor model number.
 * Bits 15-8:  Processor revision number.
 * Bits 7-0:   Number of available CPUID registers.
 *
 * As an example, my 1.3GHz Madison has
 */
unsigned long long
ia64_cpuid(int reg);

