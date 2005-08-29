
AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"


if GMX_IA32_SSE
  IA32_SSE        = nb_kernel_ia32_sse
  IA32_SSE_OBJ    = nb_kernel_ia32_sse/libnb_kernel_ia32_sse.la
endif

if GMX_IA32_3DNOW
  IA32_3DNOW      = nb_kernel_ia32_3dnow
  IA32_3DNOW_OBJ  = nb_kernel_ia32_3dnow/libnb_kernel_ia32_3dnow.la
endif

if GMX_IA32_SSE2
  IA32_SSE2       = nb_kernel_ia32_sse2
  IA32_SSE2_OBJ   = nb_kernel_ia32_sse2/libnb_kernel_ia32_sse2.la
endif

if GMX_X86_64_SSE
  X86_64_SSE       = nb_kernel_x86_64_sse
  X86_64_SSE_OBJ   = nb_kernel_x86_64_sse/libnb_kernel_x86_64_sse.la
endif

if GMX_X86_64_SSE2
  X86_64_SSE2      = nb_kernel_x86_64_sse2
  X86_64_SSE2_OBJ  = nb_kernel_x86_64_sse2/libnb_kernel_x86_64_sse2.la
endif

if GMX_PPC_ALTIVEC
  PPC_ALTIVEC     = nb_kernel_ppc_altivec
  PPC_ALTIVEC_OBJ = nb_kernel_ppc_altivec/libnb_kernel_ppc_altivec.la
endif

if GMX_IA64_ASM_SINGLE
  IA64_SINGLE     = nb_kernel_ia64_single
  IA64_SINGLE_OBJ = nb_kernel_ia64_single/libnb_kernel_ia64_single.la
endif

if GMX_IA64_ASM_DOUBLE
  IA64_DOUBLE     = nb_kernel_ia64_double
  IA64_DOUBLE_OBJ = nb_kernel_ia64_double/libnb_kernel_ia64_double.la
endif



SUBDIRS =       $(IA32_SSE)     $(IA32_SSE2)    $(IA32_3DNOW)           \
                $(X86_64_SSE)   $(X86_64_SSE2)  $(PPC_ALTIVEC)          \
                $(IA64_SINGLE)  $(IA64_DOUBLE)	nb_kernel


# Convenience library for nonbonded interactions - not installed.

noinst_LTLIBRARIES = libnonbonded.la

# Subdirectories will contain libtool convenience libraries.
# These are not installed, but added to the main library.
# We also add conditional objects directly to the main library.
libnonbonded_la_LIBADD = \
        $(IA32_SSE_OBJ)         $(IA32_SSE2_OBJ)        $(IA32_3DNOW_OBJ)  \
        $(X86_64_SSE_OBJ)       $(X86_64_SSE2_OBJ)      $(PPC_ALTIVEC_OBJ) \
        $(IA64_SINGLE_OBJ)      $(IA64_DOUBLE_OBJ)	\
	nb_kernel/libnb_kernel.la			


libnonbonded_la_SOURCES = \
	nb_kerneltype.h			nonbonded.c	\
	nb_free_energy.c		nb_free_energy.h


