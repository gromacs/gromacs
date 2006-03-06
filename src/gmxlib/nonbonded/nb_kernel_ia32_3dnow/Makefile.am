
AM_CPPFLAGS= -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

noinst_LTLIBRARIES = libnb_kernel_ia32_3dnow.la

libnb_kernel_ia32_3dnow_la_SOURCES = \
	nb_kernel010_ia32_3dnow.s	nb_kernel010_ia32_3dnow.h	\
	nb_kernel030_ia32_3dnow.s	nb_kernel030_ia32_3dnow.h	\
	nb_kernel100_ia32_3dnow.s	nb_kernel100_ia32_3dnow.h	\
	nb_kernel101_ia32_3dnow.s	nb_kernel101_ia32_3dnow.h	\
	nb_kernel102_ia32_3dnow.s	nb_kernel102_ia32_3dnow.h	\
	nb_kernel110_ia32_3dnow.s	nb_kernel110_ia32_3dnow.h	\
	nb_kernel111_ia32_3dnow.s	nb_kernel111_ia32_3dnow.h	\
	nb_kernel112_ia32_3dnow.s	nb_kernel112_ia32_3dnow.h	\
	nb_kernel300_ia32_3dnow.s	nb_kernel300_ia32_3dnow.h	\
	nb_kernel301_ia32_3dnow.s	nb_kernel301_ia32_3dnow.h	\
	nb_kernel302_ia32_3dnow.s	nb_kernel302_ia32_3dnow.h	\
	nb_kernel310_ia32_3dnow.s	nb_kernel310_ia32_3dnow.h	\
	nb_kernel311_ia32_3dnow.s	nb_kernel311_ia32_3dnow.h	\
	nb_kernel312_ia32_3dnow.s	nb_kernel312_ia32_3dnow.h	\
	nb_kernel330_ia32_3dnow.s	nb_kernel330_ia32_3dnow.h	\
	nb_kernel331_ia32_3dnow.s	nb_kernel331_ia32_3dnow.h	\
	nb_kernel332_ia32_3dnow.s	nb_kernel332_ia32_3dnow.h	\
	nb_kernel_ia32_3dnow_test_asm.s	nb_kernel_ia32_3dnow_test_asm.h	\
	nb_kernel_ia32_3dnow.c		nb_kernel_ia32_3dnow.h

#
# Kernels formatted in Intel assembly syntax (instead of AT&T).
# Currently, these are only necessary for NASM under Windows, which
# does not use the automake/autoconf build system.
#
EXTRA_DIST = \
	nb_kernel010_ia32_3dnow.intel_syntax.s	\
	nb_kernel030_ia32_3dnow.intel_syntax.s	\
	nb_kernel100_ia32_3dnow.intel_syntax.s	\
	nb_kernel101_ia32_3dnow.intel_syntax.s	\
	nb_kernel102_ia32_3dnow.intel_syntax.s	\
	nb_kernel110_ia32_3dnow.intel_syntax.s	\
	nb_kernel111_ia32_3dnow.intel_syntax.s	\
	nb_kernel112_ia32_3dnow.intel_syntax.s	\
	nb_kernel300_ia32_3dnow.intel_syntax.s	\
	nb_kernel301_ia32_3dnow.intel_syntax.s	\
	nb_kernel302_ia32_3dnow.intel_syntax.s	\
	nb_kernel310_ia32_3dnow.intel_syntax.s	\
	nb_kernel311_ia32_3dnow.intel_syntax.s	\
	nb_kernel312_ia32_3dnow.intel_syntax.s	\
	nb_kernel330_ia32_3dnow.intel_syntax.s	\
	nb_kernel331_ia32_3dnow.intel_syntax.s	\
	nb_kernel332_ia32_3dnow.intel_syntax.s	\
	nb_kernel_ia32_3dnow_test_asm.intel_syntax.s

