;#
;# $Id$
;# 
;#                This source code is part of
;# 
;#                 G   R   O   M   A   C   S
;# 
;#          GROningen MAchine for Chemical Simulations
;# 
;#                        VERSION 3.1
;# Copyright (c) 1991-2001, University of Groningen, The Netherlands
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;# 
;# If you want to redistribute modifications, please consider that
;# scientific software is very special. Version control is crucial -
;# bugs must be traceable. We will be happy to consider code for
;# inclusion in the official distribution, but derived work must not
;# be called official GROMACS. Details are found in the README & COPYING
;# files - if they are missing, get the official version at www.gromacs.org.
;# 
;# To help us fund GROMACS development, we humbly ask that you cite
;# the papers on the package - you can find them in the top README file.
;# 
;# For more info, check our website at http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
 

;#
;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).

;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as (NASM doesn't understant this
;# if syntax), but NASM doesn't see it, so the code inside is only 
;# read by NASM (NASM doesn't understand .if):

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
; .endif  # End of NASM-specific block

; .intel_syntax noprefix   # Line only read by gnu as

		
.section .text
	
.globl x86_cpuid	;# issues the cpuid instruction with supplied args 
.globl _x86_cpuid
x86_cpuid:	
_x86_cpuid:	
	push ebp
	mov  ebp,esp
	push edi
	push ebx
	push ecx
	push edx
	mov  eax, [ebp+8]	
	cpuid
	mov  edi, [ebp+12]
	mov  [edi],eax
	mov  edi, [ebp+16]
	mov  [edi],ebx
	mov  edi, [ebp+20]
	mov  [edi],ecx
	mov  edi, [ebp+24]
	mov  [edi],edx
	pop edx
	pop ecx
	pop ebx
	pop edi
	mov esp, ebp
	pop ebp
	ret










