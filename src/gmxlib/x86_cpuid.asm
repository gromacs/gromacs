;;
;;                 This source code is part of
;;
;;                  G   R   O   M   A   C   S
;;
;;           GROningen MAchine for Chemical Simulations
;;
;;                         VERSION 3.0
;;
;;  Copyright (c) 1991-2001
;;  BIOSON Research Institute, Dept. of Biophysical Chemistry
;;  University of Groningen, The Netherlands
;;
;;  This program is free software; you can redistribute it and/or
;;  modify it under the terms of the GNU General Public License
;;  as published by the Free Software Foundation; either version 2
;;  of the License, or (at your option) any later version.
;;
;;  If you want to redistribute modifications, please consider that
;;  scientific software is very special. Version control is crucial -
;;  bugs must be traceable. We will be happy to consider code for
;;  inclusion in the official distribution, but derived work must not
;;  be called official GROMACS. Details are found in the README & COPYING
;;  files - if they are missing, get the official version at www.gromacs.org.
;;
;;  To help us fund GROMACS development, we humbly ask that you cite
;;  the papers on the package - you can find them in the top README file.
;;
;;  Do check out http:	//www.gromacs.org , or mail us at gromacs@gromacs.org .
;;
;;  And Hey:
;;  GROup of MAchos and Cynical Suckers

;; this file must be processed with a version
;; of nasm that supports the extended 3dnow instructions.
;; you can find a binary of such a version on the
;; gromacs homepage.

segment .text
global x86_cpuid		;  issues the cpuid instruction with supplied args
x86_cpuid:	
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
