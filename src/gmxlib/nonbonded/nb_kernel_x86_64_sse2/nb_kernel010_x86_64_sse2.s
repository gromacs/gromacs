;#
;# $Id$
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
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
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
.equiv          .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as



.globl nb_kernel010_x86_64_sse2	 
.globl _nb_kerne010_x86_64_sse2	 
nb_kernel010_x86_64_sse2:	
_nb_kernel010_x86_64_sse2:	
.equiv          nb010_fshift,           16
.equiv          nb010_gid,              24
.equiv          nb010_pos,              32
.equiv          nb010_faction,          40
.equiv          nb010_charge,           48
.equiv          nb010_p_facel,          56
.equiv          nb010_argkrf,           64
.equiv          nb010_argcrf,           72
.equiv          nb010_Vc,               80
.equiv          nb010_type,             88
.equiv          nb010_p_ntype,          96
.equiv          nb010_vdwparam,         104
.equiv          nb010_Vvdw,             112
.equiv          nb010_p_tabscale,       120
.equiv          nb010_VFtab,            128
.equiv          nb010_invsqrta,         136
.equiv          nb010_dvda,             144
.equiv          nb010_p_gbtabscale,     152
.equiv          nb010_GBtab,            160
.equiv          nb010_p_nthreads,       168
.equiv          nb010_count,            176
.equiv          nb010_mtx,              184
.equiv          nb010_outeriter,        192
.equiv          nb010_inneriter,        200
.equiv          nb010_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb010_ix,               0
.equiv          nb010_iy,               16
.equiv          nb010_iz,               32
.equiv          nb010_dx,               48
.equiv          nb010_dy,               64
.equiv          nb010_dz,               80
.equiv          nb010_two,              96
.equiv          nb010_c6,               112
.equiv          nb010_c12,              128
.equiv          nb010_six,              144
.equiv          nb010_twelve,           160
.equiv          nb010_Vvdwtot,          176
.equiv          nb010_fix,              192
.equiv          nb010_fiy,              208
.equiv          nb010_fiz,              224
.equiv          nb010_half,             240
.equiv          nb010_three,            256
.equiv          nb010_nri,              272
.equiv          nb010_iinr,             280
.equiv          nb010_jindex,           288
.equiv          nb010_jjnr,             296
.equiv          nb010_shift,            304
.equiv          nb010_shiftvec,         312
.equiv          nb010_facel,            320
.equiv          nb010_innerjjnr,        328
.equiv          nb010_is3,              336
.equiv          nb010_ii3,              340
.equiv          nb010_ntia,             344
.equiv          nb010_innerk,           348
.equiv          nb010_n,                352
.equiv          nb010_nn1,              356
.equiv          nb010_ntype,            360
.equiv          nb010_nouter,           364
.equiv          nb010_ninner,           368
	push rbp
	mov  rbp, rsp
	push rbx
	
	femms
	sub rsp, 392		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb010_nouter], eax
	mov [rsp + nb010_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb010_nri], edi
	mov [rsp + nb010_iinr], rsi
	mov [rsp + nb010_jindex], rdx
	mov [rsp + nb010_jjnr], rcx
	mov [rsp + nb010_shift], r8
	mov [rsp + nb010_shiftvec], r9
	mov rdi, [rbp + nb010_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb010_ntype], edi

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double two IEEE (hex)
	mov ebx, 0x40000000
	mov [rsp + nb010_two], eax
	mov [rsp + nb010_two + 4], ebx
	movsd xmm1, [rsp + nb010_two]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm2, xmm1
	addpd  xmm2, xmm1       ;# 4.0
	addpd  xmm2, xmm1       ;# six
	movapd xmm3, xmm2
	addpd  xmm3, xmm3       ;# twelve
	movapd [rsp + nb010_two], xmm1
	movapd [rsp + nb010_six],  xmm2
	movapd [rsp + nb010_twelve], xmm3

.nb010_threadloop:
        mov   rsi, [rbp + nb010_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb010_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb010_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb010_n], eax
        mov [rsp + nb010_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb010_outerstart
        jmp .nb010_end

.nb010_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb010_nouter]
	mov [rsp + nb010_nouter], ebx

.nb010_outer:
	mov   rax, [rsp + nb010_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 

	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb010_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb010_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb010_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

    	mov   rdx, [rbp + nb010_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb010_ntype]
    	shl   edx, 1
    	mov   [rsp + nb010_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb010_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb010_ix], xmm0
	movapd [rsp + nb010_iy], xmm1
	movapd [rsp + nb010_iz], xmm2

	mov   [rsp + nb010_ii3], ebx
	
	;# clear Vvdwtot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb010_Vvdwtot], xmm4
	movapd [rsp + nb010_fix], xmm4
	movapd [rsp + nb010_fiy], xmm4
	movapd [rsp + nb010_fiz], xmm4
	
	mov   rax, [rsp + nb010_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb010_pos]
	mov   rdi, [rbp + nb010_faction]	
	mov   rax, [rsp + nb010_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb010_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb010_ninner]
	mov   [rsp + nb010_ninner], ecx
	add   edx, 0
	mov   [rsp + nb010_innerk], edx    ;# number of innerloop atoms 
	
	jge   .nb010_unroll_loop
	jmp   .nb010_checksingle
.nb010_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add   qword ptr [rsp + nb010_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb010_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb010_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb010_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	movd  ebx, mm1	
	
	movapd [rsp + nb010_c6], xmm4
	movapd [rsp + nb010_c12], xmm6

	mov rsi, [rbp + nb010_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	

	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]

	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb010_ix]
	movapd xmm5, [rsp + nb010_iy]
	movapd xmm6, [rsp + nb010_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [rsp + nb010_dx], xmm4 
	movapd [rsp + nb010_dy], xmm5
	movapd [rsp + nb010_dz], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6	;# rsq in xmm4 

	cvtpd2ps xmm6, xmm4	
	rcpps xmm6, xmm6
	cvtps2pd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb010_two]
	movapd xmm5, xmm4
	mulpd xmm4, xmm6	;# lu*rsq 
	subpd xmm0, xmm4	;# 2-lu*rsq 
	mulpd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb010_two]
	mulpd xmm5, xmm6	;# lu*rsq 
	subpd xmm0, xmm5	;# 2-lu*rsq 
	mulpd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulpd  xmm1, [rsp + nb010_c6]
	mulpd  xmm2, [rsp + nb010_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb010_Vvdwtot]
	mulpd  xmm1, [rsp + nb010_six]
	mulpd  xmm2, [rsp + nb010_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm0, xmm2	;# xmm4=total fscal 
	movapd xmm4, xmm0
	
	movapd xmm0, [rsp + nb010_dx]
	movapd xmm1, [rsp + nb010_dy]
	movapd xmm2, [rsp + nb010_dz]

	movapd [rsp + nb010_Vvdwtot], xmm5

	mov    rdi, [rbp + nb010_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [rsp + nb010_fix]
	movapd xmm4, [rsp + nb010_fiy]
	movapd xmm5, [rsp + nb010_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [rsp + nb010_fix], xmm3
	movapd [rsp + nb010_fiy], xmm4
	movapd [rsp + nb010_fiz], xmm5
	
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5
	
	;# should we do one more iteration? 
	sub   dword ptr [rsp + nb010_innerk],  2
	jl    .nb010_checksingle
	jmp   .nb010_unroll_loop
.nb010_checksingle:				
	mov   edx, [rsp + nb010_innerk]
	and   edx, 1
	jnz    .nb010_dosingle
	jmp    .nb010_updateouterdata
.nb010_dosingle:
	mov rdi, [rbp + nb010_pos]
	mov   rcx, [rsp + nb010_innerjjnr]
	mov   eax, [rcx]		

	movd  mm0, eax		;# use mmx registers as temp storage 	
	mov rsi, [rbp + nb010_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb010_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb010_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	
	movapd [rsp + nb010_c6], xmm4
	movapd [rsp + nb010_c12], xmm6
	
	mov rsi, [rbp + nb010_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	

	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb010_ix]
	movapd xmm5, [rsp + nb010_iy]
	movapd xmm6, [rsp + nb010_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [rsp + nb010_dx], xmm4 
	movapd [rsp + nb010_dy], xmm5
	movapd [rsp + nb010_dz], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6	;# rsq in xmm4 

	cvtsd2ss xmm6, xmm4	
	rcpss xmm6, xmm6
	cvtss2sd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb010_two]
	movapd xmm5, xmm4
	mulsd xmm4, xmm6	;# lu*rsq 
	subsd xmm0, xmm4	;# 2-lu*rsq 
	mulsd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb010_two]
	mulsd xmm5, xmm6	;# lu*rsq 
	subsd xmm0, xmm5	;# 2-lu*rsq 
	mulsd xmm0, xmm6	;# xmm0=rinvsq 
	movapd xmm4, xmm0
	
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulsd  xmm1, [rsp + nb010_c6]
	mulsd  xmm2, [rsp + nb010_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb010_Vvdwtot]
	mulsd  xmm1, [rsp + nb010_six]
	mulsd  xmm2, [rsp + nb010_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm4, xmm2	;# xmm4=total fscal 

	movapd xmm0, [rsp + nb010_dx]
	movapd xmm1, [rsp + nb010_dy]
	movapd xmm2, [rsp + nb010_dz]

	movlpd [rsp + nb010_Vvdwtot], xmm5

	mov    rdi, [rbp + nb010_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [rsp + nb010_fix]
	movlpd xmm4, [rsp + nb010_fiy]
	movlpd xmm5, [rsp + nb010_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [rsp + nb010_fix], xmm3
	movlpd [rsp + nb010_fiy], xmm4
	movlpd [rsp + nb010_fiz], xmm5
	
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb010_updateouterdata:
	mov   ecx, [rsp + nb010_ii3]
	mov   rdi, [rbp + nb010_faction]
	mov   rsi, [rbp + nb010_fshift]
	mov   edx, [rsp + nb010_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb010_fix]
	movapd xmm1, [rsp + nb010_fiy]
	movapd xmm2, [rsp + nb010_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [rsi + rdx*8]
	movsd  xmm4, [rsi + rdx*8 + 8]
	movsd  xmm5, [rsi + rdx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb010_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb010_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb010_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb010_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 

        ;# finish if last 
        mov ecx, [rsp + nb010_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb010_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb010_n], esi
        jmp .nb010_outer
.nb010_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb010_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb010_end
        ;# non-zero, do one more workunit
        jmp   .nb010_threadloop
.nb010_end:
	
	mov eax, [rsp + nb010_nouter]
	mov ebx, [rsp + nb010_ninner]
	mov rcx, [rbp + nb010_outeriter]
	mov rdx, [rbp + nb010_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 392
	femms

	pop rbx
	pop	rbp
	ret



	


	
.globl nb_kernel010nf_x86_64_sse2
.globl _nb_kernel010nf_x86_64_sse2
nb_kernel010nf_x86_64_sse2:	
_nb_kernel010nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb010nf_fshift,         16
.equiv          nb010nf_gid,            24
.equiv          nb010nf_pos,            32
.equiv          nb010nf_faction,        40
.equiv          nb010nf_charge,         48
.equiv          nb010nf_p_facel,        56
.equiv          nb010nf_argkrf,         64
.equiv          nb010nf_argcrf,         72
.equiv          nb010nf_Vc,             80
.equiv          nb010nf_type,           88
.equiv          nb010nf_p_ntype,        96
.equiv          nb010nf_vdwparam,       104
.equiv          nb010nf_Vvdw,           112
.equiv          nb010nf_p_tabscale,     120
.equiv          nb010nf_VFtab,          128
.equiv          nb010nf_invsqrta,       136
.equiv          nb010nf_dvda,           144
.equiv          nb010nf_p_gbtabscale,   152
.equiv          nb010nf_GBtab,          160
.equiv          nb010nf_p_nthreads,     168
.equiv          nb010nf_count,          176
.equiv          nb010nf_mtx,            184
.equiv          nb010nf_outeriter,      192
.equiv          nb010nf_inneriter,      200
.equiv          nb010nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb010nf_ix,             0
.equiv          nb010nf_iy,             16
.equiv          nb010nf_iz,             32
.equiv          nb010nf_two,            48
.equiv          nb010nf_c6,             64
.equiv          nb010nf_c12,            80
.equiv          nb010nf_Vvdwtot,        96
.equiv          nb010nf_half,           112
.equiv          nb010nf_three,          128
.equiv          nb010nf_nri,            144
.equiv          nb010nf_iinr,           152
.equiv          nb010nf_jindex,         160
.equiv          nb010nf_jjnr,           168
.equiv          nb010nf_shift,          176
.equiv          nb010nf_shiftvec,       184
.equiv          nb010nf_innerjjnr,      192
.equiv          nb010nf_is3,            200
.equiv          nb010nf_ii3,            204
.equiv          nb010nf_ntia,           208
.equiv          nb010nf_innerk,         212
.equiv          nb010nf_n,              216
.equiv          nb010nf_nn1,            220
.equiv          nb010nf_ntype,          224
.equiv          nb010nf_nouter,         228
.equiv          nb010nf_ninner,         232
	push rbp
	mov  rbp, rsp
	push rbx

	femms
	sub rsp, 248		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb010nf_nouter], eax
	mov [rsp + nb010nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb010nf_nri], edi
	mov [rsp + nb010nf_iinr], rsi
	mov [rsp + nb010nf_jindex], rdx
	mov [rsp + nb010nf_jjnr], rcx
	mov [rsp + nb010nf_shift], r8
	mov [rsp + nb010nf_shiftvec], r9
	mov rdi, [rbp + nb010nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb010nf_ntype], edi

        ;# create constant floating-point factors on stack
        mov eax, 0x00000000     ;# lower half of double two IEEE (hex)
        mov ebx, 0x40000000
        mov [rsp + nb010nf_two], eax
        mov [rsp + nb010nf_two + 4], ebx
        movsd xmm1, [rsp + nb010nf_two]
        shufpd xmm1, xmm1, 0    ;# splat to all elements
        movapd [rsp + nb010nf_two], xmm1

.nb010nf_threadloop:
        mov   rsi, [rbp + nb010nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb010nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb010nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb010nf_n], eax
        mov [rsp + nb010nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb010nf_outerstart
        jmp .nb010nf_end

.nb010nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb010nf_nouter]
	mov [rsp + nb010nf_nouter], ebx

.nb010nf_outer:
	mov   rax, [rsp + nb010nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 

	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb010nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb010nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

    	mov   rdx, [rbp + nb010nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb010nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb010nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb010nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb010nf_ix], xmm0
	movapd [rsp + nb010nf_iy], xmm1
	movapd [rsp + nb010nf_iz], xmm2

	mov   [rsp + nb010nf_ii3], ebx
	
	;# clear Vvdwtot 
	xorpd xmm4, xmm4
	movapd [rsp + nb010nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb010nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb010nf_pos]
	mov   rax, [rsp + nb010nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb010nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb010nf_ninner]
	mov   [rsp + nb010nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb010nf_innerk], edx    ;# number of innerloop atoms 
	
	jge   .nb010nf_unroll_loop
	jmp   .nb010nf_checksingle
.nb010nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb010nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add   qword ptr [rsp + nb010nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb010nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb010nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb010nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	movd  ebx, mm1	
	
	movapd [rsp + nb010nf_c6], xmm4
	movapd [rsp + nb010nf_c12], xmm6

	mov rsi, [rbp + nb010nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	

	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]

	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb010nf_ix]
	movapd xmm5, [rsp + nb010nf_iy]
	movapd xmm6, [rsp + nb010nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6	;# rsq in xmm4 

	cvtpd2ps xmm6, xmm4	
	rcpps xmm6, xmm6
	cvtps2pd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb010nf_two]
	movapd xmm5, xmm4
	mulpd xmm4, xmm6	;# lu*rsq 
	subpd xmm0, xmm4	;# 2-lu*rsq 
	mulpd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb010nf_two]
	mulpd xmm5, xmm6	;# lu*rsq 
	subpd xmm0, xmm5	;# 2-lu*rsq 
	mulpd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulpd  xmm1, [rsp + nb010nf_c6]
	mulpd  xmm2, [rsp + nb010nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb010nf_Vvdwtot]
	movapd [rsp + nb010nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub   dword ptr [rsp + nb010nf_innerk],  2
	jl    .nb010nf_checksingle
	jmp   .nb010nf_unroll_loop
.nb010nf_checksingle:				
	mov   edx, [rsp + nb010nf_innerk]
	and   edx, 1
	jnz    .nb010nf_dosingle
	jmp    .nb010nf_updateouterdata
.nb010nf_dosingle:
	mov rdi, [rbp + nb010nf_pos]
	mov   rcx, [rsp + nb010nf_innerjjnr]
	mov   eax, [rcx]		

	movd  mm0, eax		;# use mmx registers as temp storage 	
	mov rsi, [rbp + nb010nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb010nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb010nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	
	movapd [rsp + nb010nf_c6], xmm4
	movapd [rsp + nb010nf_c12], xmm6
	
	mov rsi, [rbp + nb010nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	

	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb010nf_ix]
	movapd xmm5, [rsp + nb010nf_iy]
	movapd xmm6, [rsp + nb010nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6	;# rsq in xmm4 

	cvtsd2ss xmm6, xmm4	
	rcpss xmm6, xmm6
	cvtss2sd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [rsp + nb010nf_two]
	movapd xmm5, xmm4
	mulsd xmm4, xmm6	;# lu*rsq 
	subsd xmm0, xmm4	;# 2-lu*rsq 
	mulsd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [rsp + nb010nf_two]
	mulsd xmm5, xmm6	;# lu*rsq 
	subsd xmm0, xmm5	;# 2-lu*rsq 
	mulsd xmm0, xmm6	;# xmm0=rinvsq 
	movapd xmm4, xmm0
	
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulsd  xmm1, [rsp + nb010nf_c6]
	mulsd  xmm2, [rsp + nb010nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb010nf_Vvdwtot]
	movlpd [rsp + nb010nf_Vvdwtot], xmm5

.nb010nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb010nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb010nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb010nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb010nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 

        ;# finish if last 
        mov ecx, [rsp + nb010nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb010nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb010nf_n], esi
        jmp .nb010nf_outer
.nb010nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb010nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb010nf_end
        ;# non-zero, do one more workunit
        jmp   .nb010nf_threadloop
.nb010nf_end:
	mov eax, [rsp + nb010nf_nouter]
	mov ebx, [rsp + nb010nf_ninner]
	mov rcx, [rbp + nb010nf_outeriter]
	mov rdx, [rbp + nb010nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 248
	femms

	pop rbx
	pop	rbp
	ret
