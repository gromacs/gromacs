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





.globl nb_kernel310_x86_64_sse2
.globl _nb_kernel310_x86_64_sse2
nb_kernel310_x86_64_sse2:	
_nb_kernel310_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb310_fshift,           16
.equiv          nb310_gid,              24
.equiv          nb310_pos,              32
.equiv          nb310_faction,          40
.equiv          nb310_charge,           48
.equiv          nb310_p_facel,          56
.equiv          nb310_argkrf,           64
.equiv          nb310_argcrf,           72
.equiv          nb310_Vc,               80
.equiv          nb310_type,             88
.equiv          nb310_p_ntype,          96
.equiv          nb310_vdwparam,         104
.equiv          nb310_Vvdw,             112
.equiv          nb310_p_tabscale,       120
.equiv          nb310_VFtab,            128
.equiv          nb310_invsqrta,         136
.equiv          nb310_dvda,             144
.equiv          nb310_p_gbtabscale,     152
.equiv          nb310_GBtab,            160
.equiv          nb310_p_nthreads,       168
.equiv          nb310_count,            176
.equiv          nb310_mtx,              184
.equiv          nb310_outeriter,        192
.equiv          nb310_inneriter,        200
.equiv          nb310_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb310_ix,               0
.equiv          nb310_iy,               16
.equiv          nb310_iz,               32
.equiv          nb310_iq,               48
.equiv          nb310_dx,               64
.equiv          nb310_dy,               80
.equiv          nb310_dz,               96
.equiv          nb310_two,              112
.equiv          nb310_six,              128
.equiv          nb310_twelve,           144
.equiv          nb310_tsc,              160
.equiv          nb310_qq,               176
.equiv          nb310_c6,               192
.equiv          nb310_c12,              208
.equiv          nb310_fscal,            224
.equiv          nb310_vctot,            240
.equiv          nb310_Vvdwtot,          256
.equiv          nb310_fix,              272
.equiv          nb310_fiy,              288
.equiv          nb310_fiz,              304
.equiv          nb310_half,             320
.equiv          nb310_three,            336
.equiv          nb310_nri,              352
.equiv          nb310_iinr,             360
.equiv          nb310_jindex,           368
.equiv          nb310_jjnr,             376
.equiv          nb310_shift,            384
.equiv          nb310_shiftvec,         392
.equiv          nb310_facel,            400
.equiv          nb310_innerjjnr,        408
.equiv          nb310_is3,              416
.equiv          nb310_ii3,              420
.equiv          nb310_ntia,             424
.equiv          nb310_innerk,           428
.equiv          nb310_n,                432
.equiv          nb310_nn1,              436
.equiv          nb310_ntype,            440
.equiv          nb310_nouter,           444
.equiv          nb310_ninner,           448
	push rbp
	mov  rbp, rsp
	push rbx	
	femms
	sub rsp, 472		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb310_nouter], eax
	mov [rsp + nb310_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb310_nri], edi
	mov [rsp + nb310_iinr], rsi
	mov [rsp + nb310_jindex], rdx
	mov [rsp + nb310_jjnr], rcx
	mov [rsp + nb310_shift], r8
	mov [rsp + nb310_shiftvec], r9
	mov rdi, [rbp + nb310_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb310_ntype], edi
	mov rsi, [rbp + nb310_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb310_facel], xmm0

	mov rax, [rbp + nb310_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb310_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb310_half], eax
	mov [rsp + nb310_half + 4], ebx
	movsd xmm1, [rsp + nb310_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb310_half], xmm1
	movapd [rsp + nb310_two], xmm2
	movapd [rsp + nb310_three], xmm3
	movapd [rsp + nb310_six], xmm4
	movapd [rsp + nb310_twelve], xmm5

.nb310_threadloop:
        mov   rsi, [rbp + nb310_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb310_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb310_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb310_n], eax
        mov [rsp + nb310_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310_outerstart
        jmp .nb310_end

.nb310_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb310_nouter]
	mov [rsp + nb310_nouter], ebx

.nb310_outer:
	mov   rax, [rsp + nb310_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb310_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb310_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb310_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb310_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb310_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb310_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb310_ntype]
    	shl   edx, 1
    	mov   [rsp + nb310_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb310_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb310_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb310_ix], xmm0
	movapd [rsp + nb310_iy], xmm1
	movapd [rsp + nb310_iz], xmm2

	mov   [rsp + nb310_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb310_vctot], xmm4
	movapd [rsp + nb310_Vvdwtot], xmm4
	movapd [rsp + nb310_fix], xmm4
	movapd [rsp + nb310_fiy], xmm4
	movapd [rsp + nb310_fiz], xmm4
	
	mov   rax, [rsp + nb310_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb310_pos]
	mov   rdi, [rbp + nb310_faction]	
	mov   rax, [rsp + nb310_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb310_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb310_ninner]
	mov   [rsp + nb310_ninner], ecx
	add   edx, 0
	mov   [rsp + nb310_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310_unroll_loop
	jmp   .nb310_checksingle
.nb310_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb310_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb310_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm2, [rsp + nb310_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb310_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb310_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb310_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb310_ntia]
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
	movapd [rsp + nb310_c6], xmm4
	movapd [rsp + nb310_c12], xmm6
	
	mov rsi, [rbp + nb310_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb310_ix]
	movapd xmm5, [rsp + nb310_iy]
	movapd xmm6, [rsp + nb310_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [rsp + nb310_dx], xmm4
	movapd [rsp + nb310_dy], xmm5
	movapd [rsp + nb310_dz], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb310_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb310_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb310_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  rsi, [rbp + nb310_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [rsp + nb310_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb310_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	
	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [rsp + nb310_vctot]

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4

	movapd [rsp + nb310_vctot], xmm5 

	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [rsp + nb310_c6]
	mulpd  xmm4, [rsp + nb310_c12]
	movapd xmm7, [rsp + nb310_Vvdwtot]
	addpd  xmm7, xmm4
	mulpd  xmm4, [rsp + nb310_twelve]
	subpd  xmm7, xmm6
	mulpd  xmm3, [rsp + nb310_tsc]
	mulpd  xmm6, [rsp + nb310_six]
	movapd [rsp + nb310_Vvdwtot], xmm7
	subpd  xmm4, xmm6
	mulpd  xmm4, xmm0
	subpd  xmm4, xmm3
	mulpd  xmm4, xmm0

	movapd xmm0, [rsp + nb310_dx]
	movapd xmm1, [rsp + nb310_dy]
	movapd xmm2, [rsp + nb310_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    rdi, [rbp + nb310_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [rsp + nb310_fix]
	movapd xmm4, [rsp + nb310_fiy]
	movapd xmm5, [rsp + nb310_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [rsp + nb310_fix], xmm3
	movapd [rsp + nb310_fiy], xmm4
	movapd [rsp + nb310_fiz], xmm5
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
	sub dword ptr [rsp + nb310_innerk],  2
	jl    .nb310_checksingle
	jmp   .nb310_unroll_loop
.nb310_checksingle:
	mov   edx, [rsp + nb310_innerk]
	and   edx, 1
	jnz    .nb310_dosingle
	jmp    .nb310_updateouterdata
.nb310_dosingle:
	mov rsi, [rbp + nb310_charge]
	mov rdi, [rbp + nb310_pos]
	mov   rcx, [rsp + nb310_innerjjnr]
	mov   eax, [rcx]
	
	xorpd  xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm2, [rsp + nb310_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb310_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov rsi, [rbp + nb310_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb310_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb310_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [rsp + nb310_c6], xmm4
	movapd [rsp + nb310_c12], xmm6
	
	mov rsi, [rbp + nb310_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb310_ix]
	movapd xmm5, [rsp + nb310_iy]
	movapd xmm6, [rsp + nb310_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [rsp + nb310_dx], xmm4
	movapd [rsp + nb310_dy], xmm5
	movapd [rsp + nb310_dz], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb310_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb310_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb310_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb310_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb310_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb310_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	
	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [rsp + nb310_vctot]

	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	movlpd [rsp + nb310_vctot], xmm5 

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [rsp + nb310_c6]
	mulsd  xmm4, [rsp + nb310_c12]
	movapd xmm7, [rsp + nb310_Vvdwtot]
	addsd  xmm7, xmm4
	mulsd  xmm4, [rsp + nb310_twelve]
	subsd  xmm7, xmm6
	mulsd  xmm3, [rsp + nb310_tsc]
	mulsd  xmm6, [rsp + nb310_six]
	movlpd [rsp + nb310_Vvdwtot], xmm7
	subsd  xmm4, xmm6
	mulsd  xmm4, xmm0
	subsd  xmm4, xmm3
	mulsd  xmm4, xmm0

	movapd xmm0, [rsp + nb310_dx]
	movapd xmm1, [rsp + nb310_dy]
	movapd xmm2, [rsp + nb310_dz]

	movd eax, mm0	

	mov    rdi, [rbp + nb310_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [rsp + nb310_fix]
	movapd xmm4, [rsp + nb310_fiy]
	movapd xmm5, [rsp + nb310_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [rsp + nb310_fix], xmm3
	movlpd [rsp + nb310_fiy], xmm4
	movlpd [rsp + nb310_fiz], xmm5
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	
.nb310_updateouterdata:
	mov   ecx, [rsp + nb310_ii3]
	mov   rdi, [rbp + nb310_faction]
	mov   rsi, [rbp + nb310_fshift]
	mov   edx, [rsp + nb310_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb310_fix]
	movapd xmm1, [rsp + nb310_fiy]
	movapd xmm2, [rsp + nb310_fiz]

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
	mov esi, [rsp + nb310_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb310_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb310_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb310_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb310_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb310_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb310_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb310_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb310_n], esi
        jmp .nb310_outer
.nb310_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb310_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb310_end
        ;# non-zero, do one more workunit
        jmp   .nb310_threadloop
.nb310_end:
	mov eax, [rsp + nb310_nouter]
	mov ebx, [rsp + nb310_ninner]
	mov rcx, [rbp + nb310_outeriter]
	mov rdx, [rbp + nb310_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 472
	femms

	pop rbx
	pop	rbp
	ret





.globl nb_kernel310nf_x86_64_sse2
.globl _nb_kernel310nf_x86_64_sse2
nb_kernel310nf_x86_64_sse2:	
_nb_kernel310nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb310nf_fshift,         16
.equiv          nb310nf_gid,            24
.equiv          nb310nf_pos,            32
.equiv          nb310nf_faction,        40
.equiv          nb310nf_charge,         48
.equiv          nb310nf_p_facel,        56
.equiv          nb310nf_argkrf,         64
.equiv          nb310nf_argcrf,         72
.equiv          nb310nf_Vc,             80
.equiv          nb310nf_type,           88
.equiv          nb310nf_p_ntype,        96
.equiv          nb310nf_vdwparam,       104
.equiv          nb310nf_Vvdw,           112
.equiv          nb310nf_p_tabscale,     120
.equiv          nb310nf_VFtab,          128
.equiv          nb310nf_invsqrta,       136
.equiv          nb310nf_dvda,           144
.equiv          nb310nf_p_gbtabscale,   152
.equiv          nb310nf_GBtab,          160
.equiv          nb310nf_p_nthreads,     168
.equiv          nb310nf_count,          176
.equiv          nb310nf_mtx,            184
.equiv          nb310nf_outeriter,      192
.equiv          nb310nf_inneriter,      200
.equiv          nb310nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb310nf_ix,             0
.equiv          nb310nf_iy,             16
.equiv          nb310nf_iz,             32
.equiv          nb310nf_iq,             48
.equiv          nb310nf_tsc,            64
.equiv          nb310nf_qq,             80
.equiv          nb310nf_c6,             96
.equiv          nb310nf_c12,            112
.equiv          nb310nf_vctot,          128
.equiv          nb310nf_Vvdwtot,        144
.equiv          nb310nf_half,           160
.equiv          nb310nf_three,          176
.equiv          nb310nf_nri,            192
.equiv          nb310nf_iinr,           200
.equiv          nb310nf_jindex,         208
.equiv          nb310nf_jjnr,           216
.equiv          nb310nf_shift,          224
.equiv          nb310nf_shiftvec,       232
.equiv          nb310nf_facel,          240
.equiv          nb310nf_innerjjnr,      248
.equiv          nb310nf_is3,            256
.equiv          nb310nf_ii3,            260
.equiv          nb310nf_ntia,           264
.equiv          nb310nf_innerk,         268
.equiv          nb310nf_n,              272
.equiv          nb310nf_nn1,            276
.equiv          nb310nf_ntype,          280
.equiv          nb310nf_nouter,         284
.equiv          nb310nf_ninner,         288
	push rbp
	mov  rbp, rsp
	push rbx
	femms
	sub rsp, 312		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb310nf_nouter], eax
	mov [rsp + nb310nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb310nf_nri], edi
	mov [rsp + nb310nf_iinr], rsi
	mov [rsp + nb310nf_jindex], rdx
	mov [rsp + nb310nf_jjnr], rcx
	mov [rsp + nb310nf_shift], r8
	mov [rsp + nb310nf_shiftvec], r9
	mov rdi, [rbp + nb310nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb310nf_ntype], edi
	mov rsi, [rbp + nb310nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb310nf_facel], xmm0

	mov rax, [rbp + nb310nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb310nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb310nf_half], eax
	mov [rsp + nb310nf_half + 4], ebx
	movsd xmm1, [rsp + nb310nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb310nf_half], xmm1
	movapd [rsp + nb310nf_three], xmm3

.nb310nf_threadloop:
        mov   rsi, [rbp + nb310nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb310nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb310nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb310nf_n], eax
        mov [rsp + nb310nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310nf_outerstart
        jmp .nb310nf_end

.nb310nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb310nf_nouter]
	mov [rsp + nb310nf_nouter], ebx

.nb310nf_outer:
	mov   rax, [rsp + nb310nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb310nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb310nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb310nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb310nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb310nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb310nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb310nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb310nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb310nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb310nf_ix], xmm0
	movapd [rsp + nb310nf_iy], xmm1
	movapd [rsp + nb310nf_iz], xmm2

	mov   [rsp + nb310nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb310nf_vctot], xmm4
	movapd [rsp + nb310nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb310nf_jindex]
	mov   ecx, [rax+rsi*4]	     		;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     	;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb310nf_pos]
	mov   rax, [rsp + nb310nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb310nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb310nf_ninner]
	mov   [rsp + nb310nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb310nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310nf_unroll_loop
	jmp   .nb310nf_checksingle
.nb310nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb310nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb310nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb310nf_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm2, [rsp + nb310nf_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb310nf_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb310nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb310nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb310nf_ntia]
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
	movapd [rsp + nb310nf_c6], xmm4
	movapd [rsp + nb310nf_c12], xmm6
	
	mov rsi, [rbp + nb310nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb310nf_ix]
	movapd xmm5, [rsp + nb310nf_iy]
	movapd xmm6, [rsp + nb310nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb310nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb310nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb310nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb310nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb310nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	
	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [rsp + nb310nf_vctot]

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4

	movapd [rsp + nb310nf_vctot], xmm5 

	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [rsp + nb310nf_c6]
	mulpd  xmm4, [rsp + nb310nf_c12]
	movapd xmm7, [rsp + nb310nf_Vvdwtot]
	addpd  xmm7, xmm4
	subpd  xmm7, xmm6
	movapd [rsp + nb310nf_Vvdwtot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb310nf_innerk],  2
	jl    .nb310nf_checksingle
	jmp   .nb310nf_unroll_loop
.nb310nf_checksingle:
	mov   edx, [rsp + nb310nf_innerk]
	and   edx, 1
	jnz    .nb310nf_dosingle
	jmp    .nb310nf_updateouterdata
.nb310nf_dosingle:
	mov rsi, [rbp + nb310nf_charge]
	mov rdi, [rbp + nb310nf_pos]
	mov   rcx, [rsp + nb310nf_innerjjnr]
	mov   eax, [rcx]
	
	xorpd  xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm2, [rsp + nb310nf_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb310nf_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov rsi, [rbp + nb310nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb310nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb310nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [rsp + nb310nf_c6], xmm4
	movapd [rsp + nb310nf_c12], xmm6
	
	mov rsi, [rbp + nb310nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb310nf_ix]
	movapd xmm5, [rsp + nb310nf_iy]
	movapd xmm6, [rsp + nb310nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb310nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb310nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb310nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb310nf_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb310nf_VFtab]

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb310nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	
	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 

	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [rsp + nb310nf_vctot]

	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	movlpd [rsp + nb310nf_vctot], xmm5 

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [rsp + nb310nf_c6]
	mulsd  xmm4, [rsp + nb310nf_c12]
	movapd xmm7, [rsp + nb310nf_Vvdwtot]
	addsd  xmm7, xmm4
	subsd  xmm7, xmm6
	movlpd [rsp + nb310nf_Vvdwtot], xmm7
	
.nb310nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [rsp + nb310nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb310nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb310nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb310nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb310nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb310nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb310nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb310nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb310nf_n], esi
        jmp .nb310nf_outer
.nb310nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb310nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb310nf_end
        ;# non-zero, do one more workunit
        jmp   .nb310nf_threadloop
.nb310nf_end:
	mov eax, [rsp + nb310nf_nouter]
	mov ebx, [rsp + nb310nf_ninner]
	mov rcx, [rbp + nb310nf_outeriter]
	mov rdx, [rbp + nb310nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 312
	femms

	pop rbx
	pop	rbp
	ret
