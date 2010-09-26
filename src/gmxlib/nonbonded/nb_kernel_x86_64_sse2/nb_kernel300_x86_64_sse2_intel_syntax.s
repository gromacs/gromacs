;#
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
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as


.section .text


.globl nb_kernel300_x86_64_sse2
.globl _nb_kernel300_x86_64_sse2
nb_kernel300_x86_64_sse2:	
_nb_kernel300_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb300_fshift,           16
.equiv          nb300_gid,              24
.equiv          nb300_pos,              32
.equiv          nb300_faction,          40
.equiv          nb300_charge,           48
.equiv          nb300_p_facel,          56
.equiv          nb300_argkrf,           64
.equiv          nb300_argcrf,           72
.equiv          nb300_Vc,               80
.equiv          nb300_type,             88
.equiv          nb300_p_ntype,          96
.equiv          nb300_vdwparam,         104
.equiv          nb300_Vvdw,             112
.equiv          nb300_p_tabscale,       120
.equiv          nb300_VFtab,            128
.equiv          nb300_invsqrta,         136
.equiv          nb300_dvda,             144
.equiv          nb300_p_gbtabscale,     152
.equiv          nb300_GBtab,            160
.equiv          nb300_p_nthreads,       168
.equiv          nb300_count,            176
.equiv          nb300_mtx,              184
.equiv          nb300_outeriter,        192
.equiv          nb300_inneriter,        200
.equiv          nb300_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb300_ix,               0
.equiv          nb300_iy,               16
.equiv          nb300_iz,               32
.equiv          nb300_iq,               48
.equiv          nb300_dx,               64
.equiv          nb300_dy,               80
.equiv          nb300_dz,               96
.equiv          nb300_two,              112
.equiv          nb300_tsc,              128
.equiv          nb300_qq,               144
.equiv          nb300_fs,               160
.equiv          nb300_vctot,            176
.equiv          nb300_fix,              192
.equiv          nb300_fiy,              208
.equiv          nb300_fiz,              224
.equiv          nb300_half,             240
.equiv          nb300_three,            256
.equiv          nb300_is3,              272
.equiv          nb300_ii3,              276
.equiv          nb300_nri,              280
.equiv          nb300_iinr,             288
.equiv          nb300_jindex,           296
.equiv          nb300_jjnr,             304
.equiv          nb300_shift,            312
.equiv          nb300_shiftvec,         320
.equiv          nb300_facel,            328
.equiv          nb300_innerjjnr,        336
.equiv          nb300_innerk,           344
.equiv          nb300_n,                348
.equiv          nb300_nn1,              352
.equiv          nb300_nouter,           356
.equiv          nb300_ninner,           360

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 368		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb300_nouter], eax
	mov [rsp + nb300_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb300_nri], edi
	mov [rsp + nb300_iinr], rsi
	mov [rsp + nb300_jindex], rdx
	mov [rsp + nb300_jjnr], rcx
	mov [rsp + nb300_shift], r8
	mov [rsp + nb300_shiftvec], r9
	mov rsi, [rbp + nb300_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb300_facel], xmm0

	mov rax, [rbp + nb300_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb300_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb300_half], eax
	mov [rsp + nb300_half+4], ebx
	movsd xmm1, [rsp + nb300_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb300_half], xmm1
	movapd [rsp + nb300_two], xmm2
	movapd [rsp + nb300_three], xmm3

.nb300_threadloop:
        mov   rsi, [rbp + nb300_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb300_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb300_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb300_n], eax
        mov [rsp + nb300_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300_outerstart
        jmp .nb300_end

.nb300_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb300_nouter]
	mov [rsp + nb300_nouter], ebx

.nb300_outer:
	mov   rax, [rsp + nb300_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb300_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb300_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb300_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb300_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb300_facel]
	shufpd xmm3, xmm3, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb300_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb300_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb300_ix], xmm0
	movapd [rsp + nb300_iy], xmm1
	movapd [rsp + nb300_iz], xmm2

	mov   [rsp + nb300_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorpd xmm12, xmm12
	movapd xmm13, xmm12
	movapd xmm14, xmm12
	movapd xmm15, xmm12
	
	mov   rax, [rsp + nb300_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb300_pos]
	mov   rdi, [rbp + nb300_faction]	
	mov   rax, [rsp + nb300_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb300_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb300_ninner]
	mov   [rsp + nb300_ninner], ecx
	add   edx, 0
	mov   [rsp + nb300_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300_unroll_loop
	jmp   .nb300_checksingle
.nb300_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb300_innerjjnr]   ;# pointer to jjnr[k] 
	mov   r10d, [rdx]
	mov   r11d, [rdx + 4]
	add qword ptr [rsp + nb300_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb300_pos]		;# base of pos[] 

	lea   rax, [r10 + r10*2]     ;# replace jnr with j3 
	lea   rbx, [r11 + r11*2]	

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		

	mov    rdi, [rbp + nb300_faction]
	
	;# calc dr 
	subpd xmm4, [rsp + nb300_ix]
	subpd xmm5, [rsp + nb300_iy]
	subpd xmm6, [rsp + nb300_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	mov rsi, [rbp + nb300_charge]    ;# base of charge[] 

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	movlpd xmm3, [rsi + r10*8]
	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movhpd xmm3, [rsi + r11*8]


	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb300_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	mulpd  xmm3, [rsp + nb300_iq]

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb300_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb300_tsc]
	movapd [rsp + nb300_qq], xmm3	

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb300_VFtab]
	movd r8d, mm6
	psrlq mm6, 32
	movd r9d, mm6		;# indices in eax/ebx 

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 	
	movapd xmm3, [rsi + r9*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + r9*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	addpd  xmm7, xmm7	;# two*Heps2 
	movapd xmm3, [rsp + nb300_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm12, xmm5

	xorpd  xmm4, xmm4

	mulpd xmm3, [rsp + nb300_tsc]
	mulpd xmm3, xmm0
	subpd  xmm4, xmm3


	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4

	mov    rdi, [rbp + nb300_faction]
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]

	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11

	addpd xmm3, xmm9
	addpd xmm4, xmm10
	addpd xmm5, xmm11
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb300_innerk],  2
	jl    .nb300_checksingle
	jmp   .nb300_unroll_loop
.nb300_checksingle:
	mov   edx, [rsp + nb300_innerk]
	and   edx, 1
	jnz    .nb300_dosingle
	jmp    .nb300_updateouterdata
.nb300_dosingle:
	mov rsi, [rbp + nb300_charge]
	mov rdi, [rbp + nb300_pos]
	mov   rcx, [rsp + nb300_innerjjnr]
	mov   eax, [rcx]	

	mov rsi, [rbp + nb300_charge]    ;# base of charge[] 

	movsd xmm3, [rsi + rax*8]

	movapd xmm2, [rsp + nb300_iq]
	mulsd  xmm3, xmm2
	movapd [rsp + nb300_qq], xmm3	
	
	mov rsi, [rbp + nb300_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm4-xmm6
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]

	mov    rdi, [rbp + nb300_faction]
	
	;# calc dr 
	subsd xmm4, [rsp + nb300_ix]
	subsd xmm5, [rsp + nb300_iy]
	subsd xmm6, [rsp + nb300_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

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
	movapd xmm1, [rsp + nb300_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb300_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb300_tsc]

	cvttsd2si r8d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r8d
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb300_VFtab]

    movsd  xmm4, [rsi + r8*8]
    movsd  xmm5, [rsi + r8*8 + 8]
    movsd  xmm6, [rsi + r8*8 + 16]
    movsd  xmm7, [rsi + r8*8 + 24]
	;# coulomb table ready, in xmm4-xmm7  		

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	addsd  xmm7, xmm7	;# two*Heps2 
	movapd xmm3, [rsp + nb300_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm12, xmm5

	xorpd  xmm4, xmm4

	mulsd xmm3, [rsp + nb300_tsc]
	mulsd xmm3, xmm0
	subsd  xmm4, xmm3

	mov    rdi, [rbp + nb300_faction]
	mulsd  xmm9, xmm4
	mulsd  xmm10, xmm4
	mulsd  xmm11, xmm4

	;# now update f_i 
	addsd  xmm13, xmm9
	addsd  xmm14, xmm10
	addsd  xmm15, xmm11

	;# the fj's - start by accumulating forces from memory 
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
	
.nb300_updateouterdata:
	mov   ecx, [rsp + nb300_ii3]
	mov   rdi, [rbp + nb300_faction]
	mov   rsi, [rbp + nb300_fshift]
	mov   edx, [rsp + nb300_is3]

	;# accumulate i forces in xmm13, xmm14, xmm15
	movhlps xmm3, xmm13
	movhlps xmm4, xmm14
	movhlps xmm5, xmm15
	addsd  xmm13, xmm3
	addsd  xmm14, xmm4
	addsd  xmm15, xmm5 ;# sum is in low xmm13-xmm15

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [rsi + rdx*8]
	movsd  xmm4, [rsi + rdx*8 + 8]
	movsd  xmm5, [rsi + rdx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb300_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb300_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total coulomb energy and update it 
	movhlps xmm6, xmm12
	addsd  xmm12, xmm6	;# low xmm12 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb300_Vc]
	addsd xmm12, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb300_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb300_n], esi
        jmp .nb300_outer
.nb300_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb300_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300_end
        ;# non-zero, do one more workunit
        jmp   .nb300_threadloop
.nb300_end:
	mov eax, [rsp + nb300_nouter]
	mov ebx, [rsp + nb300_ninner]
	mov rcx, [rbp + nb300_outeriter]
	mov rdx, [rbp + nb300_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 368
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret



	

.globl nb_kernel300nf_x86_64_sse2
.globl _nb_kernel300nf_x86_64_sse2
nb_kernel300nf_x86_64_sse2:	
_nb_kernel300nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb300nf_fshift,         16
.equiv          nb300nf_gid,            24
.equiv          nb300nf_pos,            32
.equiv          nb300nf_faction,        40
.equiv          nb300nf_charge,         48
.equiv          nb300nf_p_facel,        56
.equiv          nb300nf_argkrf,         64
.equiv          nb300nf_argcrf,         72
.equiv          nb300nf_Vc,             80
.equiv          nb300nf_type,           88
.equiv          nb300nf_p_ntype,        96
.equiv          nb300nf_vdwparam,       104
.equiv          nb300nf_Vvdw,           112
.equiv          nb300nf_p_tabscale,     120
.equiv          nb300nf_VFtab,          128
.equiv          nb300nf_invsqrta,       136
.equiv          nb300nf_dvda,           144
.equiv          nb300nf_p_gbtabscale,   152
.equiv          nb300nf_GBtab,          160
.equiv          nb300nf_p_nthreads,     168
.equiv          nb300nf_count,          176
.equiv          nb300nf_mtx,            184
.equiv          nb300nf_outeriter,      192
.equiv          nb300nf_inneriter,      200
.equiv          nb300nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb300nf_ix,             0
.equiv          nb300nf_iy,             16
.equiv          nb300nf_iz,             32
.equiv          nb300nf_iq,             48
.equiv          nb300nf_tsc,            64
.equiv          nb300nf_qq,             80
.equiv          nb300nf_vctot,          96
.equiv          nb300nf_half,           112
.equiv          nb300nf_three,          128
.equiv          nb300nf_is3,            144
.equiv          nb300nf_ii3,            148
.equiv          nb300nf_nri,            152
.equiv          nb300nf_iinr,           160
.equiv          nb300nf_jindex,         168
.equiv          nb300nf_jjnr,           176
.equiv          nb300nf_shift,          184
.equiv          nb300nf_shiftvec,       192
.equiv          nb300nf_facel,          200
.equiv          nb300nf_innerjjnr,      208
.equiv          nb300nf_innerk,         216
.equiv          nb300nf_n,              220
.equiv          nb300nf_nn1,            224
.equiv          nb300nf_nouter,         228
.equiv          nb300nf_ninner,         232

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 240		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb300nf_nouter], eax
	mov [rsp + nb300nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb300nf_nri], edi
	mov [rsp + nb300nf_iinr], rsi
	mov [rsp + nb300nf_jindex], rdx
	mov [rsp + nb300nf_jjnr], rcx
	mov [rsp + nb300nf_shift], r8
	mov [rsp + nb300nf_shiftvec], r9
	mov rsi, [rbp + nb300nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb300nf_facel], xmm0

	mov rax, [rbp + nb300nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb300nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb300nf_half], eax
	mov [rsp + nb300nf_half+4], ebx
	movsd xmm1, [rsp + nb300nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb300nf_half], xmm1
	movapd [rsp + nb300nf_three], xmm3

.nb300nf_threadloop:
        mov   rsi, [rbp + nb300nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb300nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb300nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb300nf_n], eax
        mov [rsp + nb300nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300nf_outerstart
        jmp .nb300nf_end

.nb300nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb300nf_nouter]
	mov [rsp + nb300nf_nouter], ebx

.nb300nf_outer:
	mov   rax, [rsp + nb300nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb300nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb300nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb300nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb300nf_facel]
	shufpd xmm3, xmm3, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb300nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb300nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb300nf_ix], xmm0
	movapd [rsp + nb300nf_iy], xmm1
	movapd [rsp + nb300nf_iz], xmm2

	mov   [rsp + nb300nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb300nf_vctot], xmm4
	
	mov   rax, [rsp + nb300nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb300nf_pos]
	mov   rax, [rsp + nb300nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb300nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb300nf_ninner]
	mov   [rsp + nb300nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb300nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300nf_unroll_loop
	jmp   .nb300nf_checksingle
.nb300nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb300nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb300nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb300nf_charge]    ;# base of charge[] 

	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm2, [rsp + nb300nf_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb300nf_qq], xmm3	
	
	mov rsi, [rbp + nb300nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move nb300nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb300nf_ix]
	movapd xmm5, [rsp + nb300nf_iy]
	movapd xmm6, [rsp + nb300nf_iz]

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
	movapd xmm1, [rsp + nb300nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb300nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb300nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb300nf_VFtab]
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
	movapd xmm3, [rsp + nb300nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul  
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [rsp + nb300nf_vctot]
	movapd [rsp + nb300nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb300nf_innerk],  2
	jl    .nb300nf_checksingle
	jmp   .nb300nf_unroll_loop
.nb300nf_checksingle:
	mov   edx, [rsp + nb300nf_innerk]
	and   edx, 1
	jnz    .nb300nf_dosingle
	jmp    .nb300nf_updateouterdata
.nb300nf_dosingle:
	mov rsi, [rbp + nb300nf_charge]
	mov rdi, [rbp + nb300nf_pos]
	mov   rcx, [rsp + nb300nf_innerjjnr]
	mov   eax, [rcx]	
	xorpd  xmm6, xmm6
	movlpd xmm6, [rsi + rax*8]	;# xmm6(0) has the charge 	
	mulsd  xmm6, [rsp + nb300nf_iq]
	movapd [rsp + nb300nf_qq], xmm6
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]

	;# move nb300nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb300nf_ix]
	movapd xmm5, [rsp + nb300nf_iy]
	movapd xmm6, [rsp + nb300nf_iz]

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
	movapd xmm1, [rsp + nb300nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb300nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb300nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb300nf_tsc]
	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb300nf_VFtab]

	;# Coulomb 
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
	;# table ready in xmm4-xmm7 

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb300nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [rsp + nb300nf_vctot]
	movsd [rsp + nb300nf_vctot], xmm5 

.nb300nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [rsp + nb300nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb300nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb300nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb300nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb300nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb300nf_n], esi
        jmp .nb300nf_outer
.nb300nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb300nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300nf_end
        ;# non-zero, do one more workunit
        jmp   .nb300nf_threadloop
.nb300nf_end:
	mov eax, [rsp + nb300nf_nouter]
	mov ebx, [rsp + nb300nf_ninner]
	mov rcx, [rbp + nb300nf_outeriter]
	mov rdx, [rbp + nb300nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 240
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
