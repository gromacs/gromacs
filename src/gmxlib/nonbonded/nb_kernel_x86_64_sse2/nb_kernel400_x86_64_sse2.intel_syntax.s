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




.globl nb_kernel400_x86_64_sse2
.globl _nb_kernel400_x86_64_sse2
nb_kernel400_x86_64_sse2:	
_nb_kernel400_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb400_fshift,           16
.equiv          nb400_gid,              24
.equiv          nb400_pos,              32
.equiv          nb400_faction,          40
.equiv          nb400_charge,           48
.equiv          nb400_p_facel,          56
.equiv          nb400_argkrf,           64
.equiv          nb400_argcrf,           72
.equiv          nb400_Vc,               80
.equiv          nb400_type,             88
.equiv          nb400_p_ntype,          96
.equiv          nb400_vdwparam,         104
.equiv          nb400_Vvdw,             112
.equiv          nb400_p_tabscale,       120
.equiv          nb400_VFtab,            128
.equiv          nb400_invsqrta,         136
.equiv          nb400_dvda,             144
.equiv          nb400_p_gbtabscale,     152
.equiv          nb400_GBtab,            160
.equiv          nb400_p_nthreads,       168
.equiv          nb400_count,            176
.equiv          nb400_mtx,              184
.equiv          nb400_outeriter,        192
.equiv          nb400_inneriter,        200
.equiv          nb400_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb400_ix,               0
.equiv          nb400_iy,               16
.equiv          nb400_iz,               32
.equiv          nb400_iq,               48
.equiv          nb400_dx,               64
.equiv          nb400_dy,               80
.equiv          nb400_dz,               96
.equiv          nb400_two,              112
.equiv          nb400_gbtsc,            128
.equiv          nb400_qq,               144
.equiv          nb400_r,                160
.equiv          nb400_vctot,            176
.equiv          nb400_fix,              192
.equiv          nb400_fiy,              208
.equiv          nb400_fiz,              224
.equiv          nb400_half,             240
.equiv          nb400_three,            256
.equiv          nb400_isai,             272
.equiv          nb400_isaprod,          288
.equiv          nb400_dvdasum,          304
.equiv          nb400_gbscale,          320
.equiv          nb400_nri,              336
.equiv          nb400_iinr,             344
.equiv          nb400_jindex,           352
.equiv          nb400_jjnr,             360
.equiv          nb400_shift,            368
.equiv          nb400_shiftvec,         376
.equiv          nb400_facel,            384
.equiv          nb400_innerjjnr,        392
.equiv          nb400_is3,              400
.equiv          nb400_ii3,              404
.equiv          nb400_ii,               408
.equiv          nb400_innerk,           412
.equiv          nb400_n,                416
.equiv          nb400_nn1,              420
.equiv          nb400_nouter,           424
.equiv          nb400_ninner,           428
	push rbp
	mov  rbp, rsp
	push rbx

	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 440		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb400_nouter], eax
	mov [rsp + nb400_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb400_nri], edi
	mov [rsp + nb400_iinr], rsi
	mov [rsp + nb400_jindex], rdx
	mov [rsp + nb400_jjnr], rcx
	mov [rsp + nb400_shift], r8
	mov [rsp + nb400_shiftvec], r9
	mov rsi, [rbp + nb400_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb400_facel], xmm0

	mov rbx, [rbp + nb400_p_gbtabscale]
	movsd xmm4, [rbx]
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb400_gbtsc],  xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb400_half], eax
	mov [rsp + nb400_half+4], ebx
	movsd xmm1, [rsp + nb400_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb400_half], xmm1
	movapd [rsp + nb400_two], xmm2
	movapd [rsp + nb400_three], xmm3

.nb400_threadloop:
        mov   rsi, [rbp + nb400_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb400_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb400_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb400_n], eax
        mov [rsp + nb400_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400_outerstart
        jmp .nb400_end

.nb400_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb400_nouter]
	mov [rsp + nb400_nouter], ebx

.nb400_outer:
	mov   rax, [rsp + nb400_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb400_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb400_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb400_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 
	mov   [rsp + nb400_ii], ebx
	
	mov   rdx, [rbp + nb400_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb400_facel]
	shufpd xmm3, xmm3, 0

	mov   rdx, [rbp + nb400_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [rdx + rbx*8]
	shufpd xmm4, xmm4, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb400_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb400_iq], xmm3
	movapd [rsp + nb400_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb400_ix], xmm0
	movapd [rsp + nb400_iy], xmm1
	movapd [rsp + nb400_iz], xmm2

	mov   [rsp + nb400_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd xmm8, xmm4
	movapd xmm12, xmm4
	movapd xmm13, xmm4
	movapd xmm14, xmm4
	movapd xmm15, xmm4
	
	mov   rax, [rsp + nb400_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb400_pos]
	mov   rdi, [rbp + nb400_faction]	
	mov   rax, [rsp + nb400_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb400_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb400_ninner]
	mov   [rsp + nb400_ninner], ecx
	add   edx, 0
	mov   [rsp + nb400_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400_unroll_loop
	jmp   .nb400_checksingle
.nb400_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb400_innerjjnr]   ;# pointer to jjnr[k] 
	mov   r12d, [rdx]
	mov   r13d, [rdx + 4]
	add qword ptr [rsp + nb400_innerjjnr], 8	;# advance pointer (unrolled 2) 
	
	mov rsi, [rbp + nb400_pos]		;# base of pos[] 

	lea   r8, [r12 + r12*2]     ;# j3 
	lea   r9, [r13 + r13*2]	

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + r8*8]
	movlpd xmm5, [rsi + r8*8 + 8]
	movlpd xmm6, [rsi + r8*8 + 16]
	movhpd xmm4, [rsi + r9*8]
	movhpd xmm5, [rsi + r9*8 + 8]
	movhpd xmm6, [rsi + r9*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb400_ix]
	subpd xmm5, [rsp + nb400_iy]
	subpd xmm6, [rsp + nb400_iz]


	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	mov rsi, [rbp + nb400_invsqrta]
	movlpd xmm3, [rsi + r12*8]

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movhpd xmm3, [rsi + r13*8]
	mulpd  xmm3, [rsp + nb400_isai]
	movapd [rsp + nb400_isaprod], xmm3	
    movapd xmm6, xmm3
	mulpd xmm3, [rsp + nb400_gbtsc]
	movapd [rsp + nb400_gbscale], xmm3

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb400_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	mov rsi, [rbp + nb400_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + r12*8]

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb400_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
    
    mulpd  xmm6, [rsp + nb400_iq]
	movhpd xmm3, [rsi + r13*8]
	mulpd  xmm3, xmm6
	movapd [rsp + nb400_qq], xmm3	


	movapd [rsp + nb400_r], xmm4
	mulpd xmm4, [rsp + nb400_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	
	pslld mm6, 2		;# idx *= 4 

	mov  rsi, [rbp + nb400_GBtab]
	movd r10d, mm6
	psrlq mm6, 32
	movd r11d, mm6		;# indices in r10/r11

	movapd xmm4, [rsi + r10*8]	;# Y1 F1 	
	movapd xmm3, [rsi + r11*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + r10*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + r11*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		

	mulpd  xmm7, xmm1	;# xmm7=Heps
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm1	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	addpd  xmm7, xmm7	;# two*Heps2 
	movapd xmm3, [rsp + nb400_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 

	mov rsi, [rbp + nb400_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulpd xmm3, [rsp + nb400_gbscale]
	movapd xmm6, xmm3
	mulpd  xmm6, [rsp + nb400_r]
	addpd  xmm6, xmm5

    ;# update vctot
	addpd  xmm12, xmm5

	;# xmm6=(vcoul+fijC*r)
	subpd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum
	addpd  xmm8, xmm7

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	addsd  xmm6, [rsi + r12*8]
	addsd  xmm7, [rsi + r13*8]
	movsd  [rsi + r12*8], xmm6
	movsd  [rsi + r13*8], xmm7
	
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb400_faction]
	movlpd xmm5, [rdi + r8*8]
	movlpd xmm6, [rdi + r8*8 + 8]
	movlpd xmm7, [rdi + r8*8 + 16]
	movhpd xmm5, [rdi + r9*8]
	movhpd xmm6, [rdi + r9*8 + 8]
	movhpd xmm7, [rdi + r9*8 + 16]

	xorpd  xmm4, xmm4

	mulpd xmm3, xmm0
	subpd  xmm4, xmm3

	mov    rdi, [rbp + nb400_faction]
	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4
    
	addpd xmm5, xmm9
	addpd xmm6, xmm10
	addpd xmm7, xmm11

	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11

	movlpd [rdi + r8*8], xmm5
	movlpd [rdi + r8*8 + 8], xmm6
	movlpd [rdi + r8*8 + 16], xmm7
	movhpd [rdi + r9*8], xmm5
	movhpd [rdi + r9*8 + 8], xmm6
	movhpd [rdi + r9*8 + 16], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb400_innerk],  2
	jl    .nb400_checksingle
	jmp   .nb400_unroll_loop
.nb400_checksingle:
	mov   edx, [rsp + nb400_innerk]
	and   edx, 1
	jnz    .nb400_dosingle
	jmp    .nb400_updateouterdata
.nb400_dosingle:
	mov rsi, [rbp + nb400_charge]
	mov rdx, [rbp + nb400_invsqrta]
	mov rdi, [rbp + nb400_pos]
	mov   rcx, [rsp + nb400_innerjjnr]
	mov   eax, [rcx]	

	;# load isaj
	mov rsi, [rbp + nb400_invsqrta]
	movsd xmm2, [rsi + rax*8]
	mulsd  xmm2, [rsp + nb400_isai]
	movapd [rsp + nb400_isaprod], xmm2	
	movapd xmm1, xmm2
	mulsd xmm1, [rsp + nb400_gbtsc]
	movapd [rsp + nb400_gbscale], xmm1

    mulsd xmm2, [rsp + nb400_iq]
	mov rsi, [rbp + nb400_charge]    ;# base of charge[] 
	movsd xmm3, [rsi + rax*8]
	mulsd  xmm3, xmm2
	movapd [rsp + nb400_qq], xmm3	

	mov rsi, [rbp + nb400_pos]		;# base of pos[] 

	lea   r8, [rax + rax*2]     ;# j3 

	;# move coordinate to xmm4-xmm6
	movsd xmm4, [rsi + r8*8]
	movsd xmm5, [rsi + r8*8 + 8]
	movsd xmm6, [rsi + r8*8 + 16]

	mov    rdi, [rbp + nb400_faction]
	
	;# calc dr 
	subsd xmm4, [rsp + nb400_ix]
	subsd xmm5, [rsp + nb400_iy]
	subsd xmm6, [rsp + nb400_iz]

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
	movapd xmm1, [rsp + nb400_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb400_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
    
	movapd [rsp + nb400_r], xmm4
	mulsd xmm4, [rsp + nb400_gbscale]

	cvttsd2si r10d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r10d
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	
	shl r10d, 2		;# idx *= 4 

	mov  rsi, [rbp + nb400_GBtab]

	movapd xmm4, [rsi + r10*8]	;# Y1 F1 	
	movhlps xmm5, xmm4
	movapd xmm6, [rsi + r10*8 + 16]	;# G1 H1 	
    movhlps xmm7, xmm6
	;# coulomb table ready, in xmm4-xmm7  		

	mulsd  xmm7, xmm1	;# xmm7=Heps
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm1	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	addsd  xmm7, xmm7	;# two*Heps2 
	movapd xmm3, [rsp + nb400_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 

	mov rsi, [rbp + nb400_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulsd xmm3, [rsp + nb400_gbscale]
	movapd xmm6, xmm3
	mulsd  xmm6, [rsp + nb400_r]
	addsd  xmm6, xmm5

    ;# update vctot
	addsd  xmm12, xmm5

	;# xmm6=(vcoul+fijC*r)
	subsd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum
	addsd  xmm8, xmm7

	;# update j atoms dvdaj
	addsd  xmm6, [rsi + rax*8]
	movsd  [rsi + rax*8], xmm6
	
	xorpd  xmm4, xmm4

	mulsd xmm3, xmm0
	subsd  xmm4, xmm3

	mov    rdi, [rbp + nb400_faction]
	mulsd  xmm9, xmm4
	mulsd  xmm10, xmm4
	mulsd  xmm11, xmm4
    
	;# now update f_i 
	addsd  xmm13, xmm9
	addsd  xmm14, xmm10
	addsd  xmm15, xmm11

	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb400_faction]
	addsd xmm9,  [rdi + r8*8]
	addsd xmm10, [rdi + r8*8 + 8]
	addsd xmm11, [rdi + r8*8 + 16]
	movsd [rdi + r8*8], xmm9
	movsd [rdi + r8*8 + 8], xmm10
	movsd [rdi + r8*8 + 16], xmm11
	
.nb400_updateouterdata:
	mov   ecx, [rsp + nb400_ii3]
	mov   rdi, [rbp + nb400_faction]
	mov   rsi, [rbp + nb400_fshift]
	mov   edx, [rsp + nb400_is3]

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
	mov esi, [rsp + nb400_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb400_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total coulomb energy and update it 
	movhlps xmm6, xmm12
	addsd  xmm12, xmm6	;# low xmm12 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb400_Vc]
	addsd xmm12, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm12
	
	;# accumulate dVda and update it 
	movhlps xmm6, xmm8
	addsd  xmm8, xmm6	;# low xmm8 has the sum now 
	
	mov edx, [rsp + nb400_ii]
	mov rax, [rbp + nb400_dvda]
	addsd xmm8, [rax + rdx*8]
	movsd [rax + rdx*8], xmm8
	
        ;# finish if last 
        mov ecx, [rsp + nb400_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb400_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb400_n], esi
        jmp .nb400_outer
.nb400_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb400_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb400_end
        ;# non-zero, do one more workunit
        jmp   .nb400_threadloop
.nb400_end:
	mov eax, [rsp + nb400_nouter]
	mov ebx, [rsp + nb400_ninner]
	mov rcx, [rbp + nb400_outeriter]
	mov rdx, [rbp + nb400_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 440
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret







.globl nb_kernel400nf_x86_64_sse2
.globl _nb_kernel400nf_x86_64_sse2
nb_kernel400nf_x86_64_sse2:	
_nb_kernel400nf_x86_64_sse2:	
.equiv          nb400nf_fshift,         16
.equiv          nb400nf_gid,            24
.equiv          nb400nf_pos,            32
.equiv          nb400nf_faction,        40
.equiv          nb400nf_charge,         48
.equiv          nb400nf_p_facel,        56
.equiv          nb400nf_argkrf,         64
.equiv          nb400nf_argcrf,         72
.equiv          nb400nf_Vc,             80
.equiv          nb400nf_type,           88
.equiv          nb400nf_p_ntype,        96
.equiv          nb400nf_vdwparam,       104
.equiv          nb400nf_Vvdw,           112
.equiv          nb400nf_p_tabscale,     120
.equiv          nb400nf_VFtab,          128
.equiv          nb400nf_invsqrta,       136
.equiv          nb400nf_dvda,           144
.equiv          nb400nf_p_gbtabscale,   152
.equiv          nb400nf_GBtab,          160
.equiv          nb400nf_p_nthreads,     168
.equiv          nb400nf_count,          176
.equiv          nb400nf_mtx,            184
.equiv          nb400nf_outeriter,      192
.equiv          nb400nf_inneriter,      200
.equiv          nb400nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb400nf_ix,             0
.equiv          nb400nf_iy,             16
.equiv          nb400nf_iz,             32
.equiv          nb400nf_iq,             48
.equiv          nb400nf_gbtsc,          64
.equiv          nb400nf_qq,             80
.equiv          nb400nf_vctot,          96
.equiv          nb400nf_half,           112
.equiv          nb400nf_three,          128
.equiv          nb400nf_isai,           144
.equiv          nb400nf_isaprod,        160
.equiv          nb400nf_gbscale,        176
.equiv          nb400nf_nri,            192
.equiv          nb400nf_iinr,           200
.equiv          nb400nf_jindex,         208
.equiv          nb400nf_jjnr,           216
.equiv          nb400nf_shift,          224
.equiv          nb400nf_shiftvec,       232
.equiv          nb400nf_facel,          240
.equiv          nb400nf_innerjjnr,      248
.equiv          nb400nf_is3,            256
.equiv          nb400nf_ii3,            260
.equiv          nb400nf_innerk,         264
.equiv          nb400nf_n,              268
.equiv          nb400nf_nn1,            272
.equiv          nb400nf_nouter,         276
.equiv          nb400nf_ninner,         280
	push rbp
	mov  rbp, rsp
	push rbx

	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 296		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb400nf_nouter], eax
	mov [rsp + nb400nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb400nf_nri], edi
	mov [rsp + nb400nf_iinr], rsi
	mov [rsp + nb400nf_jindex], rdx
	mov [rsp + nb400nf_jjnr], rcx
	mov [rsp + nb400nf_shift], r8
	mov [rsp + nb400nf_shiftvec], r9
	mov rsi, [rbp + nb400nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb400nf_facel], xmm0

	mov rbx, [rbp + nb400nf_p_gbtabscale]
	movsd xmm4, [rbx]
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb400nf_gbtsc],  xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb400nf_half], eax
	mov [rsp + nb400nf_half+4], ebx
	movsd xmm1, [rsp + nb400nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb400nf_half], xmm1
	movapd [rsp + nb400nf_three], xmm3

.nb400nf_threadloop:
        mov   rsi, [rbp + nb400nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb400nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb400nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb400nf_n], eax
        mov [rsp + nb400nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400nf_outerstart
        jmp .nb400nf_end

.nb400nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb400nf_nouter]
	mov [rsp + nb400nf_nouter], ebx

.nb400nf_outer:
	mov   rax, [rsp + nb400nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb400nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb400nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb400nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb400nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb400nf_facel]
	shufpd xmm3, xmm3, 0

	mov   rdx, [rbp + nb400nf_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [rdx + rbx*8]
	shufpd xmm4, xmm4, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb400nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb400nf_iq], xmm3
	movapd [rsp + nb400nf_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb400nf_ix], xmm0
	movapd [rsp + nb400nf_iy], xmm1
	movapd [rsp + nb400nf_iz], xmm2

	mov   [rsp + nb400nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb400nf_vctot], xmm4
	
	mov   rax, [rsp + nb400nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb400nf_pos]
	mov   rdi, [rbp + nb400nf_faction]	
	mov   rax, [rsp + nb400nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb400nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb400nf_ninner]
	mov   [rsp + nb400nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb400nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400nf_unroll_loop
	jmp   .nb400nf_checksingle
.nb400nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb400nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb400nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isa2
	mov rsi, [rbp + nb400nf_invsqrta]
	movlpd xmm2, [rsi + rax*8]
	movhpd xmm2, [rsi + rbx*8]
	mulpd  xmm2, [rsp + nb400nf_isai]
	movapd [rsp + nb400nf_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [rsp + nb400nf_gbtsc]
	movapd [rsp + nb400nf_gbscale], xmm1
	
	mov rsi, [rbp + nb400nf_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	mulpd xmm2, [rsp + nb400nf_iq]
    mulpd xmm3, xmm2
	movapd [rsp + nb400nf_qq], xmm3	
	
	mov rsi, [rbp + nb400nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	mov    rdi, [rbp + nb400nf_faction]
	
	;# move nb400nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb400nf_ix]
	movapd xmm5, [rsp + nb400nf_iy]
	movapd xmm6, [rsp + nb400nf_iz]

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
	movapd xmm1, [rsp + nb400nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb400nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb400nf_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  rsi, [rbp + nb400nf_GBtab]
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
	movapd xmm3, [rsp + nb400nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	addpd  xmm5, [rsp + nb400nf_vctot]
	movapd [rsp + nb400nf_vctot], xmm5  
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb400nf_innerk],  2
	jl    .nb400nf_checksingle
	jmp   .nb400nf_unroll_loop
.nb400nf_checksingle:
	mov   edx, [rsp + nb400nf_innerk]
	and   edx, 1
	jnz    .nb400nf_dosingle
	jmp    .nb400nf_updateouterdata
.nb400nf_dosingle:
	mov rsi, [rbp + nb400nf_charge]
	mov rdx, [rbp + nb400nf_invsqrta]
	mov rdi, [rbp + nb400nf_pos]
	mov   rcx, [rsp + nb400nf_innerjjnr]
	mov   eax, [rcx]	
	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [rdx + rax*8]
	movlpd xmm6, [rsi + rax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [rsp + nb400nf_isai]
	movapd [rsp + nb400nf_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [rsp + nb400nf_gbtsc]
	movapd [rsp + nb400nf_gbscale], xmm1
	
	mulsd  xmm7, [rsp + nb400nf_iq]
	mulsd  xmm6, xmm7
	movapd [rsp + nb400nf_qq], xmm6
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]

	;# move nb400nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb400nf_ix]
	movapd xmm5, [rsp + nb400nf_iy]
	movapd xmm6, [rsp + nb400nf_iz]

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
	movapd xmm1, [rsp + nb400nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb400nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb400nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb400nf_gbscale]
	
	movd mm0, eax	

	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  rsi, [rbp + nb400nf_GBtab]

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
	movapd xmm3, [rsp + nb400nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	addsd  xmm5, [rsp + nb400nf_vctot]
	movsd [rsp + nb400nf_vctot], xmm5
	
.nb400nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb400nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb400nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb400nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb400nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb400nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb400nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb400nf_n], esi
        jmp .nb400nf_outer
.nb400nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb400nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb400nf_end
        ;# non-zero, do one more workunit
        jmp   .nb400nf_threadloop
.nb400nf_end:

	mov eax, [rsp + nb400nf_nouter]
	mov ebx, [rsp + nb400nf_ninner]
	mov rcx, [rbp + nb400nf_outeriter]
	mov rdx, [rbp + nb400nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 296
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret




