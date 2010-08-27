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


	

.globl nb_kernel201_x86_64_sse2
.globl _nb_kernel201_x86_64_sse2
nb_kernel201_x86_64_sse2:	
_nb_kernel201_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb201_fshift,           16
.equiv          nb201_gid,              24
.equiv          nb201_pos,              32
.equiv          nb201_faction,          40
.equiv          nb201_charge,           48
.equiv          nb201_p_facel,          56
.equiv          nb201_argkrf,           64
.equiv          nb201_argcrf,           72
.equiv          nb201_Vc,               80
.equiv          nb201_type,             88
.equiv          nb201_p_ntype,          96
.equiv          nb201_vdwparam,         104
.equiv          nb201_Vvdw,             112
.equiv          nb201_p_tabscale,       120
.equiv          nb201_VFtab,            128
.equiv          nb201_invsqrta,         136
.equiv          nb201_dvda,             144
.equiv          nb201_p_gbtabscale,     152
.equiv          nb201_GBtab,            160
.equiv          nb201_p_nthreads,       168
.equiv          nb201_count,            176
.equiv          nb201_mtx,              184
.equiv          nb201_outeriter,        192
.equiv          nb201_inneriter,        200
.equiv          nb201_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb201_ixO,              0
.equiv          nb201_iyO,              16
.equiv          nb201_izO,              32
.equiv          nb201_ixH1,             48
.equiv          nb201_iyH1,             64
.equiv          nb201_izH1,             80
.equiv          nb201_ixH2,             96
.equiv          nb201_iyH2,             112
.equiv          nb201_izH2,             128
.equiv          nb201_iqO,              144
.equiv          nb201_iqH,              160
.equiv          nb201_dxO,              176
.equiv          nb201_dyO,              192
.equiv          nb201_dzO,              208
.equiv          nb201_dxH1,             224
.equiv          nb201_dyH1,             240
.equiv          nb201_dzH1,             256
.equiv          nb201_dxH2,             272
.equiv          nb201_dyH2,             288
.equiv          nb201_dzH2,             304
.equiv          nb201_qqO,              320
.equiv          nb201_qqH,              336
.equiv          nb201_vctot,            352
.equiv          nb201_fixO,             384
.equiv          nb201_fiyO,             400
.equiv          nb201_fizO,             416
.equiv          nb201_fixH1,            432
.equiv          nb201_fiyH1,            448
.equiv          nb201_fizH1,            464
.equiv          nb201_fixH2,            480
.equiv          nb201_fiyH2,            496
.equiv          nb201_fizH2,            512
.equiv          nb201_fjx,              528
.equiv          nb201_fjy,              544
.equiv          nb201_fjz,              560
.equiv          nb201_half,             576
.equiv          nb201_three,            592
.equiv          nb201_two,              608
.equiv          nb201_krf,              624
.equiv          nb201_crf,              640
.equiv          nb201_krsqO,            656
.equiv          nb201_krsqH1,           672
.equiv          nb201_krsqH2,           688
.equiv          nb201_nri,              704
.equiv          nb201_iinr,             712
.equiv          nb201_jindex,           720
.equiv          nb201_jjnr,             728
.equiv          nb201_shift,            736
.equiv          nb201_shiftvec,         744
.equiv          nb201_facel,            752
.equiv          nb201_innerjjnr,        760
.equiv          nb201_is3,              768
.equiv          nb201_ii3,              772
.equiv          nb201_innerk,           776
.equiv          nb201_n,                780
.equiv          nb201_nn1,              784
.equiv          nb201_nouter,           788
.equiv          nb201_ninner,           792

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
	sub rsp, 800		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb201_nouter], eax
	mov [rsp + nb201_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb201_nri], edi
	mov [rsp + nb201_iinr], rsi
	mov [rsp + nb201_jindex], rdx
	mov [rsp + nb201_jjnr], rcx
	mov [rsp + nb201_shift], r8
	mov [rsp + nb201_shiftvec], r9
	mov rsi, [rbp + nb201_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb201_facel], xmm0

	mov rsi, [rbp + nb201_argkrf]
	mov rdi, [rbp + nb201_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb201_krf], xmm1
	movapd [rsp + nb201_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb201_half], eax
	mov [rsp + nb201_half+4], ebx
	movsd xmm1, [rsp + nb201_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb201_half], xmm1
	movapd [rsp + nb201_two], xmm2
	movapd [rsp + nb201_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb201_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb201_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	

	movsd xmm5, [rsp + nb201_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb201_iqO], xmm3
	movapd [rsp + nb201_iqH], xmm4
	
.nb201_threadloop:
        mov   rsi, [rbp + nb201_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb201_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb201_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb201_n], eax
        mov [rsp + nb201_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201_outerstart
        jmp .nb201_end

.nb201_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb201_nouter]
	mov [rsp + nb201_nouter], ebx

.nb201_outer:
	mov   rax, [rsp + nb201_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb201_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb201_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb201_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb201_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb201_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb201_ixO], xmm3
	movapd [rsp + nb201_iyO], xmm4
	movapd [rsp + nb201_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb201_ixH1], xmm0
	movapd [rsp + nb201_iyH1], xmm1
	movapd [rsp + nb201_izH1], xmm2
	movapd [rsp + nb201_ixH2], xmm3
	movapd [rsp + nb201_iyH2], xmm4
	movapd [rsp + nb201_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb201_vctot], xmm4
	movapd [rsp + nb201_fixO], xmm4
	movapd [rsp + nb201_fiyO], xmm4
	movapd [rsp + nb201_fizO], xmm4
	movapd [rsp + nb201_fixH1], xmm4
	movapd [rsp + nb201_fiyH1], xmm4
	movapd [rsp + nb201_fizH1], xmm4
	movapd [rsp + nb201_fixH2], xmm4
	movapd [rsp + nb201_fiyH2], xmm4
	movapd [rsp + nb201_fizH2], xmm4
	
	mov   rax, [rsp + nb201_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb201_pos]
	mov   rdi, [rbp + nb201_faction]	
	mov   rax, [rsp + nb201_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb201_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb201_ninner]
	mov   [rsp + nb201_ninner], ecx
	add   edx, 0
	mov   [rsp + nb201_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201_unroll_loop
	jmp   .nb201_checksingle
.nb201_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb201_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb201_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb201_iqO]
	mulpd  xmm4, [rsp + nb201_iqH]
	movapd  [rsp + nb201_qqO], xmm3
	movapd  [rsp + nb201_qqH], xmm4	

	mov rsi, [rbp + nb201_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move j coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb201_ixO]
    subpd xmm1, [rsp + nb201_iyO]
    subpd xmm2, [rsp + nb201_izO]
    subpd xmm3, [rsp + nb201_ixH1]
    subpd xmm4, [rsp + nb201_iyH1]
    subpd xmm5, [rsp + nb201_izH1]
    subpd xmm6, [rsp + nb201_ixH2]
    subpd xmm7, [rsp + nb201_iyH2]
    subpd xmm8, [rsp + nb201_izH2]
    
	movapd [rsp + nb201_dxO], xmm0
	movapd [rsp + nb201_dyO], xmm1
	movapd [rsp + nb201_dzO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb201_dxH1], xmm3
	movapd [rsp + nb201_dyH1], xmm4
	movapd [rsp + nb201_dzH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb201_dxH2], xmm6
	movapd [rsp + nb201_dyH2], xmm7
	movapd [rsp + nb201_dzH2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for j atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb201_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb201_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvO
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb201_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb201_half]
	mulpd   xmm9, xmm15  ;#  rinvO
	mulpd   xmm10, xmm15 ;#   rinvH1
    mulpd   xmm11, xmm15 ;#   rinvH2
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb201_krf]    
    mulpd  xmm9, xmm9   ;# rinvsq
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, xmm2  ;# k*rsq
    mulpd  xmm3, xmm2
    mulpd  xmm6, xmm2
    movapd xmm2, xmm0 ;# copy of k*rsq
    movapd xmm5, xmm3
    movapd xmm8, xmm6
    addpd  xmm2, xmm1  ;# rinv+krsq
    addpd  xmm5, xmm4
    addpd  xmm8, xmm7
    movapd xmm14, [rsp + nb201_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb201_qqO]
    movapd xmm13, [rsp + nb201_qqH]    
    mulpd  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, xmm13
    mulpd  xmm7, xmm13
    addpd  xmm2, [rsp + nb201_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd [rsp + nb201_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j forces to xmm0-xmm2
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb201_dxO]
	mulpd xmm8, [rsp + nb201_dyO]
	mulpd xmm9, [rsp + nb201_dzO]
	mulpd xmm10, [rsp + nb201_dxH1]
	mulpd xmm11, [rsp + nb201_dyH1]
	mulpd xmm12, [rsp + nb201_dzH1]
	mulpd xmm13, [rsp + nb201_dxH2]
	mulpd xmm14, [rsp + nb201_dyH2]
	mulpd xmm15, [rsp + nb201_dzH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb201_fixO]
    addpd xmm8, [rsp + nb201_fiyO]
    addpd xmm9, [rsp + nb201_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb201_fixH1]
    addpd xmm11, [rsp + nb201_fiyH1]
    addpd xmm12, [rsp + nb201_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb201_fixH2]
    addpd xmm14, [rsp + nb201_fiyH2]
    addpd xmm15, [rsp + nb201_fizH2]

    movapd [rsp + nb201_fixO], xmm7
    movapd [rsp + nb201_fiyO], xmm8
    movapd [rsp + nb201_fizO], xmm9
    movapd [rsp + nb201_fixH1], xmm10
    movapd [rsp + nb201_fiyH1], xmm11
    movapd [rsp + nb201_fizH1], xmm12
    movapd [rsp + nb201_fixH2], xmm13
    movapd [rsp + nb201_fiyH2], xmm14
    movapd [rsp + nb201_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb201_innerk],  2
	jl    .nb201_checksingle
	jmp   .nb201_unroll_loop
.nb201_checksingle:	
	mov   edx, [rsp + nb201_innerk]
	and   edx, 1
	jnz   .nb201_dosingle
	jmp   .nb201_updateouterdata
.nb201_dosingle:
	mov   rdx, [rsp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb201_innerjjnr],  4	

	mov rsi, [rbp + nb201_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb201_iqO]
	mulpd  xmm4, [rsp + nb201_iqH]
	movapd  [rsp + nb201_qqO], xmm3
	movapd  [rsp + nb201_qqH], xmm4
	
	mov rsi, [rbp + nb201_pos]       ;# base of pos[] 
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb201_ixO]
	subsd xmm5, [rsp + nb201_iyO]
	subsd xmm6, [rsp + nb201_izO]

	;# store dr 
	movapd [rsp + nb201_dxO], xmm4
	movapd [rsp + nb201_dyO], xmm5
	movapd [rsp + nb201_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb201_ixH1]
	subsd xmm5, [rsp + nb201_iyH1]
	subsd xmm6, [rsp + nb201_izH1]

	;# store dr 
	movapd [rsp + nb201_dxH1], xmm4
	movapd [rsp + nb201_dyH1], xmm5
	movapd [rsp + nb201_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move j coords to xmm3-xmm5
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	;# calc dr 
	subsd xmm3, [rsp + nb201_ixH2]
	subsd xmm4, [rsp + nb201_iyH2]
	subsd xmm5, [rsp + nb201_izH2]

	;# store dr 
	movapd [rsp + nb201_dxH2], xmm3
	movapd [rsp + nb201_dyH2], xmm4
	movapd [rsp + nb201_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	
	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [rsp + nb201_krf]	
	mulsd  xmm1, [rsp + nb201_krf]	
	mulsd  xmm2, [rsp + nb201_krf]	

	movapd [rsp + nb201_krsqH2], xmm0
	movapd [rsp + nb201_krsqH1], xmm1
	movapd [rsp + nb201_krsqO], xmm2
	
	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm3, xmm7
	movapd  xmm0, [rsp + nb201_krsqO]
	addsd   xmm7, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [rsp + nb201_two]
	subsd   xmm7, [rsp + nb201_crf]
	subsd   xmm3, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm7, [rsp + nb201_qqO] ;# vcoul 
	mulsd   xmm3, [rsp + nb201_qqO]
	mulsd  xmm4, xmm3	;# total fsH1 in xmm4 
	
	addsd  xmm7, [rsp + nb201_vctot]

	movapd xmm0, [rsp + nb201_dxO]
	movapd xmm1, [rsp + nb201_dyO]
	movapd xmm2, [rsp + nb201_dzO]
	movlpd [rsp + nb201_vctot], xmm7
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [rsp + nb201_fixO]
	movapd xmm4, [rsp + nb201_fiyO]
	movapd xmm7, [rsp + nb201_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb201_fixO], xmm3
	movlpd [rsp + nb201_fiyO], xmm4
	movlpd [rsp + nb201_fizO], xmm7
	;# update j forces with water O 
	movlpd [rsp + nb201_fjx], xmm0
	movlpd [rsp + nb201_fjy], xmm1
	movlpd [rsp + nb201_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [rsp + nb201_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [rsp + nb201_two]
	subsd   xmm6, [rsp + nb201_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [rsp + nb201_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb201_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb201_vctot]

	movapd xmm0, [rsp + nb201_dxH1]
	movapd xmm1, [rsp + nb201_dyH1]
	movapd xmm2, [rsp + nb201_dzH1]
	movlpd [rsp + nb201_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb201_fixH1]
	movapd xmm4, [rsp + nb201_fiyH1]
	movapd xmm7, [rsp + nb201_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb201_fixH1], xmm3
	movlpd [rsp + nb201_fiyH1], xmm4
	movlpd [rsp + nb201_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb201_fjx]
	addsd  xmm1, [rsp + nb201_fjy]
	addsd  xmm2, [rsp + nb201_fjz]
	movlpd [rsp + nb201_fjx], xmm0
	movlpd [rsp + nb201_fjy], xmm1
	movlpd [rsp + nb201_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [rsp + nb201_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb201_two]
	subsd   xmm5, [rsp + nb201_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb201_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb201_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb201_vctot]

	movapd xmm0, [rsp + nb201_dxH2]
	movapd xmm1, [rsp + nb201_dyH2]
	movapd xmm2, [rsp + nb201_dzH2]
	movlpd [rsp + nb201_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb201_fixH2]
	movapd xmm4, [rsp + nb201_fiyH2]
	movapd xmm7, [rsp + nb201_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb201_fixH2], xmm3
	movlpd [rsp + nb201_fiyH2], xmm4
	movlpd [rsp + nb201_fizH2], xmm7

	mov rdi, [rbp + nb201_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb201_fjx]
	addsd  xmm1, [rsp + nb201_fjy]
	addsd  xmm2, [rsp + nb201_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5	

.nb201_updateouterdata:
	mov   ecx, [rsp + nb201_ii3]
	mov   rdi, [rbp + nb201_faction]
	mov   rsi, [rbp + nb201_fshift]
	mov   edx, [rsp + nb201_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb201_fixO]
	movapd xmm1, [rsp + nb201_fiyO]
	movapd xmm2, [rsp + nb201_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb201_fixH1]
	movapd xmm1, [rsp + nb201_fiyH1]
	movapd xmm2, [rsp + nb201_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 24]
	movsd  xmm4, [rdi + rcx*8 + 32]
	movsd  xmm5, [rdi + rcx*8 + 40]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 24], xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb201_fixH2]
	movapd xmm1, [rsp + nb201_fiyH2]
	movapd xmm2, [rsp + nb201_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 48]
	movsd  xmm4, [rdi + rcx*8 + 56]
	movsd  xmm5, [rdi + rcx*8 + 64]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 48], xmm3
	movsd  [rdi + rcx*8 + 56], xmm4
	movsd  [rdi + rcx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [rsi + rdx*8]
	movhpd xmm3, [rsi + rdx*8 + 8]
	movsd  xmm4, [rsi + rdx*8 + 16]
	subpd  xmm3, xmm6
	subsd  xmm4, xmm7
	movlpd [rsi + rdx*8],      xmm3
	movhpd [rsi + rdx*8 + 8],  xmm3
	movsd  [rsi + rdx*8 + 16], xmm4

	;# get n from stack
	mov esi, [rsp + nb201_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb201_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb201_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb201_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb201_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb201_n], esi
        jmp .nb201_outer
.nb201_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb201_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201_end
        ;# non-zero, do one more workunit
        jmp   .nb201_threadloop
.nb201_end:
	mov eax, [rsp + nb201_nouter]
	mov ebx, [rsp + nb201_ninner]
	mov rcx, [rbp + nb201_outeriter]
	mov rdx, [rbp + nb201_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 800
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





.globl nb_kernel201nf_x86_64_sse2
.globl _nb_kernel201nf_x86_64_sse2
nb_kernel201nf_x86_64_sse2:	
_nb_kernel201nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb201nf_fshift,         16
.equiv          nb201nf_gid,            24
.equiv          nb201nf_pos,            32
.equiv          nb201nf_faction,        40
.equiv          nb201nf_charge,         48
.equiv          nb201nf_p_facel,        56
.equiv          nb201nf_argkrf,         64
.equiv          nb201nf_argcrf,         72
.equiv          nb201nf_Vc,             80
.equiv          nb201nf_type,           88
.equiv          nb201nf_p_ntype,        96
.equiv          nb201nf_vdwparam,       104
.equiv          nb201nf_Vvdw,           112
.equiv          nb201nf_p_tabscale,     120
.equiv          nb201nf_VFtab,          128
.equiv          nb201nf_invsqrta,       136
.equiv          nb201nf_dvda,           144
.equiv          nb201nf_p_gbtabscale,   152
.equiv          nb201nf_GBtab,          160
.equiv          nb201nf_p_nthreads,     168
.equiv          nb201nf_count,          176
.equiv          nb201nf_mtx,            184
.equiv          nb201nf_outeriter,      192
.equiv          nb201nf_inneriter,      200
.equiv          nb201nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb201nf_ixO,            0
.equiv          nb201nf_iyO,            16
.equiv          nb201nf_izO,            32
.equiv          nb201nf_ixH1,           48
.equiv          nb201nf_iyH1,           64
.equiv          nb201nf_izH1,           80
.equiv          nb201nf_ixH2,           96
.equiv          nb201nf_iyH2,           112
.equiv          nb201nf_izH2,           128
.equiv          nb201nf_iqO,            144
.equiv          nb201nf_iqH,            160
.equiv          nb201nf_qqO,            176
.equiv          nb201nf_qqH,            192
.equiv          nb201nf_vctot,          208
.equiv          nb201nf_half,           224
.equiv          nb201nf_three,          240
.equiv          nb201nf_krf,            256
.equiv          nb201nf_crf,            272
.equiv          nb201nf_krsqO,          288
.equiv          nb201nf_krsqH1,         304
.equiv          nb201nf_krsqH2,         320
.equiv          nb201nf_nri,            336
.equiv          nb201nf_iinr,           344
.equiv          nb201nf_jindex,         352
.equiv          nb201nf_jjnr,           360
.equiv          nb201nf_shift,          368
.equiv          nb201nf_shiftvec,       376
.equiv          nb201nf_facel,          384
.equiv          nb201nf_innerjjnr,      392
.equiv          nb201nf_is3,            400
.equiv          nb201nf_ii3,            404
.equiv          nb201nf_innerk,         408
.equiv          nb201nf_n,              412
.equiv          nb201nf_nn1,            416
.equiv          nb201nf_nouter,         420
.equiv          nb201nf_ninner,         424

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
	sub rsp, 432		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb201nf_nouter], eax
	mov [rsp + nb201nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb201nf_nri], edi
	mov [rsp + nb201nf_iinr], rsi
	mov [rsp + nb201nf_jindex], rdx
	mov [rsp + nb201nf_jjnr], rcx
	mov [rsp + nb201nf_shift], r8
	mov [rsp + nb201nf_shiftvec], r9
	mov rsi, [rbp + nb201nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb201nf_facel], xmm0

	mov rsi, [rbp + nb201nf_argkrf]
	mov rdi, [rbp + nb201nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb201nf_krf], xmm1
	movapd [rsp + nb201nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb201nf_half], eax
	mov [rsp + nb201nf_half+4], ebx
	movsd xmm1, [rsp + nb201nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb201nf_half], xmm1
	movapd [rsp + nb201nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb201nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb201nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	

	movsd xmm5, [rsp + nb201nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb201nf_iqO], xmm3
	movapd [rsp + nb201nf_iqH], xmm4
	
.nb201nf_threadloop:
        mov   rsi, [rbp + nb201nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb201nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb201nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb201nf_n], eax
        mov [rsp + nb201nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201nf_outerstart
        jmp .nb201nf_end

.nb201nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb201nf_nouter]
	mov [rsp + nb201nf_nouter], ebx

.nb201nf_outer:
	mov   rax, [rsp + nb201nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb201nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb201nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb201nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb201nf_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb201nf_ixO], xmm3
	movapd [rsp + nb201nf_iyO], xmm4
	movapd [rsp + nb201nf_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb201nf_ixH1], xmm0
	movapd [rsp + nb201nf_iyH1], xmm1
	movapd [rsp + nb201nf_izH1], xmm2
	movapd [rsp + nb201nf_ixH2], xmm3
	movapd [rsp + nb201nf_iyH2], xmm4
	movapd [rsp + nb201nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb201nf_vctot], xmm4
	
	mov   rax, [rsp + nb201nf_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb201nf_pos]
	mov   rax, [rsp + nb201nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb201nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb201nf_ninner]
	mov   [rsp + nb201nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb201nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201nf_unroll_loop
	jmp   .nb201nf_checksingle
.nb201nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb201nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb201nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb201nf_iqO]
	mulpd  xmm4, [rsp + nb201nf_iqH]
	movapd  [rsp + nb201nf_qqO], xmm3
	movapd  [rsp + nb201nf_qqH], xmm4	

	mov rsi, [rbp + nb201nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb201nf_ixO]
	movapd xmm5, [rsp + nb201nf_iyO]
	movapd xmm6, [rsp + nb201nf_izO]

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
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb201nf_ixH1]
	movapd xmm5, [rsp + nb201nf_iyH1]
	movapd xmm6, [rsp + nb201nf_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb201nf_ixH2]
	movapd xmm4, [rsp + nb201nf_iyH2]
	movapd xmm5, [rsp + nb201nf_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulpd  xmm0, [rsp + nb201nf_krf]	
	mulpd  xmm1, [rsp + nb201nf_krf]	
	mulpd  xmm2, [rsp + nb201nf_krf]	

	movapd [rsp + nb201nf_krsqH2], xmm0
	movapd [rsp + nb201nf_krsqH1], xmm1
	movapd [rsp + nb201nf_krsqO], xmm2
	
	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm0, [rsp + nb201nf_krsqO]
	addpd   xmm7, xmm0	;# xmm7=rinv+ krsq 
	subpd   xmm7, [rsp + nb201nf_crf]
	mulpd   xmm7, [rsp + nb201nf_qqO] ;# vcoul 	
	addpd  xmm7, [rsp + nb201nf_vctot]

	;# H1 interactions 
	movapd  xmm0, [rsp + nb201nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [rsp + nb201nf_crf]
	mulpd   xmm6, [rsp + nb201nf_qqH] ;# vcoul 
	addpd  xmm6, xmm7

	;# H2 interactions 
	movapd  xmm0, [rsp + nb201nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [rsp + nb201nf_crf]
	mulpd   xmm5, [rsp + nb201nf_qqH] ;# vcoul 
	addpd  xmm5, xmm6
	movapd [rsp + nb201nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb201nf_innerk],  2
	jl    .nb201nf_checksingle
	jmp   .nb201nf_unroll_loop
.nb201nf_checksingle:	
	mov   edx, [rsp + nb201nf_innerk]
	and   edx, 1
	jnz   .nb201nf_dosingle
	jmp   .nb201nf_updateouterdata
.nb201nf_dosingle:
	mov   rdx, [rsp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb201nf_innerjjnr],  4	

	mov rsi, [rbp + nb201nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb201nf_iqO]
	mulpd  xmm4, [rsp + nb201nf_iqH]
	movapd  [rsp + nb201nf_qqO], xmm3
	movapd  [rsp + nb201nf_qqH], xmm4
	
	mov rsi, [rbp + nb201nf_pos]       ;# base of pos[] 
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb201nf_ixO]
	movapd xmm5, [rsp + nb201nf_iyO]
	movapd xmm6, [rsp + nb201nf_izO]

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
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb201nf_ixH1]
	movapd xmm5, [rsp + nb201nf_iyH1]
	movapd xmm6, [rsp + nb201nf_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb201nf_ixH2]
	movapd xmm4, [rsp + nb201nf_iyH2]
	movapd xmm5, [rsp + nb201nf_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	
	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [rsp + nb201nf_krf]	
	mulsd  xmm1, [rsp + nb201nf_krf]	
	mulsd  xmm2, [rsp + nb201nf_krf]	

	movapd [rsp + nb201nf_krsqH2], xmm0
	movapd [rsp + nb201nf_krsqH1], xmm1
	movapd [rsp + nb201nf_krsqO], xmm2
	
	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb201nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb201nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb201nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb201nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm0, [rsp + nb201nf_krsqO]
	addsd   xmm7, xmm0	;# xmm7=rinv+ krsq 
	subsd   xmm7, [rsp + nb201nf_crf]
	mulsd   xmm7, [rsp + nb201nf_qqO] ;# vcoul 	
	addsd  xmm7, [rsp + nb201nf_vctot]

	;# H1 interactions 
	movapd  xmm0, [rsp + nb201nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [rsp + nb201nf_crf]
	mulsd   xmm6, [rsp + nb201nf_qqH] ;# vcoul 
	addsd  xmm6, xmm7

	;# H2 interactions 
	movapd  xmm0, [rsp + nb201nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [rsp + nb201nf_crf]
	mulsd   xmm5, [rsp + nb201nf_qqH] ;# vcoul 
	addsd  xmm5, xmm6
	movlpd [rsp + nb201nf_vctot], xmm5
	
.nb201nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb201nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb201nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb201nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb201nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb201nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb201nf_n], esi
        jmp .nb201nf_outer
.nb201nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb201nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201nf_end
        ;# non-zero, do one more workunit
        jmp   .nb201nf_threadloop
.nb201nf_end:
	mov eax, [rsp + nb201nf_nouter]
	mov ebx, [rsp + nb201nf_ninner]
	mov rcx, [rbp + nb201nf_outeriter]
	mov rdx, [rbp + nb201nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 432
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
