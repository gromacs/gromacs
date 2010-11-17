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

	
	

.globl nb_kernel101_x86_64_sse2
.globl _nb_kernel101_x86_64_sse2
nb_kernel101_x86_64_sse2:	
_nb_kernel101_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb101_fshift,           16
.equiv          nb101_gid,              24
.equiv          nb101_pos,              32
.equiv          nb101_faction,          40
.equiv          nb101_charge,           48
.equiv          nb101_p_facel,          56
.equiv          nb101_argkrf,           64
.equiv          nb101_argcrf,           72
.equiv          nb101_Vc,               80
.equiv          nb101_type,             88
.equiv          nb101_p_ntype,          96
.equiv          nb101_vdwparam,         104
.equiv          nb101_Vvdw,             112
.equiv          nb101_p_tabscale,       120
.equiv          nb101_VFtab,            128
.equiv          nb101_invsqrta,         136
.equiv          nb101_dvda,             144
.equiv          nb101_p_gbtabscale,     152
.equiv          nb101_GBtab,            160
.equiv          nb101_p_nthreads,       168
.equiv          nb101_count,            176
.equiv          nb101_mtx,              184
.equiv          nb101_outeriter,        192
.equiv          nb101_inneriter,        200
.equiv          nb101_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb101_ixO,              0
.equiv          nb101_iyO,              16
.equiv          nb101_izO,              32
.equiv          nb101_ixH1,             48
.equiv          nb101_iyH1,             64
.equiv          nb101_izH1,             80
.equiv          nb101_ixH2,             96
.equiv          nb101_iyH2,             112
.equiv          nb101_izH2,             128
.equiv          nb101_iqO,              144
.equiv          nb101_iqH,              160
.equiv          nb101_dxO,              176
.equiv          nb101_dyO,              192
.equiv          nb101_dzO,              208
.equiv          nb101_dxH1,             224
.equiv          nb101_dyH1,             240
.equiv          nb101_dzH1,             256
.equiv          nb101_dxH2,             272
.equiv          nb101_dyH2,             288
.equiv          nb101_dzH2,             304
.equiv          nb101_qqO,              320
.equiv          nb101_qqH,              336
.equiv          nb101_vctot,            352
.equiv          nb101_fixO,             368
.equiv          nb101_fiyO,             384
.equiv          nb101_fizO,             400
.equiv          nb101_fixH1,            416
.equiv          nb101_fiyH1,            432
.equiv          nb101_fizH1,            448
.equiv          nb101_fixH2,            464
.equiv          nb101_fiyH2,            480
.equiv          nb101_fizH2,            496
.equiv          nb101_fjx,              512
.equiv          nb101_fjy,              528
.equiv          nb101_fjz,              544
.equiv          nb101_half,             560
.equiv          nb101_three,            576
.equiv          nb101_nri,              592
.equiv          nb101_innerjjnr,        600
.equiv          nb101_iinr,             608
.equiv          nb101_jindex,           616
.equiv          nb101_jjnr,             624
.equiv          nb101_shift,            632
.equiv          nb101_shiftvec,         640
.equiv          nb101_facel,            648
.equiv          nb101_is3,              656
.equiv          nb101_ii3,              660
.equiv          nb101_innerk,           664
.equiv          nb101_n,                668
.equiv          nb101_nn1,              672
.equiv          nb101_nouter,           676
.equiv          nb101_ninner,           680

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
	sub rsp, 688		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb101_nouter], eax
	mov [rsp + nb101_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb101_nri], edi
	mov [rsp + nb101_iinr], rsi
	mov [rsp + nb101_jindex], rdx
	mov [rsp + nb101_jjnr], rcx
	mov [rsp + nb101_shift], r8
	mov [rsp + nb101_shiftvec], r9
	mov rsi, [rbp + nb101_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb101_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb101_half], eax
	mov [rsp + nb101_half+4], ebx
	movsd xmm1, [rsp + nb101_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb101_half], xmm1
	movapd [rsp + nb101_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb101_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb101_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb101_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb101_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb101_iqO], xmm3
	movapd [rsp + nb101_iqH], xmm4
	
.nb101_threadloop:
        mov   rsi, [rbp + nb101_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb101_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb101_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb101_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb101_n], eax
        mov [rsp + nb101_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb101_outerstart
        jmp .nb101_end

.nb101_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb101_nouter]
	mov [rsp + nb101_nouter], ebx

.nb101_outer:
	mov   rax, [rsp + nb101_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb101_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb101_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb101_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb101_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb101_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb101_ixO], xmm3
	movapd [rsp + nb101_iyO], xmm4
	movapd [rsp + nb101_izO], xmm5

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
	movapd [rsp + nb101_ixH1], xmm0
	movapd [rsp + nb101_iyH1], xmm1
	movapd [rsp + nb101_izH1], xmm2
	movapd [rsp + nb101_ixH2], xmm3
	movapd [rsp + nb101_iyH2], xmm4
	movapd [rsp + nb101_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb101_vctot], xmm4
	movapd [rsp + nb101_fixO], xmm4
	movapd [rsp + nb101_fiyO], xmm4
	movapd [rsp + nb101_fizO], xmm4
	movapd [rsp + nb101_fixH1], xmm4
	movapd [rsp + nb101_fiyH1], xmm4
	movapd [rsp + nb101_fizH1], xmm4
	movapd [rsp + nb101_fixH2], xmm4
	movapd [rsp + nb101_fiyH2], xmm4
	movapd [rsp + nb101_fizH2], xmm4
	
	mov   rax, [rsp + nb101_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb101_pos]
	mov   rdi, [rbp + nb101_faction]	
	mov   rax, [rsp + nb101_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb101_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb101_ninner]
	mov   [rsp + nb101_ninner], ecx
	add   edx, 0
	mov   [rsp + nb101_innerk], edx    ;# number of innerloop atoms 
	jge   .nb101_unroll_loop
	jmp   .nb101_checksingle
.nb101_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb101_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              

	add qword ptr [rsp + nb101_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb101_charge]    ;# base of charge[] 
	
	movlpd xmm6, [rsi + rax*8]	;# jq A 
	movhpd xmm6, [rsi + rbx*8]	;# jq B 
	movapd xmm3, [rsp + nb101_iqO]
	movapd xmm4, [rsp + nb101_iqH]
	mulpd xmm3, xmm6		;# qqO 
	mulpd xmm4, xmm6		;# qqH 
	
	movapd  [rsp + nb101_qqO], xmm3
	movapd  [rsp + nb101_qqH], xmm4	

	mov rsi, [rbp + nb101_pos]       ;# base of pos[] 

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
    
    subpd xmm0, [rsp + nb101_ixO]
    subpd xmm1, [rsp + nb101_iyO]
    subpd xmm2, [rsp + nb101_izO]
    subpd xmm3, [rsp + nb101_ixH1]
    subpd xmm4, [rsp + nb101_iyH1]
    subpd xmm5, [rsp + nb101_izH1]
    subpd xmm6, [rsp + nb101_ixH2]
    subpd xmm7, [rsp + nb101_iyH2]
    subpd xmm8, [rsp + nb101_izH2]
    
	movapd [rsp + nb101_dxO], xmm0
	movapd [rsp + nb101_dyO], xmm1
	movapd [rsp + nb101_dzO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb101_dxH1], xmm3
	movapd [rsp + nb101_dyH1], xmm4
	movapd [rsp + nb101_dzH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb101_dxH2], xmm6
	movapd [rsp + nb101_dyH2], xmm7
	movapd [rsp + nb101_dzH2], xmm8
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
		
	movapd  xmm9, [rsp + nb101_three]
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

	movapd  xmm15, [rsp + nb101_half]
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
		
	movapd  xmm1, [rsp + nb101_three]
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

	movapd  xmm15, [rsp + nb101_half]
	mulpd   xmm9, xmm15  ;#  rinvO 
	mulpd   xmm10, xmm15 ;#   rinvH1
    mulpd   xmm11, xmm15 ;#   rinvH2
	
	;# interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb101_qqO] 
    mulpd  xmm1, [rsp + nb101_qqH] 
    mulpd  xmm2, [rsp + nb101_qqH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb101_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb101_vctot], xmm0
    
    ;# move j forces to xmm0-xmm2
	mov   rdi, [rbp + nb101_faction]	
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

	mulpd xmm7, [rsp + nb101_dxO]
	mulpd xmm8, [rsp + nb101_dyO]
	mulpd xmm9, [rsp + nb101_dzO]
	mulpd xmm10, [rsp + nb101_dxH1]
	mulpd xmm11, [rsp + nb101_dyH1]
	mulpd xmm12, [rsp + nb101_dzH1]
	mulpd xmm13, [rsp + nb101_dxH2]
	mulpd xmm14, [rsp + nb101_dyH2]
	mulpd xmm15, [rsp + nb101_dzH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb101_fixO]
    addpd xmm8, [rsp + nb101_fiyO]
    addpd xmm9, [rsp + nb101_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb101_fixH1]
    addpd xmm11, [rsp + nb101_fiyH1]
    addpd xmm12, [rsp + nb101_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb101_fixH2]
    addpd xmm14, [rsp + nb101_fiyH2]
    addpd xmm15, [rsp + nb101_fizH2]

    movapd [rsp + nb101_fixO], xmm7
    movapd [rsp + nb101_fiyO], xmm8
    movapd [rsp + nb101_fizO], xmm9
    movapd [rsp + nb101_fixH1], xmm10
    movapd [rsp + nb101_fiyH1], xmm11
    movapd [rsp + nb101_fizH1], xmm12
    movapd [rsp + nb101_fixH2], xmm13
    movapd [rsp + nb101_fiyH2], xmm14
    movapd [rsp + nb101_fizH2], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb101_innerk],  2
	jl    .nb101_checksingle
	jmp   .nb101_unroll_loop
.nb101_checksingle:				
	mov   edx, [rsp + nb101_innerk]
	and   edx, 1
	jnz    .nb101_dosingle
	jmp    .nb101_updateouterdata
.nb101_dosingle:
	mov   rdx, [rsp + nb101_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb101_charge]    ;# base of charge[] 
	xorpd xmm6, xmm6
	movlpd xmm6, [rsi + rax*8]	;# jq A 
	
	movapd xmm3, [rsp + nb101_iqO]
	movapd xmm4, [rsp + nb101_iqH]
	mulsd xmm3, xmm6		;# qqO 
	mulsd xmm4, xmm6		;# qqH 
	
	movapd  [rsp + nb101_qqO], xmm3
	movapd  [rsp + nb101_qqH], xmm4	

	mov rsi, [rbp + nb101_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb101_ixO]
	subsd xmm5, [rsp + nb101_iyO]
	subsd xmm6, [rsp + nb101_izO]
    

	;# store dr 
	movapd [rsp + nb101_dxO], xmm4
	movapd [rsp + nb101_dyO], xmm5
	movapd [rsp + nb101_dzO], xmm6
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
	subsd xmm4, [rsp + nb101_ixH1]
	subsd xmm5, [rsp + nb101_iyH1]
	subsd xmm6, [rsp + nb101_izH1]

	;# store dr 
	movapd [rsp + nb101_dxH1], xmm4
	movapd [rsp + nb101_dyH1], xmm5
	movapd [rsp + nb101_dzH1], xmm6
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
	subsd xmm3, [rsp + nb101_ixH2]
	subsd xmm4, [rsp + nb101_iyH2]
	subsd xmm5, [rsp + nb101_izH2]

	;# store dr 
	movapd [rsp + nb101_dxH2], xmm3
	movapd [rsp + nb101_dyH2], xmm4
	movapd [rsp + nb101_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	mulsd  xmm7, [rsp + nb101_qqO]	;# xmm7=vcoul 
	
	mulsd  xmm4, xmm7	;# total fsO in xmm4 

	addsd  xmm7, [rsp + nb101_vctot]
	
	movlpd [rsp + nb101_vctot], xmm7

	movapd xmm0, [rsp + nb101_dxO]
	movapd xmm1, [rsp + nb101_dyO]
	movapd xmm2, [rsp + nb101_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [rsp + nb101_fixO]
	movapd xmm4, [rsp + nb101_fiyO]
	movapd xmm7, [rsp + nb101_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb101_fixO], xmm3
	movlpd [rsp + nb101_fiyO], xmm4
	movlpd [rsp + nb101_fizO], xmm7
	;# update j forces with water O 
	movlpd [rsp + nb101_fjx], xmm0
	movlpd [rsp + nb101_fjy], xmm1
	movlpd [rsp + nb101_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [rsp + nb101_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb101_vctot]

	movapd xmm0, [rsp + nb101_dxH1]
	movapd xmm1, [rsp + nb101_dyH1]
	movapd xmm2, [rsp + nb101_dzH1]
	movlpd [rsp + nb101_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb101_fixH1]
	movapd xmm4, [rsp + nb101_fiyH1]
	movapd xmm7, [rsp + nb101_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb101_fixH1], xmm3
	movlpd [rsp + nb101_fiyH1], xmm4
	movlpd [rsp + nb101_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb101_fjx]
	addsd  xmm1, [rsp + nb101_fjy]
	addsd  xmm2, [rsp + nb101_fjz]
	movsd [rsp + nb101_fjx], xmm0
	movsd [rsp + nb101_fjy], xmm1
	movsd [rsp + nb101_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [rsp + nb101_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [rsp + nb101_vctot]

	movapd xmm0, [rsp + nb101_dxH2]
	movapd xmm1, [rsp + nb101_dyH2]
	movapd xmm2, [rsp + nb101_dzH2]
	movlpd [rsp + nb101_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb101_fixH2]
	movapd xmm4, [rsp + nb101_fiyH2]
	movapd xmm7, [rsp + nb101_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb101_fixH2], xmm3
	movlpd [rsp + nb101_fiyH2], xmm4
	movlpd [rsp + nb101_fizH2], xmm7

	mov rdi, [rbp + nb101_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb101_fjx]
	addsd  xmm1, [rsp + nb101_fjy]
	addsd  xmm2, [rsp + nb101_fjz]

	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb101_updateouterdata:
	mov   ecx, [rsp + nb101_ii3]
	mov   rdi, [rbp + nb101_faction]
	mov   rsi, [rbp + nb101_fshift]
	mov   edx, [rsp + nb101_is3]

	;# accumulate Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb101_fixO]
	movapd xmm1, [rsp + nb101_fiyO]
	movapd xmm2, [rsp + nb101_fizO]

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
	movapd xmm0, [rsp + nb101_fixH1]
	movapd xmm1, [rsp + nb101_fiyH1]
	movapd xmm2, [rsp + nb101_fizH1]

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
	movapd xmm0, [rsp + nb101_fixH2]
	movapd xmm1, [rsp + nb101_fiyH2]
	movapd xmm2, [rsp + nb101_fizH2]

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
	mov esi, [rsp + nb101_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb101_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb101_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   rax, [rbp + nb101_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb101_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb101_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb101_n], esi
        jmp .nb101_outer
.nb101_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb101_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb101_end
        ;# non-zero, do one more workunit
        jmp   .nb101_threadloop
.nb101_end:

	mov eax, [rsp + nb101_nouter]
	mov ebx, [rsp + nb101_ninner]
	mov rcx, [rbp + nb101_outeriter]
	mov rdx, [rbp + nb101_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

    add rsp, 688
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




.globl nb_kernel101nf_x86_64_sse2
.globl _nb_kernel101nf_x86_64_sse2
nb_kernel101nf_x86_64_sse2:	
_nb_kernel101nf_x86_64_sse2:	
.equiv          nb101nf_fshift,         16
.equiv          nb101nf_gid,            24
.equiv          nb101nf_pos,            32
.equiv          nb101nf_faction,        40
.equiv          nb101nf_charge,         48
.equiv          nb101nf_p_facel,        56
.equiv          nb101nf_argkrf,         64
.equiv          nb101nf_argcrf,         72
.equiv          nb101nf_Vc,             80
.equiv          nb101nf_type,           88
.equiv          nb101nf_p_ntype,        96
.equiv          nb101nf_vdwparam,       104
.equiv          nb101nf_Vvdw,           112
.equiv          nb101nf_p_tabscale,     120
.equiv          nb101nf_VFtab,          128
.equiv          nb101nf_invsqrta,       136
.equiv          nb101nf_dvda,           144
.equiv          nb101nf_p_gbtabscale,   152
.equiv          nb101nf_GBtab,          160
.equiv          nb101nf_p_nthreads,     168
.equiv          nb101nf_count,          176
.equiv          nb101nf_mtx,            184
.equiv          nb101nf_outeriter,      192
.equiv          nb101nf_inneriter,      200
.equiv          nb101nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb101nf_ixO,            0
.equiv          nb101nf_iyO,            16
.equiv          nb101nf_izO,            32
.equiv          nb101nf_ixH1,           48
.equiv          nb101nf_iyH1,           64
.equiv          nb101nf_izH1,           80
.equiv          nb101nf_ixH2,           96
.equiv          nb101nf_iyH2,           112
.equiv          nb101nf_izH2,           128
.equiv          nb101nf_iqO,            144
.equiv          nb101nf_iqH,            160
.equiv          nb101nf_qqO,            176
.equiv          nb101nf_qqH,            192
.equiv          nb101nf_vctot,          208
.equiv          nb101nf_half,           224
.equiv          nb101nf_three,          240
.equiv          nb101nf_is3,            256
.equiv          nb101nf_ii3,            260
.equiv          nb101nf_nri,            264
.equiv          nb101nf_iinr,           272
.equiv          nb101nf_jindex,         280
.equiv          nb101nf_jjnr,           288
.equiv          nb101nf_shift,          296
.equiv          nb101nf_shiftvec,       304
.equiv          nb101nf_facel,          312
.equiv          nb101nf_innerjjnr,      320
.equiv          nb101nf_innerk,         328
.equiv          nb101nf_n,              332
.equiv          nb101nf_nn1,            336
.equiv          nb101nf_nouter,         340
.equiv          nb101nf_ninner,         344

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
	sub rsp, 352		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb101nf_nouter], eax
	mov [rsp + nb101nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb101nf_nri], edi
	mov [rsp + nb101nf_iinr], rsi
	mov [rsp + nb101nf_jindex], rdx
	mov [rsp + nb101nf_jjnr], rcx
	mov [rsp + nb101nf_shift], r8
	mov [rsp + nb101nf_shiftvec], r9
	mov rsi, [rbp + nb101nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb101nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb101nf_half], eax
	mov [rsp + nb101nf_half+4], ebx
	movsd xmm1, [rsp + nb101nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb101nf_half], xmm1
	movapd [rsp + nb101nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb101nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb101nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb101nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb101nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb101nf_iqO], xmm3
	movapd [rsp + nb101nf_iqH], xmm4
	
.nb101nf_threadloop:
        mov   rsi, [rbp + nb101nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb101nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb101nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb101nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb101nf_n], eax
        mov [rsp + nb101nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb101nf_outerstart
        jmp .nb101nf_end

.nb101nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb101nf_nouter]
	mov [rsp + nb101nf_nouter], ebx

.nb101nf_outer:
	mov   rax, [rsp + nb101nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb101nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb101nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb101nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb101nf_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb101nf_ixO], xmm3
	movapd [rsp + nb101nf_iyO], xmm4
	movapd [rsp + nb101nf_izO], xmm5

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
	movapd [rsp + nb101nf_ixH1], xmm0
	movapd [rsp + nb101nf_iyH1], xmm1
	movapd [rsp + nb101nf_izH1], xmm2
	movapd [rsp + nb101nf_ixH2], xmm3
	movapd [rsp + nb101nf_iyH2], xmm4
	movapd [rsp + nb101nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb101nf_vctot], xmm4
	
	mov   rax, [rsp + nb101nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb101nf_pos]
	mov   rax, [rsp + nb101nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb101nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb101nf_ninner]
	mov   [rsp + nb101nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb101nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb101nf_unroll_loop
	jmp   .nb101nf_checksingle
.nb101nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb101nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              

	add qword ptr [rsp + nb101nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb101nf_charge]    ;# base of charge[] 
	
	
	movlpd xmm6, [rsi + rax*8]	;# jq A 
	movhpd xmm6, [rsi + rbx*8]	;# jq B 
	movapd xmm3, [rsp + nb101nf_iqO]
	movapd xmm4, [rsp + nb101nf_iqH]
	mulpd xmm3, xmm6		;# qqO 
	mulpd xmm4, xmm6		;# qqH 
	
	movapd  [rsp + nb101nf_qqO], xmm3
	movapd  [rsp + nb101nf_qqH], xmm4	

	mov rsi, [rbp + nb101nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [rsp + nb101nf_ixO]
	movapd xmm5, [rsp + nb101nf_iyO]
	movapd xmm6, [rsp + nb101nf_izO]

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
	movapd xmm4, [rsp + nb101nf_ixH1]
	movapd xmm5, [rsp + nb101nf_iyH1]
	movapd xmm6, [rsp + nb101nf_izH1]

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
	movapd xmm3, [rsp + nb101nf_ixH2]
	movapd xmm4, [rsp + nb101nf_iyH2]
	movapd xmm5, [rsp + nb101nf_izH2]

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

	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	mulpd  xmm7, [rsp + nb101nf_qqO]	;# xmm7=vcoul 
	addpd  xmm7, [rsp + nb101nf_vctot]
	movapd [rsp + nb101nf_vctot], xmm7

	;# H1 interactions 
	mulpd  xmm6, [rsp + nb101nf_qqH]	;# xmm6=vcoul 
	addpd  xmm6, [rsp + nb101nf_vctot]
	movapd [rsp + nb101nf_vctot], xmm6

	;# H2 interactions 
	mulpd  xmm5, [rsp + nb101nf_qqH]	;# xmm5=vcoul 
	addpd  xmm5, [rsp + nb101nf_vctot]
	movapd [rsp + nb101nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb101nf_innerk],  2
	jl    .nb101nf_checksingle
	jmp   .nb101nf_unroll_loop
.nb101nf_checksingle:				
	mov   edx, [rsp + nb101nf_innerk]
	and   edx, 1
	jnz   .nb101nf_dosingle
	jmp   .nb101nf_updateouterdata
.nb101nf_dosingle:
	mov   rdx, [rsp + nb101nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb101nf_charge]    ;# base of charge[] 
	xorpd xmm6, xmm6
	movlpd xmm6, [rsi + rax*8]	;# jq A 
	
	movapd xmm3, [rsp + nb101nf_iqO]
	movapd xmm4, [rsp + nb101nf_iqH]
	mulsd xmm3, xmm6		;# qqO 
	mulsd xmm4, xmm6		;# qqH 
	
	movapd  [rsp + nb101nf_qqO], xmm3
	movapd  [rsp + nb101nf_qqH], xmm4	

	mov rsi, [rbp + nb101nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb101nf_ixO]
	movapd xmm5, [rsp + nb101nf_iyO]
	movapd xmm6, [rsp + nb101nf_izO]

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
	movapd xmm4, [rsp + nb101nf_ixH1]
	movapd xmm5, [rsp + nb101nf_iyH1]
	movapd xmm6, [rsp + nb101nf_izH1]

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
	movapd xmm3, [rsp + nb101nf_ixH2]
	movapd xmm4, [rsp + nb101nf_iyH2]
	movapd xmm5, [rsp + nb101nf_izH2]

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

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb101nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb101nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb101nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb101nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	mulsd  xmm7, [rsp + nb101nf_qqO]	;# xmm7=vcoul 
	addsd  xmm7, [rsp + nb101nf_vctot]
	movlpd [rsp + nb101nf_vctot], xmm7

	;# H1 interactions 
	mulsd  xmm6, [rsp + nb101nf_qqH]	;# xmm6=vcoul 
	addsd  xmm6, [rsp + nb101nf_vctot]
	movlpd [rsp + nb101nf_vctot], xmm6

	;# H2 interactions 
	mulsd  xmm5, [rsp + nb101nf_qqH]	;# xmm5=vcoul 
	addsd  xmm5, [rsp + nb101nf_vctot]
	movlpd [rsp + nb101nf_vctot], xmm5

.nb101nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb101nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb101nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	movapd xmm7, [rsp + nb101nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   rax, [rbp + nb101nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb101nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb101nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb101nf_n], esi
        jmp .nb101nf_outer
.nb101nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb101nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb101nf_end
        ;# non-zero, do one more workunit
        jmp   .nb101nf_threadloop
.nb101nf_end:

	mov eax, [rsp + nb101nf_nouter]
	mov ebx, [rsp + nb101nf_ninner]
	mov rcx, [rbp + nb101nf_outeriter]
	mov rdx, [rbp + nb101nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 352
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
