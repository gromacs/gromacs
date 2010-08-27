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


	
.globl nb_kernel202_x86_64_sse2
.globl _nb_kernel202_x86_64_sse2
nb_kernel202_x86_64_sse2:	
_nb_kernel202_x86_64_sse2:	
.equiv          nb202_fshift,           16
.equiv          nb202_gid,              24
.equiv          nb202_pos,              32
.equiv          nb202_faction,          40
.equiv          nb202_charge,           48
.equiv          nb202_p_facel,          56
.equiv          nb202_argkrf,           64
.equiv          nb202_argcrf,           72
.equiv          nb202_Vc,               80
.equiv          nb202_type,             88
.equiv          nb202_p_ntype,          96
.equiv          nb202_vdwparam,         104
.equiv          nb202_Vvdw,             112
.equiv          nb202_p_tabscale,       120
.equiv          nb202_VFtab,            128
.equiv          nb202_invsqrta,         136
.equiv          nb202_dvda,             144
.equiv          nb202_p_gbtabscale,     152
.equiv          nb202_GBtab,            160
.equiv          nb202_p_nthreads,       168
.equiv          nb202_count,            176
.equiv          nb202_mtx,              184
.equiv          nb202_outeriter,        192
.equiv          nb202_inneriter,        200
.equiv          nb202_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb202_ixO,              0
.equiv          nb202_iyO,              16
.equiv          nb202_izO,              32
.equiv          nb202_ixH1,             48
.equiv          nb202_iyH1,             64
.equiv          nb202_izH1,             80
.equiv          nb202_ixH2,             96
.equiv          nb202_iyH2,             112
.equiv          nb202_izH2,             128
.equiv          nb202_jxO,              144
.equiv          nb202_jyO,              160
.equiv          nb202_jzO,              176
.equiv          nb202_jxH1,             192
.equiv          nb202_jyH1,             208
.equiv          nb202_jzH1,             224
.equiv          nb202_jxH2,             240
.equiv          nb202_jyH2,             256
.equiv          nb202_jzH2,             272
.equiv          nb202_dxOO,             288
.equiv          nb202_dyOO,             304
.equiv          nb202_dzOO,             320
.equiv          nb202_dxOH1,            336
.equiv          nb202_dyOH1,            352
.equiv          nb202_dzOH1,            368
.equiv          nb202_dxOH2,            384
.equiv          nb202_dyOH2,            400
.equiv          nb202_dzOH2,            416
.equiv          nb202_dxH1O,            432
.equiv          nb202_dyH1O,            448
.equiv          nb202_dzH1O,            464
.equiv          nb202_dxH1H1,           480
.equiv          nb202_dyH1H1,           496
.equiv          nb202_dzH1H1,           512
.equiv          nb202_dxH1H2,           528
.equiv          nb202_dyH1H2,           544
.equiv          nb202_dzH1H2,           560
.equiv          nb202_dxH2O,            576
.equiv          nb202_dyH2O,            592
.equiv          nb202_dzH2O,            608
.equiv          nb202_dxH2H1,           624
.equiv          nb202_dyH2H1,           640
.equiv          nb202_dzH2H1,           656
.equiv          nb202_dxH2H2,           672
.equiv          nb202_dyH2H2,           688
.equiv          nb202_dzH2H2,           704
.equiv          nb202_qqOO,             720
.equiv          nb202_qqOH,             736
.equiv          nb202_qqHH,             752
.equiv          nb202_vctot,            768
.equiv          nb202_fixO,             784
.equiv          nb202_fiyO,             800
.equiv          nb202_fizO,             816
.equiv          nb202_fixH1,            832
.equiv          nb202_fiyH1,            848
.equiv          nb202_fizH1,            864
.equiv          nb202_fixH2,            880
.equiv          nb202_fiyH2,            896
.equiv          nb202_fizH2,            912
.equiv          nb202_fjxO,             928
.equiv          nb202_fjyO,             944
.equiv          nb202_fjzO,             960
.equiv          nb202_fjxH1,            976
.equiv          nb202_fjyH1,            992
.equiv          nb202_fjzH1,            1008
.equiv          nb202_fjxH2,            1024
.equiv          nb202_fjyH2,            1040
.equiv          nb202_fjzH2,            1056
.equiv          nb202_half,             1072
.equiv          nb202_three,            1088
.equiv          nb202_rsqOO,            1104
.equiv          nb202_rsqOH1,           1120
.equiv          nb202_rsqOH2,           1136
.equiv          nb202_rsqH1O,           1152
.equiv          nb202_rsqH1H1,          1168
.equiv          nb202_rsqH1H2,          1184
.equiv          nb202_rsqH2O,           1200
.equiv          nb202_rsqH2H1,          1216
.equiv          nb202_rsqH2H2,          1232
.equiv          nb202_rinvOO,           1248
.equiv          nb202_rinvOH1,          1264
.equiv          nb202_rinvOH2,          1280
.equiv          nb202_rinvH1O,          1296
.equiv          nb202_rinvH1H1,         1312
.equiv          nb202_rinvH1H2,         1328
.equiv          nb202_rinvH2O,          1344
.equiv          nb202_rinvH2H1,         1360
.equiv          nb202_rinvH2H2,         1376
.equiv          nb202_two,              1392
.equiv          nb202_krf,              1408
.equiv          nb202_crf,              1424
.equiv          nb202_is3,              1440
.equiv          nb202_ii3,              1444
.equiv          nb202_nri,              1448
.equiv          nb202_iinr,             1456
.equiv          nb202_jindex,           1464
.equiv          nb202_jjnr,             1472
.equiv          nb202_shift,            1480
.equiv          nb202_shiftvec,         1488
.equiv          nb202_facel,            1496
.equiv          nb202_innerjjnr,        1504
.equiv          nb202_innerk,           1512
.equiv          nb202_n,                1516
.equiv          nb202_nn1,              1520
.equiv          nb202_nouter,           1524
.equiv          nb202_ninner,           1528

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
	sub rsp, 1536		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb202_nouter], eax
	mov [rsp + nb202_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb202_nri], edi
	mov [rsp + nb202_iinr], rsi
	mov [rsp + nb202_jindex], rdx
	mov [rsp + nb202_jjnr], rcx
	mov [rsp + nb202_shift], r8
	mov [rsp + nb202_shiftvec], r9
	mov rsi, [rbp + nb202_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb202_facel], xmm0

	mov rsi, [rbp + nb202_argkrf]
	mov rdi, [rbp + nb202_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb202_krf], xmm1
	movapd [rsp + nb202_crf], xmm2


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb202_half], eax
	mov [rsp + nb202_half+4], ebx
	movsd xmm1, [rsp + nb202_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb202_half], xmm1
	movapd [rsp + nb202_two], xmm2
	movapd [rsp + nb202_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb202_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb202_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]

	movsd xmm6, [rsp + nb202_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb202_qqOO], xmm3
	movapd [rsp + nb202_qqOH], xmm4
	movapd [rsp + nb202_qqHH], xmm5
	
.nb202_threadloop:
        mov   rsi, [rbp + nb202_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb202_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb202_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb202_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb202_n], eax
        mov [rsp + nb202_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb202_outerstart
        jmp .nb202_end

.nb202_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb202_nouter]
	mov [rsp + nb202_nouter], ebx

.nb202_outer:
	mov   rax, [rsp + nb202_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb202_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb202_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb202_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb202_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb202_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb202_ixO], xmm3
	movapd [rsp + nb202_iyO], xmm4
	movapd [rsp + nb202_izO], xmm5

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
	movapd [rsp + nb202_ixH1], xmm0
	movapd [rsp + nb202_iyH1], xmm1
	movapd [rsp + nb202_izH1], xmm2
	movapd [rsp + nb202_ixH2], xmm3
	movapd [rsp + nb202_iyH2], xmm4
	movapd [rsp + nb202_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb202_vctot], xmm4
	movapd [rsp + nb202_fixO], xmm4
	movapd [rsp + nb202_fiyO], xmm4
	movapd [rsp + nb202_fizO], xmm4
	movapd [rsp + nb202_fixH1], xmm4
	movapd [rsp + nb202_fiyH1], xmm4
	movapd [rsp + nb202_fizH1], xmm4
	movapd [rsp + nb202_fixH2], xmm4
	movapd [rsp + nb202_fiyH2], xmm4
	movapd [rsp + nb202_fizH2], xmm4
	
	mov   rax, [rsp + nb202_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb202_pos]
	mov   rdi, [rbp + nb202_faction]	
	mov   rax, [rsp + nb202_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb202_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb202_ninner]
	mov   [rsp + nb202_ninner], ecx
	add   edx, 0
	mov   [rsp + nb202_innerk], edx    ;# number of innerloop atoms 
	jge   .nb202_unroll_loop
	jmp   .nb202_checksingle
.nb202_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb202_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb202_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb202_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j O coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb202_ixO]
    subpd xmm1, [rsp + nb202_iyO]
    subpd xmm2, [rsp + nb202_izO]
    subpd xmm3, [rsp + nb202_ixH1]
    subpd xmm4, [rsp + nb202_iyH1]
    subpd xmm5, [rsp + nb202_izH1]
    subpd xmm6, [rsp + nb202_ixH2]
    subpd xmm7, [rsp + nb202_iyH2]
    subpd xmm8, [rsp + nb202_izH2]
    
	movapd [rsp + nb202_dxOO], xmm0
	movapd [rsp + nb202_dyOO], xmm1
	movapd [rsp + nb202_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb202_dxH1O], xmm3
	movapd [rsp + nb202_dyH1O], xmm4
	movapd [rsp + nb202_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb202_dxH2O], xmm6
	movapd [rsp + nb202_dyH2O], xmm7
	movapd [rsp + nb202_dzH2O], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jO atoms
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
		
	movapd  xmm9, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOO 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1O
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2O	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;#  rinvOO 
	mulpd   xmm10, xmm15 ;#   rinvH1O
    mulpd   xmm11, xmm15 ;#   rinvH2O
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb202_krf]    
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
    movapd xmm14, [rsp + nb202_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb202_qqOO]
    movapd xmm13, [rsp + nb202_qqOH]    
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
    addpd  xmm2, [rsp + nb202_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd [rsp + nb202_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j O forces to xmm0-xmm2
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

	mulpd xmm7, [rsp + nb202_dxOO]
	mulpd xmm8, [rsp + nb202_dyOO]
	mulpd xmm9, [rsp + nb202_dzOO]
	mulpd xmm10, [rsp + nb202_dxH1O]
	mulpd xmm11, [rsp + nb202_dyH1O]
	mulpd xmm12, [rsp + nb202_dzH1O]
	mulpd xmm13, [rsp + nb202_dxH2O]
	mulpd xmm14, [rsp + nb202_dyH2O]
	mulpd xmm15, [rsp + nb202_dzH2O]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb202_fixO]
    addpd xmm8, [rsp + nb202_fiyO]
    addpd xmm9, [rsp + nb202_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb202_fixH1]
    addpd xmm11, [rsp + nb202_fiyH1]
    addpd xmm12, [rsp + nb202_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb202_fixH2]
    addpd xmm14, [rsp + nb202_fiyH2]
    addpd xmm15, [rsp + nb202_fizH2]

    movapd [rsp + nb202_fixO], xmm7
    movapd [rsp + nb202_fiyO], xmm8
    movapd [rsp + nb202_fizO], xmm9
    movapd [rsp + nb202_fixH1], xmm10
    movapd [rsp + nb202_fiyH1], xmm11
    movapd [rsp + nb202_fizH1], xmm12
    movapd [rsp + nb202_fixH2], xmm13
    movapd [rsp + nb202_fiyH2], xmm14
    movapd [rsp + nb202_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# move j H1 coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8 + 24] 
    movlpd xmm1, [rsi + rax*8 + 32] 
    movlpd xmm2, [rsi + rax*8 + 40] 
    movhpd xmm0, [rsi + rbx*8 + 24] 
    movhpd xmm1, [rsi + rbx*8 + 32] 
    movhpd xmm2, [rsi + rbx*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb202_ixO]
    subpd xmm1, [rsp + nb202_iyO]
    subpd xmm2, [rsp + nb202_izO]
    subpd xmm3, [rsp + nb202_ixH1]
    subpd xmm4, [rsp + nb202_iyH1]
    subpd xmm5, [rsp + nb202_izH1]
    subpd xmm6, [rsp + nb202_ixH2]
    subpd xmm7, [rsp + nb202_iyH2]
    subpd xmm8, [rsp + nb202_izH2]
    
	movapd [rsp + nb202_dxOH1], xmm0
	movapd [rsp + nb202_dyOH1], xmm1
	movapd [rsp + nb202_dzOH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb202_dxH1H1], xmm3
	movapd [rsp + nb202_dyH1H1], xmm4
	movapd [rsp + nb202_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb202_dxH2H1], xmm6
	movapd [rsp + nb202_dyH2H1], xmm7
	movapd [rsp + nb202_dzH2H1], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
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
		
	movapd  xmm9, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOH1 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1H1
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2OH1

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;#  rinvOH1
	mulpd   xmm10, xmm15 ;#   rinvH1H1
    mulpd   xmm11, xmm15 ;#   rinvH2H1
	
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb202_krf]    
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
    movapd xmm14, [rsp + nb202_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb202_qqOH]
    movapd xmm13, [rsp + nb202_qqHH]    
    mulpd  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, xmm13
    mulpd  xmm7, xmm13
    addpd  xmm2, [rsp + nb202_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd  [rsp + nb202_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j H1 forces to xmm0-xmm2
	movlpd xmm0, [rdi + rax*8 + 24]
	movlpd xmm1, [rdi + rax*8 + 32]
	movlpd xmm2, [rdi + rax*8 + 40]
	movhpd xmm0, [rdi + rbx*8 + 24]
	movhpd xmm1, [rdi + rbx*8 + 32]
	movhpd xmm2, [rdi + rbx*8 + 40]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb202_dxOH1]
	mulpd xmm8, [rsp + nb202_dyOH1]
	mulpd xmm9, [rsp + nb202_dzOH1]
	mulpd xmm10, [rsp + nb202_dxH1H1]
	mulpd xmm11, [rsp + nb202_dyH1H1]
	mulpd xmm12, [rsp + nb202_dzH1H1]
	mulpd xmm13, [rsp + nb202_dxH2H1]
	mulpd xmm14, [rsp + nb202_dyH2H1]
	mulpd xmm15, [rsp + nb202_dzH2H1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb202_fixO]
    addpd xmm8, [rsp + nb202_fiyO]
    addpd xmm9, [rsp + nb202_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb202_fixH1]
    addpd xmm11, [rsp + nb202_fiyH1]
    addpd xmm12, [rsp + nb202_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb202_fixH2]
    addpd xmm14, [rsp + nb202_fiyH2]
    addpd xmm15, [rsp + nb202_fizH2]

    movapd [rsp + nb202_fixO], xmm7
    movapd [rsp + nb202_fiyO], xmm8
    movapd [rsp + nb202_fizO], xmm9
    movapd [rsp + nb202_fixH1], xmm10
    movapd [rsp + nb202_fiyH1], xmm11
    movapd [rsp + nb202_fizH1], xmm12
    movapd [rsp + nb202_fixH2], xmm13
    movapd [rsp + nb202_fiyH2], xmm14
    movapd [rsp + nb202_fizH2], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8 + 48] 
    movlpd xmm1, [rsi + rax*8 + 56] 
    movlpd xmm2, [rsi + rax*8 + 64] 
    movhpd xmm0, [rsi + rbx*8 + 48] 
    movhpd xmm1, [rsi + rbx*8 + 56] 
    movhpd xmm2, [rsi + rbx*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb202_ixO]
    subpd xmm1, [rsp + nb202_iyO]
    subpd xmm2, [rsp + nb202_izO]
    subpd xmm3, [rsp + nb202_ixH1]
    subpd xmm4, [rsp + nb202_iyH1]
    subpd xmm5, [rsp + nb202_izH1]
    subpd xmm6, [rsp + nb202_ixH2]
    subpd xmm7, [rsp + nb202_iyH2]
    subpd xmm8, [rsp + nb202_izH2]
    
	movapd [rsp + nb202_dxOH2], xmm0
	movapd [rsp + nb202_dyOH2], xmm1
	movapd [rsp + nb202_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb202_dxH1H2], xmm3
	movapd [rsp + nb202_dyH1H2], xmm4
	movapd [rsp + nb202_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb202_dxH2H2], xmm6
	movapd [rsp + nb202_dyH2H2], xmm7
	movapd [rsp + nb202_dzH2H2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
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
		
	movapd  xmm9, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOH2 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1H2
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2H2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb202_three]
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

	movapd  xmm15, [rsp + nb202_half]
	mulpd   xmm9, xmm15  ;#  rinvOH2
	mulpd   xmm10, xmm15 ;#   rinvH1H2
    mulpd   xmm11, xmm15 ;#   rinvH2H2
	
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb202_krf]    
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
    movapd xmm14, [rsp + nb202_crf]
    subpd  xmm2, xmm14   ;# rinv+krsq-crf
    subpd  xmm5, xmm14
    subpd  xmm8, xmm14
    movapd xmm12, [rsp + nb202_qqOH]
    movapd xmm13, [rsp + nb202_qqHH]    
    mulpd  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, xmm13
    mulpd  xmm7, xmm13
    addpd  xmm2, [rsp + nb202_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd  [rsp + nb202_vctot], xmm2
    
    mulpd  xmm9, xmm1   ;# fscal
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7

    ;# move j H2 forces to xmm0-xmm2
	movlpd xmm0, [rdi + rax*8 + 48]
	movlpd xmm1, [rdi + rax*8 + 56]
	movlpd xmm2, [rdi + rax*8 + 64]
	movhpd xmm0, [rdi + rbx*8 + 48]
	movhpd xmm1, [rdi + rbx*8 + 56]
	movhpd xmm2, [rdi + rbx*8 + 64]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb202_dxOH2]
	mulpd xmm8, [rsp + nb202_dyOH2]
	mulpd xmm9, [rsp + nb202_dzOH2]
	mulpd xmm10, [rsp + nb202_dxH1H2]
	mulpd xmm11, [rsp + nb202_dyH1H2]
	mulpd xmm12, [rsp + nb202_dzH1H2]
	mulpd xmm13, [rsp + nb202_dxH2H2]
	mulpd xmm14, [rsp + nb202_dyH2H2]
	mulpd xmm15, [rsp + nb202_dzH2H2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb202_fixO]
    addpd xmm8, [rsp + nb202_fiyO]
    addpd xmm9, [rsp + nb202_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb202_fixH1]
    addpd xmm11, [rsp + nb202_fiyH1]
    addpd xmm12, [rsp + nb202_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb202_fixH2]
    addpd xmm14, [rsp + nb202_fiyH2]
    addpd xmm15, [rsp + nb202_fizH2]

    movapd [rsp + nb202_fixO], xmm7
    movapd [rsp + nb202_fiyO], xmm8
    movapd [rsp + nb202_fizO], xmm9
    movapd [rsp + nb202_fixH1], xmm10
    movapd [rsp + nb202_fiyH1], xmm11
    movapd [rsp + nb202_fizH1], xmm12
    movapd [rsp + nb202_fixH2], xmm13
    movapd [rsp + nb202_fiyH2], xmm14
    movapd [rsp + nb202_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb202_innerk],  2
	jl    .nb202_checksingle
	jmp   .nb202_unroll_loop
.nb202_checksingle:
	mov   edx, [rsp + nb202_innerk]
	and   edx, 1
	jnz   .nb202_dosingle
	jmp   .nb202_updateouterdata
.nb202_dosingle:
	mov   rdx, [rsp + nb202_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	
	mov rsi, [rbp + nb202_pos]
	lea   rax, [rax + rax*2]  

	
	;# move j O coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8] 
    movsd xmm1, [rsi + rax*8 + 8] 
    movsd xmm2, [rsi + rax*8 + 16] 

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb202_ixO]
    subsd xmm1, [rsp + nb202_iyO]
    subsd xmm2, [rsp + nb202_izO]
    subsd xmm3, [rsp + nb202_ixH1]
    subsd xmm4, [rsp + nb202_iyH1]
    subsd xmm5, [rsp + nb202_izH1]
    subsd xmm6, [rsp + nb202_ixH2]
    subsd xmm7, [rsp + nb202_iyH2]
    subsd xmm8, [rsp + nb202_izH2]
    
	movsd [rsp + nb202_dxOO], xmm0
	movsd [rsp + nb202_dyOO], xmm1
	movsd [rsp + nb202_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb202_dxH1O], xmm3
	movsd [rsp + nb202_dyH1O], xmm4
	movsd [rsp + nb202_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb202_dxH2O], xmm6
	movsd [rsp + nb202_dyH2O], xmm7
	movsd [rsp + nb202_dzH2O], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jO atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb202_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOO 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1O
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2O	

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb202_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;#  rinvOO 
	mulsd   xmm10, xmm15 ;#   rinvH1O
    mulsd   xmm11, xmm15 ;#   rinvH2O
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd xmm1, xmm9 ;# copy of rinv
    movsd xmm4, xmm10
    movsd xmm7, xmm11
    movsd xmm2, [rsp + nb202_krf]    
    mulsd  xmm9, xmm9   ;# rinvsq
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, xmm2  ;# k*rsq
    mulsd  xmm3, xmm2
    mulsd  xmm6, xmm2
    movsd  xmm2, xmm0 ;# copy of k*rsq
    movsd  xmm5, xmm3
    movsd  xmm8, xmm6
    addsd  xmm2, xmm1  ;# rinv+krsq
    addsd  xmm5, xmm4
    addsd  xmm8, xmm7
    movsd  xmm14, [rsp + nb202_crf]
    subsd  xmm2, xmm14   ;# rinv+krsq-crf
    subsd  xmm5, xmm14
    subsd  xmm8, xmm14
    movsd  xmm12, [rsp + nb202_qqOO]
    movsd  xmm13, [rsp + nb202_qqOH]    
    mulsd  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulsd  xmm5, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    mulsd  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addsd  xmm0, xmm0 ;# 2*krsq
    addsd  xmm3, xmm3 
    addsd  xmm6, xmm6 
    subsd  xmm1, xmm0 ;# rinv-2*krsq
    subsd  xmm4, xmm3
    subsd  xmm7, xmm6
    mulsd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulsd  xmm4, xmm13
    mulsd  xmm7, xmm13
    addsd  xmm2, [rsp + nb202_vctot]
    addsd  xmm5, xmm8
    addsd  xmm2, xmm5
    movsd  [rsp + nb202_vctot], xmm2
    
    mulsd  xmm9, xmm1   ;# fscal
    mulsd  xmm10, xmm4
    mulsd  xmm11, xmm7

    ;# move j O forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8]
	movsd xmm1, [rdi + rax*8 + 8]
	movsd xmm2, [rdi + rax*8 + 16]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb202_dxOO]
	mulsd xmm8, [rsp + nb202_dyOO]
	mulsd xmm9, [rsp + nb202_dzOO]
	mulsd xmm10, [rsp + nb202_dxH1O]
	mulsd xmm11, [rsp + nb202_dyH1O]
	mulsd xmm12, [rsp + nb202_dzH1O]
	mulsd xmm13, [rsp + nb202_dxH2O]
	mulsd xmm14, [rsp + nb202_dyH2O]
	mulsd xmm15, [rsp + nb202_dzH2O]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb202_fixO]
    addsd xmm8, [rsp + nb202_fiyO]
    addsd xmm9, [rsp + nb202_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb202_fixH1]
    addsd xmm11, [rsp + nb202_fiyH1]
    addsd xmm12, [rsp + nb202_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb202_fixH2]
    addsd xmm14, [rsp + nb202_fiyH2]
    addsd xmm15, [rsp + nb202_fizH2]

    movsd [rsp + nb202_fixO], xmm7
    movsd [rsp + nb202_fiyO], xmm8
    movsd [rsp + nb202_fizO], xmm9
    movsd [rsp + nb202_fixH1], xmm10
    movsd [rsp + nb202_fiyH1], xmm11
    movsd [rsp + nb202_fizH1], xmm12
    movsd [rsp + nb202_fixH2], xmm13
    movsd [rsp + nb202_fiyH2], xmm14
    movsd [rsp + nb202_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	movsd [rdi + rax*8],      xmm0
	movsd [rdi + rax*8 + 8],  xmm1
	movsd [rdi + rax*8 + 16], xmm2

	;# move j H1 coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8 + 24] 
    movsd xmm1, [rsi + rax*8 + 32] 
    movsd xmm2, [rsi + rax*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb202_ixO]
    subsd xmm1, [rsp + nb202_iyO]
    subsd xmm2, [rsp + nb202_izO]
    subsd xmm3, [rsp + nb202_ixH1]
    subsd xmm4, [rsp + nb202_iyH1]
    subsd xmm5, [rsp + nb202_izH1]
    subsd xmm6, [rsp + nb202_ixH2]
    subsd xmm7, [rsp + nb202_iyH2]
    subsd xmm8, [rsp + nb202_izH2]
    
	movsd [rsp + nb202_dxOH1], xmm0
	movsd [rsp + nb202_dyOH1], xmm1
	movsd [rsp + nb202_dzOH1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb202_dxH1H1], xmm3
	movsd [rsp + nb202_dyH1H1], xmm4
	movsd [rsp + nb202_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb202_dxH2H1], xmm6
	movsd [rsp + nb202_dyH2H1], xmm7
	movsd [rsp + nb202_dzH2H1], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb202_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2OH1

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb202_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;#  rinvOH1
	mulsd   xmm10, xmm15 ;#   rinvH1H1
    mulsd   xmm11, xmm15 ;#   rinvH2H1
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd xmm1, xmm9 ;# copy of rinv
    movsd xmm4, xmm10
    movsd xmm7, xmm11
    movsd xmm2, [rsp + nb202_krf]    
    mulsd  xmm9, xmm9   ;# rinvsq
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, xmm2  ;# k*rsq
    mulsd  xmm3, xmm2
    mulsd  xmm6, xmm2
    movsd xmm2, xmm0 ;# copy of k*rsq
    movsd xmm5, xmm3
    movsd xmm8, xmm6
    addsd  xmm2, xmm1  ;# rinv+krsq
    addsd  xmm5, xmm4
    addsd  xmm8, xmm7
    movsd xmm14, [rsp + nb202_crf]
    subsd  xmm2, xmm14   ;# rinv+krsq-crf
    subsd  xmm5, xmm14
    subsd  xmm8, xmm14
    movsd xmm12, [rsp + nb202_qqOH]
    movsd xmm13, [rsp + nb202_qqHH]    
    mulsd  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  xmm0, xmm0 ;# 2*krsq
    addsd  xmm3, xmm3 
    addsd  xmm6, xmm6 
    subsd  xmm1, xmm0 ;# rinv-2*krsq
    subsd  xmm4, xmm3
    subsd  xmm7, xmm6
    mulsd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulsd  xmm4, xmm13
    mulsd  xmm7, xmm13
    addsd  xmm2, [rsp + nb202_vctot]
    addsd  xmm5, xmm8
    addsd  xmm2, xmm5
    movsd  [rsp + nb202_vctot], xmm2
    
    mulsd  xmm9, xmm1   ;# fscal
    mulsd  xmm10, xmm4
    mulsd  xmm11, xmm7

    ;# move j H1 forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 24]
	movsd xmm1, [rdi + rax*8 + 32]
	movsd xmm2, [rdi + rax*8 + 40]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb202_dxOH1]
	mulsd xmm8, [rsp + nb202_dyOH1]
	mulsd xmm9, [rsp + nb202_dzOH1]
	mulsd xmm10, [rsp + nb202_dxH1H1]
	mulsd xmm11, [rsp + nb202_dyH1H1]
	mulsd xmm12, [rsp + nb202_dzH1H1]
	mulsd xmm13, [rsp + nb202_dxH2H1]
	mulsd xmm14, [rsp + nb202_dyH2H1]
	mulsd xmm15, [rsp + nb202_dzH2H1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb202_fixO]
    addsd xmm8, [rsp + nb202_fiyO]
    addsd xmm9, [rsp + nb202_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb202_fixH1]
    addsd xmm11, [rsp + nb202_fiyH1]
    addsd xmm12, [rsp + nb202_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb202_fixH2]
    addsd xmm14, [rsp + nb202_fiyH2]
    addsd xmm15, [rsp + nb202_fizH2]

    movsd [rsp + nb202_fixO], xmm7
    movsd [rsp + nb202_fiyO], xmm8
    movsd [rsp + nb202_fizO], xmm9
    movsd [rsp + nb202_fixH1], xmm10
    movsd [rsp + nb202_fiyH1], xmm11
    movsd [rsp + nb202_fizH1], xmm12
    movsd [rsp + nb202_fixH2], xmm13
    movsd [rsp + nb202_fiyH2], xmm14
    movsd [rsp + nb202_fizH2], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8 + 48] 
    movsd xmm1, [rsi + rax*8 + 56] 
    movsd xmm2, [rsi + rax*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb202_ixO]
    subsd xmm1, [rsp + nb202_iyO]
    subsd xmm2, [rsp + nb202_izO]
    subsd xmm3, [rsp + nb202_ixH1]
    subsd xmm4, [rsp + nb202_iyH1]
    subsd xmm5, [rsp + nb202_izH1]
    subsd xmm6, [rsp + nb202_ixH2]
    subsd xmm7, [rsp + nb202_iyH2]
    subsd xmm8, [rsp + nb202_izH2]
    
	movsd [rsp + nb202_dxOH2], xmm0
	movsd [rsp + nb202_dyOH2], xmm1
	movsd [rsp + nb202_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb202_dxH1H2], xmm3
	movsd [rsp + nb202_dyH1H2], xmm4
	movsd [rsp + nb202_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb202_dxH2H2], xmm6
	movsd [rsp + nb202_dyH2H2], xmm7
	movsd [rsp + nb202_dzH2H2], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb202_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2H2

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb202_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb202_half]
	mulsd   xmm9, xmm15  ;#  rinvOH2
	mulsd   xmm10, xmm15 ;#   rinvH1H2
    mulsd   xmm11, xmm15 ;#   rinvH2H2
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd xmm1, xmm9 ;# copy of rinv
    movsd xmm4, xmm10
    movsd xmm7, xmm11
    movsd xmm2, [rsp + nb202_krf]    
    mulsd  xmm9, xmm9   ;# rinvsq
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, xmm2  ;# k*rsq
    mulsd  xmm3, xmm2
    mulsd  xmm6, xmm2
    movsd xmm2, xmm0 ;# copy of k*rsq
    movsd xmm5, xmm3
    movsd xmm8, xmm6
    addsd  xmm2, xmm1  ;# rinv+krsq
    addsd  xmm5, xmm4
    addsd  xmm8, xmm7
    movsd xmm14, [rsp + nb202_crf]
    subsd  xmm2, xmm14   ;# rinv+krsq-crf
    subsd  xmm5, xmm14
    subsd  xmm8, xmm14
    movsd xmm12, [rsp + nb202_qqOH]
    movsd xmm13, [rsp + nb202_qqHH]    
    mulsd  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  xmm5, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  xmm0, xmm0 ;# 2*krsq
    addsd  xmm3, xmm3 
    addsd  xmm6, xmm6 
    subsd  xmm1, xmm0 ;# rinv-2*krsq
    subsd  xmm4, xmm3
    subsd  xmm7, xmm6
    mulsd  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulsd  xmm4, xmm13
    mulsd  xmm7, xmm13
    addsd  xmm2, [rsp + nb202_vctot]
    addsd  xmm5, xmm8
    addsd  xmm2, xmm5
    movsd  [rsp + nb202_vctot], xmm2
    
    mulsd  xmm9, xmm1   ;# fscal
    mulsd  xmm10, xmm4
    mulsd  xmm11, xmm7

    ;# move j H2 forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 48]
	movsd xmm1, [rdi + rax*8 + 56]
	movsd xmm2, [rdi + rax*8 + 64]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb202_dxOH2]
	mulsd xmm8, [rsp + nb202_dyOH2]
	mulsd xmm9, [rsp + nb202_dzOH2]
	mulsd xmm10, [rsp + nb202_dxH1H2]
	mulsd xmm11, [rsp + nb202_dyH1H2]
	mulsd xmm12, [rsp + nb202_dzH1H2]
	mulsd xmm13, [rsp + nb202_dxH2H2]
	mulsd xmm14, [rsp + nb202_dyH2H2]
	mulsd xmm15, [rsp + nb202_dzH2H2]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb202_fixO]
    addsd xmm8, [rsp + nb202_fiyO]
    addsd xmm9, [rsp + nb202_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb202_fixH1]
    addsd xmm11, [rsp + nb202_fiyH1]
    addsd xmm12, [rsp + nb202_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb202_fixH2]
    addsd xmm14, [rsp + nb202_fiyH2]
    addsd xmm15, [rsp + nb202_fizH2]

    movsd [rsp + nb202_fixO], xmm7
    movsd [rsp + nb202_fiyO], xmm8
    movsd [rsp + nb202_fizO], xmm9
    movsd [rsp + nb202_fixH1], xmm10
    movsd [rsp + nb202_fiyH1], xmm11
    movsd [rsp + nb202_fizH1], xmm12
    movsd [rsp + nb202_fixH2], xmm13
    movsd [rsp + nb202_fiyH2], xmm14
    movsd [rsp + nb202_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
	
.nb202_updateouterdata:
	mov   ecx, [rsp + nb202_ii3]
	mov   rdi, [rbp + nb202_faction]
	mov   rsi, [rbp + nb202_fshift]
	mov   edx, [rsp + nb202_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb202_fixO]
	movapd xmm1, [rsp + nb202_fiyO]
	movapd xmm2, [rsp + nb202_fizO]

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
	movapd xmm0, [rsp + nb202_fixH1]
	movapd xmm1, [rsp + nb202_fiyH1]
	movapd xmm2, [rsp + nb202_fizH1]

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
	movapd xmm0, [rsp + nb202_fixH2]
	movapd xmm1, [rsp + nb202_fiyH2]
	movapd xmm2, [rsp + nb202_fizH2]

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
	mov esi, [rsp + nb202_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb202_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb202_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb202_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb202_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb202_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb202_n], esi
        jmp .nb202_outer
.nb202_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb202_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb202_end
        ;# non-zero, do one more workunit
        jmp   .nb202_threadloop
.nb202_end:
	mov eax, [rsp + nb202_nouter]
	mov ebx, [rsp + nb202_ninner]
	mov rcx, [rbp + nb202_outeriter]
	mov rdx, [rbp + nb202_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1536
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


	

	
.globl nb_kernel202nf_x86_64_sse2
.globl _nb_kernel202nf_x86_64_sse2
nb_kernel202nf_x86_64_sse2:	
_nb_kernel202nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb202nf_fshift,         16
.equiv          nb202nf_gid,            24
.equiv          nb202nf_pos,            32
.equiv          nb202nf_faction,        40
.equiv          nb202nf_charge,         48
.equiv          nb202nf_p_facel,        56
.equiv          nb202nf_argkrf,         64
.equiv          nb202nf_argcrf,         72
.equiv          nb202nf_Vc,             80
.equiv          nb202nf_type,           88
.equiv          nb202nf_p_ntype,        96
.equiv          nb202nf_vdwparam,       104
.equiv          nb202nf_Vvdw,           112
.equiv          nb202nf_p_tabscale,     120
.equiv          nb202nf_VFtab,          128
.equiv          nb202nf_invsqrta,       136
.equiv          nb202nf_dvda,           144
.equiv          nb202nf_p_gbtabscale,   152
.equiv          nb202nf_GBtab,          160
.equiv          nb202nf_p_nthreads,     168
.equiv          nb202nf_count,          176
.equiv          nb202nf_mtx,            184
.equiv          nb202nf_outeriter,      192
.equiv          nb202nf_inneriter,      200
.equiv          nb202nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb202nf_ixO,            0
.equiv          nb202nf_iyO,            16
.equiv          nb202nf_izO,            32
.equiv          nb202nf_ixH1,           48
.equiv          nb202nf_iyH1,           64
.equiv          nb202nf_izH1,           80
.equiv          nb202nf_ixH2,           96
.equiv          nb202nf_iyH2,           112
.equiv          nb202nf_izH2,           128
.equiv          nb202nf_jxO,            144
.equiv          nb202nf_jyO,            160
.equiv          nb202nf_jzO,            176
.equiv          nb202nf_jxH1,           192
.equiv          nb202nf_jyH1,           208
.equiv          nb202nf_jzH1,           224
.equiv          nb202nf_jxH2,           240
.equiv          nb202nf_jyH2,           256
.equiv          nb202nf_jzH2,           272
.equiv          nb202nf_qqOO,           288
.equiv          nb202nf_qqOH,           304
.equiv          nb202nf_qqHH,           320
.equiv          nb202nf_vctot,          336
.equiv          nb202nf_half,           352
.equiv          nb202nf_three,          368
.equiv          nb202nf_rsqOO,          384
.equiv          nb202nf_rsqOH1,         400
.equiv          nb202nf_rsqOH2,         416
.equiv          nb202nf_rsqH1O,         432
.equiv          nb202nf_rsqH1H1,        448
.equiv          nb202nf_rsqH1H2,        464
.equiv          nb202nf_rsqH2O,         480
.equiv          nb202nf_rsqH2H1,        496
.equiv          nb202nf_rsqH2H2,        512
.equiv          nb202nf_rinvOO,         528
.equiv          nb202nf_rinvOH1,        544
.equiv          nb202nf_rinvOH2,        560
.equiv          nb202nf_rinvH1O,        576
.equiv          nb202nf_rinvH1H1,       592
.equiv          nb202nf_rinvH1H2,       608
.equiv          nb202nf_rinvH2O,        624
.equiv          nb202nf_rinvH2H1,       640
.equiv          nb202nf_rinvH2H2,       656
.equiv          nb202nf_krf,            672
.equiv          nb202nf_crf,            688
.equiv          nb202nf_is3,            704
.equiv          nb202nf_ii3,            708
.equiv          nb202nf_nri,            712
.equiv          nb202nf_iinr,           720
.equiv          nb202nf_jindex,         728
.equiv          nb202nf_jjnr,           736
.equiv          nb202nf_shift,          744
.equiv          nb202nf_shiftvec,       752
.equiv          nb202nf_facel,          760
.equiv          nb202nf_innerjjnr,      768
.equiv          nb202nf_innerk,         776
.equiv          nb202nf_n,              780
.equiv          nb202nf_nn1,            784
.equiv          nb202nf_nouter,         788
.equiv          nb202nf_ninner,         792

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
	mov [rsp + nb202nf_nouter], eax
	mov [rsp + nb202nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb202nf_nri], edi
	mov [rsp + nb202nf_iinr], rsi
	mov [rsp + nb202nf_jindex], rdx
	mov [rsp + nb202nf_jjnr], rcx
	mov [rsp + nb202nf_shift], r8
	mov [rsp + nb202nf_shiftvec], r9
	mov rsi, [rbp + nb202nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb202nf_facel], xmm0

	mov rsi, [rbp + nb202nf_argkrf]
	mov rdi, [rbp + nb202nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb202nf_krf], xmm1
	movapd [rsp + nb202nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb202nf_half], eax
	mov [rsp + nb202nf_half+4], ebx
	movsd xmm1, [rsp + nb202nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb202nf_half], xmm1
	movapd [rsp + nb202nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb202nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb202nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]
	movsd xmm6, [rsp + nb202nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb202nf_qqOO], xmm3
	movapd [rsp + nb202nf_qqOH], xmm4
	movapd [rsp + nb202nf_qqHH], xmm5
	
.nb202nf_threadloop:
        mov   rsi, [rbp + nb202nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb202nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb202nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb202nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb202nf_n], eax
        mov [rsp + nb202nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb202nf_outerstart
        jmp .nb202nf_end

.nb202nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb202nf_nouter]
	mov [rsp + nb202nf_nouter], ebx

.nb202nf_outer:
	mov   rax, [rsp + nb202nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb202nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb202nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb202nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb202nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb202nf_ixO], xmm3
	movapd [rsp + nb202nf_iyO], xmm4
	movapd [rsp + nb202nf_izO], xmm5

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
	movapd [rsp + nb202nf_ixH1], xmm0
	movapd [rsp + nb202nf_iyH1], xmm1
	movapd [rsp + nb202nf_izH1], xmm2
	movapd [rsp + nb202nf_ixH2], xmm3
	movapd [rsp + nb202nf_iyH2], xmm4
	movapd [rsp + nb202nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb202nf_vctot], xmm4
	
	mov   rax, [rsp + nb202nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb202nf_pos]
	mov   rax, [rsp + nb202nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb202nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb202nf_ninner]
	mov   [rsp + nb202nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb202nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb202nf_unroll_loop
	jmp   .nb202nf_checksingle
.nb202nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb202nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb202nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb202nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movhpd xmm2, [rsi + rbx*8]
	movhpd xmm3, [rsi + rbx*8 + 8]
	movhpd xmm4, [rsi + rbx*8 + 16]
	movhpd xmm5, [rsi + rbx*8 + 24]
	movhpd xmm6, [rsi + rbx*8 + 32]
	movhpd xmm7, [rsi + rbx*8 + 40]
	movapd 	[rsp + nb202nf_jxO], xmm2
	movapd 	[rsp + nb202nf_jyO], xmm3
	movapd 	[rsp + nb202nf_jzO], xmm4
	movapd 	[rsp + nb202nf_jxH1], xmm5
	movapd 	[rsp + nb202nf_jyH1], xmm6
	movapd 	[rsp + nb202nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm2, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rbx*8 + 56]
	movhpd xmm4, [rsi + rbx*8 + 64]
	movapd 	[rsp + nb202nf_jxH2], xmm2
	movapd 	[rsp + nb202nf_jyH2], xmm3
	movapd 	[rsp + nb202nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb202nf_ixO]
	movapd xmm1, [rsp + nb202nf_iyO]
	movapd xmm2, [rsp + nb202nf_izO]
	movapd xmm3, [rsp + nb202nf_ixO]
	movapd xmm4, [rsp + nb202nf_iyO]
	movapd xmm5, [rsp + nb202nf_izO]
	subpd  xmm0, [rsp + nb202nf_jxO]
	subpd  xmm1, [rsp + nb202nf_jyO]
	subpd  xmm2, [rsp + nb202nf_jzO]
	subpd  xmm3, [rsp + nb202nf_jxH1]
	subpd  xmm4, [rsp + nb202nf_jyH1]
	subpd  xmm5, [rsp + nb202nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqOO], xmm0
	movapd [rsp + nb202nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb202nf_ixO]
	movapd xmm1, [rsp + nb202nf_iyO]
	movapd xmm2, [rsp + nb202nf_izO]
	movapd xmm3, [rsp + nb202nf_ixH1]
	movapd xmm4, [rsp + nb202nf_iyH1]
	movapd xmm5, [rsp + nb202nf_izH1]
	subpd  xmm0, [rsp + nb202nf_jxH2]
	subpd  xmm1, [rsp + nb202nf_jyH2]
	subpd  xmm2, [rsp + nb202nf_jzH2]
	subpd  xmm3, [rsp + nb202nf_jxO]
	subpd  xmm4, [rsp + nb202nf_jyO]
	subpd  xmm5, [rsp + nb202nf_jzO]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqOH2], xmm0
	movapd [rsp + nb202nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb202nf_ixH1]
	movapd xmm1, [rsp + nb202nf_iyH1]
	movapd xmm2, [rsp + nb202nf_izH1]
	movapd xmm3, [rsp + nb202nf_ixH1]
	movapd xmm4, [rsp + nb202nf_iyH1]
	movapd xmm5, [rsp + nb202nf_izH1]
	subpd  xmm0, [rsp + nb202nf_jxH1]
	subpd  xmm1, [rsp + nb202nf_jyH1]
	subpd  xmm2, [rsp + nb202nf_jzH1]
	subpd  xmm3, [rsp + nb202nf_jxH2]
	subpd  xmm4, [rsp + nb202nf_jyH2]
	subpd  xmm5, [rsp + nb202nf_jzH2]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqH1H1], xmm0
	movapd [rsp + nb202nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb202nf_ixH2]
	movapd xmm1, [rsp + nb202nf_iyH2]
	movapd xmm2, [rsp + nb202nf_izH2]
	movapd xmm3, [rsp + nb202nf_ixH2]
	movapd xmm4, [rsp + nb202nf_iyH2]
	movapd xmm5, [rsp + nb202nf_izH2]
	subpd  xmm0, [rsp + nb202nf_jxO]
	subpd  xmm1, [rsp + nb202nf_jyO]
	subpd  xmm2, [rsp + nb202nf_jzO]
	subpd  xmm3, [rsp + nb202nf_jxH1]
	subpd  xmm4, [rsp + nb202nf_jyH1]
	subpd  xmm5, [rsp + nb202nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [rsp + nb202nf_rsqH2O], xmm0
	movapd [rsp + nb202nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb202nf_ixH2]
	movapd xmm1, [rsp + nb202nf_iyH2]
	movapd xmm2, [rsp + nb202nf_izH2]
	subpd  xmm0, [rsp + nb202nf_jxH2]
	subpd  xmm1, [rsp + nb202nf_jyH2]
	subpd  xmm2, [rsp + nb202nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [rsp + nb202nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb202nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH2H2], xmm1
	movapd [rsp + nb202nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb202nf_rsqOO]
	movapd xmm4, [rsp + nb202nf_rsqOH1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb202nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb202nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb202nf_half] ;# rinv
	movapd [rsp + nb202nf_rinvOO], xmm1
	movapd [rsp + nb202nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb202nf_rsqOH2]
	movapd xmm4, [rsp + nb202nf_rsqH1O]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb202nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvOH2], xmm1
	movapd [rsp + nb202nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb202nf_rsqH1H1]
	movapd xmm4, [rsp + nb202nf_rsqH1H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb202nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb202nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH1H1], xmm1
	movapd [rsp + nb202nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb202nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb202nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb202nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm6, [rsp + nb202nf_krf]
	mulpd  xmm6, [rsp + nb202nf_rsqOO]	;# xmm5=krsq 
	addpd  xmm6, [rsp + nb202nf_rinvOO]	;# xmm6=rinv+ krsq 
	subpd  xmm6, [rsp + nb202nf_crf]
	
	mulpd  xmm6, [rsp + nb202nf_qqOO] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, [rsp + nb202nf_vctot] ;# local vctot summation variable 

	;# O-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulpd  xmm5, [rsp + nb202nf_rsqOH1]	;# xmm5=krsq 
	addpd  xmm5, [rsp + nb202nf_rinvOH1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [rsp + nb202nf_crf]
	
	mulpd  xmm5, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# O-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulpd  xmm7, [rsp + nb202nf_rsqOH2]	;# xmm5=krsq 
	addpd  xmm7, [rsp + nb202nf_rinvOH2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [rsp + nb202nf_crf]
	
	mulpd  xmm7, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 

	;# H1-O interaction 
	movapd xmm4, [rsp + nb202nf_krf]
	mulpd  xmm4, [rsp + nb202nf_rsqH1O]	;# xmm5=krsq 
	addpd  xmm4, [rsp + nb202nf_rinvH1O]	;# xmm6=rinv+ krsq 
	subpd  xmm4, [rsp + nb202nf_crf]
	
	mulpd  xmm4, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H1-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulpd  xmm5, [rsp + nb202nf_rsqH1H1]	;# xmm5=krsq 
	addpd  xmm5, [rsp + nb202nf_rinvH1H1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [rsp + nb202nf_crf]
	
	mulpd  xmm5, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# H1-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulpd  xmm7, [rsp + nb202nf_rsqH1H2]	;# xmm5=krsq 
	addpd  xmm7, [rsp + nb202nf_rinvH1H2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [rsp + nb202nf_crf]
	
	mulpd  xmm7, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 

	;# H2-O interaction 
	movapd xmm4, [rsp + nb202nf_krf]
	mulpd  xmm4, [rsp + nb202nf_rsqH2O]	;# xmm5=krsq 
	addpd  xmm4, [rsp + nb202nf_rinvH2O]	;# xmm6=rinv+ krsq 
	subpd  xmm4, [rsp + nb202nf_crf]
	
	mulpd  xmm4, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H2-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulpd  xmm5, [rsp + nb202nf_rsqH2H1]	;# xmm5=krsq 
	addpd  xmm5, [rsp + nb202nf_rinvH2H1]	;# xmm6=rinv+ krsq 
	subpd  xmm5, [rsp + nb202nf_crf]
	
	mulpd  xmm5, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm5 ;# local vctot summation variable 

	;# H2-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulpd  xmm7, [rsp + nb202nf_rsqH2H2]	;# xmm5=krsq 
	addpd  xmm7, [rsp + nb202nf_rinvH2H2]	;# xmm6=rinv+ krsq 
	subpd  xmm7, [rsp + nb202nf_crf]
	
	mulpd  xmm7, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addpd  xmm6, xmm7 ;# local vctot summation variable 
	movapd [rsp + nb202nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb202nf_innerk],  2
	jl    .nb202nf_checksingle
	jmp   .nb202nf_unroll_loop
.nb202nf_checksingle:
	mov   edx, [rsp + nb202nf_innerk]
	and   edx, 1
	jnz   .nb202nf_dosingle
	jmp   .nb202nf_updateouterdata
.nb202nf_dosingle:
	mov   rdx, [rsp + nb202nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	
	mov rsi, [rbp + nb202nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movapd 	[rsp + nb202nf_jxO], xmm2
	movapd 	[rsp + nb202nf_jyO], xmm3
	movapd 	[rsp + nb202nf_jzO], xmm4
	movapd 	[rsp + nb202nf_jxH1], xmm5
	movapd 	[rsp + nb202nf_jyH1], xmm6
	movapd 	[rsp + nb202nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movapd 	[rsp + nb202nf_jxH2], xmm2
	movapd 	[rsp + nb202nf_jyH2], xmm3
	movapd 	[rsp + nb202nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb202nf_ixO]
	movapd xmm1, [rsp + nb202nf_iyO]
	movapd xmm2, [rsp + nb202nf_izO]
	movapd xmm3, [rsp + nb202nf_ixO]
	movapd xmm4, [rsp + nb202nf_iyO]
	movapd xmm5, [rsp + nb202nf_izO]
	subsd  xmm0, [rsp + nb202nf_jxO]
	subsd  xmm1, [rsp + nb202nf_jyO]
	subsd  xmm2, [rsp + nb202nf_jzO]
	subsd  xmm3, [rsp + nb202nf_jxH1]
	subsd  xmm4, [rsp + nb202nf_jyH1]
	subsd  xmm5, [rsp + nb202nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqOO], xmm0
	movapd [rsp + nb202nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb202nf_ixO]
	movapd xmm1, [rsp + nb202nf_iyO]
	movapd xmm2, [rsp + nb202nf_izO]
	movapd xmm3, [rsp + nb202nf_ixH1]
	movapd xmm4, [rsp + nb202nf_iyH1]
	movapd xmm5, [rsp + nb202nf_izH1]
	subsd  xmm0, [rsp + nb202nf_jxH2]
	subsd  xmm1, [rsp + nb202nf_jyH2]
	subsd  xmm2, [rsp + nb202nf_jzH2]
	subsd  xmm3, [rsp + nb202nf_jxO]
	subsd  xmm4, [rsp + nb202nf_jyO]
	subsd  xmm5, [rsp + nb202nf_jzO]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqOH2], xmm0
	movapd [rsp + nb202nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb202nf_ixH1]
	movapd xmm1, [rsp + nb202nf_iyH1]
	movapd xmm2, [rsp + nb202nf_izH1]
	movapd xmm3, [rsp + nb202nf_ixH1]
	movapd xmm4, [rsp + nb202nf_iyH1]
	movapd xmm5, [rsp + nb202nf_izH1]
	subsd  xmm0, [rsp + nb202nf_jxH1]
	subsd  xmm1, [rsp + nb202nf_jyH1]
	subsd  xmm2, [rsp + nb202nf_jzH1]
	subsd  xmm3, [rsp + nb202nf_jxH2]
	subsd  xmm4, [rsp + nb202nf_jyH2]
	subsd  xmm5, [rsp + nb202nf_jzH2]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb202nf_rsqH1H1], xmm0
	movapd [rsp + nb202nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb202nf_ixH2]
	movapd xmm1, [rsp + nb202nf_iyH2]
	movapd xmm2, [rsp + nb202nf_izH2]
	movapd xmm3, [rsp + nb202nf_ixH2]
	movapd xmm4, [rsp + nb202nf_iyH2]
	movapd xmm5, [rsp + nb202nf_izH2]
	subsd  xmm0, [rsp + nb202nf_jxO]
	subsd  xmm1, [rsp + nb202nf_jyO]
	subsd  xmm2, [rsp + nb202nf_jzO]
	subsd  xmm3, [rsp + nb202nf_jxH1]
	subsd  xmm4, [rsp + nb202nf_jyH1]
	subsd  xmm5, [rsp + nb202nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [rsp + nb202nf_rsqH2O], xmm0
	movapd [rsp + nb202nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb202nf_ixH2]
	movapd xmm1, [rsp + nb202nf_iyH2]
	movapd xmm2, [rsp + nb202nf_izH2]
	subsd  xmm0, [rsp + nb202nf_jxH2]
	subsd  xmm1, [rsp + nb202nf_jyH2]
	subsd  xmm2, [rsp + nb202nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [rsp + nb202nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb202nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH2H2], xmm1
	movapd [rsp + nb202nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb202nf_rsqOO]
	movapd xmm4, [rsp + nb202nf_rsqOH1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb202nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb202nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb202nf_half] ;# rinv
	movapd [rsp + nb202nf_rinvOO], xmm1
	movapd [rsp + nb202nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb202nf_rsqOH2]
	movapd xmm4, [rsp + nb202nf_rsqH1O]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb202nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvOH2], xmm1
	movapd [rsp + nb202nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb202nf_rsqH1H1]
	movapd xmm4, [rsp + nb202nf_rsqH1H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb202nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb202nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb202nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb202nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb202nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH1H1], xmm1
	movapd [rsp + nb202nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb202nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb202nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [rsp + nb202nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb202nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [rsp + nb202nf_half] ;# rinv 
	movapd [rsp + nb202nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm6, [rsp + nb202nf_krf]
	mulsd  xmm6, [rsp + nb202nf_rsqOO]	;# xmm5=krsq 
	addsd  xmm6, [rsp + nb202nf_rinvOO]	;# xmm6=rinv+ krsq 
	subsd  xmm6, [rsp + nb202nf_crf]
	
	mulsd  xmm6, [rsp + nb202nf_qqOO] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, [rsp + nb202nf_vctot] ;# local vctot summation variable 

	;# O-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulsd  xmm5, [rsp + nb202nf_rsqOH1]	;# xmm5=krsq 
	addsd  xmm5, [rsp + nb202nf_rinvOH1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [rsp + nb202nf_crf]
	
	mulsd  xmm5, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# O-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulsd  xmm7, [rsp + nb202nf_rsqOH2]	;# xmm5=krsq 
	addsd  xmm7, [rsp + nb202nf_rinvOH2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [rsp + nb202nf_crf]
	
	mulsd  xmm7, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 

	;# H1-O interaction 
	movapd xmm4, [rsp + nb202nf_krf]
	mulsd  xmm4, [rsp + nb202nf_rsqH1O]	;# xmm5=krsq 
	addsd  xmm4, [rsp + nb202nf_rinvH1O]	;# xmm6=rinv+ krsq 
	subsd  xmm4, [rsp + nb202nf_crf]
	
	mulsd  xmm4, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H1-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulsd  xmm5, [rsp + nb202nf_rsqH1H1]	;# xmm5=krsq 
	addsd  xmm5, [rsp + nb202nf_rinvH1H1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [rsp + nb202nf_crf]
	
	mulsd  xmm5, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# H1-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulsd  xmm7, [rsp + nb202nf_rsqH1H2]	;# xmm5=krsq 
	addsd  xmm7, [rsp + nb202nf_rinvH1H2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [rsp + nb202nf_crf]
	
	mulsd  xmm7, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 

	;# H2-O interaction 
	movapd xmm4, [rsp + nb202nf_krf]
	mulsd  xmm4, [rsp + nb202nf_rsqH2O]	;# xmm5=krsq 
	addsd  xmm4, [rsp + nb202nf_rinvH2O]	;# xmm6=rinv+ krsq 
	subsd  xmm4, [rsp + nb202nf_crf]
	
	mulsd  xmm4, [rsp + nb202nf_qqOH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm4 ;# local vctot summation variable 
	
	;# H2-H1 interaction 
	movapd xmm5, [rsp + nb202nf_krf]
	mulsd  xmm5, [rsp + nb202nf_rsqH2H1]	;# xmm5=krsq 
	addsd  xmm5, [rsp + nb202nf_rinvH2H1]	;# xmm6=rinv+ krsq 
	subsd  xmm5, [rsp + nb202nf_crf]
	
	mulsd  xmm5, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm5 ;# local vctot summation variable 

	;# H2-H2 interaction 
	movapd xmm7, [rsp + nb202nf_krf]
	mulsd  xmm7, [rsp + nb202nf_rsqH2H2]	;# xmm5=krsq 
	addsd  xmm7, [rsp + nb202nf_rinvH2H2]	;# xmm6=rinv+ krsq 
	subsd  xmm7, [rsp + nb202nf_crf]
	
	mulsd  xmm7, [rsp + nb202nf_qqHH] ;# xmm6=voul=qq*(rinv+ krsq-crf) 
	addsd  xmm6, xmm7 ;# local vctot summation variable 
	movlpd [rsp + nb202nf_vctot], xmm6
	
.nb202nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb202nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb202nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb202nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb202nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb202nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb202nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb202nf_n], esi
        jmp .nb202nf_outer
.nb202nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb202nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb202nf_end
        ;# non-zero, do one more workunit
        jmp   .nb202nf_threadloop
.nb202nf_end:
	mov eax, [rsp + nb202nf_nouter]
	mov ebx, [rsp + nb202nf_ninner]
	mov rcx, [rbp + nb202nf_outeriter]
	mov rdx, [rbp + nb202nf_inneriter]
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
