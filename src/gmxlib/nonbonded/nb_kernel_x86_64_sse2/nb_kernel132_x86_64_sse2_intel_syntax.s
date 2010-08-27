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

.globl nb_kernel132_x86_64_sse2
.globl _nb_kernel132_x86_64_sse2
nb_kernel132_x86_64_sse2:	
_nb_kernel132_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb132_fshift,           16
.equiv          nb132_gid,              24
.equiv          nb132_pos,              32
.equiv          nb132_faction,          40
.equiv          nb132_charge,           48
.equiv          nb132_p_facel,          56
.equiv          nb132_argkrf,           64
.equiv          nb132_argcrf,           72
.equiv          nb132_Vc,               80
.equiv          nb132_type,             88
.equiv          nb132_p_ntype,          96
.equiv          nb132_vdwparam,         104
.equiv          nb132_Vvdw,             112
.equiv          nb132_p_tabscale,       120
.equiv          nb132_VFtab,            128
.equiv          nb132_invsqrta,         136
.equiv          nb132_dvda,             144
.equiv          nb132_p_gbtabscale,     152
.equiv          nb132_GBtab,            160
.equiv          nb132_p_nthreads,       168
.equiv          nb132_count,            176
.equiv          nb132_mtx,              184
.equiv          nb132_outeriter,        192
.equiv          nb132_inneriter,        200
.equiv          nb132_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb132_ixO,              0
.equiv          nb132_iyO,              16
.equiv          nb132_izO,              32
.equiv          nb132_ixH1,             48
.equiv          nb132_iyH1,             64
.equiv          nb132_izH1,             80
.equiv          nb132_ixH2,             96
.equiv          nb132_iyH2,             112
.equiv          nb132_izH2,             128
.equiv          nb132_jxO,              144
.equiv          nb132_jyO,              160
.equiv          nb132_jzO,              176
.equiv          nb132_jxH1,             192
.equiv          nb132_jyH1,             208
.equiv          nb132_jzH1,             224
.equiv          nb132_jxH2,             240
.equiv          nb132_jyH2,             256
.equiv          nb132_jzH2,             272
.equiv          nb132_dxOO,             288
.equiv          nb132_dyOO,             304
.equiv          nb132_dzOO,             320
.equiv          nb132_dxOH1,            336
.equiv          nb132_dyOH1,            352
.equiv          nb132_dzOH1,            368
.equiv          nb132_dxOH2,            384
.equiv          nb132_dyOH2,            400
.equiv          nb132_dzOH2,            416
.equiv          nb132_dxH1O,            432
.equiv          nb132_dyH1O,            448
.equiv          nb132_dzH1O,            464
.equiv          nb132_dxH1H1,           480
.equiv          nb132_dyH1H1,           496
.equiv          nb132_dzH1H1,           512
.equiv          nb132_dxH1H2,           528
.equiv          nb132_dyH1H2,           544
.equiv          nb132_dzH1H2,           560
.equiv          nb132_dxH2O,            576
.equiv          nb132_dyH2O,            592
.equiv          nb132_dzH2O,            608
.equiv          nb132_dxH2H1,           624
.equiv          nb132_dyH2H1,           640
.equiv          nb132_dzH2H1,           656
.equiv          nb132_dxH2H2,           672
.equiv          nb132_dyH2H2,           688
.equiv          nb132_dzH2H2,           704
.equiv          nb132_qqOO,             720
.equiv          nb132_qqOH,             736
.equiv          nb132_qqHH,             752
.equiv          nb132_c6,               768
.equiv          nb132_c12,              784
.equiv          nb132_tsc,              800
.equiv          nb132_fstmp,            816
.equiv          nb132_vctot,            832
.equiv          nb132_Vvdwtot,          848
.equiv          nb132_fixO,             864
.equiv          nb132_fiyO,             880
.equiv          nb132_fizO,             896
.equiv          nb132_fixH1,            912
.equiv          nb132_fiyH1,            928
.equiv          nb132_fizH1,            944
.equiv          nb132_fixH2,            960
.equiv          nb132_fiyH2,            976
.equiv          nb132_fizH2,            992
.equiv          nb132_fjxO,             1008
.equiv          nb132_fjyO,             1024
.equiv          nb132_fjzO,             1040
.equiv          nb132_fjxH1,            1056
.equiv          nb132_fjyH1,            1072
.equiv          nb132_fjzH1,            1088
.equiv          nb132_fjxH2,            1104
.equiv          nb132_fjyH2,            1120
.equiv          nb132_fjzH2,            1136
.equiv          nb132_half,             1152
.equiv          nb132_three,            1168
.equiv          nb132_rsqOO,            1184
.equiv          nb132_rsqOH1,           1200
.equiv          nb132_rsqOH2,           1216
.equiv          nb132_rsqH1O,           1232
.equiv          nb132_rsqH1H1,          1248
.equiv          nb132_rsqH1H2,          1264
.equiv          nb132_rsqH2O,           1280
.equiv          nb132_rsqH2H1,          1296
.equiv          nb132_rsqH2H2,          1312
.equiv          nb132_rinvOO,           1328
.equiv          nb132_rinvOH1,          1344
.equiv          nb132_rinvOH2,          1360
.equiv          nb132_rinvH1O,          1376
.equiv          nb132_rinvH1H1,         1392
.equiv          nb132_rinvH1H2,         1408
.equiv          nb132_rinvH2O,          1424
.equiv          nb132_rinvH2H1,         1440
.equiv          nb132_rinvH2H2,         1456
.equiv          nb132_two,              1472
.equiv          nb132_krf,              1488
.equiv          nb132_crf,              1504
.equiv          nb132_nri,              1520
.equiv          nb132_iinr,             1528
.equiv          nb132_jindex,           1536
.equiv          nb132_jjnr,             1544
.equiv          nb132_shift,            1552
.equiv          nb132_shiftvec,         1560
.equiv          nb132_facel,            1568
.equiv          nb132_innerjjnr,        1576
.equiv          nb132_is3,              1584
.equiv          nb132_ii3,              1588
.equiv          nb132_innerk,           1592
.equiv          nb132_n,                1596
.equiv          nb132_nn1,              1600
.equiv          nb132_nouter,           1604
.equiv          nb132_ninner,           1608

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
	sub rsp, 1616		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb132_nouter], eax
	mov [rsp + nb132_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb132_nri], edi
	mov [rsp + nb132_iinr], rsi
	mov [rsp + nb132_jindex], rdx
	mov [rsp + nb132_jjnr], rcx
	mov [rsp + nb132_shift], r8
	mov [rsp + nb132_shiftvec], r9
	mov rsi, [rbp + nb132_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb132_facel], xmm0

	mov rax, [rbp + nb132_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb132_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb132_half], eax
	mov [rsp + nb132_half+4], ebx
	movsd xmm1, [rsp + nb132_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb132_half], xmm1
	movapd [rsp + nb132_two], xmm2
	movapd [rsp + nb132_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb132_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb132_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	

	movsd xmm6, [rsp + nb132_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb132_qqOO], xmm3
	movapd [rsp + nb132_qqOH], xmm4
	movapd [rsp + nb132_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb132_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb132_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb132_vdwparam]
	movlpd xmm0, [rax + rdx*8] 
	movlpd xmm1, [rax + rdx*8 + 8] 
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb132_c6], xmm0
	movapd [rsp + nb132_c12], xmm1

.nb132_threadloop:
        mov   rsi, [rbp + nb132_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb132_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb132_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb132_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb132_n], eax
        mov [rsp + nb132_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb132_outerstart
        jmp .nb132_end

.nb132_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb132_nouter]
	mov [rsp + nb132_nouter], ebx

.nb132_outer:
	mov   rax, [rsp + nb132_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb132_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb132_shiftvec]   ;# eax = base of shiftvec[] 

	movlpd xmm0, [rax + rbx*8]
	movlpd xmm1, [rax + rbx*8 + 8]
	movlpd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb132_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb132_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb132_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb132_ixO], xmm3
	movapd [rsp + nb132_iyO], xmm4
	movapd [rsp + nb132_izO], xmm5

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
	movapd [rsp + nb132_ixH1], xmm0
	movapd [rsp + nb132_iyH1], xmm1
	movapd [rsp + nb132_izH1], xmm2
	movapd [rsp + nb132_ixH2], xmm3
	movapd [rsp + nb132_iyH2], xmm4
	movapd [rsp + nb132_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb132_vctot], xmm4
	movapd [rsp + nb132_Vvdwtot], xmm4
	movapd [rsp + nb132_fixO], xmm4
	movapd [rsp + nb132_fiyO], xmm4
	movapd [rsp + nb132_fizO], xmm4
	movapd [rsp + nb132_fixH1], xmm4
	movapd [rsp + nb132_fiyH1], xmm4
	movapd [rsp + nb132_fizH1], xmm4
	movapd [rsp + nb132_fixH2], xmm4
	movapd [rsp + nb132_fiyH2], xmm4
	movapd [rsp + nb132_fizH2], xmm4
	
	mov   rax, [rsp + nb132_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb132_pos]
	mov   rdi, [rbp + nb132_faction]	
	mov   rax, [rsp + nb132_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb132_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb132_ninner]
	mov   [rsp + nb132_ninner], ecx
	add   edx, 0
	mov   [rsp + nb132_innerk], edx    ;# number of innerloop atoms 
	jge   .nb132_unroll_loop
	jmp   .nb132_checksingle
.nb132_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb132_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb132_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb132_pos]       ;# base of pos[] 

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
    
    subpd xmm0, [rsp + nb132_ixO]
    subpd xmm1, [rsp + nb132_iyO]
    subpd xmm2, [rsp + nb132_izO]
    subpd xmm3, [rsp + nb132_ixH1]
    subpd xmm4, [rsp + nb132_iyH1]
    subpd xmm5, [rsp + nb132_izH1]
    subpd xmm6, [rsp + nb132_ixH2]
    subpd xmm7, [rsp + nb132_iyH2]
    subpd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOO], xmm0
	movapd [rsp + nb132_dyOO], xmm1
	movapd [rsp + nb132_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb132_dxH1O], xmm3
	movapd [rsp + nb132_dyH1O], xmm4
	movapd [rsp + nb132_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb132_dxH2O], xmm6
	movapd [rsp + nb132_dyH2O], xmm7
	movapd [rsp + nb132_dzH2O], xmm8
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
		
	movapd  xmm9, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
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
		
	movapd  xmm1, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
	mulpd   xmm9, xmm15  ;#  rinvOO 
	mulpd   xmm10, xmm15 ;#   rinvH1O
    mulpd   xmm11, xmm15 ;#   rinvH2O
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movapd [rsp + nb132_rsqOO], xmm0
    movapd [rsp + nb132_rsqOH1], xmm3
    movapd [rsp + nb132_rsqOH2], xmm6
    movapd [rsp + nb132_rinvOO], xmm9
    movapd [rsp + nb132_rinvOH1], xmm10
    movapd [rsp + nb132_rinvOH2], xmm11

    ;# table LJ interaction
    mulpd  xmm0, xmm9
    mulpd  xmm0, [rsp + nb132_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
         
    ;# multiply by 8
    pslld   xmm1, 3

    ;# move to integer registers
    pshufd  xmm13, xmm1, 1
    movd    r8d, xmm1
    movd    r10d, xmm13
    
    ;# calculate eps
    subpd     xmm0, xmm2
    mov  rsi, [rbp + nb132_VFtab]
            
    movlpd xmm4, [rsi + r8*8]
   	movlpd xmm5, [rsi + r8*8 + 8]
	movlpd xmm6, [rsi + r8*8 + 16]
	movlpd xmm7, [rsi + r8*8 + 24]
    movlpd xmm8, [rsi + r8*8 + 32]
   	movlpd xmm9, [rsi + r8*8 + 40]
	movlpd xmm10, [rsi + r8*8 + 48]
	movlpd xmm11, [rsi + r8*8 + 56]
    
    movhpd xmm4, [rsi + r10*8]
   	movhpd xmm5, [rsi + r10*8 + 8]
	movhpd xmm6, [rsi + r10*8 + 16]
	movhpd xmm7, [rsi + r10*8 + 24]
    movhpd xmm8, [rsi + r10*8 + 32]
   	movhpd xmm9, [rsi + r10*8 + 40]
	movhpd xmm10, [rsi + r10*8 + 48]
	movhpd xmm11, [rsi + r10*8 + 56]
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  xmm7, xmm0    ;# Heps
    mulpd  xmm11, xmm0 
    mulpd  xmm6, xmm0   ;# Geps
    mulpd  xmm10, xmm0 
    mulpd  xmm7, xmm0   ;# Heps2
    mulpd  xmm11, xmm0 
    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm0  ;# eps*Fp
    mulpd  xmm9, xmm0
    movapd xmm12, [rsp + nb132_c6]
    movapd xmm13, [rsp + nb132_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb132_Vvdwtot]
    movapd [rsp + nb132_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    mulpd  xmm7, [rsp + nb132_tsc]

    movapd xmm9, [rsp + nb132_rinvOO]
    movapd xmm10, [rsp + nb132_rinvOH1]
    movapd xmm11, [rsp + nb132_rinvOH2]

    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    
    mulpd  xmm0, [rsp + nb132_qqOO] 
    mulpd  xmm1, [rsp + nb132_qqOH] 
    mulpd  xmm2, [rsp + nb132_qqOH] 
    
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    subpd  xmm9, xmm7
    mulpd  xmm9, [rsp + nb132_rinvOO]
    
    addpd xmm0, [rsp + nb132_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb132_vctot], xmm0
    
    ;# move j O forces to xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
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

	mulpd xmm7, [rsp + nb132_dxOO]
	mulpd xmm8, [rsp + nb132_dyOO]
	mulpd xmm9, [rsp + nb132_dzOO]
	mulpd xmm10, [rsp + nb132_dxH1O]
	mulpd xmm11, [rsp + nb132_dyH1O]
	mulpd xmm12, [rsp + nb132_dzH1O]
	mulpd xmm13, [rsp + nb132_dxH2O]
	mulpd xmm14, [rsp + nb132_dyH2O]
	mulpd xmm15, [rsp + nb132_dzH2O]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb132_fixO]
    addpd xmm8, [rsp + nb132_fiyO]
    addpd xmm9, [rsp + nb132_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb132_fixH1]
    addpd xmm11, [rsp + nb132_fiyH1]
    addpd xmm12, [rsp + nb132_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb132_fixH2]
    addpd xmm14, [rsp + nb132_fiyH2]
    addpd xmm15, [rsp + nb132_fizH2]

    movapd [rsp + nb132_fixO], xmm7
    movapd [rsp + nb132_fiyO], xmm8
    movapd [rsp + nb132_fizO], xmm9
    movapd [rsp + nb132_fixH1], xmm10
    movapd [rsp + nb132_fiyH1], xmm11
    movapd [rsp + nb132_fizH1], xmm12
    movapd [rsp + nb132_fixH2], xmm13
    movapd [rsp + nb132_fiyH2], xmm14
    movapd [rsp + nb132_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb132_pos]
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
    
    subpd xmm0, [rsp + nb132_ixO]
    subpd xmm1, [rsp + nb132_iyO]
    subpd xmm2, [rsp + nb132_izO]
    subpd xmm3, [rsp + nb132_ixH1]
    subpd xmm4, [rsp + nb132_iyH1]
    subpd xmm5, [rsp + nb132_izH1]
    subpd xmm6, [rsp + nb132_ixH2]
    subpd xmm7, [rsp + nb132_iyH2]
    subpd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOH1], xmm0
	movapd [rsp + nb132_dyOH1], xmm1
	movapd [rsp + nb132_dzOH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb132_dxH1H1], xmm3
	movapd [rsp + nb132_dyH1H1], xmm4
	movapd [rsp + nb132_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb132_dxH2H1], xmm6
	movapd [rsp + nb132_dyH2H1], xmm7
	movapd [rsp + nb132_dzH2H1], xmm8
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
		
	movapd  xmm9, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
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
		
	movapd  xmm1, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
	mulpd   xmm9, xmm15  ;#  rinvOH1
	mulpd   xmm10, xmm15 ;#   rinvH1H1
    mulpd   xmm11, xmm15 ;#   rinvH2H1
	
	;# H1 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb132_qqOH] 
    mulpd  xmm1, [rsp + nb132_qqHH] 
    mulpd  xmm2, [rsp + nb132_qqHH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb132_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb132_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
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

	mulpd xmm7, [rsp + nb132_dxOH1]
	mulpd xmm8, [rsp + nb132_dyOH1]
	mulpd xmm9, [rsp + nb132_dzOH1]
	mulpd xmm10, [rsp + nb132_dxH1H1]
	mulpd xmm11, [rsp + nb132_dyH1H1]
	mulpd xmm12, [rsp + nb132_dzH1H1]
	mulpd xmm13, [rsp + nb132_dxH2H1]
	mulpd xmm14, [rsp + nb132_dyH2H1]
	mulpd xmm15, [rsp + nb132_dzH2H1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb132_fixO]
    addpd xmm8, [rsp + nb132_fiyO]
    addpd xmm9, [rsp + nb132_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb132_fixH1]
    addpd xmm11, [rsp + nb132_fiyH1]
    addpd xmm12, [rsp + nb132_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb132_fixH2]
    addpd xmm14, [rsp + nb132_fiyH2]
    addpd xmm15, [rsp + nb132_fizH2]

    movapd [rsp + nb132_fixO], xmm7
    movapd [rsp + nb132_fiyO], xmm8
    movapd [rsp + nb132_fizO], xmm9
    movapd [rsp + nb132_fixH1], xmm10
    movapd [rsp + nb132_fiyH1], xmm11
    movapd [rsp + nb132_fizH1], xmm12
    movapd [rsp + nb132_fixH2], xmm13
    movapd [rsp + nb132_fiyH2], xmm14
    movapd [rsp + nb132_fizH2], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb132_pos]
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
    
    subpd xmm0, [rsp + nb132_ixO]
    subpd xmm1, [rsp + nb132_iyO]
    subpd xmm2, [rsp + nb132_izO]
    subpd xmm3, [rsp + nb132_ixH1]
    subpd xmm4, [rsp + nb132_iyH1]
    subpd xmm5, [rsp + nb132_izH1]
    subpd xmm6, [rsp + nb132_ixH2]
    subpd xmm7, [rsp + nb132_iyH2]
    subpd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOH2], xmm0
	movapd [rsp + nb132_dyOH2], xmm1
	movapd [rsp + nb132_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb132_dxH1H2], xmm3
	movapd [rsp + nb132_dyH1H2], xmm4
	movapd [rsp + nb132_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb132_dxH2H2], xmm6
	movapd [rsp + nb132_dyH2H2], xmm7
	movapd [rsp + nb132_dzH2H2], xmm8
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
		
	movapd  xmm9, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
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
		
	movapd  xmm1, [rsp + nb132_three]
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

	movapd  xmm15, [rsp + nb132_half]
	mulpd   xmm9, xmm15  ;#  rinvOH2
	mulpd   xmm10, xmm15 ;#   rinvH1H2
    mulpd   xmm11, xmm15 ;#   rinvH2H2
	
	;# H2 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb132_qqOH] 
    mulpd  xmm1, [rsp + nb132_qqHH] 
    mulpd  xmm2, [rsp + nb132_qqHH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb132_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb132_vctot], xmm0
    
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

	mulpd xmm7, [rsp + nb132_dxOH2]
	mulpd xmm8, [rsp + nb132_dyOH2]
	mulpd xmm9, [rsp + nb132_dzOH2]
	mulpd xmm10, [rsp + nb132_dxH1H2]
	mulpd xmm11, [rsp + nb132_dyH1H2]
	mulpd xmm12, [rsp + nb132_dzH1H2]
	mulpd xmm13, [rsp + nb132_dxH2H2]
	mulpd xmm14, [rsp + nb132_dyH2H2]
	mulpd xmm15, [rsp + nb132_dzH2H2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb132_fixO]
    addpd xmm8, [rsp + nb132_fiyO]
    addpd xmm9, [rsp + nb132_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb132_fixH1]
    addpd xmm11, [rsp + nb132_fiyH1]
    addpd xmm12, [rsp + nb132_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb132_fixH2]
    addpd xmm14, [rsp + nb132_fiyH2]
    addpd xmm15, [rsp + nb132_fizH2]

    movapd [rsp + nb132_fixO], xmm7
    movapd [rsp + nb132_fiyO], xmm8
    movapd [rsp + nb132_fizO], xmm9
    movapd [rsp + nb132_fixH1], xmm10
    movapd [rsp + nb132_fiyH1], xmm11
    movapd [rsp + nb132_fizH1], xmm12
    movapd [rsp + nb132_fixH2], xmm13
    movapd [rsp + nb132_fiyH2], xmm14
    movapd [rsp + nb132_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb132_innerk],  2
	jl    .nb132_checksingle
	jmp   .nb132_unroll_loop
.nb132_checksingle:
	mov   edx, [rsp + nb132_innerk]
	and   edx, 1
	jnz   .nb132_dosingle
	jmp   .nb132_updateouterdata
.nb132_dosingle:
	mov   rdx, [rsp + nb132_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	
	add qword ptr [rsp + nb132_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb132_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move j O coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8] 
    movsd xmm1, [rsi + rax*8 + 8] 
    movsd xmm2, [rsi + rax*8 + 16] 

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subsd xmm0, [rsp + nb132_ixO]
    subsd xmm1, [rsp + nb132_iyO]
    subsd xmm2, [rsp + nb132_izO]
    subsd xmm3, [rsp + nb132_ixH1]
    subsd xmm4, [rsp + nb132_iyH1]
    subsd xmm5, [rsp + nb132_izH1]
    subsd xmm6, [rsp + nb132_ixH2]
    subsd xmm7, [rsp + nb132_iyH2]
    subsd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOO], xmm0
	movapd [rsp + nb132_dyOO], xmm1
	movapd [rsp + nb132_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [rsp + nb132_dxH1O], xmm3
	movapd [rsp + nb132_dyH1O], xmm4
	movapd [rsp + nb132_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movapd [rsp + nb132_dxH2O], xmm6
	movapd [rsp + nb132_dyH2O], xmm7
	movapd [rsp + nb132_dzH2O], xmm8
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
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb132_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOO 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1O
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2O	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb132_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;#  rinvOO 
	mulsd   xmm10, xmm15 ;#   rinvH1O
    mulsd   xmm11, xmm15 ;#   rinvH2O
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movapd [rsp + nb132_rsqOO], xmm0
    movapd [rsp + nb132_rsqOH1], xmm3
    movapd [rsp + nb132_rsqOH2], xmm6
    movapd [rsp + nb132_rinvOO], xmm9
    movapd [rsp + nb132_rinvOH1], xmm10
    movapd [rsp + nb132_rinvOH2], xmm11

    ;# table LJ interaction
    mulsd  xmm0, xmm9
    mulsd  xmm0, [rsp + nb132_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
         
    ;# multiply by 8
    shl   r8d, 3

    ;# calculate eps
    subsd     xmm0, xmm2
    mov  rsi, [rbp + nb132_VFtab]
            
    movsd xmm4, [rsi + r8*8]
   	movsd xmm5, [rsi + r8*8 + 8]
	movsd xmm6, [rsi + r8*8 + 16]
	movsd xmm7, [rsi + r8*8 + 24]
    movsd xmm8, [rsi + r8*8 + 32]
   	movsd xmm9, [rsi + r8*8 + 40]
	movsd xmm10, [rsi + r8*8 + 48]
	movsd xmm11, [rsi + r8*8 + 56]
    
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulsd  xmm7, xmm0    ;# Heps
    mulsd  xmm11, xmm0 
    mulsd  xmm6, xmm0   ;# Geps
    mulsd  xmm10, xmm0 
    mulsd  xmm7, xmm0   ;# Heps2
    mulsd  xmm11, xmm0 
    addsd  xmm5, xmm6  ;# F+Geps
    addsd  xmm9, xmm10 
    addsd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addsd  xmm9, xmm11 
    addsd  xmm7, xmm7    ;# 2*Heps2
    addsd  xmm11, xmm11
    addsd  xmm7, xmm6   ;# 2*Heps2+Geps
    addsd  xmm11, xmm10
    
    addsd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm11, xmm9
    mulsd  xmm5, xmm0  ;# eps*Fp
    mulsd  xmm9, xmm0
    movapd xmm12, [rsp + nb132_c6]
    movapd xmm13, [rsp + nb132_c12]
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulsd  xmm9, xmm13  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb132_Vvdwtot]
    movsd [rsp + nb132_Vvdwtot], xmm5
        
    mulsd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulsd  xmm11, xmm13   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    mulsd  xmm7, [rsp + nb132_tsc]

    movapd xmm9, [rsp + nb132_rinvOO]
    movapd xmm10, [rsp + nb132_rinvOH1]
    movapd xmm11, [rsp + nb132_rinvOH2]

    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    
    mulsd  xmm0, [rsp + nb132_qqOO] 
    mulsd  xmm1, [rsp + nb132_qqOH] 
    mulsd  xmm2, [rsp + nb132_qqOH] 
    
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    subsd  xmm9, xmm7
    mulsd  xmm9, [rsp + nb132_rinvOO]
    
    addsd xmm0, [rsp + nb132_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb132_vctot], xmm0
    
    ;# move j O forces to xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movsd xmm0, [rdi + rax*8]
	movsd xmm1, [rdi + rax*8 + 8]
	movsd xmm2, [rdi + rax*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulsd xmm7, [rsp + nb132_dxOO]
	mulsd xmm8, [rsp + nb132_dyOO]
	mulsd xmm9, [rsp + nb132_dzOO]
	mulsd xmm10, [rsp + nb132_dxH1O]
	mulsd xmm11, [rsp + nb132_dyH1O]
	mulsd xmm12, [rsp + nb132_dzH1O]
	mulsd xmm13, [rsp + nb132_dxH2O]
	mulsd xmm14, [rsp + nb132_dyH2O]
	mulsd xmm15, [rsp + nb132_dzH2O]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb132_fixO]
    addsd xmm8, [rsp + nb132_fiyO]
    addsd xmm9, [rsp + nb132_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb132_fixH1]
    addsd xmm11, [rsp + nb132_fiyH1]
    addsd xmm12, [rsp + nb132_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb132_fixH2]
    addsd xmm14, [rsp + nb132_fiyH2]
    addsd xmm15, [rsp + nb132_fizH2]

    movsd [rsp + nb132_fixO], xmm7
    movsd [rsp + nb132_fiyO], xmm8
    movsd [rsp + nb132_fizO], xmm9
    movsd [rsp + nb132_fixH1], xmm10
    movsd [rsp + nb132_fiyH1], xmm11
    movsd [rsp + nb132_fizH1], xmm12
    movsd [rsp + nb132_fixH2], xmm13
    movsd [rsp + nb132_fiyH2], xmm14
    movsd [rsp + nb132_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movsd [rdi + rax*8],      xmm0
	movsd [rdi + rax*8 + 8],  xmm1
	movsd [rdi + rax*8 + 16], xmm2

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb132_pos]
    movsd xmm0, [rsi + rax*8 + 24] 
    movsd xmm1, [rsi + rax*8 + 32] 
    movsd xmm2, [rsi + rax*8 + 40] 
   
    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subsd xmm0, [rsp + nb132_ixO]
    subsd xmm1, [rsp + nb132_iyO]
    subsd xmm2, [rsp + nb132_izO]
    subsd xmm3, [rsp + nb132_ixH1]
    subsd xmm4, [rsp + nb132_iyH1]
    subsd xmm5, [rsp + nb132_izH1]
    subsd xmm6, [rsp + nb132_ixH2]
    subsd xmm7, [rsp + nb132_iyH2]
    subsd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOH1], xmm0
	movapd [rsp + nb132_dyOH1], xmm1
	movapd [rsp + nb132_dzOH1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [rsp + nb132_dxH1H1], xmm3
	movapd [rsp + nb132_dyH1H1], xmm4
	movapd [rsp + nb132_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movapd [rsp + nb132_dxH2H1], xmm6
	movapd [rsp + nb132_dyH2H1], xmm7
	movapd [rsp + nb132_dzH2H1], xmm8
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
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb132_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2OH1

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb132_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;#  rinvOH1
	mulsd   xmm10, xmm15 ;#   rinvH1H1
    mulsd   xmm11, xmm15 ;#   rinvH2H1
	
	;# H1 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb132_qqOH] 
    mulsd  xmm1, [rsp + nb132_qqHH] 
    mulsd  xmm2, [rsp + nb132_qqHH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb132_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb132_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movsd xmm0, [rdi + rax*8 + 24]
	movsd xmm1, [rdi + rax*8 + 32]
	movsd xmm2, [rdi + rax*8 + 40]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulsd xmm7, [rsp + nb132_dxOH1]
	mulsd xmm8, [rsp + nb132_dyOH1]
	mulsd xmm9, [rsp + nb132_dzOH1]
	mulsd xmm10, [rsp + nb132_dxH1H1]
	mulsd xmm11, [rsp + nb132_dyH1H1]
	mulsd xmm12, [rsp + nb132_dzH1H1]
	mulsd xmm13, [rsp + nb132_dxH2H1]
	mulsd xmm14, [rsp + nb132_dyH2H1]
	mulsd xmm15, [rsp + nb132_dzH2H1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb132_fixO]
    addsd xmm8, [rsp + nb132_fiyO]
    addsd xmm9, [rsp + nb132_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb132_fixH1]
    addsd xmm11, [rsp + nb132_fiyH1]
    addsd xmm12, [rsp + nb132_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb132_fixH2]
    addsd xmm14, [rsp + nb132_fiyH2]
    addsd xmm15, [rsp + nb132_fizH2]

    movsd [rsp + nb132_fixO], xmm7
    movsd [rsp + nb132_fiyO], xmm8
    movsd [rsp + nb132_fizO], xmm9
    movsd [rsp + nb132_fixH1], xmm10
    movsd [rsp + nb132_fiyH1], xmm11
    movsd [rsp + nb132_fizH1], xmm12
    movsd [rsp + nb132_fixH2], xmm13
    movsd [rsp + nb132_fiyH2], xmm14
    movsd [rsp + nb132_fizH2], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	mov   rdi, [rbp + nb132_faction]	
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb132_pos]
    movsd xmm0, [rsi + rax*8 + 48] 
    movsd xmm1, [rsi + rax*8 + 56] 
    movsd xmm2, [rsi + rax*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subsd xmm0, [rsp + nb132_ixO]
    subsd xmm1, [rsp + nb132_iyO]
    subsd xmm2, [rsp + nb132_izO]
    subsd xmm3, [rsp + nb132_ixH1]
    subsd xmm4, [rsp + nb132_iyH1]
    subsd xmm5, [rsp + nb132_izH1]
    subsd xmm6, [rsp + nb132_ixH2]
    subsd xmm7, [rsp + nb132_iyH2]
    subsd xmm8, [rsp + nb132_izH2]
    
	movapd [rsp + nb132_dxOH2], xmm0
	movapd [rsp + nb132_dyOH2], xmm1
	movapd [rsp + nb132_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movapd [rsp + nb132_dxH1H2], xmm3
	movapd [rsp + nb132_dyH1H2], xmm4
	movapd [rsp + nb132_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movapd [rsp + nb132_dxH2H2], xmm6
	movapd [rsp + nb132_dyH2H2], xmm7
	movapd [rsp + nb132_dzH2H2], xmm8
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
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb132_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2H2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb132_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb132_half]
	mulsd   xmm9, xmm15  ;#  rinvOH2
	mulsd   xmm10, xmm15 ;#   rinvH1H2
    mulsd   xmm11, xmm15 ;#   rinvH2H2
	
	;# H2 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb132_qqOH] 
    mulsd  xmm1, [rsp + nb132_qqHH] 
    mulsd  xmm2, [rsp + nb132_qqHH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb132_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb132_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 48]
	movsd xmm1, [rdi + rax*8 + 56]
	movsd xmm2, [rdi + rax*8 + 64]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulsd xmm7, [rsp + nb132_dxOH2]
	mulsd xmm8, [rsp + nb132_dyOH2]
	mulsd xmm9, [rsp + nb132_dzOH2]
	mulsd xmm10, [rsp + nb132_dxH1H2]
	mulsd xmm11, [rsp + nb132_dyH1H2]
	mulsd xmm12, [rsp + nb132_dzH1H2]
	mulsd xmm13, [rsp + nb132_dxH2H2]
	mulsd xmm14, [rsp + nb132_dyH2H2]
	mulsd xmm15, [rsp + nb132_dzH2H2]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb132_fixO]
    addsd xmm8, [rsp + nb132_fiyO]
    addsd xmm9, [rsp + nb132_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb132_fixH1]
    addsd xmm11, [rsp + nb132_fiyH1]
    addsd xmm12, [rsp + nb132_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb132_fixH2]
    addsd xmm14, [rsp + nb132_fiyH2]
    addsd xmm15, [rsp + nb132_fizH2]

    movsd [rsp + nb132_fixO], xmm7
    movsd [rsp + nb132_fiyO], xmm8
    movsd [rsp + nb132_fizO], xmm9
    movsd [rsp + nb132_fixH1], xmm10
    movsd [rsp + nb132_fiyH1], xmm11
    movsd [rsp + nb132_fizH1], xmm12
    movsd [rsp + nb132_fixH2], xmm13
    movsd [rsp + nb132_fiyH2], xmm14
    movsd [rsp + nb132_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2

.nb132_updateouterdata:
	mov   ecx, [rsp + nb132_ii3]
	mov   rdi, [rbp + nb132_faction]
	mov   rsi, [rbp + nb132_fshift]
	mov   edx, [rsp + nb132_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb132_fixO]
	movapd xmm1, [rsp + nb132_fiyO]
	movapd xmm2, [rsp + nb132_fizO]

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
	movapd xmm0, [rsp + nb132_fixH1]
	movapd xmm1, [rsp + nb132_fiyH1]
	movapd xmm2, [rsp + nb132_fizH1]

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
	movapd xmm0, [rsp + nb132_fixH2]
	movapd xmm1, [rsp + nb132_fiyH2]
	movapd xmm2, [rsp + nb132_fizH2]

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
	mov esi, [rsp + nb132_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb132_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb132_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb132_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb132_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb132_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb132_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb132_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb132_n], esi
        jmp .nb132_outer
.nb132_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb132_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb132_end
        ;# non-zero, do one more workunit
        jmp   .nb132_threadloop
.nb132_end:
	mov eax, [rsp + nb132_nouter]
	mov ebx, [rsp + nb132_ninner]
	mov rcx, [rbp + nb132_outeriter]
	mov rdx, [rbp + nb132_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1616
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


	
.globl nb_kernel132nf_x86_64_sse2
.globl _nb_kernel132nf_x86_64_sse2
nb_kernel132nf_x86_64_sse2:	
_nb_kernel132nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb132nf_fshift,         16
.equiv          nb132nf_gid,            24
.equiv          nb132nf_pos,            32
.equiv          nb132nf_faction,        40
.equiv          nb132nf_charge,         48
.equiv          nb132nf_p_facel,        56
.equiv          nb132nf_argkrf,         64
.equiv          nb132nf_argcrf,         72
.equiv          nb132nf_Vc,             80
.equiv          nb132nf_type,           88
.equiv          nb132nf_p_ntype,        96
.equiv          nb132nf_vdwparam,       104
.equiv          nb132nf_Vvdw,           112
.equiv          nb132nf_p_tabscale,     120
.equiv          nb132nf_VFtab,          128
.equiv          nb132nf_invsqrta,       136
.equiv          nb132nf_dvda,           144
.equiv          nb132nf_p_gbtabscale,   152
.equiv          nb132nf_GBtab,          160
.equiv          nb132nf_p_nthreads,     168
.equiv          nb132nf_count,          176
.equiv          nb132nf_mtx,            184
.equiv          nb132nf_outeriter,      192
.equiv          nb132nf_inneriter,      200
.equiv          nb132nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb132nf_ixO,            0
.equiv          nb132nf_iyO,            16
.equiv          nb132nf_izO,            32
.equiv          nb132nf_ixH1,           48
.equiv          nb132nf_iyH1,           64
.equiv          nb132nf_izH1,           80
.equiv          nb132nf_ixH2,           96
.equiv          nb132nf_iyH2,           112
.equiv          nb132nf_izH2,           128
.equiv          nb132nf_jxO,            144
.equiv          nb132nf_jyO,            160
.equiv          nb132nf_jzO,            176
.equiv          nb132nf_jxH1,           192
.equiv          nb132nf_jyH1,           208
.equiv          nb132nf_jzH1,           224
.equiv          nb132nf_jxH2,           240
.equiv          nb132nf_jyH2,           256
.equiv          nb132nf_jzH2,           272
.equiv          nb132nf_qqOO,           288
.equiv          nb132nf_qqOH,           304
.equiv          nb132nf_qqHH,           320
.equiv          nb132nf_c6,             336
.equiv          nb132nf_c12,            352
.equiv          nb132nf_vctot,          368
.equiv          nb132nf_Vvdwtot,        384
.equiv          nb132nf_half,           400
.equiv          nb132nf_three,          416
.equiv          nb132nf_rsqOO,          432
.equiv          nb132nf_rsqOH1,         448
.equiv          nb132nf_rsqOH2,         464
.equiv          nb132nf_rsqH1O,         480
.equiv          nb132nf_rsqH1H1,        496
.equiv          nb132nf_rsqH1H2,        512
.equiv          nb132nf_rsqH2O,         528
.equiv          nb132nf_rsqH2H1,        544
.equiv          nb132nf_rsqH2H2,        560
.equiv          nb132nf_rinvOO,         576
.equiv          nb132nf_rinvOH1,        592
.equiv          nb132nf_rinvOH2,        608
.equiv          nb132nf_rinvH1O,        624
.equiv          nb132nf_rinvH1H1,       640
.equiv          nb132nf_rinvH1H2,       656
.equiv          nb132nf_rinvH2O,        672
.equiv          nb132nf_rinvH2H1,       688
.equiv          nb132nf_rinvH2H2,       704
.equiv          nb132nf_krf,            720
.equiv          nb132nf_crf,            736
.equiv          nb132nf_tsc,            752
.equiv          nb132nf_nri,            768
.equiv          nb132nf_iinr,           776
.equiv          nb132nf_jindex,         784
.equiv          nb132nf_jjnr,           792
.equiv          nb132nf_shift,          800
.equiv          nb132nf_shiftvec,       808
.equiv          nb132nf_facel,          816
.equiv          nb132nf_innerjjnr,      824
.equiv          nb132nf_is3,            832
.equiv          nb132nf_ii3,            836
.equiv          nb132nf_innerk,         840
.equiv          nb132nf_n,              844
.equiv          nb132nf_nn1,            848
.equiv          nb132nf_nouter,         852
.equiv          nb132nf_ninner,         856

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
	sub rsp, 864		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb132nf_nouter], eax
	mov [rsp + nb132nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb132nf_nri], edi
	mov [rsp + nb132nf_iinr], rsi
	mov [rsp + nb132nf_jindex], rdx
	mov [rsp + nb132nf_jjnr], rcx
	mov [rsp + nb132nf_shift], r8
	mov [rsp + nb132nf_shiftvec], r9
	mov rsi, [rbp + nb132nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb132nf_facel], xmm0

	mov rax, [rbp + nb132nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb132nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb132nf_half], eax
	mov [rsp + nb132nf_half+4], ebx
	movsd xmm1, [rsp + nb132nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb132nf_half], xmm1
	movapd [rsp + nb132nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb132nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb132nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	

	movsd xmm6, [rsp + nb132nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb132nf_qqOO], xmm3
	movapd [rsp + nb132nf_qqOH], xmm4
	movapd [rsp + nb132nf_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb132nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb132nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb132nf_vdwparam]
	movlpd xmm0, [rax + rdx*8] 
	movlpd xmm1, [rax + rdx*8 + 8] 
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb132nf_c6], xmm0
	movapd [rsp + nb132nf_c12], xmm1

.nb132nf_threadloop:
        mov   rsi, [rbp + nb132nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb132nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb132nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb132nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb132nf_n], eax
        mov [rsp + nb132nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb132nf_outerstart
        jmp .nb132nf_end

.nb132nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb132nf_nouter]
	mov [rsp + nb132nf_nouter], ebx

.nb132nf_outer:
	mov   rax, [rsp + nb132nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb132nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb132nf_shiftvec]   ;# eax = base of shiftvec[] 

	movlpd xmm0, [rax + rbx*8]
	movlpd xmm1, [rax + rbx*8 + 8]
	movlpd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb132nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb132nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb132nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb132nf_ixO], xmm3
	movapd [rsp + nb132nf_iyO], xmm4
	movapd [rsp + nb132nf_izO], xmm5

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
	movapd [rsp + nb132nf_ixH1], xmm0
	movapd [rsp + nb132nf_iyH1], xmm1
	movapd [rsp + nb132nf_izH1], xmm2
	movapd [rsp + nb132nf_ixH2], xmm3
	movapd [rsp + nb132nf_iyH2], xmm4
	movapd [rsp + nb132nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb132nf_vctot], xmm4
	movapd [rsp + nb132nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb132nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb132nf_pos]
	mov   rdi, [rbp + nb132nf_faction]	
	mov   rax, [rsp + nb132nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb132nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb132nf_ninner]
	mov   [rsp + nb132nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb132nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb132nf_unroll_loop
	jmp   .nb132nf_checksingle
.nb132nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb132nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb132nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb132nf_pos]       ;# base of pos[] 

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
	movapd 	[rsp + nb132nf_jxO], xmm2
	movapd 	[rsp + nb132nf_jyO], xmm3
	movapd 	[rsp + nb132nf_jzO], xmm4
	movapd 	[rsp + nb132nf_jxH1], xmm5
	movapd 	[rsp + nb132nf_jyH1], xmm6
	movapd 	[rsp + nb132nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm2, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rbx*8 + 56]
	movhpd xmm4, [rsi + rbx*8 + 64]
	movapd 	[rsp + nb132nf_jxH2], xmm2
	movapd 	[rsp + nb132nf_jyH2], xmm3
	movapd 	[rsp + nb132nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb132nf_ixO]
	movapd xmm1, [rsp + nb132nf_iyO]
	movapd xmm2, [rsp + nb132nf_izO]
	movapd xmm3, [rsp + nb132nf_ixO]
	movapd xmm4, [rsp + nb132nf_iyO]
	movapd xmm5, [rsp + nb132nf_izO]
	subpd  xmm0, [rsp + nb132nf_jxO]
	subpd  xmm1, [rsp + nb132nf_jyO]
	subpd  xmm2, [rsp + nb132nf_jzO]
	subpd  xmm3, [rsp + nb132nf_jxH1]
	subpd  xmm4, [rsp + nb132nf_jyH1]
	subpd  xmm5, [rsp + nb132nf_jzH1]
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
	movapd [rsp + nb132nf_rsqOO], xmm0
	movapd [rsp + nb132nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb132nf_ixO]
	movapd xmm1, [rsp + nb132nf_iyO]
	movapd xmm2, [rsp + nb132nf_izO]
	movapd xmm3, [rsp + nb132nf_ixH1]
	movapd xmm4, [rsp + nb132nf_iyH1]
	movapd xmm5, [rsp + nb132nf_izH1]
	subpd  xmm0, [rsp + nb132nf_jxH2]
	subpd  xmm1, [rsp + nb132nf_jyH2]
	subpd  xmm2, [rsp + nb132nf_jzH2]
	subpd  xmm3, [rsp + nb132nf_jxO]
	subpd  xmm4, [rsp + nb132nf_jyO]
	subpd  xmm5, [rsp + nb132nf_jzO]
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
	movapd [rsp + nb132nf_rsqOH2], xmm0
	movapd [rsp + nb132nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb132nf_ixH1]
	movapd xmm1, [rsp + nb132nf_iyH1]
	movapd xmm2, [rsp + nb132nf_izH1]
	movapd xmm3, [rsp + nb132nf_ixH1]
	movapd xmm4, [rsp + nb132nf_iyH1]
	movapd xmm5, [rsp + nb132nf_izH1]
	subpd  xmm0, [rsp + nb132nf_jxH1]
	subpd  xmm1, [rsp + nb132nf_jyH1]
	subpd  xmm2, [rsp + nb132nf_jzH1]
	subpd  xmm3, [rsp + nb132nf_jxH2]
	subpd  xmm4, [rsp + nb132nf_jyH2]
	subpd  xmm5, [rsp + nb132nf_jzH2]
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
	movapd [rsp + nb132nf_rsqH1H1], xmm0
	movapd [rsp + nb132nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb132nf_ixH2]
	movapd xmm1, [rsp + nb132nf_iyH2]
	movapd xmm2, [rsp + nb132nf_izH2]
	movapd xmm3, [rsp + nb132nf_ixH2]
	movapd xmm4, [rsp + nb132nf_iyH2]
	movapd xmm5, [rsp + nb132nf_izH2]
	subpd  xmm0, [rsp + nb132nf_jxO]
	subpd  xmm1, [rsp + nb132nf_jyO]
	subpd  xmm2, [rsp + nb132nf_jzO]
	subpd  xmm3, [rsp + nb132nf_jxH1]
	subpd  xmm4, [rsp + nb132nf_jyH1]
	subpd  xmm5, [rsp + nb132nf_jzH1]
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
	movapd [rsp + nb132nf_rsqH2O], xmm0
	movapd [rsp + nb132nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb132nf_ixH2]
	movapd xmm1, [rsp + nb132nf_iyH2]
	movapd xmm2, [rsp + nb132nf_izH2]
	subpd  xmm0, [rsp + nb132nf_jxH2]
	subpd  xmm1, [rsp + nb132nf_jyH2]
	subpd  xmm2, [rsp + nb132nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [rsp + nb132nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb132nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH2H2], xmm1
	movapd [rsp + nb132nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb132nf_rsqOO]
	movapd xmm4, [rsp + nb132nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb132nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb132nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb132nf_half] ;# rinv
	movapd [rsp + nb132nf_rinvOO], xmm1
	movapd [rsp + nb132nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb132nf_rsqOH2]
	movapd xmm4, [rsp + nb132nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb132nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvOH2], xmm1
	movapd [rsp + nb132nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb132nf_rsqH1H1]
	movapd xmm4, [rsp + nb132nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb132nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb132nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH1H1], xmm1
	movapd [rsp + nb132nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb132nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb132nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb132nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movapd xmm0, [rsp + nb132nf_rinvOO]
	movapd xmm7, xmm0		;# xmm7=rinv 
	mulpd  xmm7, [rsp + nb132nf_qqOO]
	
	addpd  xmm7, [rsp + nb132nf_vctot]
	movapd [rsp + nb132nf_vctot], xmm7

	;# LJ table interaction
	movapd xmm4, [rsp + nb132nf_rsqOO]
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb132nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb132nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [rsi + rax*8]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 16]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 16]	;# G2
	movhpd xmm6, [rsi + rax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 24]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [rsp + nb132nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb132nf_Vvdwtot]
	movapd [rsp + nb132nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rax*8 + 32]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8 + 32]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 48]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 48]	;# G2
	movhpd xmm6, [rsi + rax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [rsp + nb132nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb132nf_Vvdwtot]
	movapd [rsp + nb132nf_Vvdwtot], xmm5

	;# Other O--H interaction 
	movapd xmm0, [rsp + nb132nf_rinvOH1]
	movapd xmm1, [rsp + nb132nf_rinvH1H1]
	addpd  xmm0, [rsp + nb132nf_rinvOH2]
	addpd  xmm1, [rsp + nb132nf_rinvH1H2]
	addpd  xmm0, [rsp + nb132nf_rinvH1O]
	addpd  xmm1, [rsp + nb132nf_rinvH2H1]
	addpd  xmm0, [rsp + nb132nf_rinvH2O]
	addpd  xmm1, [rsp + nb132nf_rinvH2H2]
	mulpd  xmm0, [rsp + nb132nf_qqOH]
	mulpd  xmm1, [rsp + nb132nf_qqHH]
	addpd  xmm0, xmm1

	addpd  xmm0, [rsp + nb132nf_vctot]
	movapd [rsp + nb132nf_vctot], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb132nf_innerk],  2
	jl    .nb132nf_checksingle
	jmp   .nb132nf_unroll_loop
.nb132nf_checksingle:
	mov   edx, [rsp + nb132nf_innerk]
	and   edx, 1
	jnz   .nb132nf_dosingle
	jmp   .nb132nf_updateouterdata
.nb132nf_dosingle:
	mov   rdx, [rsp + nb132nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	
	mov rsi, [rbp + nb132nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movapd 	[rsp + nb132nf_jxO], xmm2
	movapd 	[rsp + nb132nf_jyO], xmm3
	movapd 	[rsp + nb132nf_jzO], xmm4
	movapd 	[rsp + nb132nf_jxH1], xmm5
	movapd 	[rsp + nb132nf_jyH1], xmm6
	movapd 	[rsp + nb132nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movapd 	[rsp + nb132nf_jxH2], xmm2
	movapd 	[rsp + nb132nf_jyH2], xmm3
	movapd 	[rsp + nb132nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb132nf_ixO]
	movapd xmm1, [rsp + nb132nf_iyO]
	movapd xmm2, [rsp + nb132nf_izO]
	movapd xmm3, [rsp + nb132nf_ixO]
	movapd xmm4, [rsp + nb132nf_iyO]
	movapd xmm5, [rsp + nb132nf_izO]
	subsd  xmm0, [rsp + nb132nf_jxO]
	subsd  xmm1, [rsp + nb132nf_jyO]
	subsd  xmm2, [rsp + nb132nf_jzO]
	subsd  xmm3, [rsp + nb132nf_jxH1]
	subsd  xmm4, [rsp + nb132nf_jyH1]
	subsd  xmm5, [rsp + nb132nf_jzH1]
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
	movapd [rsp + nb132nf_rsqOO], xmm0
	movapd [rsp + nb132nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb132nf_ixO]
	movapd xmm1, [rsp + nb132nf_iyO]
	movapd xmm2, [rsp + nb132nf_izO]
	movapd xmm3, [rsp + nb132nf_ixH1]
	movapd xmm4, [rsp + nb132nf_iyH1]
	movapd xmm5, [rsp + nb132nf_izH1]
	subsd  xmm0, [rsp + nb132nf_jxH2]
	subsd  xmm1, [rsp + nb132nf_jyH2]
	subsd  xmm2, [rsp + nb132nf_jzH2]
	subsd  xmm3, [rsp + nb132nf_jxO]
	subsd  xmm4, [rsp + nb132nf_jyO]
	subsd  xmm5, [rsp + nb132nf_jzO]
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
	movapd [rsp + nb132nf_rsqOH2], xmm0
	movapd [rsp + nb132nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb132nf_ixH1]
	movapd xmm1, [rsp + nb132nf_iyH1]
	movapd xmm2, [rsp + nb132nf_izH1]
	movapd xmm3, [rsp + nb132nf_ixH1]
	movapd xmm4, [rsp + nb132nf_iyH1]
	movapd xmm5, [rsp + nb132nf_izH1]
	subsd  xmm0, [rsp + nb132nf_jxH1]
	subsd  xmm1, [rsp + nb132nf_jyH1]
	subsd  xmm2, [rsp + nb132nf_jzH1]
	subsd  xmm3, [rsp + nb132nf_jxH2]
	subsd  xmm4, [rsp + nb132nf_jyH2]
	subsd  xmm5, [rsp + nb132nf_jzH2]
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
	movapd [rsp + nb132nf_rsqH1H1], xmm0
	movapd [rsp + nb132nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb132nf_ixH2]
	movapd xmm1, [rsp + nb132nf_iyH2]
	movapd xmm2, [rsp + nb132nf_izH2]
	movapd xmm3, [rsp + nb132nf_ixH2]
	movapd xmm4, [rsp + nb132nf_iyH2]
	movapd xmm5, [rsp + nb132nf_izH2]
	subsd  xmm0, [rsp + nb132nf_jxO]
	subsd  xmm1, [rsp + nb132nf_jyO]
	subsd  xmm2, [rsp + nb132nf_jzO]
	subsd  xmm3, [rsp + nb132nf_jxH1]
	subsd  xmm4, [rsp + nb132nf_jyH1]
	subsd  xmm5, [rsp + nb132nf_jzH1]
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
	movapd [rsp + nb132nf_rsqH2O], xmm0
	movapd [rsp + nb132nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb132nf_ixH2]
	movapd xmm1, [rsp + nb132nf_iyH2]
	movapd xmm2, [rsp + nb132nf_izH2]
	subsd  xmm0, [rsp + nb132nf_jxH2]
	subsd  xmm1, [rsp + nb132nf_jyH2]
	subsd  xmm2, [rsp + nb132nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [rsp + nb132nf_rsqH2H2], xmm0
		
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb132nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH2H2], xmm1
	movapd [rsp + nb132nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb132nf_rsqOO]
	movapd xmm4, [rsp + nb132nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb132nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb132nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb132nf_half] ;# rinv
	movapd [rsp + nb132nf_rinvOO], xmm1
	movapd [rsp + nb132nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb132nf_rsqOH2]
	movapd xmm4, [rsp + nb132nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb132nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvOH2], xmm1
	movapd [rsp + nb132nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb132nf_rsqH1H1]
	movapd xmm4, [rsp + nb132nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb132nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb132nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb132nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb132nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb132nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH1H1], xmm1
	movapd [rsp + nb132nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb132nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb132nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [rsp + nb132nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb132nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [rsp + nb132nf_half] ;# rinv 
	movapd [rsp + nb132nf_rinvH2O], xmm1
	
	;# start with OO interaction 	
	movsd xmm0, [rsp + nb132nf_rinvOO]
	movsd xmm7, xmm0		;# xmm7=rinv 
	mulsd xmm7, [rsp + nb132nf_qqOO]
	addsd  xmm7, [rsp + nb132nf_vctot]
	movsd [rsp + nb132nf_vctot], xmm7
	
	;# LJ table interaction
	movsd xmm4, [rsp + nb132nf_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb132nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  rsi, [rbp + nb132nf_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb132nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb132nf_Vvdwtot]
	movsd [rsp + nb132nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb132nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb132nf_Vvdwtot]
	movsd [rsp + nb132nf_Vvdwtot], xmm5

	;# Other O--H interaction 
	movsd xmm0, [rsp + nb132nf_rinvOH1]
	movsd xmm1, [rsp + nb132nf_rinvH1H1]
	addsd  xmm0, [rsp + nb132nf_rinvOH2]
	addsd  xmm1, [rsp + nb132nf_rinvH1H2]
	addsd  xmm0, [rsp + nb132nf_rinvH1O]
	addsd  xmm1, [rsp + nb132nf_rinvH2H1]
	addsd  xmm0, [rsp + nb132nf_rinvH2O]
	addsd  xmm1, [rsp + nb132nf_rinvH2H2]
	mulsd  xmm0, [rsp + nb132nf_qqOH]
	mulsd  xmm1, [rsp + nb132nf_qqHH]
	addsd  xmm0, xmm1

	addsd  xmm0, [rsp + nb132nf_vctot]
	movsd [rsp + nb132nf_vctot], xmm0
		
.nb132nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb132nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb132nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb132nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb132nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb132nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb132nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb132nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb132nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb132nf_n], esi
        jmp .nb132nf_outer
.nb132nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb132nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb132nf_end
        ;# non-zero, do one more workunit
        jmp   .nb132nf_threadloop
.nb132nf_end:
	mov eax, [rsp + nb132nf_nouter]
	mov ebx, [rsp + nb132nf_ninner]
	mov rcx, [rbp + nb132nf_outeriter]
	mov rdx, [rbp + nb132nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 864
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
