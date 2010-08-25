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


	
.globl nb_kernel204_x86_64_sse
.globl _nb_kernel204_x86_64_sse
nb_kernel204_x86_64_sse:	
_nb_kernel204_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb204_fshift,           16
.equiv          nb204_gid,              24
.equiv          nb204_pos,              32
.equiv          nb204_faction,          40
.equiv          nb204_charge,           48
.equiv          nb204_p_facel,          56
.equiv          nb204_argkrf,           64
.equiv          nb204_argcrf,           72
.equiv          nb204_Vc,               80
.equiv          nb204_type,             88
.equiv          nb204_p_ntype,          96
.equiv          nb204_vdwparam,         104
.equiv          nb204_Vvdw,             112
.equiv          nb204_p_tabscale,       120
.equiv          nb204_VFtab,            128
.equiv          nb204_invsqrta,         136
.equiv          nb204_dvda,             144
.equiv          nb204_p_gbtabscale,     152
.equiv          nb204_GBtab,            160
.equiv          nb204_p_nthreads,       168
.equiv          nb204_count,            176
.equiv          nb204_mtx,              184
.equiv          nb204_outeriter,        192
.equiv          nb204_inneriter,        200
.equiv          nb204_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb204_ixH1,             0
.equiv          nb204_iyH1,             16
.equiv          nb204_izH1,             32
.equiv          nb204_ixH2,             48
.equiv          nb204_iyH2,             64
.equiv          nb204_izH2,             80
.equiv          nb204_ixM,              96
.equiv          nb204_iyM,              112
.equiv          nb204_izM,              128
.equiv          nb204_jxH1,             144
.equiv          nb204_jyH1,             160
.equiv          nb204_jzH1,             176
.equiv          nb204_jxH2,             192
.equiv          nb204_jyH2,             208
.equiv          nb204_jzH2,             224
.equiv          nb204_jxM,              240
.equiv          nb204_jyM,              256
.equiv          nb204_jzM,              272
.equiv          nb204_dxH1H1,           288
.equiv          nb204_dyH1H1,           304
.equiv          nb204_dzH1H1,           320
.equiv          nb204_dxH1H2,           336
.equiv          nb204_dyH1H2,           352
.equiv          nb204_dzH1H2,           368
.equiv          nb204_dxH1M,            384
.equiv          nb204_dyH1M,            400
.equiv          nb204_dzH1M,            416
.equiv          nb204_dxH2H1,           432
.equiv          nb204_dyH2H1,           448
.equiv          nb204_dzH2H1,           464
.equiv          nb204_dxH2H2,           480
.equiv          nb204_dyH2H2,           496
.equiv          nb204_dzH2H2,           512
.equiv          nb204_dxH2M,            528
.equiv          nb204_dyH2M,            544
.equiv          nb204_dzH2M,            560
.equiv          nb204_dxMH1,            576
.equiv          nb204_dyMH1,            592
.equiv          nb204_dzMH1,            608
.equiv          nb204_dxMH2,            624
.equiv          nb204_dyMH2,            640
.equiv          nb204_dzMH2,            656
.equiv          nb204_dxMM,             672
.equiv          nb204_dyMM,             688
.equiv          nb204_dzMM,             704
.equiv          nb204_qqHH,             720
.equiv          nb204_qqMH,             736
.equiv          nb204_qqMM,             752
.equiv          nb204_vctot,            768
.equiv          nb204_fixH1,            784
.equiv          nb204_fiyH1,            800
.equiv          nb204_fizH1,            816
.equiv          nb204_fixH2,            832
.equiv          nb204_fiyH2,            848
.equiv          nb204_fizH2,            864
.equiv          nb204_fixM,             880
.equiv          nb204_fiyM,             896
.equiv          nb204_fizM,             912
.equiv          nb204_fjxH1,            928
.equiv          nb204_fjyH1,            944
.equiv          nb204_fjzH1,            960
.equiv          nb204_fjxH2,            976
.equiv          nb204_fjyH2,            992
.equiv          nb204_fjzH2,            1008
.equiv          nb204_fjxM,             1024
.equiv          nb204_fjyM,             1040
.equiv          nb204_fjzM,             1056
.equiv          nb204_half,             1072
.equiv          nb204_three,            1088
.equiv          nb204_rsqH1H1,          1104
.equiv          nb204_rsqH1H2,          1120
.equiv          nb204_rsqH1M,           1136
.equiv          nb204_rsqH2H1,          1152
.equiv          nb204_rsqH2H2,          1168
.equiv          nb204_rsqH2M,           1184
.equiv          nb204_rsqMH1,           1200
.equiv          nb204_rsqMH2,           1216
.equiv          nb204_rsqMM,            1232
.equiv          nb204_rinvH1H1,         1248
.equiv          nb204_rinvH1H2,         1264
.equiv          nb204_rinvH1M,          1280
.equiv          nb204_rinvH2H1,         1296
.equiv          nb204_rinvH2H2,         1312
.equiv          nb204_rinvH2M,          1328
.equiv          nb204_rinvMH1,          1344
.equiv          nb204_rinvMH2,          1360
.equiv          nb204_rinvMM,           1376
.equiv          nb204_two,              1392
.equiv          nb204_krf,              1408
.equiv          nb204_crf,              1424
.equiv          nb204_is3,              1440
.equiv          nb204_ii3,              1444
.equiv          nb204_innerjjnr,        1448
.equiv          nb204_nri,              1456
.equiv          nb204_iinr,             1464
.equiv          nb204_jindex,           1472
.equiv          nb204_jjnr,             1480
.equiv          nb204_shift,            1488
.equiv          nb204_shiftvec,         1496
.equiv          nb204_facel,            1504
.equiv          nb204_innerk,           1512
.equiv          nb204_n,                1516
.equiv          nb204_nn1,              1520
.equiv          nb204_nouter,           1524
.equiv          nb204_ninner,           1528

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
    
	emms
	sub rsp, 1536		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb204_nouter], eax
	mov [rsp + nb204_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb204_nri], edi
	mov [rsp + nb204_iinr], rsi
	mov [rsp + nb204_jindex], rdx
	mov [rsp + nb204_jjnr], rcx
	mov [rsp + nb204_shift], r8
	mov [rsp + nb204_shiftvec], r9
	mov rsi, [rbp + nb204_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb204_facel], xmm0


	mov rsi, [rbp + nb204_argkrf]
	mov rdi, [rbp + nb204_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb204_krf], xmm1
	movaps [rsp + nb204_crf], xmm2
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb204_half], eax
	movss xmm1, [rsp + nb204_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb204_half],  xmm1
	movaps [rsp + nb204_two],  xmm2
	movaps [rsp + nb204_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb204_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb204_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb204_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb204_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_qqHH], xmm3
	movaps [rsp + nb204_qqMH], xmm4
	movaps [rsp + nb204_qqMM], xmm5
	
.nb204_threadloop:
        mov   rsi, [rbp + nb204_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb204_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb204_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb204_n], eax
        mov [rsp + nb204_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204_outerstart
        jmp .nb204_end
	
.nb204_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb204_nouter]
	mov [rsp + nb204_nouter], ebx

.nb204_outer:
	mov   rax, [rsp + nb204_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb204_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb204_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb204_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb204_pos]	;# rax = base of pos[]  
	mov   [rsp + nb204_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_ixH1], xmm3
	movaps [rsp + nb204_iyH1], xmm4
	movaps [rsp + nb204_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204_ixH2], xmm0
	movaps [rsp + nb204_iyH2], xmm1
	movaps [rsp + nb204_izH2], xmm2
	movaps [rsp + nb204_ixM], xmm3
	movaps [rsp + nb204_iyM], xmm4
	movaps [rsp + nb204_izM], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb204_vctot], xmm4
	movaps [rsp + nb204_fixH1], xmm4
	movaps [rsp + nb204_fiyH1], xmm4
	movaps [rsp + nb204_fizH1], xmm4
	movaps [rsp + nb204_fixH2], xmm4
	movaps [rsp + nb204_fiyH2], xmm4
	movaps [rsp + nb204_fizH2], xmm4
	movaps [rsp + nb204_fixM], xmm4
	movaps [rsp + nb204_fiyM], xmm4
	movaps [rsp + nb204_fizM], xmm4
	
	mov   rax, [rsp + nb204_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb204_pos]
	mov   rdi, [rbp + nb204_faction]	
	mov   rax, [rsp + nb204_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb204_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb204_ninner]
	mov   [rsp + nb204_ninner], ecx
	add   edx, 0
	mov   [rsp + nb204_innerk], edx	;# number of innerloop atoms 
	jge   .nb204_unroll_loop
	jmp   .nb204_single_check
.nb204_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb204_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb204_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb204_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j H1 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb204_ixH1]
    subps xmm1, [rsp + nb204_iyH1]
    subps xmm2, [rsp + nb204_izH1]
    subps xmm3, [rsp + nb204_ixH2]
    subps xmm4, [rsp + nb204_iyH2]
    subps xmm5, [rsp + nb204_izH2]
    subps xmm6, [rsp + nb204_ixM]
    subps xmm7, [rsp + nb204_iyM]
    subps xmm8, [rsp + nb204_izM]
    
	movaps [rsp + nb204_dxH1H1], xmm0
	movaps [rsp + nb204_dyH1H1], xmm1
	movaps [rsp + nb204_dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH2H1], xmm3
	movaps [rsp + nb204_dyH2H1], xmm4
	movaps [rsp + nb204_dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb204_dxMH1], xmm6
	movaps [rsp + nb204_dyMH1], xmm7
	movaps [rsp + nb204_dzMH1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb204_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb204_half]
	mulps   xmm9, xmm4  ;# rinvH1H1 
	mulps   xmm10, xmm4 ;# rinvH2H1
    mulps   xmm11, xmm4 ;# rinvMH1
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb204_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb204_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb204_qqHH]
    movaps xmm13, [rsp + nb204_qqMH]    
    mulps  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm2, [rsp + nb204_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps xmm15, xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j H1 forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm10, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm9, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm10, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm11, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm12, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm11, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm12, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm9: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm10: jxH1c jyH1c jxH1d jyH1d
    ;# xmm11: jzH1a jzH1b jzH1c jzH1d

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb204_dxH1H1]
	mulps xmm1, [rsp + nb204_dyH1H1]
	mulps xmm2, [rsp + nb204_dzH1H1]
	mulps xmm3, [rsp + nb204_dxH2H1]
	mulps xmm4, [rsp + nb204_dyH2H1]
	mulps xmm5, [rsp + nb204_dzH2H1]
	mulps xmm6, [rsp + nb204_dxMH1]
	mulps xmm7, [rsp + nb204_dyMH1]
	mulps xmm8, [rsp + nb204_dzMH1]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb204_fixH1]
    addps xmm1, [rsp + nb204_fiyH1]
    addps xmm2, [rsp + nb204_fizH1]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb204_fixH2]
    addps xmm4, [rsp + nb204_fiyH2]
    addps xmm5, [rsp + nb204_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb204_fixM]
    addps xmm7, [rsp + nb204_fiyM]
    addps xmm8, [rsp + nb204_fizM]

    movaps [rsp + nb204_fixH1], xmm0
    movaps [rsp + nb204_fiyH1], xmm1
    movaps [rsp + nb204_fizH1], xmm2
    movaps [rsp + nb204_fixH2], xmm3
    movaps [rsp + nb204_fiyH2], xmm4
    movaps [rsp + nb204_fizH2], xmm5
    movaps [rsp + nb204_fixM], xmm6
    movaps [rsp + nb204_fiyM], xmm7
    movaps [rsp + nb204_fizM], xmm8
    
    ;# xmm9 = fH1x
    ;# xmm10 = fH1y
    ;# xmm11 = fH1z
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fH1zc fH1zd
    
    movlps [rdi + rax*4 + 12], xmm9
    movhps [rdi + rbx*4 + 12], xmm9
    movlps [rdi + rcx*4 + 12], xmm10
    movhps [rdi + rdx*4 + 12], xmm10
    movss  [rdi + rax*4 + 20], xmm11
    movss  [rdi + rcx*4 + 20], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 20], xmm11
    movss  [rdi + rdx*4 + 20], xmm12

	;# move j H2 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm3, [rsi + rdx*4 + 32] ;# jzH2c  -  jzH2d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    
    subps xmm0, [rsp + nb204_ixH1]
    subps xmm1, [rsp + nb204_iyH1]
    subps xmm2, [rsp + nb204_izH1]
    subps xmm3, [rsp + nb204_ixH2]
    subps xmm4, [rsp + nb204_iyH2]
    subps xmm5, [rsp + nb204_izH2]
    subps xmm6, [rsp + nb204_ixM]
    subps xmm7, [rsp + nb204_iyM]
    subps xmm8, [rsp + nb204_izM]
    
	movaps [rsp + nb204_dxH1H2], xmm0
	movaps [rsp + nb204_dyH1H2], xmm1
	movaps [rsp + nb204_dzH1H2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH2H2], xmm3
	movaps [rsp + nb204_dyH2H2], xmm4
	movaps [rsp + nb204_dzH2H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb204_dxMH2], xmm6
	movaps [rsp + nb204_dyMH2], xmm7
	movaps [rsp + nb204_dzMH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb204_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb204_half]
	mulps   xmm9, xmm4  ;# rinvH1H2
	mulps   xmm10, xmm4 ;# rinvH2H2
    mulps   xmm11, xmm4 ;# rinvMH2
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb204_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb204_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb204_qqHH]
    movaps xmm13, [rsp + nb204_qqMH]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm15, xmm2
    addps  xmm5, xmm8
    addps  xmm15, xmm5
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j H2 forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm10, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm9, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm10, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm11, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm12, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movhps xmm11, [rdi + rbx*4 + 32] ;# jzH2a  -  jzH2b  -
    movhps xmm12, [rdi + rdx*4 + 32] ;# jzH2c  -  jzH2d -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm9: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm10: jxH2c jyH2c jxH2d jyH2d
    ;# xmm11: jzH2a jzH2b jzH2c jzH2d

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb204_dxH1H2]
	mulps xmm1, [rsp + nb204_dyH1H2]
	mulps xmm2, [rsp + nb204_dzH1H2]
	mulps xmm3, [rsp + nb204_dxH2H2]
	mulps xmm4, [rsp + nb204_dyH2H2]
	mulps xmm5, [rsp + nb204_dzH2H2]
	mulps xmm6, [rsp + nb204_dxMH2]
	mulps xmm7, [rsp + nb204_dyMH2]
	mulps xmm8, [rsp + nb204_dzMH2]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb204_fixH1]
    addps xmm1, [rsp + nb204_fiyH1]
    addps xmm2, [rsp + nb204_fizH1]

    addps xmm13, xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb204_fixH2]
    addps xmm4, [rsp + nb204_fiyH2]
    addps xmm5, [rsp + nb204_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb204_fixM]
    addps xmm7, [rsp + nb204_fiyM]
    addps xmm8, [rsp + nb204_fizM]

    movaps [rsp + nb204_fixH1], xmm0
    movaps [rsp + nb204_fiyH1], xmm1
    movaps [rsp + nb204_fizH1], xmm2
    movaps [rsp + nb204_fixH2], xmm3
    movaps [rsp + nb204_fiyH2], xmm4
    movaps [rsp + nb204_fizH2], xmm5
    movaps [rsp + nb204_fixM], xmm6
    movaps [rsp + nb204_fiyM], xmm7
    movaps [rsp + nb204_fizM], xmm8
    
    ;# xmm9  = fH2x
    ;# xmm10 = fH2y
    ;# xmm11 = fH2z
    movaps xmm0, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm0, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm0

    movhlps  xmm12, xmm11 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm9
    movhps [rdi + rbx*4 + 24], xmm9
    movlps [rdi + rcx*4 + 24], xmm10
    movhps [rdi + rdx*4 + 24], xmm10
    movss  [rdi + rax*4 + 32], xmm11
    movss  [rdi + rcx*4 + 32], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 32], xmm11
    movss  [rdi + rdx*4 + 32], xmm12

	;# move j M coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm1, [rsi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm0, [rsi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm1, [rsi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm2, [rsi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm5, [rsi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm6, [rsi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm2, xmm5 ;# jzMa  -  jzMb  -
    movlhps xmm3, xmm6 ;# jzMc  -  jzMd -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxMa jxMc jyMa jyMc        
    unpckhps xmm4, xmm1  ;# jxMb jxMd jyMb jyMd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb204_ixH1]
    subps xmm1, [rsp + nb204_iyH1]
    subps xmm2, [rsp + nb204_izH1]
    subps xmm3, [rsp + nb204_ixH2]
    subps xmm4, [rsp + nb204_iyH2]
    subps xmm5, [rsp + nb204_izH2]
    subps xmm6, [rsp + nb204_ixM]
    subps xmm7, [rsp + nb204_iyM]
    subps xmm8, [rsp + nb204_izM]
    
	movaps [rsp + nb204_dxH1M], xmm0
	movaps [rsp + nb204_dyH1M], xmm1
	movaps [rsp + nb204_dzH1M], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb204_dxH2M], xmm3
	movaps [rsp + nb204_dyH2M], xmm4
	movaps [rsp + nb204_dzH2M], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb204_dxMM], xmm6
	movaps [rsp + nb204_dyMM], xmm7
	movaps [rsp + nb204_dzMM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jM atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb204_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb204_half]
	mulps   xmm9, xmm4  ;# rinvH1M
	mulps   xmm10, xmm4 ;# rinvH2M
    mulps   xmm11, xmm4 ;# rinvMM
	
	;# M interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb204_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb204_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb204_qqMH]
    movaps xmm13, [rsp + nb204_qqMM]    
    mulps  xmm2, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# xmm6=voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm5, xmm8
    addps  xmm2, xmm15
    addps  xmm2, xmm5
    movaps [rsp + nb204_vctot], xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j M forces to local temp variables 
    movlps xmm9, [rdi + rax*4 + 36] ;# jxMa jyMa  -   -
    movlps xmm10, [rdi + rcx*4 + 36] ;# jxMc jyMc  -   -
    movhps xmm9, [rdi + rbx*4 + 36] ;# jxMa jyMa jxMb jyMb 
    movhps xmm10, [rdi + rdx*4 + 36] ;# jxMc jyMc jxMd jyMd 

    movss  xmm11, [rdi + rax*4 + 44] ;# jzMa  -  -  -
    movss  xmm12, [rdi + rcx*4 + 44] ;# jzMc  -  -  -
    movss  xmm2,  [rdi + rbx*4 + 44] ;# jzMb  -  -  -
    movss  xmm3,  [rdi + rdx*4 + 44] ;# jzMd  -  -  -
    movlhps xmm11, xmm2 ;# jzMa  -  jzMb  -
    movlhps xmm12, xmm3 ;# jzMc  -  jzMd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jzMa jzMb jzMc jzMd

    ;# xmm9: jxMa jyMa jxMb jyMb 
    ;# xmm10: jxMc jyMc jxMd jyMd
    ;# xmm11: jzMa jzMb jzMc jzMd

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb204_dxH1M]
	mulps xmm1, [rsp + nb204_dyH1M]
	mulps xmm2, [rsp + nb204_dzH1M]
	mulps xmm3, [rsp + nb204_dxH2M]
	mulps xmm4, [rsp + nb204_dyH2M]
	mulps xmm5, [rsp + nb204_dzH2M]
	mulps xmm6, [rsp + nb204_dxMM]
	mulps xmm7, [rsp + nb204_dyMM]
	mulps xmm8, [rsp + nb204_dzMM]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb204_fixH1]
    addps xmm1, [rsp + nb204_fiyH1]
    addps xmm2, [rsp + nb204_fizH1]

    addps xmm13, xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb204_fixH2]
    addps xmm4, [rsp + nb204_fiyH2]
    addps xmm5, [rsp + nb204_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb204_fixM]
    addps xmm7, [rsp + nb204_fiyM]
    addps xmm8, [rsp + nb204_fizM]

    movaps [rsp + nb204_fixH1], xmm0
    movaps [rsp + nb204_fiyH1], xmm1
    movaps [rsp + nb204_fizH1], xmm2
    movaps [rsp + nb204_fixH2], xmm3
    movaps [rsp + nb204_fiyH2], xmm4
    movaps [rsp + nb204_fizH2], xmm5
    movaps [rsp + nb204_fixM], xmm6
    movaps [rsp + nb204_fiyM], xmm7
    movaps [rsp + nb204_fizM], xmm8
    
    ;# xmm0 = fMx
    ;# xmm1 = fMy
    ;# xmm2 = fMz
    movaps xmm15, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm15, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm15

    movhlps  xmm12, xmm11 ;# fMzc fMzd
    
    movlps [rdi + rax*4 + 36], xmm9
    movhps [rdi + rbx*4 + 36], xmm9
    movlps [rdi + rcx*4 + 36], xmm10
    movhps [rdi + rdx*4 + 36], xmm10
    movss  [rdi + rax*4 + 44], xmm11
    movss  [rdi + rcx*4 + 44], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 44], xmm11
    movss  [rdi + rdx*4 + 44], xmm12
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb204_innerk],  4
	jl    .nb204_single_check
	jmp   .nb204_unroll_loop
.nb204_single_check:
	add dword ptr [rsp + nb204_innerk],  4
	jnz   .nb204_single_loop
	jmp   .nb204_updateouterdata
.nb204_single_loop:
	mov   rdx, [rsp + nb204_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb204_innerjjnr],  4	

	mov rsi, [rbp + nb204_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	movss xmm0, [rsi + rax*4 + 36]		;# jxM  -  -  -
	movss xmm1, [rsi + rax*4 + 40]		;# jyM  -  -  -
	movss xmm2, [rsi + rax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -

	movlhps xmm0, xmm6			;# xmm0 = jxM   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm1 = jyM   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm2 = jzM   0   jzH1 jzH2

	;# store all j coordinates in jM 
	movaps [rsp + nb204_jxM], xmm0
	movaps [rsp + nb204_jyM], xmm1
	movaps [rsp + nb204_jzM], xmm2
	subps  xmm0, [rsp + nb204_ixM]
	subps  xmm1, [rsp + nb204_iyM]
	subps  xmm2, [rsp + nb204_izM]
	movaps [rsp + nb204_dxMM], xmm0
	movaps [rsp + nb204_dyMM], xmm1
	movaps [rsp + nb204_dzMM], xmm2

	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [rsp + nb204_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6 	;# xmm7=krsq 
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	mulps   xmm7, [rsp + nb204_two]
	subps   xmm6, [rsp + nb204_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	;# xmm3=rinv-2*krsq 
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb204_qqMM]
	movhps  xmm4, [rsp + nb204_qqMH]

	mulps xmm6, xmm4	;# vcoul  
	mulps xmm3, xmm4	;# coul part of fs  


	addps   xmm6, [rsp + nb204_vctot]
	mulps   xmm0, xmm3	;# total fscal 
	movaps  [rsp + nb204_vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb204_dxMM]
	mulps   xmm1, [rsp + nb204_dyMM]
	mulps   xmm2, [rsp + nb204_dzMM]
	
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixM]
	addps   xmm1, [rsp + nb204_fiyM]
	addps   xmm2, [rsp + nb204_fizM]
	movaps  [rsp + nb204_fixM], xmm0
	movaps  [rsp + nb204_fiyM], xmm1
	movaps  [rsp + nb204_fizM], xmm2

	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb204_jxM]
    movaps  xmm1, [rsp + nb204_jyM]
    movaps  xmm2, [rsp + nb204_jzM]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
 
	subps   xmm0, [rsp + nb204_ixH1]
	subps   xmm1, [rsp + nb204_iyH1]
	subps   xmm2, [rsp + nb204_izH1]
	subps   xmm3, [rsp + nb204_ixH2] 
	subps   xmm4, [rsp + nb204_iyH2]
	subps   xmm5, [rsp + nb204_izH2] 

	movaps [rsp + nb204_dxH1M], xmm0
	movaps [rsp + nb204_dyH1M], xmm1
	movaps [rsp + nb204_dzH1M], xmm2
	movaps [rsp + nb204_dxH2M], xmm3
	movaps [rsp + nb204_dyH2M], xmm4
	movaps [rsp + nb204_dzH2M], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb204_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb204_krf] ;# krsq 
	mulps xmm4, [rsp + nb204_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb204_qqMH]
	movhps  xmm6, [rsp + nb204_qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [rsp + nb204_crf]
	subps   xmm4, [rsp + nb204_crf]
	mulps   xmm1, [rsp + nb204_two]
	mulps   xmm5, [rsp + nb204_two]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb204_vctot]
	movaps  [rsp + nb204_vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;# fscal 
	mulps   xmm7, xmm4	;# fscal 
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb204_dxH1M]
	mulps   xmm1, [rsp + nb204_dyH1M]
	mulps   xmm2, [rsp + nb204_dzH1M]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixH1]
	addps   xmm1, [rsp + nb204_fiyH1]
	addps   xmm2, [rsp + nb204_fizH1]
	movaps  [rsp + nb204_fixH1], xmm0
	movaps  [rsp + nb204_fiyH1], xmm1
	movaps  [rsp + nb204_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [rsp + nb204_dxH2M]
	mulps   xmm1, [rsp + nb204_dyH2M]
	mulps   xmm2, [rsp + nb204_dzH2M]
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb204_faction]
	movaps  [rsp + nb204_fjxM], xmm3
	movaps  [rsp + nb204_fjyM], xmm4
	movaps  [rsp + nb204_fjzM], xmm5
	addps   xmm0, [rsp + nb204_fixH2]
	addps   xmm1, [rsp + nb204_fiyH2]
	addps   xmm2, [rsp + nb204_fizH2]
	movaps  [rsp + nb204_fixH2], xmm0
	movaps  [rsp + nb204_fiyH2], xmm1
	movaps  [rsp + nb204_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4 + 36]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb204_fjxM]
	movaps  xmm4, [rsp + nb204_fjyM]
	movaps  xmm5, [rsp + nb204_fjzM]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, [rsi + rax*4 + 44]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 44], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4 + 36], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + nb204_innerk]
	jz    .nb204_updateouterdata
	jmp   .nb204_single_loop
.nb204_updateouterdata:
	mov   ecx, [rsp + nb204_ii3]
	mov   rdi, [rbp + nb204_faction]
	mov   rsi, [rbp + nb204_fshift]
	mov   edx, [rsp + nb204_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixH1]
	movaps xmm1, [rsp + nb204_fiyH1] 
	movaps xmm2, [rsp + nb204_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	
 
	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixH2]
	movaps xmm1, [rsp + nb204_fiyH2]
	movaps xmm2, [rsp + nb204_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate M i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb204_fixM]
	movaps xmm1, [rsp + nb204_fiyM]
	movaps xmm2, [rsp + nb204_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 36]
	movss  xmm4, [rdi + rcx*4 + 40]
	movss  xmm5, [rdi + rcx*4 + 44]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 36], xmm3
	movss  [rdi + rcx*4 + 40], xmm4
	movss  [rdi + rcx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb204_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb204_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb204_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb204_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb204_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb204_n], esi
        jmp .nb204_outer
.nb204_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb204_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204_end
        ;# non-zero, do one more workunit
        jmp   .nb204_threadloop
.nb204_end:

	mov eax, [rsp + nb204_nouter]
	mov ebx, [rsp + nb204_ninner]
	mov rcx, [rbp + nb204_outeriter]
	mov rdx, [rbp + nb204_inneriter]
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






	
.globl nb_kernel204nf_x86_64_sse
.globl _nb_kernel204nf_x86_64_sse
nb_kernel204nf_x86_64_sse:	
_nb_kernel204nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb204nf_fshift,         16
.equiv          nb204nf_gid,            24
.equiv          nb204nf_pos,            32
.equiv          nb204nf_faction,        40
.equiv          nb204nf_charge,         48
.equiv          nb204nf_p_facel,        56
.equiv          nb204nf_argkrf,         64
.equiv          nb204nf_argcrf,         72
.equiv          nb204nf_Vc,             80
.equiv          nb204nf_type,           88
.equiv          nb204nf_p_ntype,        96
.equiv          nb204nf_vdwparam,       104
.equiv          nb204nf_Vvdw,           112
.equiv          nb204nf_p_tabscale,     120
.equiv          nb204nf_VFtab,          128
.equiv          nb204nf_invsqrta,       136
.equiv          nb204nf_dvda,           144
.equiv          nb204nf_p_gbtabscale,   152
.equiv          nb204nf_GBtab,          160
.equiv          nb204nf_p_nthreads,     168
.equiv          nb204nf_count,          176
.equiv          nb204nf_mtx,            184
.equiv          nb204nf_outeriter,      192
.equiv          nb204nf_inneriter,      200
.equiv          nb204nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb204nf_ixH1,           0
.equiv          nb204nf_iyH1,           16
.equiv          nb204nf_izH1,           32
.equiv          nb204nf_ixH2,           48
.equiv          nb204nf_iyH2,           64
.equiv          nb204nf_izH2,           80
.equiv          nb204nf_ixM,            96
.equiv          nb204nf_iyM,            112
.equiv          nb204nf_izM,            128
.equiv          nb204nf_jxH1,           144
.equiv          nb204nf_jyH1,           160
.equiv          nb204nf_jzH1,           176
.equiv          nb204nf_jxH2,           192
.equiv          nb204nf_jyH2,           208
.equiv          nb204nf_jzH2,           224
.equiv          nb204nf_jxM,            240
.equiv          nb204nf_jyM,            256
.equiv          nb204nf_jzM,            272
.equiv          nb204nf_qqHH,           288
.equiv          nb204nf_qqMH,           304
.equiv          nb204nf_qqMM,           320
.equiv          nb204nf_vctot,          336
.equiv          nb204nf_half,           352
.equiv          nb204nf_three,          368
.equiv          nb204nf_rsqH1H1,        384
.equiv          nb204nf_rsqH1H2,        400
.equiv          nb204nf_rsqH1M,         416
.equiv          nb204nf_rsqH2H1,        432
.equiv          nb204nf_rsqH2H2,        448
.equiv          nb204nf_rsqH2M,         464
.equiv          nb204nf_rsqMH1,         480
.equiv          nb204nf_rsqMH2,         496
.equiv          nb204nf_rsqMM,          512
.equiv          nb204nf_rinvH1H1,       528
.equiv          nb204nf_rinvH1H2,       544
.equiv          nb204nf_rinvH1M,        560
.equiv          nb204nf_rinvH2H1,       576
.equiv          nb204nf_rinvH2H2,       592
.equiv          nb204nf_rinvH2M,        608
.equiv          nb204nf_rinvMH1,        624
.equiv          nb204nf_rinvMH2,        640
.equiv          nb204nf_rinvMM,         656
.equiv          nb204nf_krf,            672
.equiv          nb204nf_crf,            688
.equiv          nb204nf_is3,            704
.equiv          nb204nf_ii3,            708
.equiv          nb204nf_innerjjnr,      712
.equiv          nb204nf_nri,            720
.equiv          nb204nf_iinr,           728
.equiv          nb204nf_jindex,         736
.equiv          nb204nf_jjnr,           744
.equiv          nb204nf_shift,          752
.equiv          nb204nf_shiftvec,       760
.equiv          nb204nf_facel,          768
.equiv          nb204nf_innerk,         776
.equiv          nb204nf_n,              780
.equiv          nb204nf_nn1,            784
.equiv          nb204nf_nouter,         788
.equiv          nb204nf_ninner,         792

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
    
	emms
	sub rsp, 800		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb204nf_nouter], eax
	mov [rsp + nb204nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb204nf_nri], edi
	mov [rsp + nb204nf_iinr], rsi
	mov [rsp + nb204nf_jindex], rdx
	mov [rsp + nb204nf_jjnr], rcx
	mov [rsp + nb204nf_shift], r8
	mov [rsp + nb204nf_shiftvec], r9
	mov rsi, [rbp + nb204nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb204nf_facel], xmm0

	mov rsi, [rbp + nb204nf_argkrf]
	mov rdi, [rbp + nb204nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb204nf_krf], xmm1
	movaps [rsp + nb204nf_crf], xmm2

	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb204nf_half], eax
	movss xmm1, [rsp + nb204nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb204nf_half],  xmm1
	movaps [rsp + nb204nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb204nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb204nf_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb204nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb204nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_qqHH], xmm3
	movaps [rsp + nb204nf_qqMH], xmm4
	movaps [rsp + nb204nf_qqMM], xmm5
	
.nb204nf_threadloop:
        mov   rsi, [rbp + nb204nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb204nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb204nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb204nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb204nf_n], eax
        mov [rsp + nb204nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb204nf_outerstart
        jmp .nb204nf_end
	
.nb204nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb204nf_nouter]
	mov [rsp + nb204nf_nouter], ebx

.nb204nf_outer:
	mov   rax, [rsp + nb204nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb204nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb204nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb204nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb204nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb204nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_ixH1], xmm3
	movaps [rsp + nb204nf_iyH1], xmm4
	movaps [rsp + nb204nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb204nf_ixH2], xmm0
	movaps [rsp + nb204nf_iyH2], xmm1
	movaps [rsp + nb204nf_izH2], xmm2
	movaps [rsp + nb204nf_ixM], xmm3
	movaps [rsp + nb204nf_iyM], xmm4
	movaps [rsp + nb204nf_izM], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb204nf_vctot], xmm4
	
	mov   rax, [rsp + nb204nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb204nf_pos]
	mov   rdi, [rbp + nb204nf_faction]	
	mov   rax, [rsp + nb204nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb204nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb204nf_ninner]
	mov   [rsp + nb204nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb204nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb204nf_unroll_loop
	jmp   .nb204nf_single_check
.nb204nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb204nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb204nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4 + 12]
	movlps xmm3, [rsi + rax*4 + 24]
	movlps xmm4, [rsi + rax*4 + 36]

	movlps xmm5, [rsi + rbx*4 + 12]
	movlps xmm6, [rsi + rbx*4 + 24]
	movlps xmm7, [rsi + rbx*4 + 36]

	movhps xmm2, [rsi + rcx*4 + 12]
	movhps xmm3, [rsi + rcx*4 + 24]
	movhps xmm4, [rsi + rcx*4 + 36]

	movhps xmm5, [rsi + rdx*4 + 12]
	movhps xmm6, [rsi + rdx*4 + 24]
	movhps xmm7, [rsi + rdx*4 + 36]
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5
	unpcklps xmm1, xmm6
	unpckhps xmm2, xmm5
	unpckhps xmm3, xmm6
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7
	unpckhps xmm5, xmm7
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2
	movaps [rsp + nb204nf_jxH1], xmm0
	movhlps  xmm2, xmm6
	movaps [rsp + nb204nf_jyH1], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb204nf_jxH2], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb204nf_jyH2], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb204nf_jxM], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb204nf_jyM], xmm5

	movss  xmm0, [rsi + rax*4 + 20]
	movss  xmm1, [rsi + rax*4 + 32]
	movss  xmm2, [rsi + rax*4 + 44]

	movss  xmm3, [rsi + rcx*4 + 20]
	movss  xmm4, [rsi + rcx*4 + 32]
	movss  xmm5, [rsi + rcx*4 + 44]

	movhps xmm0, [rsi + rbx*4 + 16]
	movhps xmm1, [rsi + rbx*4 + 28]
	movhps xmm2, [rsi + rbx*4 + 40]
	
	movhps xmm3, [rsi + rdx*4 + 16]
	movhps xmm4, [rsi + rdx*4 + 28]
	movhps xmm5, [rsi + rdx*4 + 40]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb204nf_jzH1],  xmm0
	movaps [rsp + nb204nf_jzH2],  xmm1
	movaps [rsp + nb204nf_jzM],  xmm2

	movaps xmm0, [rsp + nb204nf_ixH1]
	movaps xmm1, [rsp + nb204nf_iyH1]
	movaps xmm2, [rsp + nb204nf_izH1]
	movaps xmm3, [rsp + nb204nf_ixH1]
	movaps xmm4, [rsp + nb204nf_iyH1]
	movaps xmm5, [rsp + nb204nf_izH1]
	subps  xmm0, [rsp + nb204nf_jxH1]
	subps  xmm1, [rsp + nb204nf_jyH1]
	subps  xmm2, [rsp + nb204nf_jzH1]
	subps  xmm3, [rsp + nb204nf_jxH2]
	subps  xmm4, [rsp + nb204nf_jyH2]
	subps  xmm5, [rsp + nb204nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204nf_rsqH1H1], xmm0
	movaps [rsp + nb204nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb204nf_ixH1]
	movaps xmm1, [rsp + nb204nf_iyH1]
	movaps xmm2, [rsp + nb204nf_izH1]
	movaps xmm3, [rsp + nb204nf_ixH2]
	movaps xmm4, [rsp + nb204nf_iyH2]
	movaps xmm5, [rsp + nb204nf_izH2]
	subps  xmm0, [rsp + nb204nf_jxM]
	subps  xmm1, [rsp + nb204nf_jyM]
	subps  xmm2, [rsp + nb204nf_jzM]
	subps  xmm3, [rsp + nb204nf_jxH1]
	subps  xmm4, [rsp + nb204nf_jyH1]
	subps  xmm5, [rsp + nb204nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204nf_rsqH1M], xmm0
	movaps [rsp + nb204nf_rsqH2H1], xmm3

	movaps xmm0, [rsp + nb204nf_ixH2]
	movaps xmm1, [rsp + nb204nf_iyH2]
	movaps xmm2, [rsp + nb204nf_izH2]
	movaps xmm3, [rsp + nb204nf_ixH2]
	movaps xmm4, [rsp + nb204nf_iyH2]
	movaps xmm5, [rsp + nb204nf_izH2]
	subps  xmm0, [rsp + nb204nf_jxH2]
	subps  xmm1, [rsp + nb204nf_jyH2]
	subps  xmm2, [rsp + nb204nf_jzH2]
	subps  xmm3, [rsp + nb204nf_jxM]
	subps  xmm4, [rsp + nb204nf_jyM]
	subps  xmm5, [rsp + nb204nf_jzM]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb204nf_rsqH2H2], xmm0
	movaps [rsp + nb204nf_rsqH2M], xmm3

	movaps xmm0, [rsp + nb204nf_ixM]
	movaps xmm1, [rsp + nb204nf_iyM]
	movaps xmm2, [rsp + nb204nf_izM]
	movaps xmm3, [rsp + nb204nf_ixM]
	movaps xmm4, [rsp + nb204nf_iyM]
	movaps xmm5, [rsp + nb204nf_izM]
	subps  xmm0, [rsp + nb204nf_jxH1]
	subps  xmm1, [rsp + nb204nf_jyH1]
	subps  xmm2, [rsp + nb204nf_jzH1]
	subps  xmm3, [rsp + nb204nf_jxH2]
	subps  xmm4, [rsp + nb204nf_jyH2]
	subps  xmm5, [rsp + nb204nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb204nf_rsqMH1], xmm0
	movaps [rsp + nb204nf_rsqMH2], xmm4

	movaps xmm0, [rsp + nb204nf_ixM]
	movaps xmm1, [rsp + nb204nf_iyM]
	movaps xmm2, [rsp + nb204nf_izM]
	subps  xmm0, [rsp + nb204nf_jxM]
	subps  xmm1, [rsp + nb204nf_jyM]
	subps  xmm2, [rsp + nb204nf_jzM]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb204nf_rsqMM], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb204nf_half] ;# rinvH2H1 
	movaps  [rsp + nb204nf_rinvMM], xmm3
	movaps  [rsp + nb204nf_rinvMH2], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb204nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH1H1]
	mulps   xmm5, [rsp + nb204nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH1H1], xmm3
	movaps  [rsp + nb204nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH1M]
	rsqrtps xmm5, [rsp + nb204nf_rsqH2H1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH1M]
	mulps   xmm5, [rsp + nb204nf_rsqH2H1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH1M], xmm3
	movaps  [rsp + nb204nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqH2H2]
	rsqrtps xmm5, [rsp + nb204nf_rsqH2M]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb204nf_rsqH2H2]
	mulps   xmm5, [rsp + nb204nf_rsqH2M]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] 
	mulps   xmm7, [rsp + nb204nf_half]
	movaps  [rsp + nb204nf_rinvH2H2], xmm3
	movaps  [rsp + nb204nf_rinvH2M], xmm7
	
	rsqrtps xmm1, [rsp + nb204nf_rsqMH1]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204nf_three]
	mulps   xmm1, [rsp + nb204nf_rsqMH1]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204nf_half] 
	movaps  [rsp + nb204nf_rinvMH1], xmm3

	;# Coulomb interactions 
	;# add all H-H rsq in xmm2, H-M rsq xmm4
	;# H-H rinv in xmm0, H-M in xmm1
	movaps xmm0, [rsp + nb204nf_rinvH1H1]
	movaps xmm1, [rsp + nb204nf_rinvH1M]
	movaps xmm2, [rsp + nb204nf_rsqH1H1]
	movaps xmm4, [rsp + nb204nf_rsqH1M]
	addps  xmm0, [rsp + nb204nf_rinvH1H2]
	addps  xmm1, [rsp + nb204nf_rinvH2M]
	addps  xmm2, [rsp + nb204nf_rsqH1H2]
	addps  xmm4, [rsp + nb204nf_rsqH2M]
	addps  xmm0, [rsp + nb204nf_rinvH2H1]
	addps  xmm1, [rsp + nb204nf_rinvMH1]
	addps  xmm2, [rsp + nb204nf_rsqH2H1]
	addps  xmm4, [rsp + nb204nf_rsqMH1]
	addps  xmm0, [rsp + nb204nf_rinvH2H2]
	addps  xmm1, [rsp + nb204nf_rinvMH2]
	addps  xmm2, [rsp + nb204nf_rsqH2H2]
	addps  xmm4, [rsp + nb204nf_rsqMH2]
	movaps xmm5, [rsp + nb204nf_krf]
	movaps xmm6, [rsp + nb204nf_crf]

	;# calc 4*crf in xmm7
	movaps xmm7, xmm6
	addps  xmm7, xmm7
	addps  xmm7, xmm7
	mulps  xmm2, xmm5 ;# H-H krsq
	mulps  xmm4, xmm5 ;# H-M krsq
	addps  xmm0, xmm2 ;# H-H rinv+krsq
	addps  xmm1, xmm4 ;# H-M rinv+krsq
	subps  xmm0, xmm7 ;# H-H rinv+krsq-crf
	subps  xmm1, xmm7 ;# H-M rinv+krsq-crf
	mulps  xmm0, [rsp + nb204nf_qqHH]
	mulps  xmm1, [rsp + nb204nf_qqMH]
	addps  xmm0, xmm1 
	addps  xmm0, [rsp + nb204nf_vctot]
	;# M-M interaction
	movaps xmm4, [rsp + nb204nf_rinvMM]
	mulps  xmm5, [rsp + nb204nf_rsqMM] ;# krsq
	addps  xmm5, xmm4                  ;# rinv+krsq
	subps xmm5, [rsp + nb204nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [rsp + nb204nf_qqMM]
	addps xmm5, xmm0
	movaps [rsp + nb204nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb204nf_innerk],  4
	jl    .nb204nf_single_check
	jmp   .nb204nf_unroll_loop
.nb204nf_single_check:
	add dword ptr [rsp + nb204nf_innerk],  4
	jnz   .nb204nf_single_loop
	jmp   .nb204nf_updateouterdata
.nb204nf_single_loop:
	mov   rdx, [rsp + nb204nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb204nf_innerjjnr],  4	

	mov rsi, [rbp + nb204nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [rsi + rax*4 + 36]		;# jxM  -  -  -
	movss xmm4, [rsi + rax*4 + 40]		;# jyM  -  -  -
	movss xmm5, [rsi + rax*4 + 44]		;# jzM  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb204nf_ixM]     
	movaps  xmm1, [rsp + nb204nf_iyM]
	movaps  xmm2, [rsp + nb204nf_izM]	
	movlhps xmm3, xmm6			;# xmm3 = jxM   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyM   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzM   0   jzH1 jzH2

	;# store all j coordinates in jM
	movaps [rsp + nb204nf_jxM], xmm3
	movaps [rsp + nb204nf_jyM], xmm4
	movaps [rsp + nb204nf_jzM], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 

	movaps xmm6, xmm0
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	mulps   xmm6, [rsp + nb204nf_krf] ;# xmm6=krsq 
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb204nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb204nf_half] ;# rinv iO - j water 
	
	addps   xmm6, xmm3	;# xmm6=rinv+ krsq 
	subps   xmm6, [rsp + nb204nf_crf] ;# xmm6=rinv+ krsq-crf 
	
	xorps   xmm4, xmm4
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb204nf_qqMM]
	movhps  xmm4, [rsp + nb204nf_qqMH]
	mulps xmm6, xmm4	;# vcoul  
	addps   xmm6, [rsp + nb204nf_vctot]
	movaps  [rsp + nb204nf_vctot], xmm6	
	
	;# done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb204nf_ixH1]
	movaps  xmm1, [rsp + nb204nf_iyH1]
	movaps  xmm2, [rsp + nb204nf_izH1]	
	movaps  xmm3, [rsp + nb204nf_ixH2] 
	movaps  xmm4, [rsp + nb204nf_iyH2] 
	movaps  xmm5, [rsp + nb204nf_izH2] 
	subps   xmm0, [rsp + nb204nf_jxM]
	subps   xmm1, [rsp + nb204nf_jyM]
	subps   xmm2, [rsp + nb204nf_jzM]
	subps   xmm3, [rsp + nb204nf_jxM]
	subps   xmm4, [rsp + nb204nf_jyM]
	subps   xmm5, [rsp + nb204nf_jzM]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb204nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb204nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb204nf_half] ;# rinv H2 - j water  

	mulps xmm0, [rsp + nb204nf_krf] ;# krsq 
	mulps xmm4, [rsp + nb204nf_krf] ;# krsq  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	movss   xmm6, [rsp + nb204nf_qqMH]
	movhps  xmm6, [rsp + nb204nf_qqHH]
	addps   xmm0, xmm3	;# krsq+ rinv 
	addps   xmm4, xmm7	;# krsq+ rinv 
	subps   xmm0, [rsp + nb204nf_crf]
	subps   xmm4, [rsp + nb204nf_crf]
	mulps   xmm0, xmm6	;# vcoul 
	mulps   xmm4, xmm6	;# vcoul 
	addps   xmm4, xmm0		
	addps   xmm4, [rsp + nb204nf_vctot]
	movaps  [rsp + nb204nf_vctot], xmm4
	dec dword ptr [rsp + nb204nf_innerk]
	jz    .nb204nf_updateouterdata
	jmp   .nb204nf_single_loop
.nb204nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb204nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb204nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb204nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb204nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb204nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb204nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb204nf_n], esi
        jmp .nb204nf_outer
.nb204nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb204nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb204nf_end
        ;# non-zero, do one more workunit
        jmp   .nb204nf_threadloop
.nb204nf_end:

	mov eax, [rsp + nb204nf_nouter]
	mov ebx, [rsp + nb204nf_ninner]
	mov rcx, [rbp + nb204nf_outeriter]
	mov rdx, [rbp + nb204nf_inneriter]
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
