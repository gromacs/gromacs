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

	
.globl nb_kernel134_x86_64_sse2
.globl _nb_kernel134_x86_64_sse2
nb_kernel134_x86_64_sse2:	
_nb_kernel134_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb134_fshift,           16
.equiv          nb134_gid,              24
.equiv          nb134_pos,              32
.equiv          nb134_faction,          40
.equiv          nb134_charge,           48
.equiv          nb134_p_facel,          56
.equiv          nb134_argkrf,           64
.equiv          nb134_argcrf,           72
.equiv          nb134_Vc,               80
.equiv          nb134_type,             88
.equiv          nb134_p_ntype,          96
.equiv          nb134_vdwparam,         104
.equiv          nb134_Vvdw,             112
.equiv          nb134_p_tabscale,       120
.equiv          nb134_VFtab,            128
.equiv          nb134_invsqrta,         136
.equiv          nb134_dvda,             144
.equiv          nb134_p_gbtabscale,     152
.equiv          nb134_GBtab,            160
.equiv          nb134_p_nthreads,       168
.equiv          nb134_count,            176
.equiv          nb134_mtx,              184
.equiv          nb134_outeriter,        192
.equiv          nb134_inneriter,        200
.equiv          nb134_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb134_ixO,              0
.equiv          nb134_iyO,              16
.equiv          nb134_izO,              32
.equiv          nb134_ixH1,             48
.equiv          nb134_iyH1,             64
.equiv          nb134_izH1,             80
.equiv          nb134_ixH2,             96
.equiv          nb134_iyH2,             112
.equiv          nb134_izH2,             128
.equiv          nb134_ixM,              144
.equiv          nb134_iyM,              160
.equiv          nb134_izM,              176
.equiv          nb134_jxO,              192
.equiv          nb134_jyO,              208
.equiv          nb134_jzO,              224
.equiv          nb134_jxH1,             240
.equiv          nb134_jyH1,             256
.equiv          nb134_jzH1,             272
.equiv          nb134_jxH2,             288
.equiv          nb134_jyH2,             304
.equiv          nb134_jzH2,             320
.equiv          nb134_jxM,              336
.equiv          nb134_jyM,              352
.equiv          nb134_jzM,              368
.equiv          nb134_dxOO,             384
.equiv          nb134_dyOO,             400
.equiv          nb134_dzOO,             416
.equiv          nb134_dxH1H1,           432
.equiv          nb134_dyH1H1,           448
.equiv          nb134_dzH1H1,           464
.equiv          nb134_dxH1H2,           480
.equiv          nb134_dyH1H2,           496
.equiv          nb134_dzH1H2,           512
.equiv          nb134_dxH1M,            528
.equiv          nb134_dyH1M,            544
.equiv          nb134_dzH1M,            560
.equiv          nb134_dxH2H1,           576
.equiv          nb134_dyH2H1,           592
.equiv          nb134_dzH2H1,           608
.equiv          nb134_dxH2H2,           624
.equiv          nb134_dyH2H2,           640
.equiv          nb134_dzH2H2,           656
.equiv          nb134_dxH2M,            672
.equiv          nb134_dyH2M,            688
.equiv          nb134_dzH2M,            704
.equiv          nb134_dxMH1,            720
.equiv          nb134_dyMH1,            736
.equiv          nb134_dzMH1,            752
.equiv          nb134_dxMH2,            768
.equiv          nb134_dyMH2,            784
.equiv          nb134_dzMH2,            800
.equiv          nb134_dxMM,             816
.equiv          nb134_dyMM,             832
.equiv          nb134_dzMM,             848
.equiv          nb134_qqMM,             864
.equiv          nb134_qqMH,             880
.equiv          nb134_qqHH,             896
.equiv          nb134_two,              912
.equiv          nb134_c6,               944
.equiv          nb134_c12,              960
.equiv          nb134_vctot,            976
.equiv          nb134_Vvdwtot,          992
.equiv          nb134_fixO,             1008
.equiv          nb134_fiyO,             1024
.equiv          nb134_fizO,             1040
.equiv          nb134_fixH1,            1056
.equiv          nb134_fiyH1,            1072
.equiv          nb134_fizH1,            1088
.equiv          nb134_fixH2,            1104
.equiv          nb134_fiyH2,            1120
.equiv          nb134_fizH2,            1136
.equiv          nb134_fixM,             1152
.equiv          nb134_fiyM,             1168
.equiv          nb134_fizM,             1184
.equiv          nb134_fjxO,             1200
.equiv          nb134_fjyO,             1216
.equiv          nb134_fjzO,             1232
.equiv          nb134_fjxH1,            1248
.equiv          nb134_fjyH1,            1264
.equiv          nb134_fjzH1,            1280
.equiv          nb134_fjxH2,            1296
.equiv          nb134_fjyH2,            1312
.equiv          nb134_fjzH2,            1328
.equiv          nb134_fjxM,             1344
.equiv          nb134_fjyM,             1360
.equiv          nb134_fjzM,             1376
.equiv          nb134_half,             1392
.equiv          nb134_three,            1408
.equiv          nb134_tsc,              1424
.equiv          nb134_fstmp,            1440
.equiv          nb134_rsqOO,            1456
.equiv          nb134_rsqH1H1,          1472
.equiv          nb134_rsqH1H2,          1488
.equiv          nb134_rsqH1M,           1504
.equiv          nb134_rsqH2H1,          1520
.equiv          nb134_rsqH2H2,          1536
.equiv          nb134_rsqH2M,           1552
.equiv          nb134_rsqMH1,           1568
.equiv          nb134_rsqMH2,           1584
.equiv          nb134_rsqMM,            1600
.equiv          nb134_rinvOO,           1616
.equiv          nb134_rinvH1H1,         1632
.equiv          nb134_rinvH1H2,         1648
.equiv          nb134_rinvH1M,          1664
.equiv          nb134_rinvH2H1,         1680
.equiv          nb134_rinvH2H2,         1696
.equiv          nb134_rinvH2M,          1712
.equiv          nb134_rinvMH1,          1728
.equiv          nb134_rinvMH2,          1744
.equiv          nb134_rinvMM,           1760
.equiv          nb134_krf,              1776
.equiv          nb134_crf,              1792
.equiv          nb134_is3,              1808
.equiv          nb134_ii3,              1812
.equiv          nb134_nri,              1816
.equiv          nb134_iinr,             1824
.equiv          nb134_jindex,           1832
.equiv          nb134_jjnr,             1840
.equiv          nb134_shift,            1848
.equiv          nb134_shiftvec,         1856
.equiv          nb134_facel,            1864
.equiv          nb134_innerjjnr,        1872
.equiv          nb134_innerk,           1880
.equiv          nb134_n,                1884
.equiv          nb134_nn1,              1888
.equiv          nb134_nouter,           1892
.equiv          nb134_ninner,           1896

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
	sub rsp, 1904		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb134_nouter], eax
	mov [rsp + nb134_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb134_nri], edi
	mov [rsp + nb134_iinr], rsi
	mov [rsp + nb134_jindex], rdx
	mov [rsp + nb134_jjnr], rcx
	mov [rsp + nb134_shift], r8
	mov [rsp + nb134_shiftvec], r9
	mov rsi, [rbp + nb134_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb134_facel], xmm0

	mov rax, [rbp + nb134_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb134_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb134_half], eax
	mov [rsp + nb134_half+4], ebx
	movsd xmm1, [rsp + nb134_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb134_half], xmm1
	movapd [rsp + nb134_two], xmm2
	movapd [rsp + nb134_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb134_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb134_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	

	movsd xmm6, [rsp + nb134_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb134_qqMM], xmm3
	movapd [rsp + nb134_qqMH], xmm4
	movapd [rsp + nb134_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb134_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb134_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb134_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb134_c6], xmm0
	movapd [rsp + nb134_c12], xmm1

.nb134_threadloop:
        mov   rsi, [rbp + nb134_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb134_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb134_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb134_n], eax
        mov [rsp + nb134_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134_outerstart
        jmp .nb134_end

.nb134_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb134_nouter]
	mov [rsp + nb134_nouter], ebx

.nb134_outer:
	mov   rax, [rsp + nb134_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb134_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb134_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb134_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb134_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb134_ii3], ebx		

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb134_ixO], xmm3
	movapd [rsp + nb134_iyO], xmm4
	movapd [rsp + nb134_izO], xmm5
	movapd [rsp + nb134_ixH1], xmm6
	movapd [rsp + nb134_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb134_izH1], xmm6
	movapd [rsp + nb134_ixH2], xmm0
	movapd [rsp + nb134_iyH2], xmm1
	movapd [rsp + nb134_izH2], xmm2
	movapd [rsp + nb134_ixM], xmm3
	movapd [rsp + nb134_iyM], xmm4
	movapd [rsp + nb134_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb134_vctot], xmm4
	movapd [rsp + nb134_Vvdwtot], xmm4
	movapd [rsp + nb134_fixO], xmm4
	movapd [rsp + nb134_fiyO], xmm4
	movapd [rsp + nb134_fizO], xmm4
	movapd [rsp + nb134_fixH1], xmm4
	movapd [rsp + nb134_fiyH1], xmm4
	movapd [rsp + nb134_fizH1], xmm4
	movapd [rsp + nb134_fixH2], xmm4
	movapd [rsp + nb134_fiyH2], xmm4
	movapd [rsp + nb134_fizH2], xmm4
	movapd [rsp + nb134_fixM], xmm4
	movapd [rsp + nb134_fiyM], xmm4
	movapd [rsp + nb134_fizM], xmm4
	
	mov   rax, [rsp + nb134_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb134_pos] 
	mov   rdi, [rbp + nb134_faction]	
	mov   rax, [rsp + nb134_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb134_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb134_ninner]
	mov   [rsp + nb134_ninner], ecx
	add   edx, 0
	mov   [rsp + nb134_innerk], edx    ;# number of innerloop atoms 
	jge   .nb134_unroll_loop
	jmp   .nb134_checksingle
.nb134_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb134_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb134_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb134_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# load j O coordinates
    movlpd xmm4, [rsi + rax*8] 
    movlpd xmm5, [rsi + rax*8 + 8] 
    movlpd xmm6, [rsi + rax*8 + 16] 
    movhpd xmm4, [rsi + rbx*8] 
    movhpd xmm5, [rsi + rbx*8 + 8] 
    movhpd xmm6, [rsi + rbx*8 + 16] 

    ;# xmm4 = Ox
    ;# xmm5 = Oy
    ;# xmm6 = Oz
        
    subpd xmm4, [rsp + nb134_ixO]
    subpd xmm5, [rsp + nb134_iyO]
    subpd xmm6, [rsp + nb134_izO]

    ;# store dx/dy/dz
    movapd xmm13, xmm4
    movapd xmm14, xmm5
    movapd xmm15, xmm6
    
    ;# square it
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	mulpd  xmm6, xmm6
    
   	addpd  xmm4, xmm5
	addpd  xmm4, xmm6
    ;# rsq in xmm4
    
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb134_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb134_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb134_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb134_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm2, xmm0	;# xmm0=iter2 of rinv (new lu) 
	
	mulpd xmm4, xmm2	;# xmm4=r 
	mulpd xmm4, [rsp + nb134_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4
    ;# xmm1=eps 
    ;# xmm2=rinv
    movapd xmm3, xmm4   ;# eps
	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb134_VFtab]
	movd r10d, mm6
	psrlq mm6, 32
	movd r11d, mm6

    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
	movapd xmm4, [rsi + r10*8]          ;# Y1d F1d	
	movapd xmm0, [rsi + r11*8]         ;# Y2d F2d 
	movapd xmm8, [rsi + r10*8 + 32]     ;# Y1r F1r 	
	movapd xmm3, [rsi + r11*8 + 32]	;# Y2r F2r 
	movapd xmm5, xmm4
	movapd xmm9, xmm8
	unpcklpd xmm4, xmm0	;# Y1d Y2d 
	unpckhpd xmm5, xmm0	;# F1d F2d 
	unpcklpd xmm8, xmm3	;# Y1r Y2r 
	unpckhpd xmm9, xmm3	;# F1r F2r 

	movapd xmm6, [rsi + r10*8 + 16]     ;# G1d H1d 	
	movapd xmm0, [rsi + r11*8 + 16]  	;# G2d H2d 
	movapd xmm10, [rsi + r10*8 + 48]	;# G1r H1r 	
	movapd xmm3, [rsi + r11*8 + 48]	    ;# G2r H2r 
	movapd xmm7, xmm6
	movapd xmm11, xmm10
	unpcklpd xmm6, xmm0	;# G1d G2d 
	unpckhpd xmm7, xmm0	;# H1d H2d 
	unpcklpd xmm10, xmm3	;# G1r G2r 
	unpckhpd xmm11, xmm3	;# H1r H2r 
	;# tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulpd  xmm7, xmm1    ;# Heps
    mulpd  xmm11, xmm1 
    mulpd  xmm6, xmm1   ;# Geps
    mulpd  xmm10, xmm1 
    mulpd  xmm7, xmm1   ;# Heps2
    mulpd  xmm11, xmm1 
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
    mulpd  xmm5, xmm1  ;# eps*Fp
    mulpd  xmm9, xmm1
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, [rsp + nb134_c6]  ;# VV*c6 = vnb6
    mulpd  xmm9, [rsp + nb134_c12]  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb134_Vvdwtot]
    movapd [rsp + nb134_Vvdwtot], xmm5
        
    mulpd  xmm7, [rsp + nb134_c6]   ;# FF*c6 = fnb6
    mulpd  xmm11, [rsp + nb134_c12]   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mulpd  xmm7, [rsp + nb134_tsc]
    mulpd  xmm7, xmm2
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7
    
    mulpd xmm13, xmm9
    mulpd xmm14, xmm9
    mulpd xmm15, xmm9
    
    movapd xmm0, [rsp + nb134_fixO]
    movapd xmm1, [rsp + nb134_fiyO]
    movapd xmm2, [rsp + nb134_fizO]
    
    ;# accumulate i forces
    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    movapd [rsp + nb134_fixO], xmm0
    movapd [rsp + nb134_fiyO], xmm1
    movapd [rsp + nb134_fizO], xmm2
    
	;# the fj's - start by accumulating forces from memory 
	mov   rdi, [rbp + nb134_faction]	
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]
	addpd xmm3, xmm13
	addpd xmm4, xmm14
	addpd xmm5, xmm15
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5
    ;# done with OO interaction
    
    ;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
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
    
    subpd xmm0, [rsp + nb134_ixH1]
    subpd xmm1, [rsp + nb134_iyH1]
    subpd xmm2, [rsp + nb134_izH1]
    subpd xmm3, [rsp + nb134_ixH2]
    subpd xmm4, [rsp + nb134_iyH2]
    subpd xmm5, [rsp + nb134_izH2]
    subpd xmm6, [rsp + nb134_ixM]
    subpd xmm7, [rsp + nb134_iyM]
    subpd xmm8, [rsp + nb134_izM]
    
	movapd [rsp + nb134_dxH1H1], xmm0
	movapd [rsp + nb134_dyH1H1], xmm1
	movapd [rsp + nb134_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb134_dxH2H1], xmm3
	movapd [rsp + nb134_dyH2H1], xmm4
	movapd [rsp + nb134_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb134_dxMH1], xmm6
	movapd [rsp + nb134_dyMH1], xmm7
	movapd [rsp + nb134_dzMH1], xmm8
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
		
	movapd  xmm9, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H1 
	mulpd   xmm10, xmm15 ;#   rinvH2H1
    mulpd   xmm11, xmm15 ;#   rinvMH1
	
	;# H1 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb134_qqHH] 
    mulpd  xmm1, [rsp + nb134_qqHH] 
    mulpd  xmm2, [rsp + nb134_qqMH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb134_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb134_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
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

	mulpd xmm7, [rsp + nb134_dxH1H1]
	mulpd xmm8, [rsp + nb134_dyH1H1]
	mulpd xmm9, [rsp + nb134_dzH1H1]
	mulpd xmm10, [rsp + nb134_dxH2H1]
	mulpd xmm11, [rsp + nb134_dyH2H1]
	mulpd xmm12, [rsp + nb134_dzH2H1]
	mulpd xmm13, [rsp + nb134_dxMH1]
	mulpd xmm14, [rsp + nb134_dyMH1]
	mulpd xmm15, [rsp + nb134_dzMH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb134_fixH1]
    addpd xmm8, [rsp + nb134_fiyH1]
    addpd xmm9, [rsp + nb134_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb134_fixH2]
    addpd xmm11, [rsp + nb134_fiyH2]
    addpd xmm12, [rsp + nb134_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb134_fixM]
    addpd xmm14, [rsp + nb134_fiyM]
    addpd xmm15, [rsp + nb134_fizM]

    movapd [rsp + nb134_fixH1], xmm7
    movapd [rsp + nb134_fiyH1], xmm8
    movapd [rsp + nb134_fizH1], xmm9
    movapd [rsp + nb134_fixH2], xmm10
    movapd [rsp + nb134_fiyH2], xmm11
    movapd [rsp + nb134_fizH2], xmm12
    movapd [rsp + nb134_fixM], xmm13
    movapd [rsp + nb134_fiyM], xmm14
    movapd [rsp + nb134_fizM], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
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
    
    subpd xmm0, [rsp + nb134_ixH1]
    subpd xmm1, [rsp + nb134_iyH1]
    subpd xmm2, [rsp + nb134_izH1]
    subpd xmm3, [rsp + nb134_ixH2]
    subpd xmm4, [rsp + nb134_iyH2]
    subpd xmm5, [rsp + nb134_izH2]
    subpd xmm6, [rsp + nb134_ixM]
    subpd xmm7, [rsp + nb134_iyM]
    subpd xmm8, [rsp + nb134_izM]
    
	movapd [rsp + nb134_dxH1H2], xmm0
	movapd [rsp + nb134_dyH1H2], xmm1
	movapd [rsp + nb134_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb134_dxH2H2], xmm3
	movapd [rsp + nb134_dyH2H2], xmm4
	movapd [rsp + nb134_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb134_dxMH2], xmm6
	movapd [rsp + nb134_dyMH2], xmm7
	movapd [rsp + nb134_dzMH2], xmm8
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
		
	movapd  xmm9, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH1H2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H2
	mulpd   xmm10, xmm15 ;#   rinvH2H2
    mulpd   xmm11, xmm15 ;#   rinvMH2
	
	;# H2 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb134_qqHH] 
    mulpd  xmm1, [rsp + nb134_qqHH] 
    mulpd  xmm2, [rsp + nb134_qqMH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb134_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb134_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
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

	mulpd xmm7, [rsp + nb134_dxH1H2]
	mulpd xmm8, [rsp + nb134_dyH1H2]
	mulpd xmm9, [rsp + nb134_dzH1H2]
	mulpd xmm10, [rsp + nb134_dxH2H2]
	mulpd xmm11, [rsp + nb134_dyH2H2]
	mulpd xmm12, [rsp + nb134_dzH2H2]
	mulpd xmm13, [rsp + nb134_dxMH2]
	mulpd xmm14, [rsp + nb134_dyMH2]
	mulpd xmm15, [rsp + nb134_dzMH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb134_fixH1]
    addpd xmm8, [rsp + nb134_fiyH1]
    addpd xmm9, [rsp + nb134_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb134_fixH2]
    addpd xmm11, [rsp + nb134_fiyH2]
    addpd xmm12, [rsp + nb134_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb134_fixM]
    addpd xmm14, [rsp + nb134_fiyM]
    addpd xmm15, [rsp + nb134_fizM]

    movapd [rsp + nb134_fixH1], xmm7
    movapd [rsp + nb134_fiyH1], xmm8
    movapd [rsp + nb134_fizH1], xmm9
    movapd [rsp + nb134_fixH2], xmm10
    movapd [rsp + nb134_fiyH2], xmm11
    movapd [rsp + nb134_fizH2], xmm12
    movapd [rsp + nb134_fixM], xmm13
    movapd [rsp + nb134_fiyM], xmm14
    movapd [rsp + nb134_fizM], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
    movlpd xmm0, [rsi + rax*8 + 72] 
    movlpd xmm1, [rsi + rax*8 + 80] 
    movlpd xmm2, [rsi + rax*8 + 88] 
    movhpd xmm0, [rsi + rbx*8 + 72] 
    movhpd xmm1, [rsi + rbx*8 + 80] 
    movhpd xmm2, [rsi + rbx*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb134_ixH1]
    subpd xmm1, [rsp + nb134_iyH1]
    subpd xmm2, [rsp + nb134_izH1]
    subpd xmm3, [rsp + nb134_ixH2]
    subpd xmm4, [rsp + nb134_iyH2]
    subpd xmm5, [rsp + nb134_izH2]
    subpd xmm6, [rsp + nb134_ixM]
    subpd xmm7, [rsp + nb134_iyM]
    subpd xmm8, [rsp + nb134_izM]
    
	movapd [rsp + nb134_dxH1M], xmm0
	movapd [rsp + nb134_dyH1M], xmm1
	movapd [rsp + nb134_dzH1M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb134_dxH2M], xmm3
	movapd [rsp + nb134_dyH2M], xmm4
	movapd [rsp + nb134_dzH2M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb134_dxMM], xmm6
	movapd [rsp + nb134_dyMM], xmm7
	movapd [rsp + nb134_dzMM], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
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
		
	movapd  xmm9, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulpd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb134_three]
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

	movapd  xmm15, [rsp + nb134_half]
	mulpd   xmm9, xmm15  ;#  rinvH1M
	mulpd   xmm10, xmm15 ;#   rinvH2M
    mulpd   xmm11, xmm15 ;#   rinvMM
	
	;# M interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb134_qqMH] 
    mulpd  xmm1, [rsp + nb134_qqMH] 
    mulpd  xmm2, [rsp + nb134_qqMM] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb134_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb134_vctot], xmm0
    
    ;# move j M forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
	movlpd xmm0, [rdi + rax*8 + 72]
	movlpd xmm1, [rdi + rax*8 + 80]
	movlpd xmm2, [rdi + rax*8 + 88]
	movhpd xmm0, [rdi + rbx*8 + 72]
	movhpd xmm1, [rdi + rbx*8 + 80]
	movhpd xmm2, [rdi + rbx*8 + 88]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb134_dxH1M]
	mulpd xmm8, [rsp + nb134_dyH1M]
	mulpd xmm9, [rsp + nb134_dzH1M]
	mulpd xmm10, [rsp + nb134_dxH2M]
	mulpd xmm11, [rsp + nb134_dyH2M]
	mulpd xmm12, [rsp + nb134_dzH2M]
	mulpd xmm13, [rsp + nb134_dxMM]
	mulpd xmm14, [rsp + nb134_dyMM]
	mulpd xmm15, [rsp + nb134_dzMM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb134_fixH1]
    addpd xmm8, [rsp + nb134_fiyH1]
    addpd xmm9, [rsp + nb134_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb134_fixH2]
    addpd xmm11, [rsp + nb134_fiyH2]
    addpd xmm12, [rsp + nb134_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb134_fixM]
    addpd xmm14, [rsp + nb134_fiyM]
    addpd xmm15, [rsp + nb134_fizM]

    movapd [rsp + nb134_fixH1], xmm7
    movapd [rsp + nb134_fiyH1], xmm8
    movapd [rsp + nb134_fizH1], xmm9
    movapd [rsp + nb134_fixH2], xmm10
    movapd [rsp + nb134_fiyH2], xmm11
    movapd [rsp + nb134_fizH2], xmm12
    movapd [rsp + nb134_fixM], xmm13
    movapd [rsp + nb134_fiyM], xmm14
    movapd [rsp + nb134_fizM], xmm15
   
    ;# store back j M forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 72], xmm0
	movlpd [rdi + rax*8 + 80], xmm1
	movlpd [rdi + rax*8 + 88], xmm2
	movhpd [rdi + rbx*8 + 72], xmm0
	movhpd [rdi + rbx*8 + 80], xmm1
	movhpd [rdi + rbx*8 + 88], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb134_innerk],  2
	jl    .nb134_checksingle
	jmp   .nb134_unroll_loop
.nb134_checksingle:
	mov   edx, [rsp + nb134_innerk]
	and   edx, 1
	jnz   .nb134_dosingle
	jmp   .nb134_updateouterdata
.nb134_dosingle:
	mov   rdx, [rsp + nb134_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb134_pos]
	lea   rax, [rax + rax*2]  

	;# load j O coordinates
    movsd xmm4, [rsi + rax*8] 
    movsd xmm5, [rsi + rax*8 + 8] 
    movsd xmm6, [rsi + rax*8 + 16] 

    ;# xmm4 = Ox
    ;# xmm5 = Oy
    ;# xmm6 = Oz
        
    subsd xmm4, [rsp + nb134_ixO]
    subsd xmm5, [rsp + nb134_iyO]
    subsd xmm6, [rsp + nb134_izO]

    ;# store dx/dy/dz
    movapd xmm13, xmm4
    movapd xmm14, xmm5
    movapd xmm15, xmm6
    
    ;# square it
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	mulsd  xmm6, xmm6
    
   	addsd  xmm4, xmm5
	addsd  xmm4, xmm6
    ;# rsq in xmm4
    
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb134_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb134_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb134_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb134_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm2, xmm0	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm2	;# xmm4=r 
	mulsd xmm4, [rsp + nb134_tsc]
	
	cvttsd2si r10d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r10d
	subsd xmm4, xmm5
	movapd xmm1, xmm4
    ;# xmm1=eps 
    ;# xmm2=rinv
    movapd xmm3, xmm4   ;# eps
	shl r10d, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb134_VFtab]


    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
	movapd xmm4, [rsi + r10*8]          ;# Y1d F1d	
	movapd xmm8, [rsi + r10*8 + 32]     ;# Y1r F1r 	
	movhlps xmm5, xmm4
	movhlps xmm9, xmm8

	movapd xmm6, [rsi + r10*8 + 16]     ;# G1d H1d 	
	movapd xmm10, [rsi + r10*8 + 48]	;# G1r H1r 	
	movhlps xmm7, xmm6
	movhlps xmm11, xmm10
	;# tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulsd  xmm7, xmm1    ;# Heps
    mulsd  xmm11, xmm1 
    mulsd  xmm6, xmm1   ;# Geps
    mulsd  xmm10, xmm1 
    mulsd  xmm7, xmm1   ;# Heps2
    mulsd  xmm11, xmm1 
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
    mulsd  xmm5, xmm1  ;# eps*Fp
    mulsd  xmm9, xmm1
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, [rsp + nb134_c6]  ;# VV*c6 = vnb6
    mulsd  xmm9, [rsp + nb134_c12]  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb134_Vvdwtot]
    movsd [rsp + nb134_Vvdwtot], xmm5
        
    mulsd  xmm7, [rsp + nb134_c6]   ;# FF*c6 = fnb6
    mulsd  xmm11, [rsp + nb134_c12]   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    
    mulsd  xmm7, [rsp + nb134_tsc]
    mulsd  xmm7, xmm2
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7
    
    mulsd xmm13, xmm9
    mulsd xmm14, xmm9
    mulsd xmm15, xmm9
    
    movapd xmm0, [rsp + nb134_fixO]
    movapd xmm1, [rsp + nb134_fiyO]
    movapd xmm2, [rsp + nb134_fizO]
    
    ;# accumulate i forces
    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    movsd [rsp + nb134_fixO], xmm0
    movsd [rsp + nb134_fiyO], xmm1
    movsd [rsp + nb134_fizO], xmm2
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb134_faction]
	addsd xmm13, [rdi + rax*8]
	addsd xmm14, [rdi + rax*8 + 8]
	addsd xmm15, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm13
	movsd [rdi + rax*8 + 8], xmm14
	movsd [rdi + rax*8 + 16], xmm15
    ;# done with OO interaction
    
	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
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
    
    subsd xmm0, [rsp + nb134_ixH1]
    subsd xmm1, [rsp + nb134_iyH1]
    subsd xmm2, [rsp + nb134_izH1]
    subsd xmm3, [rsp + nb134_ixH2]
    subsd xmm4, [rsp + nb134_iyH2]
    subsd xmm5, [rsp + nb134_izH2]
    subsd xmm6, [rsp + nb134_ixM]
    subsd xmm7, [rsp + nb134_iyM]
    subsd xmm8, [rsp + nb134_izM]
    
	movsd [rsp + nb134_dxH1H1], xmm0
	movsd [rsp + nb134_dyH1H1], xmm1
	movsd [rsp + nb134_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb134_dxH2H1], xmm3
	movsd [rsp + nb134_dyH2H1], xmm4
	movsd [rsp + nb134_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb134_dxMH1], xmm6
	movsd [rsp + nb134_dyMH1], xmm7
	movsd [rsp + nb134_dzMH1], xmm8
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
		
	movsd  xmm9, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H1 
	mulsd   xmm10, xmm15 ;#   rinvH2H1
    mulsd   xmm11, xmm15 ;#   rinvMH1
	
	;# H1 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb134_qqHH] 
    mulsd  xmm1, [rsp + nb134_qqHH] 
    mulsd  xmm2, [rsp + nb134_qqMH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb134_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb134_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
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

	mulsd xmm7, [rsp + nb134_dxH1H1]
	mulsd xmm8, [rsp + nb134_dyH1H1]
	mulsd xmm9, [rsp + nb134_dzH1H1]
	mulsd xmm10, [rsp + nb134_dxH2H1]
	mulsd xmm11, [rsp + nb134_dyH2H1]
	mulsd xmm12, [rsp + nb134_dzH2H1]
	mulsd xmm13, [rsp + nb134_dxMH1]
	mulsd xmm14, [rsp + nb134_dyMH1]
	mulsd xmm15, [rsp + nb134_dzMH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb134_fixH1]
    addsd xmm8, [rsp + nb134_fiyH1]
    addsd xmm9, [rsp + nb134_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb134_fixH2]
    addsd xmm11, [rsp + nb134_fiyH2]
    addsd xmm12, [rsp + nb134_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb134_fixM]
    addsd xmm14, [rsp + nb134_fiyM]
    addsd xmm15, [rsp + nb134_fizM]

    movsd [rsp + nb134_fixH1], xmm7
    movsd [rsp + nb134_fiyH1], xmm8
    movsd [rsp + nb134_fizH1], xmm9
    movsd [rsp + nb134_fixH2], xmm10
    movsd [rsp + nb134_fiyH2], xmm11
    movsd [rsp + nb134_fizH2], xmm12
    movsd [rsp + nb134_fixM], xmm13
    movsd [rsp + nb134_fiyM], xmm14
    movsd [rsp + nb134_fizM], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
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
    
    subsd xmm0, [rsp + nb134_ixH1]
    subsd xmm1, [rsp + nb134_iyH1]
    subsd xmm2, [rsp + nb134_izH1]
    subsd xmm3, [rsp + nb134_ixH2]
    subsd xmm4, [rsp + nb134_iyH2]
    subsd xmm5, [rsp + nb134_izH2]
    subsd xmm6, [rsp + nb134_ixM]
    subsd xmm7, [rsp + nb134_iyM]
    subsd xmm8, [rsp + nb134_izM]
    
	movsd [rsp + nb134_dxH1H2], xmm0
	movsd [rsp + nb134_dyH1H2], xmm1
	movsd [rsp + nb134_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb134_dxH2H2], xmm3
	movsd [rsp + nb134_dyH2H2], xmm4
	movsd [rsp + nb134_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb134_dxMH2], xmm6
	movsd [rsp + nb134_dyMH2], xmm7
	movsd [rsp + nb134_dzMH2], xmm8
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
		
	movsd  xmm9, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH2

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H2
	mulsd   xmm10, xmm15 ;#   rinvH2H2
    mulsd   xmm11, xmm15 ;#   rinvMH2
	
	;# H2 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb134_qqHH] 
    mulsd  xmm1, [rsp + nb134_qqHH] 
    mulsd  xmm2, [rsp + nb134_qqMH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb134_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb134_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
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

	mulsd xmm7, [rsp + nb134_dxH1H2]
	mulsd xmm8, [rsp + nb134_dyH1H2]
	mulsd xmm9, [rsp + nb134_dzH1H2]
	mulsd xmm10, [rsp + nb134_dxH2H2]
	mulsd xmm11, [rsp + nb134_dyH2H2]
	mulsd xmm12, [rsp + nb134_dzH2H2]
	mulsd xmm13, [rsp + nb134_dxMH2]
	mulsd xmm14, [rsp + nb134_dyMH2]
	mulsd xmm15, [rsp + nb134_dzMH2]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb134_fixH1]
    addsd xmm8, [rsp + nb134_fiyH1]
    addsd xmm9, [rsp + nb134_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb134_fixH2]
    addsd xmm11, [rsp + nb134_fiyH2]
    addsd xmm12, [rsp + nb134_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb134_fixM]
    addsd xmm14, [rsp + nb134_fiyM]
    addsd xmm15, [rsp + nb134_fizM]

    movsd [rsp + nb134_fixH1], xmm7
    movsd [rsp + nb134_fiyH1], xmm8
    movsd [rsp + nb134_fizH1], xmm9
    movsd [rsp + nb134_fixH2], xmm10
    movsd [rsp + nb134_fiyH2], xmm11
    movsd [rsp + nb134_fizH2], xmm12
    movsd [rsp + nb134_fixM], xmm13
    movsd [rsp + nb134_fiyM], xmm14
    movsd [rsp + nb134_fizM], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
    mov rsi, [rbp + nb134_pos]
    movsd xmm0, [rsi + rax*8 + 72] 
    movsd xmm1, [rsi + rax*8 + 80] 
    movsd xmm2, [rsi + rax*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb134_ixH1]
    subsd xmm1, [rsp + nb134_iyH1]
    subsd xmm2, [rsp + nb134_izH1]
    subsd xmm3, [rsp + nb134_ixH2]
    subsd xmm4, [rsp + nb134_iyH2]
    subsd xmm5, [rsp + nb134_izH2]
    subsd xmm6, [rsp + nb134_ixM]
    subsd xmm7, [rsp + nb134_iyM]
    subsd xmm8, [rsp + nb134_izM]
    
	movsd [rsp + nb134_dxH1M], xmm0
	movsd [rsp + nb134_dyH1M], xmm1
	movsd [rsp + nb134_dzH1M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb134_dxH2M], xmm3
	movsd [rsp + nb134_dyH2M], xmm4
	movsd [rsp + nb134_dzH2M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb134_dxMM], xmm6
	movsd [rsp + nb134_dyMM], xmm7
	movsd [rsp + nb134_dzMM], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
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
		
	movsd  xmm9, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulsd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb134_three]
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

	movsd  xmm15, [rsp + nb134_half]
	mulsd   xmm9, xmm15  ;#  rinvH1M
	mulsd   xmm10, xmm15 ;#   rinvH2M
    mulsd   xmm11, xmm15 ;#   rinvMM
	
	;# M interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb134_qqMH] 
    mulsd  xmm1, [rsp + nb134_qqMH] 
    mulsd  xmm2, [rsp + nb134_qqMM] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb134_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb134_vctot], xmm0
    
    ;# move j M forces to xmm0-xmm2
    mov rdi, [rbp + nb134_faction]
	movsd xmm0, [rdi + rax*8 + 72]
	movsd xmm1, [rdi + rax*8 + 80]
	movsd xmm2, [rdi + rax*8 + 88]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb134_dxH1M]
	mulsd xmm8, [rsp + nb134_dyH1M]
	mulsd xmm9, [rsp + nb134_dzH1M]
	mulsd xmm10, [rsp + nb134_dxH2M]
	mulsd xmm11, [rsp + nb134_dyH2M]
	mulsd xmm12, [rsp + nb134_dzH2M]
	mulsd xmm13, [rsp + nb134_dxMM]
	mulsd xmm14, [rsp + nb134_dyMM]
	mulsd xmm15, [rsp + nb134_dzMM]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb134_fixH1]
    addsd xmm8, [rsp + nb134_fiyH1]
    addsd xmm9, [rsp + nb134_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb134_fixH2]
    addsd xmm11, [rsp + nb134_fiyH2]
    addsd xmm12, [rsp + nb134_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb134_fixM]
    addsd xmm14, [rsp + nb134_fiyM]
    addsd xmm15, [rsp + nb134_fizM]

    movsd [rsp + nb134_fixH1], xmm7
    movsd [rsp + nb134_fiyH1], xmm8
    movsd [rsp + nb134_fizH1], xmm9
    movsd [rsp + nb134_fixH2], xmm10
    movsd [rsp + nb134_fiyH2], xmm11
    movsd [rsp + nb134_fizH2], xmm12
    movsd [rsp + nb134_fixM], xmm13
    movsd [rsp + nb134_fiyM], xmm14
    movsd [rsp + nb134_fizM], xmm15
   
    ;# store back j M forces from xmm0-xmm2
	movsd [rdi + rax*8 + 72], xmm0
	movsd [rdi + rax*8 + 80], xmm1
	movsd [rdi + rax*8 + 88], xmm2
	
.nb134_updateouterdata:
	mov   ecx, [rsp + nb134_ii3]
	mov   rdi, [rbp + nb134_faction]
	mov   rsi, [rbp + nb134_fshift]
	mov   edx, [rsp + nb134_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb134_fixO]
	movapd xmm1, [rsp + nb134_fiyO]
	movapd xmm2, [rsp + nb134_fizO]

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
	movapd xmm0, [rsp + nb134_fixH1]
	movapd xmm1, [rsp + nb134_fiyH1]
	movapd xmm2, [rsp + nb134_fizH1]

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
	movapd xmm0, [rsp + nb134_fixH2]
	movapd xmm1, [rsp + nb134_fiyH2]
	movapd xmm2, [rsp + nb134_fizH2]

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

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb134_fixM]
	movapd xmm1, [rsp + nb134_fiyM]
	movapd xmm2, [rsp + nb134_fizM]

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
	movsd  xmm3, [rdi + rcx*8 + 72]
	movsd  xmm4, [rdi + rcx*8 + 80]
	movsd  xmm5, [rdi + rcx*8 + 88]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 72], xmm3
	movsd  [rdi + rcx*8 + 80], xmm4
	movsd  [rdi + rcx*8 + 88], xmm5

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
	mov esi, [rsp + nb134_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb134_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb134_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb134_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb134_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb134_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb134_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb134_n], esi
        jmp .nb134_outer
.nb134_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb134_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134_end
        ;# non-zero, do one more workunit
        jmp   .nb134_threadloop
.nb134_end:
	mov eax, [rsp + nb134_nouter]
	mov ebx, [rsp + nb134_ninner]
	mov rcx, [rbp + nb134_outeriter]
	mov rdx, [rbp + nb134_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1904
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


	
.globl nb_kernel134nf_x86_64_sse2
.globl _nb_kernel134nf_x86_64_sse2
nb_kernel134nf_x86_64_sse2:	
_nb_kernel134nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb134nf_fshift,         16
.equiv          nb134nf_gid,            24
.equiv          nb134nf_pos,            32
.equiv          nb134nf_faction,        40
.equiv          nb134nf_charge,         48
.equiv          nb134nf_p_facel,        56
.equiv          nb134nf_argkrf,         64
.equiv          nb134nf_argcrf,         72
.equiv          nb134nf_Vc,             80
.equiv          nb134nf_type,           88
.equiv          nb134nf_p_ntype,        96
.equiv          nb134nf_vdwparam,       104
.equiv          nb134nf_Vvdw,           112
.equiv          nb134nf_p_tabscale,     120
.equiv          nb134nf_VFtab,          128
.equiv          nb134nf_invsqrta,       136
.equiv          nb134nf_dvda,           144
.equiv          nb134nf_p_gbtabscale,   152
.equiv          nb134nf_GBtab,          160
.equiv          nb134nf_p_nthreads,     168
.equiv          nb134nf_count,          176
.equiv          nb134nf_mtx,            184
.equiv          nb134nf_outeriter,      192
.equiv          nb134nf_inneriter,      200
.equiv          nb134nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb134nf_ixO,            0
.equiv          nb134nf_iyO,            16
.equiv          nb134nf_izO,            32
.equiv          nb134nf_ixH1,           48
.equiv          nb134nf_iyH1,           64
.equiv          nb134nf_izH1,           80
.equiv          nb134nf_ixH2,           96
.equiv          nb134nf_iyH2,           112
.equiv          nb134nf_izH2,           128
.equiv          nb134nf_ixM,            144
.equiv          nb134nf_iyM,            160
.equiv          nb134nf_izM,            176
.equiv          nb134nf_jxO,            192
.equiv          nb134nf_jyO,            208
.equiv          nb134nf_jzO,            224
.equiv          nb134nf_jxH1,           240
.equiv          nb134nf_jyH1,           256
.equiv          nb134nf_jzH1,           272
.equiv          nb134nf_jxH2,           288
.equiv          nb134nf_jyH2,           304
.equiv          nb134nf_jzH2,           320
.equiv          nb134nf_jxM,            336
.equiv          nb134nf_jyM,            352
.equiv          nb134nf_jzM,            368
.equiv          nb134nf_qqMM,           384
.equiv          nb134nf_qqMH,           400
.equiv          nb134nf_qqHH,           416
.equiv          nb134nf_tsc,            432
.equiv          nb134nf_c6,             448
.equiv          nb134nf_c12,            464
.equiv          nb134nf_vctot,          480
.equiv          nb134nf_Vvdwtot,        496
.equiv          nb134nf_half,           512
.equiv          nb134nf_three,          528
.equiv          nb134nf_rsqOO,          544
.equiv          nb134nf_rsqH1H1,        560
.equiv          nb134nf_rsqH1H2,        576
.equiv          nb134nf_rsqH1M,         592
.equiv          nb134nf_rsqH2H1,        608
.equiv          nb134nf_rsqH2H2,        624
.equiv          nb134nf_rsqH2M,         640
.equiv          nb134nf_rsqMH1,         656
.equiv          nb134nf_rsqMH2,         672
.equiv          nb134nf_rsqMM,          688
.equiv          nb134nf_rinvOO,         704
.equiv          nb134nf_rinvH1H1,       720
.equiv          nb134nf_rinvH1H2,       736
.equiv          nb134nf_rinvH1M,        752
.equiv          nb134nf_rinvH2H1,       768
.equiv          nb134nf_rinvH2H2,       784
.equiv          nb134nf_rinvH2M,        800
.equiv          nb134nf_rinvMH1,        816
.equiv          nb134nf_rinvMH2,        832
.equiv          nb134nf_rinvMM,         848
.equiv          nb134nf_krf,            864
.equiv          nb134nf_crf,            880
.equiv          nb134nf_is3,            896
.equiv          nb134nf_ii3,            900
.equiv          nb134nf_nri,            904
.equiv          nb134nf_iinr,           912
.equiv          nb134nf_jindex,         920
.equiv          nb134nf_jjnr,           928
.equiv          nb134nf_shift,          936
.equiv          nb134nf_shiftvec,       944
.equiv          nb134nf_facel,          952
.equiv          nb134nf_innerjjnr,      960
.equiv          nb134nf_innerk,         968
.equiv          nb134nf_n,              972
.equiv          nb134nf_nn1,            976
.equiv          nb134nf_nouter,         980
.equiv          nb134nf_ninner,         984

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
	sub rsp, 992		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb134nf_nouter], eax
	mov [rsp + nb134nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb134nf_nri], edi
	mov [rsp + nb134nf_iinr], rsi
	mov [rsp + nb134nf_jindex], rdx
	mov [rsp + nb134nf_jjnr], rcx
	mov [rsp + nb134nf_shift], r8
	mov [rsp + nb134nf_shiftvec], r9
	mov rsi, [rbp + nb134nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb134nf_facel], xmm0

	mov rax, [rbp + nb134nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb134nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb134nf_half], eax
	mov [rsp + nb134nf_half+4], ebx
	movsd xmm1, [rsp + nb134nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb134nf_half], xmm1
	movapd [rsp + nb134nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb134nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb134nf_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	

	movsd xmm6, [rsp + nb134nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb134nf_qqMM], xmm3
	movapd [rsp + nb134nf_qqMH], xmm4
	movapd [rsp + nb134nf_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb134nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb134nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb134nf_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb134nf_c6], xmm0
	movapd [rsp + nb134nf_c12], xmm1

.nb134nf_threadloop:
        mov   rsi, [rbp + nb134nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb134nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb134nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb134nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb134nf_n], eax
        mov [rsp + nb134nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb134nf_outerstart
        jmp .nb134nf_end

.nb134nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb134nf_nouter]
	mov [rsp + nb134nf_nouter], ebx

.nb134nf_outer:
	mov   rax, [rsp + nb134nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb134nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb134nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb134nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	movapd xmm6, xmm0
	movapd xmm7, xmm1

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb134nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb134nf_ii3], ebx		

	addsd xmm3, [rax + rbx*8] 	;# ox
	addsd xmm4, [rax + rbx*8 + 8] 	;# oy
	addsd xmm5, [rax + rbx*8 + 16]	;# oz	
	addsd xmm6, [rax + rbx*8 + 24] 	;# h1x
	addsd xmm7, [rax + rbx*8 + 32] 	;# h1y
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	shufpd xmm7, xmm7, 0
	movapd [rsp + nb134nf_ixO], xmm3
	movapd [rsp + nb134nf_iyO], xmm4
	movapd [rsp + nb134nf_izO], xmm5
	movapd [rsp + nb134nf_ixH1], xmm6
	movapd [rsp + nb134nf_iyH1], xmm7

	movsd xmm6, xmm2
	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm6, [rax + rbx*8 + 40] ;# h1z
	addsd xmm0, [rax + rbx*8 + 48] ;# h2x
	addsd xmm1, [rax + rbx*8 + 56] ;# h2y
	addsd xmm2, [rax + rbx*8 + 64] ;# h2z
	addsd xmm3, [rax + rbx*8 + 72] ;# mx
	addsd xmm4, [rax + rbx*8 + 80] ;# my
	addsd xmm5, [rax + rbx*8 + 88] ;# mz

	shufpd xmm6, xmm6, 0
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb134nf_izH1], xmm6
	movapd [rsp + nb134nf_ixH2], xmm0
	movapd [rsp + nb134nf_iyH2], xmm1
	movapd [rsp + nb134nf_izH2], xmm2
	movapd [rsp + nb134nf_ixM], xmm3
	movapd [rsp + nb134nf_iyM], xmm4
	movapd [rsp + nb134nf_izM], xmm5

	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb134nf_vctot], xmm4
	movapd [rsp + nb134nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb134nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb134nf_pos] 
	mov   rdi, [rbp + nb134nf_faction]	
	mov   rax, [rsp + nb134nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb134nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb134nf_ninner]
	mov   [rsp + nb134nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb134nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb134nf_unroll_loop
	jmp   .nb134nf_checksingle
.nb134nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb134nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb134nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb134nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm2, [rsi + rbx*8]
	movhpd xmm0, [rsi + rax*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 8]
	movlpd xmm3, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rbx*8 + 16]
	movhpd xmm3, [rsi + rax*8 + 24]
	movhpd xmm5, [rsi + rbx*8 + 24]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# ox 
	unpckhpd xmm1, xmm2 ;# oy
	unpcklpd xmm3, xmm5 ;# ox 
	unpckhpd xmm4, xmm5 ;# oy
	movapd 	[rsp + nb134nf_jxO], xmm0
	movapd 	[rsp + nb134nf_jyO], xmm1
	movapd 	[rsp + nb134nf_jzO], xmm3
	movapd 	[rsp + nb134nf_jxH1], xmm4
	
	;# load h1y, h1z, h2x, h2y 
	movlpd xmm0, [rsi + rax*8 + 32]
	movlpd xmm2, [rsi + rbx*8 + 32]
	movhpd xmm0, [rsi + rax*8 + 40]
	movhpd xmm2, [rsi + rbx*8 + 40]
	movlpd xmm3, [rsi + rax*8 + 48]
	movlpd xmm5, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rax*8 + 56]
	movhpd xmm5, [rsi + rbx*8 + 56]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h1y
	unpckhpd xmm1, xmm2 ;# h1z
	unpcklpd xmm3, xmm5 ;# h2x
	unpckhpd xmm4, xmm5 ;# h2y
	movapd 	[rsp + nb134nf_jyH1], xmm0
	movapd 	[rsp + nb134nf_jzH1], xmm1
	movapd 	[rsp + nb134nf_jxH2], xmm3
	movapd 	[rsp + nb134nf_jyH2], xmm4
	
	;# load h2z, mx, my, mz
	movlpd xmm0, [rsi + rax*8 + 64]
	movlpd xmm2, [rsi + rbx*8 + 64]
	movhpd xmm0, [rsi + rax*8 + 72]
	movhpd xmm2, [rsi + rbx*8 + 72]
	movlpd xmm3, [rsi + rax*8 + 80]
	movlpd xmm5, [rsi + rbx*8 + 80]
	movhpd xmm3, [rsi + rax*8 + 88]
	movhpd xmm5, [rsi + rbx*8 + 88]
	movapd xmm1, xmm0 
	movapd xmm4, xmm3
	unpcklpd xmm0, xmm2 ;# h2z
	unpckhpd xmm1, xmm2 ;# mx
	unpcklpd xmm3, xmm5 ;# my
	unpckhpd xmm4, xmm5 ;# mz
	movapd 	[rsp + nb134nf_jzH2], xmm0
	movapd 	[rsp + nb134nf_jxM], xmm1
	movapd 	[rsp + nb134nf_jyM], xmm3
	movapd 	[rsp + nb134nf_jzM], xmm4
	
	;# start calculating pairwise distances
	movapd xmm0, [rsp + nb134nf_ixO]
	movapd xmm1, [rsp + nb134nf_iyO]
	movapd xmm2, [rsp + nb134nf_izO]
	movapd xmm3, [rsp + nb134nf_ixH1]
	movapd xmm4, [rsp + nb134nf_iyH1]
	movapd xmm5, [rsp + nb134nf_izH1]
	subpd  xmm0, [rsp + nb134nf_jxO]
	subpd  xmm1, [rsp + nb134nf_jyO]
	subpd  xmm2, [rsp + nb134nf_jzO]
	subpd  xmm3, [rsp + nb134nf_jxH1]
	subpd  xmm4, [rsp + nb134nf_jyH1]
	subpd  xmm5, [rsp + nb134nf_jzH1]
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
	movapd [rsp + nb134nf_rsqOO], xmm0
	movapd [rsp + nb134nf_rsqH1H1], xmm3

	movapd xmm0, [rsp + nb134nf_ixH1]
	movapd xmm1, [rsp + nb134nf_iyH1]
	movapd xmm2, [rsp + nb134nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb134nf_jxH2]
	subpd  xmm1, [rsp + nb134nf_jyH2]
	subpd  xmm2, [rsp + nb134nf_jzH2]
	subpd  xmm3, [rsp + nb134nf_jxM]
	subpd  xmm4, [rsp + nb134nf_jyM]
	subpd  xmm5, [rsp + nb134nf_jzM]
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
	movapd [rsp + nb134nf_rsqH1H2], xmm0
	movapd [rsp + nb134nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb134nf_ixH2]
	movapd xmm1, [rsp + nb134nf_iyH2]
	movapd xmm2, [rsp + nb134nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb134nf_jxH1]
	subpd  xmm1, [rsp + nb134nf_jyH1]
	subpd  xmm2, [rsp + nb134nf_jzH1]
	subpd  xmm3, [rsp + nb134nf_jxH2]
	subpd  xmm4, [rsp + nb134nf_jyH2]
	subpd  xmm5, [rsp + nb134nf_jzH2]
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
	movapd [rsp + nb134nf_rsqH2H1], xmm0
	movapd [rsp + nb134nf_rsqH2H2], xmm3

	movapd xmm0, [rsp + nb134nf_ixH2]
	movapd xmm1, [rsp + nb134nf_iyH2]
	movapd xmm2, [rsp + nb134nf_izH2]
	movapd xmm3, [rsp + nb134nf_ixM]
	movapd xmm4, [rsp + nb134nf_iyM]
	movapd xmm5, [rsp + nb134nf_izM]
	subpd  xmm0, [rsp + nb134nf_jxM]
	subpd  xmm1, [rsp + nb134nf_jyM]
	subpd  xmm2, [rsp + nb134nf_jzM]
	subpd  xmm3, [rsp + nb134nf_jxH1]
	subpd  xmm4, [rsp + nb134nf_jyH1]
	subpd  xmm5, [rsp + nb134nf_jzH1]
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
	movapd [rsp + nb134nf_rsqH2M], xmm0
	movapd [rsp + nb134nf_rsqMH1], xmm4

	movapd xmm0, [rsp + nb134nf_ixM]
	movapd xmm1, [rsp + nb134nf_iyM]
	movapd xmm2, [rsp + nb134nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subpd  xmm0, [rsp + nb134nf_jxH2]
	subpd  xmm1, [rsp + nb134nf_jyH2]
	subpd  xmm2, [rsp + nb134nf_jzH2]
	subpd  xmm3, [rsp + nb134nf_jxM]
	subpd  xmm4, [rsp + nb134nf_jyM]
	subpd  xmm5, [rsp + nb134nf_jzM]
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
	movapd [rsp + nb134nf_rsqMH2], xmm0
	movapd [rsp + nb134nf_rsqMM], xmm4
	
	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtpd2ps xmm1, xmm0
	cvtpd2ps xmm5, xmm4
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1   ;# luA
	cvtps2pd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb134nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvMH2], xmm1
	movapd [rsp + nb134nf_rinvMM], xmm5

	movapd xmm0, [rsp + nb134nf_rsqOO]
	movapd xmm4, [rsp + nb134nf_rsqH1H1]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb134nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb134nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb134nf_half] ;# rinv
	movapd [rsp + nb134nf_rinvOO], xmm1
	movapd [rsp + nb134nf_rinvH1H1], xmm5

	movapd xmm0, [rsp + nb134nf_rsqH1H2]
	movapd xmm4, [rsp + nb134nf_rsqH1M]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb134nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvH1H2], xmm1
	movapd [rsp + nb134nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb134nf_rsqH2H1]
	movapd xmm4, [rsp + nb134nf_rsqH2H2]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb134nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvH2H1], xmm1
	movapd [rsp + nb134nf_rinvH2H2], xmm5

	movapd xmm0, [rsp + nb134nf_rsqMH1]
	movapd xmm4, [rsp + nb134nf_rsqH2M]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb134nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvMH1], xmm1
	movapd [rsp + nb134nf_rinvH2M], xmm5

	;# start with OO interaction 
	movapd xmm0, [rsp + nb134nf_rinvOO] 
	movapd xmm4, [rsp + nb134nf_rsqOO]
	
		mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb134nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  rsi, [rbp + nb134nf_VFtab]
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

	movapd xmm4, [rsp + nb134nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb134nf_Vvdwtot]
	movapd [rsp + nb134nf_Vvdwtot], xmm5

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
	
	movapd xmm4, [rsp + nb134nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb134nf_Vvdwtot]
	movapd [rsp + nb134nf_Vvdwtot], xmm5

	;# All Coulomb interactions
	movapd xmm0, [rsp + nb134nf_rinvH1H1]
	movapd xmm1, [rsp + nb134nf_rinvH1M]
	addpd  xmm0, [rsp + nb134nf_rinvH1H2]
	addpd  xmm1, [rsp + nb134nf_rinvH2M]
	addpd  xmm0, [rsp + nb134nf_rinvH2H1]
	addpd  xmm1, [rsp + nb134nf_rinvMH1]
	addpd  xmm0, [rsp + nb134nf_rinvH2H2]
	addpd  xmm1, [rsp + nb134nf_rinvMH2]
	movapd xmm2, [rsp + nb134nf_rinvMM]
	
	mulpd  xmm0, [rsp + nb134nf_qqHH]
	mulpd  xmm1, [rsp + nb134nf_qqMH]
	mulpd  xmm2, [rsp + nb134nf_qqMM]
	addpd  xmm0, xmm1
	addpd  xmm2, [rsp + nb134nf_vctot]
	addpd  xmm0, xmm2
	movapd [rsp + nb134nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb134nf_innerk],  2
	jl    .nb134nf_checksingle
	jmp   .nb134nf_unroll_loop
.nb134nf_checksingle:
	mov   edx, [rsp + nb134nf_innerk]
	and   edx, 1
	jnz   .nb134nf_dosingle
	jmp   .nb134nf_updateouterdata
.nb134nf_dosingle:
	mov   rdx, [rsp + nb134nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb134nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	
	;# move j coordinates to local temp variables 
	;# load ox, oy, oz, h1x
	movlpd xmm0, [rsi + rax*8]
	movhpd xmm0, [rsi + rax*8 + 8]
	movlpd xmm1, [rsi + rax*8 + 16]
	movhpd xmm1, [rsi + rax*8 + 24]
	movlpd xmm2, [rsi + rax*8 + 32]
	movhpd xmm2, [rsi + rax*8 + 40]
	movlpd xmm3, [rsi + rax*8 + 48]
	movhpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm4, [rsi + rax*8 + 72]
	movlpd xmm5, [rsi + rax*8 + 80]
	movhpd xmm5, [rsi + rax*8 + 88]
	movsd  [rsp + nb134nf_jxO], xmm0
	movsd  [rsp + nb134nf_jzO], xmm1
	movsd  [rsp + nb134nf_jyH1], xmm2
	movsd  [rsp + nb134nf_jxH2], xmm3
	movsd  [rsp + nb134nf_jzH2], xmm4
	movsd  [rsp + nb134nf_jyM], xmm5
	unpckhpd xmm0, xmm0
	unpckhpd xmm1, xmm1
	unpckhpd xmm2, xmm2
	unpckhpd xmm3, xmm3
	unpckhpd xmm4, xmm4
	unpckhpd xmm5, xmm5
	movsd  [rsp + nb134nf_jyO], xmm0
	movsd  [rsp + nb134nf_jxH1], xmm1
	movsd  [rsp + nb134nf_jzH1], xmm2
	movsd  [rsp + nb134nf_jyH2], xmm3
	movsd  [rsp + nb134nf_jxM], xmm4
	movsd  [rsp + nb134nf_jzM], xmm5

	;# start calculating pairwise distances
	movapd xmm0, [rsp + nb134nf_ixO]
	movapd xmm1, [rsp + nb134nf_iyO]
	movapd xmm2, [rsp + nb134nf_izO]
	movapd xmm3, [rsp + nb134nf_ixH1]
	movapd xmm4, [rsp + nb134nf_iyH1]
	movapd xmm5, [rsp + nb134nf_izH1]
	subsd  xmm0, [rsp + nb134nf_jxO]
	subsd  xmm1, [rsp + nb134nf_jyO]
	subsd  xmm2, [rsp + nb134nf_jzO]
	subsd  xmm3, [rsp + nb134nf_jxH1]
	subsd  xmm4, [rsp + nb134nf_jyH1]
	subsd  xmm5, [rsp + nb134nf_jzH1]
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
	movapd [rsp + nb134nf_rsqOO], xmm0
	movapd [rsp + nb134nf_rsqH1H1], xmm3

	movapd xmm0, [rsp + nb134nf_ixH1]
	movapd xmm1, [rsp + nb134nf_iyH1]
	movapd xmm2, [rsp + nb134nf_izH1]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb134nf_jxH2]
	subsd  xmm1, [rsp + nb134nf_jyH2]
	subsd  xmm2, [rsp + nb134nf_jzH2]
	subsd  xmm3, [rsp + nb134nf_jxM]
	subsd  xmm4, [rsp + nb134nf_jyM]
	subsd  xmm5, [rsp + nb134nf_jzM]
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
	movapd [rsp + nb134nf_rsqH1H2], xmm0
	movapd [rsp + nb134nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb134nf_ixH2]
	movapd xmm1, [rsp + nb134nf_iyH2]
	movapd xmm2, [rsp + nb134nf_izH2]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb134nf_jxH1]
	subsd  xmm1, [rsp + nb134nf_jyH1]
	subsd  xmm2, [rsp + nb134nf_jzH1]
	subsd  xmm3, [rsp + nb134nf_jxH2]
	subsd  xmm4, [rsp + nb134nf_jyH2]
	subsd  xmm5, [rsp + nb134nf_jzH2]
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
	movapd [rsp + nb134nf_rsqH2H1], xmm0
	movapd [rsp + nb134nf_rsqH2H2], xmm3

	movapd xmm0, [rsp + nb134nf_ixH2]
	movapd xmm1, [rsp + nb134nf_iyH2]
	movapd xmm2, [rsp + nb134nf_izH2]
	movapd xmm3, [rsp + nb134nf_ixM]
	movapd xmm4, [rsp + nb134nf_iyM]
	movapd xmm5, [rsp + nb134nf_izM]
	subsd  xmm0, [rsp + nb134nf_jxM]
	subsd  xmm1, [rsp + nb134nf_jyM]
	subsd  xmm2, [rsp + nb134nf_jzM]
	subsd  xmm3, [rsp + nb134nf_jxH1]
	subsd  xmm4, [rsp + nb134nf_jyH1]
	subsd  xmm5, [rsp + nb134nf_jzH1]
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
	movapd [rsp + nb134nf_rsqH2M], xmm0
	movapd [rsp + nb134nf_rsqMH1], xmm4

	movapd xmm0, [rsp + nb134nf_ixM]
	movapd xmm1, [rsp + nb134nf_iyM]
	movapd xmm2, [rsp + nb134nf_izM]
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	subsd  xmm0, [rsp + nb134nf_jxH2]
	subsd  xmm1, [rsp + nb134nf_jyH2]
	subsd  xmm2, [rsp + nb134nf_jzH2]
	subsd  xmm3, [rsp + nb134nf_jxM]
	subsd  xmm4, [rsp + nb134nf_jyM]
	subsd  xmm5, [rsp + nb134nf_jzM]
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
	movapd [rsp + nb134nf_rsqMH2], xmm0
	movapd [rsp + nb134nf_rsqMM], xmm4

	;# Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
	cvtsd2ss xmm1, xmm0
	cvtsd2ss xmm5, xmm4
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1   ;# luA
	cvtss2sd xmm5, xmm5   ;# luB
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb134nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvMH2], xmm1
	movapd [rsp + nb134nf_rinvMM], xmm5

	movapd xmm0, [rsp + nb134nf_rsqOO]
	movapd xmm4, [rsp + nb134nf_rsqH1H1]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb134nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb134nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb134nf_half] ;# rinv
	movapd [rsp + nb134nf_rinvOO], xmm1
	movapd [rsp + nb134nf_rinvH1H1], xmm5

	movapd xmm0, [rsp + nb134nf_rsqH1H2]
	movapd xmm4, [rsp + nb134nf_rsqH1M]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb134nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb134nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvH1H2], xmm1
	movapd [rsp + nb134nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb134nf_rsqH2H1]
	movapd xmm4, [rsp + nb134nf_rsqH2H2]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb134nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvH2H1], xmm1
	movapd [rsp + nb134nf_rinvH2H2], xmm5

	movapd xmm0, [rsp + nb134nf_rsqMH1]
	movapd xmm4, [rsp + nb134nf_rsqH2M]	
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
	movapd  xmm3, [rsp + nb134nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb134nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb134nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb134nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb134nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb134nf_half] ;# rinv 
	movapd [rsp + nb134nf_rinvMH1], xmm1
	movapd [rsp + nb134nf_rinvH2M], xmm5

	;# start with OO interaction 
	movsd xmm0, [rsp + nb134nf_rinvOO] 
	movsd xmm4, [rsp + nb134nf_rsqOO]
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb134nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  rsi, [rbp + nb134nf_VFtab]

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

	movsd xmm4, [rsp + nb134nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb134nf_Vvdwtot]
	movsd [rsp + nb134nf_Vvdwtot], xmm5

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
	
	movsd xmm4, [rsp + nb134nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb134nf_Vvdwtot]
	movsd [rsp + nb134nf_Vvdwtot], xmm5

	;# All Coulomb interactions
	movsd xmm0, [rsp + nb134nf_rinvH1H1]
	movsd xmm1, [rsp + nb134nf_rinvH1M]
	addsd  xmm0, [rsp + nb134nf_rinvH1H2]
	addsd  xmm1, [rsp + nb134nf_rinvH2M]
	addsd  xmm0, [rsp + nb134nf_rinvH2H1]
	addsd  xmm1, [rsp + nb134nf_rinvMH1]
	addsd  xmm0, [rsp + nb134nf_rinvH2H2]
	addsd  xmm1, [rsp + nb134nf_rinvMH2]
	movsd xmm2, [rsp + nb134nf_rinvMM]
	
	mulsd  xmm0, [rsp + nb134nf_qqHH]
	mulsd  xmm1, [rsp + nb134nf_qqMH]
	mulsd  xmm2, [rsp + nb134nf_qqMM]
	addsd  xmm0, xmm1
	addsd  xmm2, [rsp + nb134nf_vctot]
	addsd  xmm0, xmm2
	movsd [rsp + nb134nf_vctot], xmm0
	
.nb134nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb134nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb134nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb134nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb134nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb134nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb134nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb134nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb134nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb134nf_n], esi
        jmp .nb134nf_outer
.nb134nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb134nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb134nf_end
        ;# non-zero, do one more workunit
        jmp   .nb134nf_threadloop
.nb134nf_end:
	mov eax, [rsp + nb134nf_nouter]
	mov ebx, [rsp + nb134nf_ninner]
	mov rcx, [rbp + nb134nf_outeriter]
	mov rdx, [rbp + nb134nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 992
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
