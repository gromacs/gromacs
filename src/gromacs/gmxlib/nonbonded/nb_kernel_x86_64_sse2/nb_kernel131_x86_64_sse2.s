##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##





.globl nb_kernel131_x86_64_sse2
.globl _nb_kernel131_x86_64_sse2
nb_kernel131_x86_64_sse2:       
_nb_kernel131_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb131_fshift, 16
.set nb131_gid, 24
.set nb131_pos, 32
.set nb131_faction, 40
.set nb131_charge, 48
.set nb131_p_facel, 56
.set nb131_argkrf, 64
.set nb131_argcrf, 72
.set nb131_Vc, 80
.set nb131_type, 88
.set nb131_p_ntype, 96
.set nb131_vdwparam, 104
.set nb131_Vvdw, 112
.set nb131_p_tabscale, 120
.set nb131_VFtab, 128
.set nb131_invsqrta, 136
.set nb131_dvda, 144
.set nb131_p_gbtabscale, 152
.set nb131_GBtab, 160
.set nb131_p_nthreads, 168
.set nb131_count, 176
.set nb131_mtx, 184
.set nb131_outeriter, 192
.set nb131_inneriter, 200
.set nb131_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb131_ixO, 0
.set nb131_iyO, 16
.set nb131_izO, 32
.set nb131_ixH1, 48
.set nb131_iyH1, 64
.set nb131_izH1, 80
.set nb131_ixH2, 96
.set nb131_iyH2, 112
.set nb131_izH2, 128
.set nb131_iqO, 144
.set nb131_iqH, 160
.set nb131_dxO, 176
.set nb131_dyO, 192
.set nb131_dzO, 208
.set nb131_dxH1, 224
.set nb131_dyH1, 240
.set nb131_dzH1, 256
.set nb131_dxH2, 272
.set nb131_dyH2, 288
.set nb131_dzH2, 304
.set nb131_qqO, 320
.set nb131_qqH, 336
.set nb131_c6, 352
.set nb131_c12, 368
.set nb131_tsc, 384
.set nb131_fstmp, 400
.set nb131_vctot, 416
.set nb131_Vvdwtot, 432
.set nb131_fixO, 448
.set nb131_fiyO, 464
.set nb131_fizO, 480
.set nb131_fixH1, 496
.set nb131_fiyH1, 512
.set nb131_fizH1, 528
.set nb131_fixH2, 544
.set nb131_fiyH2, 560
.set nb131_fizH2, 576
.set nb131_rsqO, 592
.set nb131_rsqH1, 608
.set nb131_rsqH2, 624
.set nb131_half, 640
.set nb131_three, 656
.set nb131_two, 672
.set nb131_krf, 688
.set nb131_crf, 704
.set nb131_krsqO, 720
.set nb131_krsqH1, 736
.set nb131_krsqH2, 752
.set nb131_rinvO, 768
.set nb131_rinvH1, 784
.set nb131_rinvH2, 800
.set nb131_facel, 816
.set nb131_iinr, 824
.set nb131_jindex, 832
.set nb131_jjnr, 840
.set nb131_shift, 848
.set nb131_shiftvec, 856
.set nb131_innerjjnr, 864
.set nb131_nri, 872
.set nb131_is3, 876
.set nb131_ii3, 880
.set nb131_ntia, 884
.set nb131_innerk, 888
.set nb131_n, 892
.set nb131_nn1, 896
.set nb131_nouter, 900
.set nb131_ninner, 904

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $920,%rsp  # # local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb131_nouter(%rsp)
        movl %eax,nb131_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb131_nri(%rsp)
        movq %rsi,nb131_iinr(%rsp)
        movq %rdx,nb131_jindex(%rsp)
        movq %rcx,nb131_jjnr(%rsp)
        movq %r8,nb131_shift(%rsp)
        movq %r9,nb131_shiftvec(%rsp)
        movq nb131_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb131_facel(%rsp)

        movq nb131_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb131_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb131_half(%rsp)
        movl %ebx,nb131_half+4(%rsp)
        movsd nb131_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb131_half(%rsp)
        movapd %xmm2,nb131_two(%rsp)
        movapd %xmm3,nb131_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb131_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb131_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4

        movsd nb131_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb131_iqO(%rsp)
        movapd %xmm4,nb131_iqH(%rsp)

        movq  nb131_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb131_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131_ntia(%rsp)
_nb_kernel131_x86_64_sse2.nb131_threadloop: 
        movq  nb131_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel131_x86_64_sse2.nb131_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131_x86_64_sse2.nb131_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131_n(%rsp)
        movl %ebx,nb131_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131_x86_64_sse2.nb131_outerstart
        jmp _nb_kernel131_x86_64_sse2.nb131_end

_nb_kernel131_x86_64_sse2.nb131_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131_nouter(%rsp),%ebx
        movl %ebx,nb131_nouter(%rsp)

_nb_kernel131_x86_64_sse2.nb131_outer: 
        movq  nb131_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb131_is3(%rsp)      ## store is3 

        movq  nb131_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb131_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb131_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb131_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb131_ixO(%rsp)
        movapd %xmm4,nb131_iyO(%rsp)
        movapd %xmm5,nb131_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb131_ixH1(%rsp)
        movapd %xmm1,nb131_iyH1(%rsp)
        movapd %xmm2,nb131_izH1(%rsp)
        movapd %xmm3,nb131_ixH2(%rsp)
        movapd %xmm4,nb131_iyH2(%rsp)
        movapd %xmm5,nb131_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb131_vctot(%rsp)
        movapd %xmm4,nb131_Vvdwtot(%rsp)
        movapd %xmm4,nb131_fixO(%rsp)
        movapd %xmm4,nb131_fiyO(%rsp)
        movapd %xmm4,nb131_fizO(%rsp)
        movapd %xmm4,nb131_fixH1(%rsp)
        movapd %xmm4,nb131_fiyH1(%rsp)
        movapd %xmm4,nb131_fizH1(%rsp)
        movapd %xmm4,nb131_fixH2(%rsp)
        movapd %xmm4,nb131_fiyH2(%rsp)
        movapd %xmm4,nb131_fizH2(%rsp)

        movq  nb131_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb131_pos(%rbp),%rsi
        movq  nb131_faction(%rbp),%rdi
        movq  nb131_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb131_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb131_ninner(%rsp),%ecx
        movl  %ecx,nb131_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb131_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel131_x86_64_sse2.nb131_unroll_loop
        jmp   _nb_kernel131_x86_64_sse2.nb131_checksingle
_nb_kernel131_x86_64_sse2.nb131_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb131_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb131_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb131_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131_iqO(%rsp),%xmm3
        mulpd  nb131_iqH(%rsp),%xmm4

        movapd  %xmm3,nb131_qqO(%rsp)
        movapd  %xmm4,nb131_qqH(%rsp)

        movq nb131_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb131_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb131_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb131_c6(%rsp)
        movapd %xmm6,nb131_c12(%rsp)

        movq nb131_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb131_ixO(%rsp),%xmm0
    subpd nb131_iyO(%rsp),%xmm1
    subpd nb131_izO(%rsp),%xmm2
    subpd nb131_ixH1(%rsp),%xmm3
    subpd nb131_iyH1(%rsp),%xmm4
    subpd nb131_izH1(%rsp),%xmm5
    subpd nb131_ixH2(%rsp),%xmm6
    subpd nb131_iyH2(%rsp),%xmm7
    subpd nb131_izH2(%rsp),%xmm8

        movapd %xmm0,nb131_dxO(%rsp)
        movapd %xmm1,nb131_dyO(%rsp)
        movapd %xmm2,nb131_dzO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb131_dxH1(%rsp)
        movapd %xmm4,nb131_dyH1(%rsp)
        movapd %xmm5,nb131_dzH1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb131_dxH2(%rsp)
        movapd %xmm7,nb131_dyH2(%rsp)
        movapd %xmm8,nb131_dzH2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb131_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb131_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvO
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb131_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb131_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvO 
        mulpd   %xmm15,%xmm10 ##   rinvH1
    mulpd   %xmm15,%xmm11 ##   rinvH2

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm0,nb131_rsqO(%rsp)
    movapd %xmm3,nb131_rsqH1(%rsp)
    movapd %xmm6,nb131_rsqH2(%rsp)
    movapd %xmm9,nb131_rinvO(%rsp)
    movapd %xmm10,nb131_rinvH1(%rsp)
    movapd %xmm11,nb131_rinvH2(%rsp)

    ## table LJ interaction
    mulpd  %xmm9,%xmm0
    mulpd  nb131_tsc(%rsp),%xmm0   ## rtab

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2

    ## multiply by 8
    pslld   $3,%xmm1

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    movd    %xmm1,%r8d
    movd    %xmm13,%r10d

    ## calculate eps
    subpd     %xmm2,%xmm0
    movq nb131_VFtab(%rbp),%rsi

    movlpd (%rsi,%r8,8),%xmm4
        movlpd 8(%rsi,%r8,8),%xmm5
        movlpd 16(%rsi,%r8,8),%xmm6
        movlpd 24(%rsi,%r8,8),%xmm7
    movlpd 32(%rsi,%r8,8),%xmm8
        movlpd 40(%rsi,%r8,8),%xmm9
        movlpd 48(%rsi,%r8,8),%xmm10
        movlpd 56(%rsi,%r8,8),%xmm11

    movhpd (%rsi,%r10,8),%xmm4
        movhpd 8(%rsi,%r10,8),%xmm5
        movhpd 16(%rsi,%r10,8),%xmm6
        movhpd 24(%rsi,%r10,8),%xmm7
    movhpd 32(%rsi,%r10,8),%xmm8
        movhpd 40(%rsi,%r10,8),%xmm9
        movhpd 48(%rsi,%r10,8),%xmm10
        movhpd 56(%rsi,%r10,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  %xmm0,%xmm7   ## Heps
    mulpd  %xmm0,%xmm11
    mulpd  %xmm0,%xmm6  ## Geps
    mulpd  %xmm0,%xmm10
    mulpd  %xmm0,%xmm7  ## Heps2
    mulpd  %xmm0,%xmm11
    addpd  %xmm6,%xmm5 ## F+Geps
    addpd  %xmm10,%xmm9
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm11,%xmm9
    addpd  %xmm7,%xmm7   ## 2*Heps2
    addpd  %xmm11,%xmm11
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm10,%xmm11

    addpd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm9,%xmm11
    mulpd  %xmm0,%xmm5 ## eps*Fp
    mulpd  %xmm0,%xmm9
    movapd nb131_c6(%rsp),%xmm12
    movapd nb131_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb131_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb131_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7
    mulpd  nb131_tsc(%rsp),%xmm7

    movapd nb131_rinvO(%rsp),%xmm9
    movapd nb131_rinvH1(%rsp),%xmm10
    movapd nb131_rinvH2(%rsp),%xmm11

    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2

    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11

    mulpd  nb131_qqO(%rsp),%xmm0
    mulpd  nb131_qqH(%rsp),%xmm1
    mulpd  nb131_qqH(%rsp),%xmm2

    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    subpd  %xmm7,%xmm9
    mulpd  nb131_rinvO(%rsp),%xmm9

    addpd nb131_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb131_vctot(%rsp)

    ## move j forces to xmm0-xmm2
    movq nb131_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb131_dxO(%rsp),%xmm7
        mulpd nb131_dyO(%rsp),%xmm8
        mulpd nb131_dzO(%rsp),%xmm9
        mulpd nb131_dxH1(%rsp),%xmm10
        mulpd nb131_dyH1(%rsp),%xmm11
        mulpd nb131_dzH1(%rsp),%xmm12
        mulpd nb131_dxH2(%rsp),%xmm13
        mulpd nb131_dyH2(%rsp),%xmm14
        mulpd nb131_dzH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb131_fixO(%rsp),%xmm7
    addpd nb131_fiyO(%rsp),%xmm8
    addpd nb131_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb131_fixH1(%rsp),%xmm10
    addpd nb131_fiyH1(%rsp),%xmm11
    addpd nb131_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb131_fixH2(%rsp),%xmm13
    addpd nb131_fiyH2(%rsp),%xmm14
    addpd nb131_fizH2(%rsp),%xmm15

    movapd %xmm7,nb131_fixO(%rsp)
    movapd %xmm8,nb131_fiyO(%rsp)
    movapd %xmm9,nb131_fizO(%rsp)
    movapd %xmm10,nb131_fixH1(%rsp)
    movapd %xmm11,nb131_fiyH1(%rsp)
    movapd %xmm12,nb131_fizH1(%rsp)
    movapd %xmm13,nb131_fixH2(%rsp)
    movapd %xmm14,nb131_fiyH2(%rsp)
    movapd %xmm15,nb131_fizH2(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb131_innerk(%rsp)
        jl    _nb_kernel131_x86_64_sse2.nb131_checksingle
        jmp   _nb_kernel131_x86_64_sse2.nb131_unroll_loop

_nb_kernel131_x86_64_sse2.nb131_checksingle:    
        movl  nb131_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel131_x86_64_sse2.nb131_dosingle
        jmp    _nb_kernel131_x86_64_sse2.nb131_updateouterdata
_nb_kernel131_x86_64_sse2.nb131_dosingle: 
        movq  nb131_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        addq $8,nb131_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb131_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb131_iqO(%rsp),%xmm3
        mulsd  nb131_iqH(%rsp),%xmm4

        movapd  %xmm3,nb131_qqO(%rsp)
        movapd  %xmm4,nb131_qqH(%rsp)

        movq nb131_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb131_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb131_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm6        ## c6a
        movsd 8(%rsi,%r8,8),%xmm7       ## c12a
        movapd %xmm6,nb131_c6(%rsp)
        movapd %xmm7,nb131_c12(%rsp)

        movq nb131_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
    movsd (%rsi,%rax,8),%xmm0
    movsd 8(%rsi,%rax,8),%xmm1
    movsd 16(%rsi,%rax,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subsd nb131_ixO(%rsp),%xmm0
    subsd nb131_iyO(%rsp),%xmm1
    subsd nb131_izO(%rsp),%xmm2
    subsd nb131_ixH1(%rsp),%xmm3
    subsd nb131_iyH1(%rsp),%xmm4
    subsd nb131_izH1(%rsp),%xmm5
    subsd nb131_ixH2(%rsp),%xmm6
    subsd nb131_iyH2(%rsp),%xmm7
    subsd nb131_izH2(%rsp),%xmm8

        movapd %xmm0,nb131_dxO(%rsp)
        movapd %xmm1,nb131_dyO(%rsp)
        movapd %xmm2,nb131_dzO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb131_dxH1(%rsp)
        movapd %xmm4,nb131_dyH1(%rsp)
        movapd %xmm5,nb131_dzH1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movapd %xmm6,nb131_dxH2(%rsp)
        movapd %xmm7,nb131_dyH2(%rsp)
        movapd %xmm8,nb131_dzH2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movapd  nb131_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb131_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvO
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movapd  nb131_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb131_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvO 
        mulsd   %xmm15,%xmm10 ##   rinvH1
    mulsd   %xmm15,%xmm11 ##   rinvH2

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm0,nb131_rsqO(%rsp)
    movapd %xmm3,nb131_rsqH1(%rsp)
    movapd %xmm6,nb131_rsqH2(%rsp)
    movapd %xmm9,nb131_rinvO(%rsp)
    movapd %xmm10,nb131_rinvH1(%rsp)
    movapd %xmm11,nb131_rinvH2(%rsp)

    ## table LJ interaction
    mulsd  %xmm9,%xmm0
    mulsd  nb131_tsc(%rsp),%xmm0   ## rtab

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2

    ## mult. by 8
    shll $3,%r8d

    ## calculate eps
    subsd     %xmm2,%xmm0
    movq nb131_VFtab(%rbp),%rsi

    movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm5
        movsd 16(%rsi,%r8,8),%xmm6
        movsd 24(%rsi,%r8,8),%xmm7
    movsd 32(%rsi,%r8,8),%xmm8
        movsd 40(%rsi,%r8,8),%xmm9
        movsd 48(%rsi,%r8,8),%xmm10
        movsd 56(%rsi,%r8,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulsd  %xmm0,%xmm7   ## Heps
    mulsd  %xmm0,%xmm11
    mulsd  %xmm0,%xmm6  ## Geps
    mulsd  %xmm0,%xmm10
    mulsd  %xmm0,%xmm7  ## Heps2
    mulsd  %xmm0,%xmm11
    addsd  %xmm6,%xmm5 ## F+Geps
    addsd  %xmm10,%xmm9
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm11,%xmm9
    addsd  %xmm7,%xmm7   ## 2*Heps2
    addsd  %xmm11,%xmm11
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm10,%xmm11

    addsd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm9,%xmm11
    mulsd  %xmm0,%xmm5 ## eps*Fp
    mulsd  %xmm0,%xmm9
    movapd nb131_c6(%rsp),%xmm12
    movapd nb131_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb131_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb131_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7
    mulsd  nb131_tsc(%rsp),%xmm7

    movapd nb131_rinvO(%rsp),%xmm9
    movapd nb131_rinvH1(%rsp),%xmm10
    movapd nb131_rinvH2(%rsp),%xmm11

    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2

    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11

    mulsd  nb131_qqO(%rsp),%xmm0
    mulsd  nb131_qqH(%rsp),%xmm1
    mulsd  nb131_qqH(%rsp),%xmm2

    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    subsd  %xmm7,%xmm9
    mulsd  nb131_rinvO(%rsp),%xmm9

    addsd nb131_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb131_vctot(%rsp)

    ## move j forces to xmm0-xmm2
    movq nb131_faction(%rbp),%rdi
        movsd (%rdi,%rax,8),%xmm0
        movsd 8(%rdi,%rax,8),%xmm1
        movsd 16(%rdi,%rax,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulsd nb131_dxO(%rsp),%xmm7
        mulsd nb131_dyO(%rsp),%xmm8
        mulsd nb131_dzO(%rsp),%xmm9
        mulsd nb131_dxH1(%rsp),%xmm10
        mulsd nb131_dyH1(%rsp),%xmm11
        mulsd nb131_dzH1(%rsp),%xmm12
        mulsd nb131_dxH2(%rsp),%xmm13
        mulsd nb131_dyH2(%rsp),%xmm14
        mulsd nb131_dzH2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb131_fixO(%rsp),%xmm7
    addsd nb131_fiyO(%rsp),%xmm8
    addsd nb131_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb131_fixH1(%rsp),%xmm10
    addsd nb131_fiyH1(%rsp),%xmm11
    addsd nb131_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb131_fixH2(%rsp),%xmm13
    addsd nb131_fiyH2(%rsp),%xmm14
    addsd nb131_fizH2(%rsp),%xmm15

    movsd %xmm7,nb131_fixO(%rsp)
    movsd %xmm8,nb131_fiyO(%rsp)
    movsd %xmm9,nb131_fizO(%rsp)
    movsd %xmm10,nb131_fixH1(%rsp)
    movsd %xmm11,nb131_fiyH1(%rsp)
    movsd %xmm12,nb131_fizH1(%rsp)
    movsd %xmm13,nb131_fixH2(%rsp)
    movsd %xmm14,nb131_fiyH2(%rsp)
    movsd %xmm15,nb131_fizH2(%rsp)

    ## store back j forces from xmm0-xmm2
        movsd %xmm0,(%rdi,%rax,8)
        movsd %xmm1,8(%rdi,%rax,8)
        movsd %xmm2,16(%rdi,%rax,8)

_nb_kernel131_x86_64_sse2.nb131_updateouterdata: 
        movl  nb131_ii3(%rsp),%ecx
        movq  nb131_faction(%rbp),%rdi
        movq  nb131_fshift(%rbp),%rsi
        movl  nb131_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb131_fixO(%rsp),%xmm0
        movapd nb131_fiyO(%rsp),%xmm1
        movapd nb131_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb131_fixH1(%rsp),%xmm0
        movapd nb131_fiyH1(%rsp),%xmm1
        movapd nb131_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb131_fixH2(%rsp),%xmm0
        movapd nb131_fiyH2(%rsp),%xmm1
        movapd nb131_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb131_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb131_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb131_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb131_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb131_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb131_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb131_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131_x86_64_sse2.nb131_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131_n(%rsp)
        jmp _nb_kernel131_x86_64_sse2.nb131_outer
_nb_kernel131_x86_64_sse2.nb131_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131_x86_64_sse2.nb131_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131_x86_64_sse2.nb131_threadloop
_nb_kernel131_x86_64_sse2.nb131_end: 
        movl nb131_nouter(%rsp),%eax
        movl nb131_ninner(%rsp),%ebx
        movq nb131_outeriter(%rbp),%rcx
        movq nb131_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $920,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel131nf_x86_64_sse2
.globl _nb_kernel131nf_x86_64_sse2
nb_kernel131nf_x86_64_sse2:     
_nb_kernel131nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb131nf_fshift, 16
.set nb131nf_gid, 24
.set nb131nf_pos, 32
.set nb131nf_faction, 40
.set nb131nf_charge, 48
.set nb131nf_p_facel, 56
.set nb131nf_argkrf, 64
.set nb131nf_argcrf, 72
.set nb131nf_Vc, 80
.set nb131nf_type, 88
.set nb131nf_p_ntype, 96
.set nb131nf_vdwparam, 104
.set nb131nf_Vvdw, 112
.set nb131nf_p_tabscale, 120
.set nb131nf_VFtab, 128
.set nb131nf_invsqrta, 136
.set nb131nf_dvda, 144
.set nb131nf_p_gbtabscale, 152
.set nb131nf_GBtab, 160
.set nb131nf_p_nthreads, 168
.set nb131nf_count, 176
.set nb131nf_mtx, 184
.set nb131nf_outeriter, 192
.set nb131nf_inneriter, 200
.set nb131nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb131nf_ixO, 0
.set nb131nf_iyO, 16
.set nb131nf_izO, 32
.set nb131nf_ixH1, 48
.set nb131nf_iyH1, 64
.set nb131nf_izH1, 80
.set nb131nf_ixH2, 96
.set nb131nf_iyH2, 112
.set nb131nf_izH2, 128
.set nb131nf_iqO, 144
.set nb131nf_iqH, 160
.set nb131nf_qqO, 176
.set nb131nf_qqH, 192
.set nb131nf_c6, 208
.set nb131nf_c12, 224
.set nb131nf_vctot, 240
.set nb131nf_Vvdwtot, 256
.set nb131nf_half, 272
.set nb131nf_three, 288
.set nb131nf_rsqO, 304
.set nb131nf_rinvH1, 320
.set nb131nf_rinvH2, 336
.set nb131nf_krsqO, 352
.set nb131nf_krsqH1, 368
.set nb131nf_krsqH2, 384
.set nb131nf_tsc, 400
.set nb131nf_facel, 416
.set nb131nf_iinr, 424
.set nb131nf_jindex, 432
.set nb131nf_jjnr, 440
.set nb131nf_shift, 448
.set nb131nf_shiftvec, 456
.set nb131nf_innerjjnr, 464
.set nb131nf_nri, 472
.set nb131nf_ntia, 476
.set nb131nf_is3, 480
.set nb131nf_ii3, 484
.set nb131nf_innerk, 488
.set nb131nf_n, 492
.set nb131nf_nn1, 496
.set nb131nf_nouter, 500
.set nb131nf_ninner, 504

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $520,%rsp # # local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb131nf_nouter(%rsp)
        movl %eax,nb131nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb131nf_nri(%rsp)
        movq %rsi,nb131nf_iinr(%rsp)
        movq %rdx,nb131nf_jindex(%rsp)
        movq %rcx,nb131nf_jjnr(%rsp)
        movq %r8,nb131nf_shift(%rsp)
        movq %r9,nb131nf_shiftvec(%rsp)
        movq nb131nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb131nf_facel(%rsp)

        movq nb131nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb131nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb131nf_half(%rsp)
        movl %ebx,nb131nf_half+4(%rsp)
        movsd nb131nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb131nf_half(%rsp)
        movapd %xmm3,nb131nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb131nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb131nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movsd nb131nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb131nf_iqO(%rsp)
        movapd %xmm4,nb131nf_iqH(%rsp)

        movq  nb131nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb131nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131nf_ntia(%rsp)
_nb_kernel131nf_x86_64_sse2.nb131nf_threadloop: 
        movq  nb131nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel131nf_x86_64_sse2.nb131nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131nf_x86_64_sse2.nb131nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131nf_n(%rsp)
        movl %ebx,nb131nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131nf_x86_64_sse2.nb131nf_outerstart
        jmp _nb_kernel131nf_x86_64_sse2.nb131nf_end

_nb_kernel131nf_x86_64_sse2.nb131nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131nf_nouter(%rsp),%ebx
        movl %ebx,nb131nf_nouter(%rsp)

_nb_kernel131nf_x86_64_sse2.nb131nf_outer: 
        movq  nb131nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb131nf_is3(%rsp)            ## store is3 

        movq  nb131nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb131nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb131nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb131nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb131nf_ixO(%rsp)
        movapd %xmm4,nb131nf_iyO(%rsp)
        movapd %xmm5,nb131nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb131nf_ixH1(%rsp)
        movapd %xmm1,nb131nf_iyH1(%rsp)
        movapd %xmm2,nb131nf_izH1(%rsp)
        movapd %xmm3,nb131nf_ixH2(%rsp)
        movapd %xmm4,nb131nf_iyH2(%rsp)
        movapd %xmm5,nb131nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb131nf_vctot(%rsp)
        movapd %xmm4,nb131nf_Vvdwtot(%rsp)

        movq  nb131nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb131nf_pos(%rbp),%rsi
        movq  nb131nf_faction(%rbp),%rdi
        movq  nb131nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb131nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb131nf_ninner(%rsp),%ecx
        movl  %ecx,nb131nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb131nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel131nf_x86_64_sse2.nb131nf_unroll_loop
        jmp   _nb_kernel131nf_x86_64_sse2.nb131nf_checksingle
_nb_kernel131nf_x86_64_sse2.nb131nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb131nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb131nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb131nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131nf_iqO(%rsp),%xmm3
        mulpd  nb131nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb131nf_qqO(%rsp)
        movapd  %xmm4,nb131nf_qqH(%rsp)

        movq nb131nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb131nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb131nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb131nf_c6(%rsp)
        movapd %xmm6,nb131nf_c12(%rsp)

        movq nb131nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb131nf_ixO(%rsp),%xmm4
        movapd nb131nf_iyO(%rsp),%xmm5
        movapd nb131nf_izO(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb131nf_rsqO(%rsp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131nf_ixH1(%rsp),%xmm4
        movapd nb131nf_iyH1(%rsp),%xmm5
        movapd nb131nf_izH1(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb131nf_ixH2(%rsp),%xmm3
        movapd nb131nf_iyH2(%rsp),%xmm4
        movapd nb131nf_izH2(%rsp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movapd  %xmm6,nb131nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movapd  %xmm5,nb131nf_rinvH2(%rsp)

        ## do O interactions 
        movapd %xmm7,%xmm0
        mulpd  nb131nf_qqO(%rsp),%xmm7

        addpd  nb131nf_vctot(%rsp),%xmm7
        movapd %xmm7,nb131nf_vctot(%rsp)

        movapd nb131nf_rsqO(%rsp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb131nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movq nb131nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%rsi,%rax,8),%xmm4      ## Y1   
        movlpd (%rsi,%rbx,8),%xmm3      ## Y2 
        movhpd 8(%rsi,%rax,8),%xmm4     ## Y1 F1        
        movhpd 8(%rsi,%rbx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rax,8),%xmm6    ## G1
        movlpd 16(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 24(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 24(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb131nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ##  Update Vvdwtot directly 
        addpd  nb131nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb131nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rax,8),%xmm4    ## Y1   
        movlpd 32(%rsi,%rbx,8),%xmm3    ## Y2 
        movhpd 40(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movhpd 40(%rsi,%rbx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rax,8),%xmm6    ## G1
        movlpd 48(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 56(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 56(%rsi,%rbx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb131nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb131nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb131nf_Vvdwtot(%rsp)

        ## H1 & H2 interactions 
        movapd  nb131nf_rinvH1(%rsp),%xmm6
        addpd   nb131nf_rinvH2(%rsp),%xmm6
        mulpd   nb131nf_qqH(%rsp),%xmm6   ## vcoul 

        addpd  nb131nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb131nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb131nf_innerk(%rsp)
        jl    _nb_kernel131nf_x86_64_sse2.nb131nf_checksingle
        jmp   _nb_kernel131nf_x86_64_sse2.nb131nf_unroll_loop
_nb_kernel131nf_x86_64_sse2.nb131nf_checksingle: 
        movl  nb131nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel131nf_x86_64_sse2.nb131nf_dosingle
        jmp   _nb_kernel131nf_x86_64_sse2.nb131nf_updateouterdata
_nb_kernel131nf_x86_64_sse2.nb131nf_dosingle: 
        movq  nb131nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb131nf_innerjjnr(%rsp)

        movq nb131nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb131nf_iqO(%rsp),%xmm3
        mulpd  nb131nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb131nf_qqO(%rsp)
        movapd  %xmm4,nb131nf_qqH(%rsp)

        movq nb131nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb131nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb131nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb131nf_c6(%rsp)
        movapd %xmm6,nb131nf_c12(%rsp)

        movq nb131nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb131nf_ixO(%rsp),%xmm4
        movapd nb131nf_iyO(%rsp),%xmm5
        movapd nb131nf_izO(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6


        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb131nf_rsqO(%rsp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb131nf_ixH1(%rsp),%xmm4
        movapd nb131nf_iyH1(%rsp),%xmm5
        movapd nb131nf_izH1(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb131nf_ixH2(%rsp),%xmm3
        movapd nb131nf_iyH2(%rsp),%xmm4
        movapd nb131nf_izH2(%rsp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 
        movsd  %xmm6,nb131nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb131nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb131nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb131nf_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb131nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 
        movsd %xmm5,nb131nf_rinvH2(%rsp)

        ## do O interactions 
        movsd %xmm7,%xmm0
        mulsd  nb131nf_qqO(%rsp),%xmm7

        addsd  nb131nf_vctot(%rsp),%xmm7
        movsd %xmm7,nb131nf_vctot(%rsp)

        movsd nb131nf_rsqO(%rsp),%xmm4
        ## LJ table interaction. xmm0=rinv, xmm4=rsq

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb131nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movq nb131nf_VFtab(%rbp),%rsi

        ## dispersion 
        movlpd (%rsi,%rbx,8),%xmm4      ## Y1   
        movhpd 8(%rsi,%rbx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 24(%rsi,%rbx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb131nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb131nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb131nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rbx,8),%xmm4    ## Y1   
        movhpd 40(%rsi,%rbx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 56(%rsi,%rbx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb131nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb131nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb131nf_Vvdwtot(%rsp)


        ## H1 & H2 interactions 
        movsd  nb131nf_rinvH1(%rsp),%xmm6
        addsd   nb131nf_rinvH2(%rsp),%xmm6
        mulsd   nb131nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd   nb131nf_vctot(%rsp),%xmm6
        movsd  %xmm6,nb131nf_vctot(%rsp)

_nb_kernel131nf_x86_64_sse2.nb131nf_updateouterdata: 
        ## get n from stack
        movl nb131nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb131nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb131nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb131nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb131nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb131nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb131nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131nf_x86_64_sse2.nb131nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131nf_n(%rsp)
        jmp _nb_kernel131nf_x86_64_sse2.nb131nf_outer
_nb_kernel131nf_x86_64_sse2.nb131nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131nf_x86_64_sse2.nb131nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131nf_x86_64_sse2.nb131nf_threadloop
_nb_kernel131nf_x86_64_sse2.nb131nf_end: 
        movl nb131nf_nouter(%rsp),%eax
        movl nb131nf_ninner(%rsp),%ebx
        movq nb131nf_outeriter(%rbp),%rcx
        movq nb131nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $520,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

