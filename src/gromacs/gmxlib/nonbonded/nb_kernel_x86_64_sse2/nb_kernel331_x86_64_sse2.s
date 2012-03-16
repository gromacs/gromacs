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








.globl nb_kernel331_x86_64_sse2
.globl _nb_kernel331_x86_64_sse2
nb_kernel331_x86_64_sse2:       
_nb_kernel331_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb331_fshift, 16
.set nb331_gid, 24
.set nb331_pos, 32
.set nb331_faction, 40
.set nb331_charge, 48
.set nb331_p_facel, 56
.set nb331_argkrf, 64
.set nb331_argcrf, 72
.set nb331_Vc, 80
.set nb331_type, 88
.set nb331_p_ntype, 96
.set nb331_vdwparam, 104
.set nb331_Vvdw, 112
.set nb331_p_tabscale, 120
.set nb331_VFtab, 128
.set nb331_invsqrta, 136
.set nb331_dvda, 144
.set nb331_p_gbtabscale, 152
.set nb331_GBtab, 160
.set nb331_p_nthreads, 168
.set nb331_count, 176
.set nb331_mtx, 184
.set nb331_outeriter, 192
.set nb331_inneriter, 200
.set nb331_work, 208
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb331_ixO, 0
.set nb331_iyO, 16
.set nb331_izO, 32
.set nb331_ixH1, 48
.set nb331_iyH1, 64
.set nb331_izH1, 80
.set nb331_ixH2, 96
.set nb331_iyH2, 112
.set nb331_izH2, 128
.set nb331_iqO, 144
.set nb331_iqH, 160
.set nb331_dxO, 176
.set nb331_dyO, 192
.set nb331_dzO, 208
.set nb331_dxH1, 224
.set nb331_dyH1, 240
.set nb331_dzH1, 256
.set nb331_dxH2, 272
.set nb331_dyH2, 288
.set nb331_dzH2, 304
.set nb331_qqO, 320
.set nb331_qqH, 336
.set nb331_rinvO, 352
.set nb331_rinvH1, 368
.set nb331_rinvH2, 384
.set nb331_rO, 400
.set nb331_rH1, 416
.set nb331_rH2, 432
.set nb331_tsc, 448
.set nb331_two, 464
.set nb331_c6, 480
.set nb331_c12, 496
.set nb331_vctot, 512
.set nb331_Vvdwtot, 528
.set nb331_fixO, 544
.set nb331_fiyO, 560
.set nb331_fizO, 576
.set nb331_fixH1, 592
.set nb331_fiyH1, 608
.set nb331_fizH1, 624
.set nb331_fixH2, 640
.set nb331_fiyH2, 656
.set nb331_fizH2, 672
.set nb331_fjx, 688
.set nb331_fjy, 704
.set nb331_fjz, 720
.set nb331_half, 736
.set nb331_three, 752
.set nb331_is3, 768
.set nb331_ii3, 772
.set nb331_nri, 776
.set nb331_iinr, 784
.set nb331_jindex, 792
.set nb331_jjnr, 800
.set nb331_shift, 808
.set nb331_shiftvec, 816
.set nb331_facel, 824
.set nb331_innerjjnr, 832
.set nb331_ntia, 840
.set nb331_innerk, 844
.set nb331_n, 848
.set nb331_nn1, 852
.set nb331_nouter, 856
.set nb331_ninner, 860
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $872,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb331_nouter(%rsp)
        movl %eax,nb331_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb331_nri(%rsp)
        movq %rsi,nb331_iinr(%rsp)
        movq %rdx,nb331_jindex(%rsp)
        movq %rcx,nb331_jjnr(%rsp)
        movq %r8,nb331_shift(%rsp)
        movq %r9,nb331_shiftvec(%rsp)
        movq nb331_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb331_facel(%rsp)

        movq nb331_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb331_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb331_half(%rsp)
        movl %ebx,nb331_half+4(%rsp)
        movsd nb331_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb331_half(%rsp)
        movapd %xmm2,nb331_two(%rsp)
        movapd %xmm3,nb331_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb331_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb331_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movsd nb331_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb331_iqO(%rsp)
        movapd %xmm4,nb331_iqH(%rsp)

        movq  nb331_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb331_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331_ntia(%rsp)

_nb_kernel331_x86_64_sse2.nb331_threadloop: 
        movq  nb331_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel331_x86_64_sse2.nb331_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331_x86_64_sse2.nb331_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331_n(%rsp)
        movl %ebx,nb331_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331_x86_64_sse2.nb331_outerstart
        jmp _nb_kernel331_x86_64_sse2.nb331_end

_nb_kernel331_x86_64_sse2.nb331_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331_nouter(%rsp),%ebx
        movl %ebx,nb331_nouter(%rsp)

_nb_kernel331_x86_64_sse2.nb331_outer: 
        movq  nb331_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb331_is3(%rsp)      ## store is3 

        movq  nb331_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb331_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb331_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb331_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb331_ixO(%rsp)
        movapd %xmm4,nb331_iyO(%rsp)
        movapd %xmm5,nb331_izO(%rsp)

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
        movapd %xmm0,nb331_ixH1(%rsp)
        movapd %xmm1,nb331_iyH1(%rsp)
        movapd %xmm2,nb331_izH1(%rsp)
        movapd %xmm3,nb331_ixH2(%rsp)
        movapd %xmm4,nb331_iyH2(%rsp)
        movapd %xmm5,nb331_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb331_vctot(%rsp)
        movapd %xmm4,nb331_Vvdwtot(%rsp)
        movapd %xmm4,nb331_fixO(%rsp)
        movapd %xmm4,nb331_fiyO(%rsp)
        movapd %xmm4,nb331_fizO(%rsp)
        movapd %xmm4,nb331_fixH1(%rsp)
        movapd %xmm4,nb331_fiyH1(%rsp)
        movapd %xmm4,nb331_fizH1(%rsp)
        movapd %xmm4,nb331_fixH2(%rsp)
        movapd %xmm4,nb331_fiyH2(%rsp)
        movapd %xmm4,nb331_fizH2(%rsp)

        movq  nb331_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb331_pos(%rbp),%rsi
        movq  nb331_faction(%rbp),%rdi
        movq  nb331_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb331_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb331_ninner(%rsp),%ecx
        movl  %ecx,nb331_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb331_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel331_x86_64_sse2.nb331_unroll_loop
        jmp   _nb_kernel331_x86_64_sse2.nb331_checksingle
_nb_kernel331_x86_64_sse2.nb331_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb331_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb331_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb331_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331_iqO(%rsp),%xmm3
        mulpd  nb331_iqH(%rsp),%xmm4
        movapd  %xmm3,nb331_qqO(%rsp)
        movapd  %xmm4,nb331_qqH(%rsp)

        movq nb331_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb331_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb331_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb331_c6(%rsp)
        movapd %xmm6,nb331_c12(%rsp)

        movq nb331_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb331_ixO(%rsp),%xmm0
    subpd nb331_iyO(%rsp),%xmm1
    subpd nb331_izO(%rsp),%xmm2
    subpd nb331_ixH1(%rsp),%xmm3
    subpd nb331_iyH1(%rsp),%xmm4
    subpd nb331_izH1(%rsp),%xmm5
    subpd nb331_ixH2(%rsp),%xmm6
    subpd nb331_iyH2(%rsp),%xmm7
    subpd nb331_izH2(%rsp),%xmm8

        movapd %xmm0,nb331_dxO(%rsp)
        movapd %xmm1,nb331_dyO(%rsp)
        movapd %xmm2,nb331_dzO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb331_dxH1(%rsp)
        movapd %xmm4,nb331_dyH1(%rsp)
        movapd %xmm5,nb331_dzH1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb331_dxH2(%rsp)
        movapd %xmm7,nb331_dyH2(%rsp)
        movapd %xmm8,nb331_dzH2(%rsp)
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

        movapd  nb331_three(%rsp),%xmm9
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

        movapd  nb331_half(%rsp),%xmm15
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

        movapd  nb331_three(%rsp),%xmm1
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

        movapd  nb331_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvO
        mulpd   %xmm15,%xmm10 ##   rinvH1
    mulpd   %xmm15,%xmm11 ##   rinvH2

        movapd  %xmm9,nb331_rinvO(%rsp)
        movapd  %xmm10,nb331_rinvH1(%rsp)
        movapd  %xmm11,nb331_rinvH2(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm9,nb331_rinvO(%rsp)
    movapd nb331_tsc(%rsp),%xmm1

    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## multiply by three (copy, mult. by two, add back)
    movapd  %xmm1,%xmm10
    movapd  %xmm4,%xmm11
    movapd  %xmm7,%xmm12
    pslld   $1,%xmm1
    pslld   $1,%xmm4
    pslld   $1,%xmm7
    paddd   %xmm10,%xmm1
    paddd   %xmm11,%xmm4
    paddd   %xmm12,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb331_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,%xmm12
    movapd    %xmm3,%xmm13
    movapd    %xmm6,%xmm14

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2 ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  nb331_qqO(%rsp),%xmm1     ## VV*qq = vcoul
    mulpd  nb331_qqH(%rsp),%xmm5
    mulpd  nb331_qqH(%rsp),%xmm9
    mulpd  nb331_qqO(%rsp),%xmm3      ## FF*qq = fij
    mulpd  nb331_qqH(%rsp),%xmm7
    mulpd  nb331_qqH(%rsp),%xmm11

    ## accumulate vctot
    addpd  nb331_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb331_vctot(%rsp)

    movapd %xmm7,%xmm2
    movapd %xmm11,%xmm1

    ## fij coul in xmm3, xmm2, xmm1    

    ## calculate LJ table
    movlpd 32(%rsi,%r8,8),%xmm4
    movlpd 40(%rsi,%r8,8),%xmm5
    movlpd 48(%rsi,%r8,8),%xmm6
    movlpd 56(%rsi,%r8,8),%xmm7
    movlpd 64(%rsi,%r8,8),%xmm8
    movlpd 72(%rsi,%r8,8),%xmm9
    movlpd 80(%rsi,%r8,8),%xmm10
    movlpd 88(%rsi,%r8,8),%xmm11
    movhpd 32(%rsi,%r9,8),%xmm4
    movhpd 40(%rsi,%r9,8),%xmm5
    movhpd 48(%rsi,%r9,8),%xmm6
    movhpd 56(%rsi,%r9,8),%xmm7
    movhpd 64(%rsi,%r9,8),%xmm8
    movhpd 72(%rsi,%r9,8),%xmm9
    movhpd 80(%rsi,%r9,8),%xmm10
    movhpd 88(%rsi,%r9,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    ## xmm12 = epsO

    mulpd  %xmm12,%xmm7   ## Heps
    mulpd  %xmm12,%xmm11
    mulpd  %xmm12,%xmm6  ## Geps
    mulpd  %xmm12,%xmm10
    mulpd  %xmm12,%xmm7  ## Heps2
    mulpd  %xmm12,%xmm11
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
    mulpd  %xmm12,%xmm5 ## eps*Fp
    mulpd  %xmm12,%xmm9
    movapd nb331_c6(%rsp),%xmm12
    movapd nb331_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb331_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb331_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    addpd  %xmm7,%xmm3
    movapd nb331_tsc(%rsp),%xmm10

    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm2
    mulpd  %xmm10,%xmm1

    ## move j forces to xmm11-xmm13
    movq nb331_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm11
        movlpd 8(%rdi,%rax,8),%xmm12
        movlpd 16(%rdi,%rax,8),%xmm13
        movhpd (%rdi,%rbx,8),%xmm11
        movhpd 8(%rdi,%rbx,8),%xmm12
        movhpd 16(%rdi,%rbx,8),%xmm13

    xorpd  %xmm0,%xmm0
    xorpd  %xmm4,%xmm4
    xorpd  %xmm8,%xmm8

    subpd  %xmm3,%xmm0
    subpd  %xmm2,%xmm4
    subpd  %xmm1,%xmm8

    mulpd  nb331_rinvO(%rsp),%xmm0
    mulpd  nb331_rinvH1(%rsp),%xmm4
    mulpd  nb331_rinvH2(%rsp),%xmm8

    movapd %xmm0,%xmm1
    movapd %xmm0,%xmm2
    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm6
    movapd %xmm8,%xmm7

        mulpd nb331_dxO(%rsp),%xmm0
        mulpd nb331_dyO(%rsp),%xmm1
        mulpd nb331_dzO(%rsp),%xmm2
        mulpd nb331_dxH1(%rsp),%xmm3
        mulpd nb331_dyH1(%rsp),%xmm4
        mulpd nb331_dzH1(%rsp),%xmm5
        mulpd nb331_dxH2(%rsp),%xmm6
        mulpd nb331_dyH2(%rsp),%xmm7
        mulpd nb331_dzH2(%rsp),%xmm8

    addpd %xmm0,%xmm11
    addpd %xmm1,%xmm12
    addpd %xmm2,%xmm13
    addpd nb331_fixO(%rsp),%xmm0
    addpd nb331_fiyO(%rsp),%xmm1
    addpd nb331_fizO(%rsp),%xmm2

    addpd %xmm3,%xmm11
    addpd %xmm4,%xmm12
    addpd %xmm5,%xmm13
    addpd nb331_fixH1(%rsp),%xmm3
    addpd nb331_fiyH1(%rsp),%xmm4
    addpd nb331_fizH1(%rsp),%xmm5

    addpd %xmm6,%xmm11
    addpd %xmm7,%xmm12
    addpd %xmm8,%xmm13
    addpd nb331_fixH2(%rsp),%xmm6
    addpd nb331_fiyH2(%rsp),%xmm7
    addpd nb331_fizH2(%rsp),%xmm8

    movapd %xmm0,nb331_fixO(%rsp)
    movapd %xmm1,nb331_fiyO(%rsp)
    movapd %xmm2,nb331_fizO(%rsp)
    movapd %xmm3,nb331_fixH1(%rsp)
    movapd %xmm4,nb331_fiyH1(%rsp)
    movapd %xmm5,nb331_fizH1(%rsp)
    movapd %xmm6,nb331_fixH2(%rsp)
    movapd %xmm7,nb331_fiyH2(%rsp)
    movapd %xmm8,nb331_fizH2(%rsp)

    ## store back j forces from xmm11-xmm13
        movlpd %xmm11,(%rdi,%rax,8)
        movlpd %xmm12,8(%rdi,%rax,8)
        movlpd %xmm13,16(%rdi,%rax,8)
        movhpd %xmm11,(%rdi,%rbx,8)
        movhpd %xmm12,8(%rdi,%rbx,8)
        movhpd %xmm13,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb331_innerk(%rsp)
        jl    _nb_kernel331_x86_64_sse2.nb331_checksingle
        jmp   _nb_kernel331_x86_64_sse2.nb331_unroll_loop
_nb_kernel331_x86_64_sse2.nb331_checksingle: 
        movl  nb331_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel331_x86_64_sse2.nb331_dosingle
        jmp   _nb_kernel331_x86_64_sse2.nb331_updateouterdata
_nb_kernel331_x86_64_sse2.nb331_dosingle: 
        movq  nb331_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb331_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
    mulsd nb331_iqO(%rsp),%xmm3
    mulsd nb331_iqH(%rsp),%xmm4
        movapd  %xmm3,nb331_qqO(%rsp)
        movapd  %xmm4,nb331_qqH(%rsp)

        movq nb331_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb331_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb331_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm6             ## c6a
        movsd 8(%rsi,%r8,8),%xmm7        ## c12a

        movapd %xmm6,nb331_c6(%rsp)
        movapd %xmm7,nb331_c12(%rsp)

        movq nb331_pos(%rbp),%rsi        ## base of pos[] 

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

    subsd nb331_ixO(%rsp),%xmm0
    subsd nb331_iyO(%rsp),%xmm1
    subsd nb331_izO(%rsp),%xmm2
    subsd nb331_ixH1(%rsp),%xmm3
    subsd nb331_iyH1(%rsp),%xmm4
    subsd nb331_izH1(%rsp),%xmm5
    subsd nb331_ixH2(%rsp),%xmm6
    subsd nb331_iyH2(%rsp),%xmm7
    subsd nb331_izH2(%rsp),%xmm8

        movapd %xmm0,nb331_dxO(%rsp)
        movapd %xmm1,nb331_dyO(%rsp)
        movapd %xmm2,nb331_dzO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb331_dxH1(%rsp)
        movapd %xmm4,nb331_dyH1(%rsp)
        movapd %xmm5,nb331_dzH1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movapd %xmm6,nb331_dxH2(%rsp)
        movapd %xmm7,nb331_dyH2(%rsp)
        movapd %xmm8,nb331_dzH2(%rsp)
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

        movapd  nb331_three(%rsp),%xmm9
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

        movapd  nb331_half(%rsp),%xmm15
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

        movapd  nb331_three(%rsp),%xmm1
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

        movapd  nb331_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvO
        mulsd   %xmm15,%xmm10 ##   rinvH1
    mulsd   %xmm15,%xmm11 ##   rinvH2

        movapd  %xmm9,nb331_rinvO(%rsp)
        movapd  %xmm10,nb331_rinvH1(%rsp)
        movapd  %xmm11,nb331_rinvH2(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm9,nb331_rinvO(%rsp)
    movapd nb331_tsc(%rsp),%xmm1

    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shll  $2,%r8d
    shll  $2,%r10d
    shll  $2,%r12d

    ## multiply by 3
        lea  (%r8,%r8,2),%r8
        lea  (%r10,%r10,2),%r10
        lea  (%r12,%r12,2),%r12

    movq nb331_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movapd    %xmm0,%xmm12
    movapd    %xmm3,%xmm13
    movapd    %xmm6,%xmm14

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2 ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  nb331_qqO(%rsp),%xmm1     ## VV*qq = vcoul
    mulsd  nb331_qqH(%rsp),%xmm5
    mulsd  nb331_qqH(%rsp),%xmm9
    mulsd  nb331_qqO(%rsp),%xmm3      ## FF*qq = fij
    mulsd  nb331_qqH(%rsp),%xmm7
    mulsd  nb331_qqH(%rsp),%xmm11

    ## accumulate vctot
    addsd  nb331_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb331_vctot(%rsp)

    movapd %xmm7,%xmm2
    movapd %xmm11,%xmm1

    ## fij coul in xmm3, xmm2, xmm1    

    ## calculate LJ table
    movsd 32(%rsi,%r8,8),%xmm4
    movsd 40(%rsi,%r8,8),%xmm5
    movsd 48(%rsi,%r8,8),%xmm6
    movsd 56(%rsi,%r8,8),%xmm7
    movsd 64(%rsi,%r8,8),%xmm8
    movsd 72(%rsi,%r8,8),%xmm9
    movsd 80(%rsi,%r8,8),%xmm10
    movsd 88(%rsi,%r8,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    ## xmm12 = epsO

    mulsd  %xmm12,%xmm7   ## Heps
    mulsd  %xmm12,%xmm11
    mulsd  %xmm12,%xmm6  ## Geps
    mulsd  %xmm12,%xmm10
    mulsd  %xmm12,%xmm7  ## Heps2
    mulsd  %xmm12,%xmm11
    addpd  %xmm6,%xmm5 ## F+Geps
    addsd  %xmm10,%xmm9
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm11,%xmm9
    addsd  %xmm7,%xmm7   ## 2*Heps2
    addsd  %xmm11,%xmm11
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm10,%xmm11

    addsd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm5 ## eps*Fp
    mulsd  %xmm12,%xmm9
    movapd nb331_c6(%rsp),%xmm12
    movapd nb331_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb331_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb331_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    addsd  %xmm7,%xmm3
    movapd nb331_tsc(%rsp),%xmm10

    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm2
    mulsd  %xmm10,%xmm1

    ## move j forces to xmm11-xmm13
    movq nb331_faction(%rbp),%rdi
        movsd (%rdi,%rax,8),%xmm11
        movsd 8(%rdi,%rax,8),%xmm12
        movsd 16(%rdi,%rax,8),%xmm13

    xorpd  %xmm0,%xmm0
    xorpd  %xmm4,%xmm4
    xorpd  %xmm8,%xmm8

    subsd  %xmm3,%xmm0
    subsd  %xmm2,%xmm4
    subsd  %xmm1,%xmm8

    mulsd  nb331_rinvO(%rsp),%xmm0
    mulsd  nb331_rinvH1(%rsp),%xmm4
    mulsd  nb331_rinvH2(%rsp),%xmm8

    movapd %xmm0,%xmm1
    movapd %xmm0,%xmm2
    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm6
    movapd %xmm8,%xmm7

        mulsd nb331_dxO(%rsp),%xmm0
        mulsd nb331_dyO(%rsp),%xmm1
        mulsd nb331_dzO(%rsp),%xmm2
        mulsd nb331_dxH1(%rsp),%xmm3
        mulsd nb331_dyH1(%rsp),%xmm4
        mulsd nb331_dzH1(%rsp),%xmm5
        mulsd nb331_dxH2(%rsp),%xmm6
        mulsd nb331_dyH2(%rsp),%xmm7
        mulsd nb331_dzH2(%rsp),%xmm8

    addsd %xmm0,%xmm11
    addsd %xmm1,%xmm12
    addsd %xmm2,%xmm13
    addsd nb331_fixO(%rsp),%xmm0
    addsd nb331_fiyO(%rsp),%xmm1
    addsd nb331_fizO(%rsp),%xmm2

    addsd %xmm3,%xmm11
    addsd %xmm4,%xmm12
    addsd %xmm5,%xmm13
    addsd nb331_fixH1(%rsp),%xmm3
    addsd nb331_fiyH1(%rsp),%xmm4
    addsd nb331_fizH1(%rsp),%xmm5

    addsd %xmm6,%xmm11
    addsd %xmm7,%xmm12
    addsd %xmm8,%xmm13
    addsd nb331_fixH2(%rsp),%xmm6
    addsd nb331_fiyH2(%rsp),%xmm7
    addsd nb331_fizH2(%rsp),%xmm8

    movsd %xmm0,nb331_fixO(%rsp)
    movsd %xmm1,nb331_fiyO(%rsp)
    movsd %xmm2,nb331_fizO(%rsp)
    movsd %xmm3,nb331_fixH1(%rsp)
    movsd %xmm4,nb331_fiyH1(%rsp)
    movsd %xmm5,nb331_fizH1(%rsp)
    movsd %xmm6,nb331_fixH2(%rsp)
    movsd %xmm7,nb331_fiyH2(%rsp)
    movsd %xmm8,nb331_fizH2(%rsp)

    ## store back j forces from xmm11-xmm13
        movsd %xmm11,(%rdi,%rax,8)
        movsd %xmm12,8(%rdi,%rax,8)
        movsd %xmm13,16(%rdi,%rax,8)

_nb_kernel331_x86_64_sse2.nb331_updateouterdata: 
        movl  nb331_ii3(%rsp),%ecx
        movq  nb331_faction(%rbp),%rdi
        movq  nb331_fshift(%rbp),%rsi
        movl  nb331_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb331_fixO(%rsp),%xmm0
        movapd nb331_fiyO(%rsp),%xmm1
        movapd nb331_fizO(%rsp),%xmm2

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
        movapd nb331_fixH1(%rsp),%xmm0
        movapd nb331_fiyH1(%rsp),%xmm1
        movapd nb331_fizH1(%rsp),%xmm2

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
        movapd nb331_fixH2(%rsp),%xmm0
        movapd nb331_fiyH2(%rsp),%xmm1
        movapd nb331_fizH2(%rsp),%xmm2

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
        movl nb331_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb331_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb331_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb331_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb331_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb331_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb331_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331_x86_64_sse2.nb331_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331_n(%rsp)
        jmp _nb_kernel331_x86_64_sse2.nb331_outer
_nb_kernel331_x86_64_sse2.nb331_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331_x86_64_sse2.nb331_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331_x86_64_sse2.nb331_threadloop
_nb_kernel331_x86_64_sse2.nb331_end: 
        movl nb331_nouter(%rsp),%eax
        movl nb331_ninner(%rsp),%ebx
        movq nb331_outeriter(%rbp),%rcx
        movq nb331_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $872,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel331nf_x86_64_sse2
.globl _nb_kernel331nf_x86_64_sse2
nb_kernel331nf_x86_64_sse2:     
_nb_kernel331nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb331nf_fshift, 16
.set nb331nf_gid, 24
.set nb331nf_pos, 32
.set nb331nf_faction, 40
.set nb331nf_charge, 48
.set nb331nf_p_facel, 56
.set nb331nf_argkrf, 64
.set nb331nf_argcrf, 72
.set nb331nf_Vc, 80
.set nb331nf_type, 88
.set nb331nf_p_ntype, 96
.set nb331nf_vdwparam, 104
.set nb331nf_Vvdw, 112
.set nb331nf_p_tabscale, 120
.set nb331nf_VFtab, 128
.set nb331nf_invsqrta, 136
.set nb331nf_dvda, 144
.set nb331nf_p_gbtabscale, 152
.set nb331nf_GBtab, 160
.set nb331nf_p_nthreads, 168
.set nb331nf_count, 176
.set nb331nf_mtx, 184
.set nb331nf_outeriter, 192
.set nb331nf_inneriter, 200
.set nb331nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb331nf_ixO, 0
.set nb331nf_iyO, 16
.set nb331nf_izO, 32
.set nb331nf_ixH1, 48
.set nb331nf_iyH1, 64
.set nb331nf_izH1, 80
.set nb331nf_ixH2, 96
.set nb331nf_iyH2, 112
.set nb331nf_izH2, 128
.set nb331nf_iqO, 144
.set nb331nf_iqH, 160
.set nb331nf_qqO, 176
.set nb331nf_qqH, 192
.set nb331nf_rinvO, 208
.set nb331nf_rinvH1, 224
.set nb331nf_rinvH2, 240
.set nb331nf_rO, 256
.set nb331nf_rH1, 272
.set nb331nf_rH2, 288
.set nb331nf_tsc, 304
.set nb331nf_c6, 320
.set nb331nf_c12, 336
.set nb331nf_vctot, 352
.set nb331nf_Vvdwtot, 368
.set nb331nf_half, 384
.set nb331nf_three, 400
.set nb331nf_is3, 416
.set nb331nf_ii3, 420
.set nb331nf_nri, 424
.set nb331nf_iinr, 432
.set nb331nf_jindex, 440
.set nb331nf_jjnr, 448
.set nb331nf_shift, 456
.set nb331nf_shiftvec, 464
.set nb331nf_facel, 472
.set nb331nf_innerjjnr, 480
.set nb331nf_ntia, 488
.set nb331nf_innerk, 492
.set nb331nf_n, 496
.set nb331nf_nn1, 500
.set nb331nf_nouter, 504
.set nb331nf_ninner, 508
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $520,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb331nf_nouter(%rsp)
        movl %eax,nb331nf_ninner(%rsp)



        movl (%rdi),%edi
        movl %edi,nb331nf_nri(%rsp)
        movq %rsi,nb331nf_iinr(%rsp)
        movq %rdx,nb331nf_jindex(%rsp)
        movq %rcx,nb331nf_jjnr(%rsp)
        movq %r8,nb331nf_shift(%rsp)
        movq %r9,nb331nf_shiftvec(%rsp)
        movq nb331nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb331nf_facel(%rsp)

        movq nb331nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb331nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb331nf_half(%rsp)
        movl %ebx,nb331nf_half+4(%rsp)
        movsd nb331nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb331nf_half(%rsp)
        movapd %xmm3,nb331nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb331nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb331nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movq nb331nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb331nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb331nf_iqO(%rsp)
        movapd %xmm4,nb331nf_iqH(%rsp)

        movq  nb331nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb331nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb331nf_ntia(%rsp)
_nb_kernel331nf_x86_64_sse2.nb331nf_threadloop: 
        movq  nb331nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel331nf_x86_64_sse2.nb331nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel331nf_x86_64_sse2.nb331nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb331nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb331nf_n(%rsp)
        movl %ebx,nb331nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel331nf_x86_64_sse2.nb331nf_outerstart
        jmp _nb_kernel331nf_x86_64_sse2.nb331nf_end

_nb_kernel331nf_x86_64_sse2.nb331nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb331nf_nouter(%rsp),%ebx
        movl %ebx,nb331nf_nouter(%rsp)

_nb_kernel331nf_x86_64_sse2.nb331nf_outer: 
        movq  nb331nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb331nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb331nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb331nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb331nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb331nf_ixO(%rsp)
        movapd %xmm4,nb331nf_iyO(%rsp)
        movapd %xmm5,nb331nf_izO(%rsp)

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
        movapd %xmm0,nb331nf_ixH1(%rsp)
        movapd %xmm1,nb331nf_iyH1(%rsp)
        movapd %xmm2,nb331nf_izH1(%rsp)
        movapd %xmm3,nb331nf_ixH2(%rsp)
        movapd %xmm4,nb331nf_iyH2(%rsp)
        movapd %xmm5,nb331nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb331nf_vctot(%rsp)
        movapd %xmm4,nb331nf_Vvdwtot(%rsp)

        movq  nb331nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb331nf_pos(%rbp),%rsi
        movq  nb331nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb331nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb331nf_ninner(%rsp),%ecx
        movl  %ecx,nb331nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb331nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel331nf_x86_64_sse2.nb331nf_unroll_loop
        jmp   _nb_kernel331nf_x86_64_sse2.nb331nf_checksingle
_nb_kernel331nf_x86_64_sse2.nb331nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb331nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb331nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb331nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331nf_iqO(%rsp),%xmm3
        mulpd  nb331nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb331nf_qqO(%rsp)
        movapd  %xmm4,nb331nf_qqH(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movq nb331nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb331nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb331nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb331nf_c6(%rsp)
        movapd %xmm6,nb331nf_c12(%rsp)

        movq nb331nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb331nf_ixO(%rsp),%xmm4
        movapd nb331nf_iyO(%rsp),%xmm5
        movapd nb331nf_izO(%rsp),%xmm6

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

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331nf_ixH1(%rsp),%xmm4
        movapd nb331nf_iyH1(%rsp),%xmm5
        movapd nb331nf_izH1(%rsp),%xmm6

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
        movapd nb331nf_ixH2(%rsp),%xmm3
        movapd nb331nf_iyH2(%rsp),%xmm4
        movapd nb331nf_izH2(%rsp),%xmm5

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
        movapd  nb331nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb331nf_rinvO(%rsp)       ## rinvO in xmm4 
        mulpd   %xmm4,%xmm7
        movapd  %xmm7,nb331nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH1(%rsp)       ## rinvH1 
        mulpd  %xmm4,%xmm6
        movapd %xmm6,nb331nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb331nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subpd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH2(%rsp)   ## rinv 
        mulpd %xmm4,%xmm5
        movapd %xmm5,nb331nf_rH2(%rsp)   ## r 

        ## do O interactions 
        ## rO is still in xmm7 
        movapd nb331nf_rinvO(%rsp),%xmm0
        mulpd   nb331nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movq nb331nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now) 
        lea  (%rbx,%rbx,2),%rbx

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7  
    mulpd  %xmm1,%xmm6      ## xmm6=Geps 
    mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
    addpd  %xmm6,%xmm5
    addpd  %xmm7,%xmm5      ## xmm5=Fp        
    movapd nb331nf_qqO(%rsp),%xmm0
    mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addpd  %xmm4,%xmm5 ## xmm5=VV 
    mulpd  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb331nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb331nf_vctot(%rsp)

        ## Dispersion 
        movapd 32(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movapd 32(%rsi,%rbx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 48(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 48(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb331nf_c6(%rsp),%xmm5    ## Vvdw6 

        addpd  nb331nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb331nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 64(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movapd 64(%rsi,%rbx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 80(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 80(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb331nf_c12(%rsp),%xmm5   ## Vvdw12 
        addpd  nb331nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb331nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb331nf_rH1(%rsp),%xmm7
        mulpd nb331nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb331nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   
        lea  (%rbx,%rbx,2),%rbx

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul 
    addpd  nb331nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb331nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331nf_rH2(%rsp),%xmm7
        mulpd   nb331nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb331nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)
        lea  (%rbx,%rbx,2),%rbx

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addpd  nb331nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb331nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb331nf_innerk(%rsp)
        jl    _nb_kernel331nf_x86_64_sse2.nb331nf_checksingle
        jmp   _nb_kernel331nf_x86_64_sse2.nb331nf_unroll_loop
_nb_kernel331nf_x86_64_sse2.nb331nf_checksingle: 
        movl  nb331nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel331nf_x86_64_sse2.nb331nf_dosingle
        jmp   _nb_kernel331nf_x86_64_sse2.nb331nf_updateouterdata
_nb_kernel331nf_x86_64_sse2.nb331nf_dosingle: 
        movq  nb331nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb331nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb331nf_iqO(%rsp),%xmm3
        mulpd  nb331nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb331nf_qqO(%rsp)
        movapd  %xmm4,nb331nf_qqH(%rsp)

        movq nb331nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb331nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb331nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb331nf_c6(%rsp)
        movapd %xmm6,nb331nf_c12(%rsp)

        movq nb331nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb331nf_ixO(%rsp),%xmm4
        movapd nb331nf_iyO(%rsp),%xmm5
        movapd nb331nf_izO(%rsp),%xmm6

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

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb331nf_ixH1(%rsp),%xmm4
        movapd nb331nf_iyH1(%rsp),%xmm5
        movapd nb331nf_izH1(%rsp),%xmm6

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
        movapd nb331nf_ixH2(%rsp),%xmm3
        movapd nb331nf_iyH2(%rsp),%xmm4
        movapd nb331nf_izH2(%rsp),%xmm5

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
        movapd  nb331nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm7,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,nb331nf_rinvO(%rsp)       ## rinvO in xmm4 
        mulsd   %xmm4,%xmm7
        movapd  %xmm7,nb331nf_rO(%rsp)          ## r in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH1(%rsp)       ## rinvH1 
        mulsd  %xmm4,%xmm6
        movapd %xmm6,nb331nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb331nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb331nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu 
        movapd nb331nf_three(%rsp),%xmm4
        subsd %xmm2,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb331nf_half(%rsp),%xmm4   ## rinv 
        movapd %xmm4,nb331nf_rinvH2(%rsp)   ## rinv 
        mulsd %xmm4,%xmm5
        movapd %xmm5,nb331nf_rH2(%rsp)   ## r 

        ## do O interactions 
        movd %eax,%mm0
        ## rO is still in xmm7 
        mulsd   nb331nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm7,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax
        movq nb331nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7  
    mulsd  %xmm1,%xmm6      ## xmm6=Geps 
    mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
    addsd  %xmm6,%xmm5  ## F+Geps 
    addsd  %xmm7,%xmm5      ## xmm5=Fp=F+Geps+Heps2        
    mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addsd  %xmm4,%xmm5 ## xmm5=VV 

    movapd  nb331nf_qqO(%rsp),%xmm0
    mulsd  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb331nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb331nf_vctot(%rsp)

        ## Dispersion 
        movapd 32(%rsi,%rax,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 48(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb331nf_c6(%rsp),%xmm5    ## Vvdw6 

        addsd  nb331nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb331nf_Vvdwtot(%rsp)

        ## Repulsion 
        movapd 64(%rsi,%rax,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 80(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb331nf_c12(%rsp),%xmm5   ## Vvdw12 

        addsd  nb331nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb331nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb331nf_rH1(%rsp),%xmm7
        mulpd nb331nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb331nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addsd  nb331nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb331nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movapd nb331nf_rH2(%rsp),%xmm7
        mulsd   nb331nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb331nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax ## idx *= 3 (total *=12 now)   

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb331nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addsd  nb331nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb331nf_vctot(%rsp)

_nb_kernel331nf_x86_64_sse2.nb331nf_updateouterdata: 
        ## get n from stack
        movl nb331nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb331nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb331nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb331nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb331nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb331nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb331nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel331nf_x86_64_sse2.nb331nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb331nf_n(%rsp)
        jmp _nb_kernel331nf_x86_64_sse2.nb331nf_outer
_nb_kernel331nf_x86_64_sse2.nb331nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb331nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel331nf_x86_64_sse2.nb331nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel331nf_x86_64_sse2.nb331nf_threadloop
_nb_kernel331nf_x86_64_sse2.nb331nf_end: 
        movl nb331nf_nouter(%rsp),%eax
        movl nb331nf_ninner(%rsp),%ebx
        movq nb331nf_outeriter(%rbp),%rcx
        movq nb331nf_inneriter(%rbp),%rdx
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

