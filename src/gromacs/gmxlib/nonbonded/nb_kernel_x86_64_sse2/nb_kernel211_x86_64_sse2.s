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






.globl nb_kernel211_x86_64_sse2
.globl _nb_kernel211_x86_64_sse2
nb_kernel211_x86_64_sse2:       
_nb_kernel211_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb211_fshift, 16
.set nb211_gid, 24
.set nb211_pos, 32
.set nb211_faction, 40
.set nb211_charge, 48
.set nb211_p_facel, 56
.set nb211_argkrf, 64
.set nb211_argcrf, 72
.set nb211_Vc, 80
.set nb211_type, 88
.set nb211_p_ntype, 96
.set nb211_vdwparam, 104
.set nb211_Vvdw, 112
.set nb211_p_tabscale, 120
.set nb211_VFtab, 128
.set nb211_invsqrta, 136
.set nb211_dvda, 144
.set nb211_p_gbtabscale, 152
.set nb211_GBtab, 160
.set nb211_p_nthreads, 168
.set nb211_count, 176
.set nb211_mtx, 184
.set nb211_outeriter, 192
.set nb211_inneriter, 200
.set nb211_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb211_ixO, 0
.set nb211_iyO, 16
.set nb211_izO, 32
.set nb211_ixH1, 48
.set nb211_iyH1, 64
.set nb211_izH1, 80
.set nb211_ixH2, 96
.set nb211_iyH2, 112
.set nb211_izH2, 128
.set nb211_iqO, 144
.set nb211_iqH, 160
.set nb211_dxO, 176
.set nb211_dyO, 192
.set nb211_dzO, 208
.set nb211_dxH1, 224
.set nb211_dyH1, 240
.set nb211_dzH1, 256
.set nb211_dxH2, 272
.set nb211_dyH2, 288
.set nb211_dzH2, 304
.set nb211_qqO, 320
.set nb211_qqH, 336
.set nb211_c6, 352
.set nb211_c12, 368
.set nb211_six, 384
.set nb211_twelve, 400
.set nb211_vctot, 416
.set nb211_Vvdwtot, 432
.set nb211_fixO, 448
.set nb211_fiyO, 464
.set nb211_fizO, 480
.set nb211_fixH1, 496
.set nb211_fiyH1, 512
.set nb211_fizH1, 528
.set nb211_fixH2, 544
.set nb211_fiyH2, 560
.set nb211_fizH2, 576
.set nb211_fjx, 592
.set nb211_fjy, 608
.set nb211_fjz, 624
.set nb211_half, 640
.set nb211_three, 656
.set nb211_two, 672
.set nb211_krf, 688
.set nb211_crf, 704
.set nb211_krsqO, 720
.set nb211_krsqH1, 736
.set nb211_krsqH2, 752
.set nb211_nri, 768
.set nb211_iinr, 776
.set nb211_jindex, 784
.set nb211_jjnr, 792
.set nb211_shift, 800
.set nb211_shiftvec, 808
.set nb211_facel, 816
.set nb211_innerjjnr, 824
.set nb211_is3, 832
.set nb211_ii3, 836
.set nb211_ntia, 840
.set nb211_innerk, 844
.set nb211_n, 848
.set nb211_nn1, 852
.set nb211_nouter, 856
.set nb211_ninner, 860
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
        movl %eax,nb211_nouter(%rsp)
        movl %eax,nb211_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb211_nri(%rsp)
        movq %rsi,nb211_iinr(%rsp)
        movq %rdx,nb211_jindex(%rsp)
        movq %rcx,nb211_jjnr(%rsp)
        movq %r8,nb211_shift(%rsp)
        movq %r9,nb211_shiftvec(%rsp)
        movq nb211_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb211_facel(%rsp)

        movq nb211_argkrf(%rbp),%rsi
        movq nb211_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb211_krf(%rsp)
        movapd %xmm2,nb211_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb211_half(%rsp)
        movl %ebx,nb211_half+4(%rsp)
        movsd nb211_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb211_half(%rsp)
        movapd %xmm2,nb211_two(%rsp)
        movapd %xmm3,nb211_three(%rsp)
        movapd %xmm4,nb211_six(%rsp)
        movapd %xmm5,nb211_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb211_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb211_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4

        movsd nb211_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb211_iqO(%rsp)
        movapd %xmm4,nb211_iqH(%rsp)

        movq  nb211_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb211_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211_ntia(%rsp)
_nb_kernel211_x86_64_sse2.nb211_threadloop: 
        movq  nb211_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel211_x86_64_sse2.nb211_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211_x86_64_sse2.nb211_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211_n(%rsp)
        movl %ebx,nb211_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211_x86_64_sse2.nb211_outerstart
        jmp _nb_kernel211_x86_64_sse2.nb211_end

_nb_kernel211_x86_64_sse2.nb211_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211_nouter(%rsp),%ebx
        movl %ebx,nb211_nouter(%rsp)

_nb_kernel211_x86_64_sse2.nb211_outer: 
        movq  nb211_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb211_is3(%rsp)      ## store is3 

        movq  nb211_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb211_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb211_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb211_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb211_ixO(%rsp)
        movapd %xmm4,nb211_iyO(%rsp)
        movapd %xmm5,nb211_izO(%rsp)

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
        movapd %xmm0,nb211_ixH1(%rsp)
        movapd %xmm1,nb211_iyH1(%rsp)
        movapd %xmm2,nb211_izH1(%rsp)
        movapd %xmm3,nb211_ixH2(%rsp)
        movapd %xmm4,nb211_iyH2(%rsp)
        movapd %xmm5,nb211_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb211_vctot(%rsp)
        movapd %xmm4,nb211_Vvdwtot(%rsp)
        movapd %xmm4,nb211_fixO(%rsp)
        movapd %xmm4,nb211_fiyO(%rsp)
        movapd %xmm4,nb211_fizO(%rsp)
        movapd %xmm4,nb211_fixH1(%rsp)
        movapd %xmm4,nb211_fiyH1(%rsp)
        movapd %xmm4,nb211_fizH1(%rsp)
        movapd %xmm4,nb211_fixH2(%rsp)
        movapd %xmm4,nb211_fiyH2(%rsp)
        movapd %xmm4,nb211_fizH2(%rsp)

        movq  nb211_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb211_pos(%rbp),%rsi
        movq  nb211_faction(%rbp),%rdi
        movq  nb211_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb211_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb211_ninner(%rsp),%ecx
        movl  %ecx,nb211_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb211_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel211_x86_64_sse2.nb211_unroll_loop
        jmp   _nb_kernel211_x86_64_sse2.nb211_checksingle
_nb_kernel211_x86_64_sse2.nb211_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb211_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb211_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb211_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211_iqO(%rsp),%xmm3
        mulpd  nb211_iqH(%rsp),%xmm4

        movapd  %xmm3,nb211_qqO(%rsp)
        movapd  %xmm4,nb211_qqH(%rsp)

        movq nb211_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb211_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb211_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb211_c6(%rsp)
        movapd %xmm6,nb211_c12(%rsp)

        movq nb211_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb211_ixO(%rsp),%xmm0
    subpd nb211_iyO(%rsp),%xmm1
    subpd nb211_izO(%rsp),%xmm2
    subpd nb211_ixH1(%rsp),%xmm3
    subpd nb211_iyH1(%rsp),%xmm4
    subpd nb211_izH1(%rsp),%xmm5
    subpd nb211_ixH2(%rsp),%xmm6
    subpd nb211_iyH2(%rsp),%xmm7
    subpd nb211_izH2(%rsp),%xmm8

        movapd %xmm0,nb211_dxO(%rsp)
        movapd %xmm1,nb211_dyO(%rsp)
        movapd %xmm2,nb211_dzO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb211_dxH1(%rsp)
        movapd %xmm4,nb211_dyH1(%rsp)
        movapd %xmm5,nb211_dzH1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb211_dxH2(%rsp)
        movapd %xmm7,nb211_dyH2(%rsp)
        movapd %xmm8,nb211_dzH2(%rsp)
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

        movapd  nb211_three(%rsp),%xmm9
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

        movapd  nb211_half(%rsp),%xmm15
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

        movapd  nb211_three(%rsp),%xmm1
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

        movapd  nb211_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvO
        mulpd   %xmm15,%xmm10 ##   rinvH1
    mulpd   %xmm15,%xmm11 ##   rinvH2

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb211_krf(%rsp),%xmm2
    mulpd  %xmm9,%xmm9  ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  %xmm2,%xmm0 ## k*rsq
    mulpd  %xmm2,%xmm3
    mulpd  %xmm2,%xmm6
    movapd %xmm0,%xmm2 ## copy of k*rsq
    movapd %xmm3,%xmm5
    movapd %xmm6,%xmm8
    addpd  %xmm1,%xmm2 ## rinv+krsq
    addpd  %xmm4,%xmm5
    addpd  %xmm7,%xmm8
    movapd %xmm9,%xmm12
    mulpd  %xmm12,%xmm12 ## rinv4
    mulpd  %xmm9,%xmm12 ## rinv6
    subpd  nb211_crf(%rsp),%xmm2     ## rinv+krsq-crf
    subpd  nb211_crf(%rsp),%xmm5
    subpd  nb211_crf(%rsp),%xmm8
    mulpd  nb211_qqO(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf)
    mulpd  nb211_qqH(%rsp),%xmm5   ## voul=qq*(rinv+ krsq-crf)
    mulpd  nb211_qqH(%rsp),%xmm8   ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    movapd %xmm12,%xmm13 ## rinv6
    mulpd %xmm12,%xmm12 ## rinv12
        mulpd  nb211_c6(%rsp),%xmm13
        mulpd  nb211_c12(%rsp),%xmm12
    movapd %xmm12,%xmm14
    subpd  %xmm13,%xmm14
    mulpd  nb211_qqO(%rsp),%xmm1     ## (rinv-2*krsq)*qq
    mulpd  nb211_qqH(%rsp),%xmm4
    mulpd  nb211_qqH(%rsp),%xmm7
    addpd  nb211_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb211_vctot(%rsp)

        addpd  nb211_Vvdwtot(%rsp),%xmm14
        mulpd  nb211_six(%rsp),%xmm13
        mulpd  nb211_twelve(%rsp),%xmm12
        movapd %xmm14,nb211_Vvdwtot(%rsp)
    subpd  %xmm13,%xmm12 ## LJ fscal        

    addpd %xmm12,%xmm1

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j forces to xmm0-xmm2
    movq nb211_faction(%rbp),%rdi
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

        mulpd nb211_dxO(%rsp),%xmm7
        mulpd nb211_dyO(%rsp),%xmm8
        mulpd nb211_dzO(%rsp),%xmm9
        mulpd nb211_dxH1(%rsp),%xmm10
        mulpd nb211_dyH1(%rsp),%xmm11
        mulpd nb211_dzH1(%rsp),%xmm12
        mulpd nb211_dxH2(%rsp),%xmm13
        mulpd nb211_dyH2(%rsp),%xmm14
        mulpd nb211_dzH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb211_fixO(%rsp),%xmm7
    addpd nb211_fiyO(%rsp),%xmm8
    addpd nb211_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb211_fixH1(%rsp),%xmm10
    addpd nb211_fiyH1(%rsp),%xmm11
    addpd nb211_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb211_fixH2(%rsp),%xmm13
    addpd nb211_fiyH2(%rsp),%xmm14
    addpd nb211_fizH2(%rsp),%xmm15

    movapd %xmm7,nb211_fixO(%rsp)
    movapd %xmm8,nb211_fiyO(%rsp)
    movapd %xmm9,nb211_fizO(%rsp)
    movapd %xmm10,nb211_fixH1(%rsp)
    movapd %xmm11,nb211_fiyH1(%rsp)
    movapd %xmm12,nb211_fizH1(%rsp)
    movapd %xmm13,nb211_fixH2(%rsp)
    movapd %xmm14,nb211_fiyH2(%rsp)
    movapd %xmm15,nb211_fizH2(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb211_innerk(%rsp)
        jl    _nb_kernel211_x86_64_sse2.nb211_checksingle
        jmp   _nb_kernel211_x86_64_sse2.nb211_unroll_loop
_nb_kernel211_x86_64_sse2.nb211_checksingle: 
        movl  nb211_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel211_x86_64_sse2.nb211_dosingle
        jmp   _nb_kernel211_x86_64_sse2.nb211_updateouterdata
_nb_kernel211_x86_64_sse2.nb211_dosingle: 
        movq  nb211_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb211_innerjjnr(%rsp)

        movq nb211_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211_iqO(%rsp),%xmm3
        mulpd  nb211_iqH(%rsp),%xmm4


        movapd  %xmm3,nb211_qqO(%rsp)
        movapd  %xmm4,nb211_qqH(%rsp)

        movq nb211_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb211_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb211_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb211_c6(%rsp)
        movapd %xmm6,nb211_c12(%rsp)

        movq nb211_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb211_ixO(%rsp),%xmm4
        subsd nb211_iyO(%rsp),%xmm5
        subsd nb211_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxO(%rsp)
        movapd %xmm5,nb211_dyO(%rsp)
        movapd %xmm6,nb211_dzO(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb211_ixH1(%rsp),%xmm4
        subsd nb211_iyH1(%rsp),%xmm5
        subsd nb211_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb211_dxH1(%rsp)
        movapd %xmm5,nb211_dyH1(%rsp)
        movapd %xmm6,nb211_dzH1(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move j coords to xmm3-xmm5
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## calc dr 
        subsd nb211_ixH2(%rsp),%xmm3
        subsd nb211_iyH2(%rsp),%xmm4
        subsd nb211_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb211_dxH2(%rsp)
        movapd %xmm4,nb211_dyH2(%rsp)
        movapd %xmm5,nb211_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb211_krf(%rsp),%xmm0
        mulsd  nb211_krf(%rsp),%xmm1
        mulsd  nb211_krf(%rsp),%xmm2

        movapd %xmm0,nb211_krsqH2(%rsp)
        movapd %xmm1,nb211_krsqH1(%rsp)
        movapd %xmm2,nb211_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb211_c6(%rsp),%xmm1
        mulsd  nb211_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb211_Vvdwtot(%rsp),%xmm3
        mulsd  nb211_six(%rsp),%xmm1
        mulsd  nb211_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2      ## nb part of fs  

        movapd %xmm7,%xmm0
        movapd nb211_krsqO(%rsp),%xmm1
        addsd  %xmm1,%xmm0
        mulsd  nb211_two(%rsp),%xmm1
        subsd  nb211_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subsd  %xmm1,%xmm7
        mulsd  nb211_qqO(%rsp),%xmm0
        mulsd  nb211_qqO(%rsp),%xmm7
        addsd  %xmm7,%xmm2

        mulsd  %xmm2,%xmm4      ## total fsO in xmm4 

        addsd  nb211_vctot(%rsp),%xmm0
        movlpd %xmm3,nb211_Vvdwtot(%rsp)
        movlpd %xmm0,nb211_vctot(%rsp)

        movapd nb211_dxO(%rsp),%xmm0
        movapd nb211_dyO(%rsp),%xmm1
        movapd nb211_dzO(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb211_fixO(%rsp),%xmm3
        movapd nb211_fiyO(%rsp),%xmm4
        movapd nb211_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixO(%rsp)
        movlpd %xmm4,nb211_fiyO(%rsp)
        movlpd %xmm7,nb211_fizO(%rsp)
        ## update j forces with water O 
        movlpd %xmm0,nb211_fjx(%rsp)
        movlpd %xmm1,nb211_fjy(%rsp)
        movlpd %xmm2,nb211_fjz(%rsp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb211_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb211_two(%rsp),%xmm0
        subsd   nb211_crf(%rsp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb211_qqH(%rsp),%xmm6   ## vcoul 
        mulsd   nb211_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb211_vctot(%rsp),%xmm6

        movapd nb211_dxH1(%rsp),%xmm0
        movapd nb211_dyH1(%rsp),%xmm1
        movapd nb211_dzH1(%rsp),%xmm2
        movlpd %xmm6,nb211_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb211_fixH1(%rsp),%xmm3
        movapd nb211_fiyH1(%rsp),%xmm4
        movapd nb211_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixH1(%rsp)
        movlpd %xmm4,nb211_fiyH1(%rsp)
        movlpd %xmm7,nb211_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb211_fjx(%rsp),%xmm0
        addsd  nb211_fjy(%rsp),%xmm1
        addsd  nb211_fjz(%rsp),%xmm2
        movlpd %xmm0,nb211_fjx(%rsp)
        movlpd %xmm1,nb211_fjy(%rsp)
        movlpd %xmm2,nb211_fjz(%rsp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb211_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb211_two(%rsp),%xmm0
        subsd   nb211_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb211_qqH(%rsp),%xmm5   ## vcoul 
        mulsd   nb211_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb211_vctot(%rsp),%xmm5

        movapd nb211_dxH2(%rsp),%xmm0
        movapd nb211_dyH2(%rsp),%xmm1
        movapd nb211_dzH2(%rsp),%xmm2
        movlpd %xmm5,nb211_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb211_fixH2(%rsp),%xmm3
        movapd nb211_fiyH2(%rsp),%xmm4
        movapd nb211_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb211_fixH2(%rsp)
        movlpd %xmm4,nb211_fiyH2(%rsp)
        movlpd %xmm7,nb211_fizH2(%rsp)

        movq nb211_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb211_fjx(%rsp),%xmm0
        addsd  nb211_fjy(%rsp),%xmm1
        addsd  nb211_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel211_x86_64_sse2.nb211_updateouterdata: 
        movl  nb211_ii3(%rsp),%ecx
        movq  nb211_faction(%rbp),%rdi
        movq  nb211_fshift(%rbp),%rsi
        movl  nb211_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb211_fixO(%rsp),%xmm0
        movapd nb211_fiyO(%rsp),%xmm1
        movapd nb211_fizO(%rsp),%xmm2

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
        movapd nb211_fixH1(%rsp),%xmm0
        movapd nb211_fiyH1(%rsp),%xmm1
        movapd nb211_fizH1(%rsp),%xmm2

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
        movapd nb211_fixH2(%rsp),%xmm0
        movapd nb211_fiyH2(%rsp),%xmm1
        movapd nb211_fizH2(%rsp),%xmm2

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
        movl nb211_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb211_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb211_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb211_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb211_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb211_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb211_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211_x86_64_sse2.nb211_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211_n(%rsp)
        jmp _nb_kernel211_x86_64_sse2.nb211_outer
_nb_kernel211_x86_64_sse2.nb211_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211_x86_64_sse2.nb211_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211_x86_64_sse2.nb211_threadloop
_nb_kernel211_x86_64_sse2.nb211_end: 
        movl nb211_nouter(%rsp),%eax
        movl nb211_ninner(%rsp),%ebx
        movq nb211_outeriter(%rbp),%rcx
        movq nb211_inneriter(%rbp),%rdx
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




.globl nb_kernel211nf_x86_64_sse2
.globl _nb_kernel211nf_x86_64_sse2
nb_kernel211nf_x86_64_sse2:     
_nb_kernel211nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb211nf_fshift, 16
.set nb211nf_gid, 24
.set nb211nf_pos, 32
.set nb211nf_faction, 40
.set nb211nf_charge, 48
.set nb211nf_p_facel, 56
.set nb211nf_argkrf, 64
.set nb211nf_argcrf, 72
.set nb211nf_Vc, 80
.set nb211nf_type, 88
.set nb211nf_p_ntype, 96
.set nb211nf_vdwparam, 104
.set nb211nf_Vvdw, 112
.set nb211nf_p_tabscale, 120
.set nb211nf_VFtab, 128
.set nb211nf_invsqrta, 136
.set nb211nf_dvda, 144
.set nb211nf_p_gbtabscale, 152
.set nb211nf_GBtab, 160
.set nb211nf_p_nthreads, 168
.set nb211nf_count, 176
.set nb211nf_mtx, 184
.set nb211nf_outeriter, 192
.set nb211nf_inneriter, 200
.set nb211nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb211nf_ixO, 0
.set nb211nf_iyO, 16
.set nb211nf_izO, 32
.set nb211nf_ixH1, 48
.set nb211nf_iyH1, 64
.set nb211nf_izH1, 80
.set nb211nf_ixH2, 96
.set nb211nf_iyH2, 112
.set nb211nf_izH2, 128
.set nb211nf_iqO, 144
.set nb211nf_iqH, 160
.set nb211nf_qqO, 176
.set nb211nf_qqH, 192
.set nb211nf_c6, 208
.set nb211nf_c12, 224
.set nb211nf_vctot, 240
.set nb211nf_Vvdwtot, 256
.set nb211nf_half, 272
.set nb211nf_three, 288
.set nb211nf_krf, 304
.set nb211nf_crf, 320
.set nb211nf_krsqO, 336
.set nb211nf_krsqH1, 352
.set nb211nf_krsqH2, 368
.set nb211nf_nri, 384
.set nb211nf_iinr, 392
.set nb211nf_jindex, 400
.set nb211nf_jjnr, 408
.set nb211nf_shift, 416
.set nb211nf_shiftvec, 424
.set nb211nf_facel, 432
.set nb211nf_innerjjnr, 440
.set nb211nf_is3, 448
.set nb211nf_ii3, 452
.set nb211nf_ntia, 456
.set nb211nf_innerk, 460
.set nb211nf_n, 464
.set nb211nf_nn1, 468
.set nb211nf_nouter, 472
.set nb211nf_ninner, 476
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $488,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb211nf_nouter(%rsp)
        movl %eax,nb211nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb211nf_nri(%rsp)
        movq %rsi,nb211nf_iinr(%rsp)
        movq %rdx,nb211nf_jindex(%rsp)
        movq %rcx,nb211nf_jjnr(%rsp)
        movq %r8,nb211nf_shift(%rsp)
        movq %r9,nb211nf_shiftvec(%rsp)
        movq nb211nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb211nf_facel(%rsp)

        movq nb211nf_argkrf(%rbp),%rsi
        movq nb211nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb211nf_krf(%rsp)
        movapd %xmm2,nb211nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb211nf_half(%rsp)
        movl %ebx,nb211nf_half+4(%rsp)
        movsd nb211nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb211nf_half(%rsp)
        movapd %xmm3,nb211nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb211nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb211nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd 8(%rdx,%rbx,8),%xmm4
        movsd nb211nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb211nf_iqO(%rsp)
        movapd %xmm4,nb211nf_iqH(%rsp)

        movq  nb211nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb211nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211nf_ntia(%rsp)
_nb_kernel211nf_x86_64_sse2.nb211nf_threadloop: 
        movq  nb211nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel211nf_x86_64_sse2.nb211nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211nf_x86_64_sse2.nb211nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211nf_n(%rsp)
        movl %ebx,nb211nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211nf_x86_64_sse2.nb211nf_outerstart
        jmp _nb_kernel211nf_x86_64_sse2.nb211nf_end

_nb_kernel211nf_x86_64_sse2.nb211nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211nf_nouter(%rsp),%ebx
        movl %ebx,nb211nf_nouter(%rsp)

_nb_kernel211nf_x86_64_sse2.nb211nf_outer: 
        movq  nb211nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb211nf_is3(%rsp)            ## store is3 

        movq  nb211nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb211nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb211nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb211nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb211nf_ixO(%rsp)
        movapd %xmm4,nb211nf_iyO(%rsp)
        movapd %xmm5,nb211nf_izO(%rsp)

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
        movapd %xmm0,nb211nf_ixH1(%rsp)
        movapd %xmm1,nb211nf_iyH1(%rsp)
        movapd %xmm2,nb211nf_izH1(%rsp)
        movapd %xmm3,nb211nf_ixH2(%rsp)
        movapd %xmm4,nb211nf_iyH2(%rsp)
        movapd %xmm5,nb211nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb211nf_vctot(%rsp)
        movapd %xmm4,nb211nf_Vvdwtot(%rsp)

        movq  nb211nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb211nf_pos(%rbp),%rsi
        movq  nb211nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb211nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb211nf_ninner(%rsp),%ecx
        movl  %ecx,nb211nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb211nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel211nf_x86_64_sse2.nb211nf_unroll_loop
        jmp   _nb_kernel211nf_x86_64_sse2.nb211nf_checksingle
_nb_kernel211nf_x86_64_sse2.nb211nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb211nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb211nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb211nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211nf_iqO(%rsp),%xmm3
        mulpd  nb211nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb211nf_qqO(%rsp)
        movapd  %xmm4,nb211nf_qqH(%rsp)

        movq nb211nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb211nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb211nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb211nf_c6(%rsp)
        movapd %xmm6,nb211nf_c12(%rsp)

        movq nb211nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb211nf_ixO(%rsp),%xmm4
        movapd nb211nf_iyO(%rsp),%xmm5
        movapd nb211nf_izO(%rsp),%xmm6

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
        movapd nb211nf_ixH1(%rsp),%xmm4
        movapd nb211nf_iyH1(%rsp),%xmm5
        movapd nb211nf_izH1(%rsp),%xmm6

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
        movapd nb211nf_ixH2(%rsp),%xmm3
        movapd nb211nf_iyH2(%rsp),%xmm4
        movapd nb211nf_izH2(%rsp),%xmm5

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

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulpd  nb211nf_krf(%rsp),%xmm0
        mulpd  nb211nf_krf(%rsp),%xmm1
        mulpd  nb211nf_krf(%rsp),%xmm2

        movapd %xmm0,nb211nf_krsqH2(%rsp)
        movapd %xmm1,nb211nf_krsqH1(%rsp)
        movapd %xmm2,nb211nf_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb211nf_c6(%rsp),%xmm1
        mulpd  nb211nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb211nf_Vvdwtot(%rsp),%xmm3

        movapd %xmm7,%xmm0
        movapd nb211nf_krsqO(%rsp),%xmm1
        addpd  %xmm1,%xmm0
        subpd  nb211nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulpd  nb211nf_qqO(%rsp),%xmm0

        addpd  nb211nf_vctot(%rsp),%xmm0
        movapd %xmm3,nb211nf_Vvdwtot(%rsp)
        movapd %xmm0,nb211nf_vctot(%rsp)

        ## H1 interactions 
        movapd  nb211nf_krsqH1(%rsp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb211nf_crf(%rsp),%xmm6
        mulpd   nb211nf_qqH(%rsp),%xmm6   ## vcoul      
        addpd  nb211nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb211nf_vctot(%rsp)

        ## H2 interactions 
        movapd  nb211nf_krsqH2(%rsp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb211nf_crf(%rsp),%xmm5
        mulpd   nb211nf_qqH(%rsp),%xmm5   ## vcoul 
        addpd  nb211nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb211nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb211nf_innerk(%rsp)
        jl    _nb_kernel211nf_x86_64_sse2.nb211nf_checksingle
        jmp   _nb_kernel211nf_x86_64_sse2.nb211nf_unroll_loop
_nb_kernel211nf_x86_64_sse2.nb211nf_checksingle: 
        movl  nb211nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel211nf_x86_64_sse2.nb211nf_dosingle
        jmp   _nb_kernel211nf_x86_64_sse2.nb211nf_updateouterdata
_nb_kernel211nf_x86_64_sse2.nb211nf_dosingle: 
        movq  nb211nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb211nf_innerjjnr(%rsp)

        movq nb211nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb211nf_iqO(%rsp),%xmm3
        mulpd  nb211nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb211nf_qqO(%rsp)
        movapd  %xmm4,nb211nf_qqH(%rsp)

        movq nb211nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb211nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb211nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb211nf_c6(%rsp)
        movapd %xmm6,nb211nf_c12(%rsp)

        movq nb211nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb211nf_ixO(%rsp),%xmm4
        movapd nb211nf_iyO(%rsp),%xmm5
        movapd nb211nf_izO(%rsp),%xmm6

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
        movapd nb211nf_ixH1(%rsp),%xmm4
        movapd nb211nf_iyH1(%rsp),%xmm5
        movapd nb211nf_izH1(%rsp),%xmm6

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
        movapd nb211nf_ixH2(%rsp),%xmm3
        movapd nb211nf_iyH2(%rsp),%xmm4
        movapd nb211nf_izH2(%rsp),%xmm5

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

        movapd %xmm5,%xmm0
        movapd %xmm6,%xmm1
        movapd %xmm7,%xmm2

        mulsd  nb211nf_krf(%rsp),%xmm0
        mulsd  nb211nf_krf(%rsp),%xmm1
        mulsd  nb211nf_krf(%rsp),%xmm2

        movapd %xmm0,nb211nf_krsqH2(%rsp)
        movapd %xmm1,nb211nf_krsqH1(%rsp)
        movapd %xmm2,nb211nf_krsqO(%rsp)

        ## start with rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb211nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb211nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb211nf_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb211nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb211nf_c6(%rsp),%xmm1
        mulsd  nb211nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb211nf_Vvdwtot(%rsp),%xmm3

        movapd %xmm7,%xmm0
        movapd nb211nf_krsqO(%rsp),%xmm1
        addsd  %xmm1,%xmm0
        subsd  nb211nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulsd  nb211nf_qqO(%rsp),%xmm0
        addsd  nb211nf_vctot(%rsp),%xmm0
        movlpd %xmm3,nb211nf_Vvdwtot(%rsp)
        movlpd %xmm0,nb211nf_vctot(%rsp)

        ## H1 interactions 
        movapd  nb211nf_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb211nf_crf(%rsp),%xmm6
        mulsd   nb211nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd  nb211nf_vctot(%rsp),%xmm6
        movlpd %xmm6,nb211nf_vctot(%rsp)

        ## H2 interactions 
        movapd  nb211nf_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb211nf_crf(%rsp),%xmm5
        mulsd   nb211nf_qqH(%rsp),%xmm5   ## vcoul 
        addsd  nb211nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb211nf_vctot(%rsp)

_nb_kernel211nf_x86_64_sse2.nb211nf_updateouterdata: 
        ## get n from stack
        movl nb211nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb211nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb211nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb211nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb211nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb211nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb211nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211nf_x86_64_sse2.nb211nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211nf_n(%rsp)
        jmp _nb_kernel211nf_x86_64_sse2.nb211nf_outer
_nb_kernel211nf_x86_64_sse2.nb211nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211nf_x86_64_sse2.nb211nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211nf_x86_64_sse2.nb211nf_threadloop
_nb_kernel211nf_x86_64_sse2.nb211nf_end: 
        movl nb211nf_nouter(%rsp),%eax
        movl nb211nf_ninner(%rsp),%ebx
        movq nb211nf_outeriter(%rbp),%rcx
        movq nb211nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $488,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

