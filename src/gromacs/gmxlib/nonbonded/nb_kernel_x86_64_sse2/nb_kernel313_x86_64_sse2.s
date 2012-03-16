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








.globl nb_kernel313_x86_64_sse2
.globl _nb_kernel313_x86_64_sse2
nb_kernel313_x86_64_sse2:       
_nb_kernel313_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb313_fshift, 16
.set nb313_gid, 24
.set nb313_pos, 32
.set nb313_faction, 40
.set nb313_charge, 48
.set nb313_p_facel, 56
.set nb313_argkrf, 64
.set nb313_argcrf, 72
.set nb313_Vc, 80
.set nb313_type, 88
.set nb313_p_ntype, 96
.set nb313_vdwparam, 104
.set nb313_Vvdw, 112
.set nb313_p_tabscale, 120
.set nb313_VFtab, 128
.set nb313_invsqrta, 136
.set nb313_dvda, 144
.set nb313_p_gbtabscale, 152
.set nb313_GBtab, 160
.set nb313_p_nthreads, 168
.set nb313_count, 176
.set nb313_mtx, 184
.set nb313_outeriter, 192
.set nb313_inneriter, 200
.set nb313_work, 208
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb313_ixO, 0
.set nb313_iyO, 16
.set nb313_izO, 32
.set nb313_ixH1, 48
.set nb313_iyH1, 64
.set nb313_izH1, 80
.set nb313_ixH2, 96
.set nb313_iyH2, 112
.set nb313_izH2, 128
.set nb313_ixM, 144
.set nb313_iyM, 160
.set nb313_izM, 176
.set nb313_iqM, 192
.set nb313_iqH, 208
.set nb313_dxO, 224
.set nb313_dyO, 240
.set nb313_dzO, 256
.set nb313_dxH1, 272
.set nb313_dyH1, 288
.set nb313_dzH1, 304
.set nb313_dxH2, 320
.set nb313_dyH2, 336
.set nb313_dzH2, 352
.set nb313_dxM, 368
.set nb313_dyM, 384
.set nb313_dzM, 400
.set nb313_qqM, 416
.set nb313_qqH, 432
.set nb313_rinvsqO, 448
.set nb313_rinvH1, 464
.set nb313_rinvH2, 480
.set nb313_rinvM, 496
.set nb313_rO, 512
.set nb313_rH1, 528
.set nb313_rH2, 544
.set nb313_rM, 560
.set nb313_tsc, 576
.set nb313_two, 592
.set nb313_c6, 608
.set nb313_c12, 624
.set nb313_six, 640
.set nb313_twelve, 656
.set nb313_vctot, 672
.set nb313_Vvdwtot, 688
.set nb313_fixO, 704
.set nb313_fiyO, 720
.set nb313_fizO, 736
.set nb313_fixH1, 752
.set nb313_fiyH1, 768
.set nb313_fizH1, 784
.set nb313_fixH2, 800
.set nb313_fiyH2, 816
.set nb313_fizH2, 832
.set nb313_fixM, 848
.set nb313_fiyM, 864
.set nb313_fizM, 880
.set nb313_fjx, 896
.set nb313_fjy, 912
.set nb313_fjz, 928
.set nb313_half, 944
.set nb313_three, 960
.set nb313_is3, 976
.set nb313_ii3, 980
.set nb313_nri, 984
.set nb313_iinr, 992
.set nb313_jindex, 1000
.set nb313_jjnr, 1008
.set nb313_shift, 1016
.set nb313_shiftvec, 1024
.set nb313_facel, 1032
.set nb313_innerjjnr, 1040
.set nb313_ntia, 1048
.set nb313_innerk, 1052
.set nb313_n, 1056
.set nb313_nn1, 1060
.set nb313_nouter, 1064
.set nb313_ninner, 1068
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1080,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb313_nouter(%rsp)
        movl %eax,nb313_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb313_nri(%rsp)
        movq %rsi,nb313_iinr(%rsp)
        movq %rdx,nb313_jindex(%rsp)
        movq %rcx,nb313_jjnr(%rsp)
        movq %r8,nb313_shift(%rsp)
        movq %r9,nb313_shiftvec(%rsp)
        movq nb313_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb313_facel(%rsp)

        movq nb313_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb313_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb313_half(%rsp)
        movl %ebx,nb313_half+4(%rsp)
        movsd nb313_half(%rsp),%xmm1
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
        movapd %xmm1,nb313_half(%rsp)
        movapd %xmm2,nb313_two(%rsp)
        movapd %xmm3,nb313_three(%rsp)
        movapd %xmm4,nb313_six(%rsp)
        movapd %xmm5,nb313_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb313_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb313_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb313_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb313_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb313_iqH(%rsp)
        movapd %xmm4,nb313_iqM(%rsp)

        movq  nb313_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb313_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313_ntia(%rsp)
_nb_kernel313_x86_64_sse2.nb313_threadloop: 
        movq  nb313_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel313_x86_64_sse2.nb313_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313_x86_64_sse2.nb313_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313_n(%rsp)
        movl %ebx,nb313_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313_x86_64_sse2.nb313_outerstart
        jmp _nb_kernel313_x86_64_sse2.nb313_end

_nb_kernel313_x86_64_sse2.nb313_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313_nouter(%rsp),%ebx
        movl %ebx,nb313_nouter(%rsp)

_nb_kernel313_x86_64_sse2.nb313_outer: 
        movq  nb313_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb313_is3(%rsp)      ## store is3 

        movq  nb313_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb313_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb313_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb313_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb313_ixO(%rsp)
        movapd %xmm4,nb313_iyO(%rsp)
        movapd %xmm5,nb313_izO(%rsp)
        movapd %xmm6,nb313_ixH1(%rsp)
        movapd %xmm7,nb313_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb313_izH1(%rsp)
        movapd %xmm0,nb313_ixH2(%rsp)
        movapd %xmm1,nb313_iyH2(%rsp)
        movapd %xmm2,nb313_izH2(%rsp)
        movapd %xmm3,nb313_ixM(%rsp)
        movapd %xmm4,nb313_iyM(%rsp)
        movapd %xmm5,nb313_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb313_vctot(%rsp)
        movapd %xmm4,nb313_Vvdwtot(%rsp)
        movapd %xmm4,nb313_fixO(%rsp)
        movapd %xmm4,nb313_fiyO(%rsp)
        movapd %xmm4,nb313_fizO(%rsp)
        movapd %xmm4,nb313_fixH1(%rsp)
        movapd %xmm4,nb313_fiyH1(%rsp)
        movapd %xmm4,nb313_fizH1(%rsp)
        movapd %xmm4,nb313_fixH2(%rsp)
        movapd %xmm4,nb313_fiyH2(%rsp)
        movapd %xmm4,nb313_fizH2(%rsp)
        movapd %xmm4,nb313_fixM(%rsp)
        movapd %xmm4,nb313_fiyM(%rsp)
        movapd %xmm4,nb313_fizM(%rsp)

        movq  nb313_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb313_pos(%rbp),%rsi
        movq  nb313_faction(%rbp),%rdi
        movq  nb313_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb313_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb313_ninner(%rsp),%ecx
        movl  %ecx,nb313_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb313_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel313_x86_64_sse2.nb313_unroll_loop
        jmp   _nb_kernel313_x86_64_sse2.nb313_checksingle
_nb_kernel313_x86_64_sse2.nb313_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb313_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb313_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb313_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313_iqM(%rsp),%xmm3
        mulpd  nb313_iqH(%rsp),%xmm4
        movapd  %xmm3,nb313_qqM(%rsp)
        movapd  %xmm4,nb313_qqH(%rsp)

        movq nb313_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb313_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb313_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb313_c6(%rsp)
        movapd %xmm6,nb313_c12(%rsp)

        movq nb313_pos(%rbp),%rsi        ## base of pos[] 

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

    ## O interaction
    ## copy to xmm3-xmm5
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5

    subpd nb313_ixO(%rsp),%xmm3
    subpd nb313_iyO(%rsp),%xmm4
    subpd nb313_izO(%rsp),%xmm5

    movapd %xmm3,%xmm13
    movapd %xmm4,%xmm14
    movapd %xmm5,%xmm15

        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5

        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3

    ## calc 1/rsq
    cvtpd2ps %xmm3,%xmm6
    rcpps %xmm6,%xmm6
    cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

    ## 1/x lookup seed in xmm6 
    movapd nb313_two(%rsp),%xmm4
    movapd %xmm3,%xmm5
    mulpd %xmm6,%xmm3       ## lu*rsq 
    subpd %xmm3,%xmm4       ## 2-lu*rsq 
    mulpd %xmm4,%xmm6       ## (new lu) 

    movapd nb313_two(%rsp),%xmm4
    mulpd %xmm6,%xmm5       ## lu*rsq 
    subpd %xmm5,%xmm4       ## 2-lu*rsq 
    mulpd %xmm6,%xmm4       ## xmm4=rinvsq 

    movapd %xmm4,%xmm3      ## rinvsq
    mulpd  %xmm4,%xmm4      ## rinv4
    mulpd  %xmm3,%xmm4      ## rinv6
    movapd %xmm4,%xmm5
    mulpd  %xmm5,%xmm5      ## rinv12
    mulpd  nb313_c6(%rsp),%xmm4
    mulpd  nb313_c12(%rsp),%xmm5
    movapd %xmm5,%xmm6
    subpd  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulpd  nb313_six(%rsp),%xmm4
    mulpd  nb313_twelve(%rsp),%xmm5
    subpd  %xmm4,%xmm5
    mulpd  %xmm5,%xmm3  ## fscal

    addpd  nb313_Vvdwtot(%rsp),%xmm6
    movapd %xmm6,nb313_Vvdwtot(%rsp)

    mulpd  %xmm3,%xmm13 ## fx
    mulpd  %xmm3,%xmm14 ## fy
    mulpd  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movapd %xmm13,nb313_fjx(%rsp)
    movapd %xmm14,nb313_fjy(%rsp)
    movapd %xmm15,nb313_fjz(%rsp)

    ## increment i O force
    addpd nb313_fixO(%rsp),%xmm13
    addpd nb313_fiyO(%rsp),%xmm14
    addpd nb313_fizO(%rsp),%xmm15
    movapd %xmm13,nb313_fixO(%rsp)
    movapd %xmm14,nb313_fiyO(%rsp)
    movapd %xmm15,nb313_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb313_ixH1(%rsp),%xmm0
    subpd nb313_iyH1(%rsp),%xmm1
    subpd nb313_izH1(%rsp),%xmm2
    subpd nb313_ixH2(%rsp),%xmm3
    subpd nb313_iyH2(%rsp),%xmm4
    subpd nb313_izH2(%rsp),%xmm5
    subpd nb313_ixM(%rsp),%xmm6
    subpd nb313_iyM(%rsp),%xmm7
    subpd nb313_izM(%rsp),%xmm8

        movapd %xmm0,nb313_dxH1(%rsp)
        movapd %xmm1,nb313_dyH1(%rsp)
        movapd %xmm2,nb313_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb313_dxH2(%rsp)
        movapd %xmm4,nb313_dyH2(%rsp)
        movapd %xmm5,nb313_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb313_dxM(%rsp)
        movapd %xmm7,nb313_dyM(%rsp)
        movapd %xmm8,nb313_dzM(%rsp)
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

        movapd  nb313_three(%rsp),%xmm9
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

        movapd  nb313_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb313_three(%rsp),%xmm1
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

        movapd  nb313_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        movapd  %xmm9,nb313_rinvH1(%rsp)
        movapd  %xmm10,nb313_rinvH2(%rsp)
        movapd  %xmm11,nb313_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb313_tsc(%rsp),%xmm1
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

    movq nb313_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd %xmm0,%xmm12 ## epsH1
    movapd %xmm3,%xmm13 ## epsH2
    movapd %xmm6,%xmm14 ## epsM

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
    mulpd  %xmm12,%xmm2  ## Geps
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
    movapd nb313_qqH(%rsp),%xmm12
    movapd nb313_qqM(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb313_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb313_vctot(%rsp)

    movapd nb313_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb313_rinvH1(%rsp),%xmm4
    mulpd nb313_rinvH2(%rsp),%xmm8
    mulpd nb313_rinvM(%rsp),%xmm11

    ## move j forces to xmm0-xmm2
    movq nb313_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

    ## add forces from O interaction
    addpd nb313_fjx(%rsp),%xmm0
    addpd nb313_fjy(%rsp),%xmm1
    addpd nb313_fjz(%rsp),%xmm2

        mulpd nb313_dxH1(%rsp),%xmm3
        mulpd nb313_dyH1(%rsp),%xmm4
        mulpd nb313_dzH1(%rsp),%xmm5
        mulpd nb313_dxH2(%rsp),%xmm7
        mulpd nb313_dyH2(%rsp),%xmm8
        mulpd nb313_dzH2(%rsp),%xmm9
        mulpd nb313_dxM(%rsp),%xmm10
        mulpd nb313_dyM(%rsp),%xmm11
        mulpd nb313_dzM(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb313_fixH1(%rsp),%xmm3
    addpd nb313_fiyH1(%rsp),%xmm4
    addpd nb313_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb313_fixH2(%rsp),%xmm7
    addpd nb313_fiyH2(%rsp),%xmm8
    addpd nb313_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb313_fixM(%rsp),%xmm10
    addpd nb313_fiyM(%rsp),%xmm11
    addpd nb313_fizM(%rsp),%xmm12

    movapd %xmm3,nb313_fixH1(%rsp)
    movapd %xmm4,nb313_fiyH1(%rsp)
    movapd %xmm5,nb313_fizH1(%rsp)
    movapd %xmm7,nb313_fixH2(%rsp)
    movapd %xmm8,nb313_fiyH2(%rsp)
    movapd %xmm9,nb313_fizH2(%rsp)
    movapd %xmm10,nb313_fixM(%rsp)
    movapd %xmm11,nb313_fiyM(%rsp)
    movapd %xmm12,nb313_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb313_innerk(%rsp)
        jl    _nb_kernel313_x86_64_sse2.nb313_checksingle
        jmp   _nb_kernel313_x86_64_sse2.nb313_unroll_loop
_nb_kernel313_x86_64_sse2.nb313_checksingle: 
        movl  nb313_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel313_x86_64_sse2.nb313_dosingle
        jmp   _nb_kernel313_x86_64_sse2.nb313_updateouterdata
_nb_kernel313_x86_64_sse2.nb313_dosingle: 
        movq  nb313_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb313_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313_iqM(%rsp),%xmm3
        mulpd  nb313_iqH(%rsp),%xmm4

        movapd  %xmm3,nb313_qqM(%rsp)
        movapd  %xmm4,nb313_qqH(%rsp)

        movq nb313_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb313_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb313_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb313_c6(%rsp)
        movapd %xmm6,nb313_c12(%rsp)

        movq nb313_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2  and xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb313_ixO(%rsp),%xmm4
        subsd nb313_iyO(%rsp),%xmm5
        subsd nb313_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxO(%rsp)
        movapd %xmm5,nb313_dyO(%rsp)
        movapd %xmm6,nb313_dzO(%rsp)
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
        subsd nb313_ixH1(%rsp),%xmm4
        subsd nb313_iyH1(%rsp),%xmm5
        subsd nb313_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb313_dxH1(%rsp)
        movapd %xmm5,nb313_dyH1(%rsp)
        movapd %xmm6,nb313_dzH1(%rsp)
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
        subsd nb313_ixH2(%rsp),%xmm3
        subsd nb313_iyH2(%rsp),%xmm4
        subsd nb313_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb313_dxH2(%rsp)
        movapd %xmm4,nb313_dyH2(%rsp)
        movapd %xmm5,nb313_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5

        ## move j coords to xmm4-xmm2
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm3
    ## xmm2 already contains z

        ## calc dr 
        subsd nb313_ixM(%rsp),%xmm4
        subsd nb313_iyM(%rsp),%xmm3
        subsd nb313_izM(%rsp),%xmm2

        ## store dr 
        movapd %xmm4,nb313_dxM(%rsp)
        movapd %xmm3,nb313_dyM(%rsp)
        movapd %xmm2,nb313_dzM(%rsp)

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtsd2ss %xmm7,%xmm2
        movsd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movsd   nb313_two(%rsp),%xmm1
        movsd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm3
        subsd   %xmm3,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movsd  %xmm0,nb313_rinvsqO(%rsp)

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%rsp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH1(%rsp)         ## rinvH1 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb313_rH1(%rsp)    ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%rsp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvH2(%rsp)   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb313_rH2(%rsp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313_three(%rsp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313_rinvM(%rsp)   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb313_rM(%rsp)   ## r 

        ## do O interactions 
        movapd  nb313_rinvsqO(%rsp),%xmm0
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb313_c6(%rsp),%xmm1
        mulsd  nb313_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb313_Vvdwtot(%rsp),%xmm3
        mulsd  nb313_six(%rsp),%xmm1
        mulsd  nb313_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb313_Vvdwtot(%rsp)

        movapd nb313_dxO(%rsp),%xmm0
        movapd nb313_dyO(%rsp),%xmm1
        movapd nb313_dzO(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## tx in xmm0-xmm2 

        ## update O forces 
        movapd nb313_fixO(%rsp),%xmm3
        movapd nb313_fiyO(%rsp),%xmm4
        movapd nb313_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixO(%rsp)
        movlpd %xmm4,nb313_fiyO(%rsp)
        movlpd %xmm7,nb313_fizO(%rsp)
        ## update j forces with water O 
        movlpd %xmm0,nb313_fjx(%rsp)
        movlpd %xmm1,nb313_fjy(%rsp)
        movlpd %xmm2,nb313_fjz(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb313_rH1(%rsp),%xmm7
        mulpd nb313_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb313_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul 
        xorpd  %xmm4,%xmm4
    addsd  nb313_vctot(%rsp),%xmm5
        mulsd  nb313_rinvH1(%rsp),%xmm3
    movlpd %xmm5,nb313_vctot(%rsp)
        mulsd  nb313_tsc(%rsp),%xmm3
        subsd %xmm3,%xmm4

        movapd nb313_dxH1(%rsp),%xmm0
        movapd nb313_dyH1(%rsp),%xmm1
        movapd nb313_dzH1(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb313_fixH1(%rsp),%xmm3
        movapd nb313_fiyH1(%rsp),%xmm4
        movapd nb313_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixH1(%rsp)
        movlpd %xmm4,nb313_fiyH1(%rsp)
        movlpd %xmm7,nb313_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb313_fjx(%rsp),%xmm0
        addsd  nb313_fjy(%rsp),%xmm1
        addsd  nb313_fjz(%rsp),%xmm2
        movlpd %xmm0,nb313_fjx(%rsp)
        movlpd %xmm1,nb313_fjy(%rsp)
        movlpd %xmm2,nb313_fjz(%rsp)

        ##  H2 interactions 
        movapd nb313_rH2(%rsp),%xmm7
        mulsd   nb313_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb313_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb313_qqH(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb313_vctot(%rsp),%xmm5
        mulsd  nb313_rinvH2(%rsp),%xmm3
        movlpd %xmm5,nb313_vctot(%rsp)
        mulsd  nb313_tsc(%rsp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb313_dxH2(%rsp),%xmm0
        movapd nb313_dyH2(%rsp),%xmm1
        movapd nb313_dzH2(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb313_fixH2(%rsp),%xmm3
        movapd nb313_fiyH2(%rsp),%xmm4
        movapd nb313_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixH2(%rsp)
        movlpd %xmm4,nb313_fiyH2(%rsp)
        movlpd %xmm7,nb313_fizH2(%rsp)
        ## update j forces with water H1 
        addsd  nb313_fjx(%rsp),%xmm0
        addsd  nb313_fjy(%rsp),%xmm1
        addsd  nb313_fjz(%rsp),%xmm2
        movlpd %xmm0,nb313_fjx(%rsp)
        movlpd %xmm1,nb313_fjy(%rsp)
        movlpd %xmm2,nb313_fjz(%rsp)

        ## M interactions 
        movapd nb313_rM(%rsp),%xmm7
        mulsd   nb313_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 
        movq nb313_VFtab(%rbp),%rsi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb313_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb313_qqM(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 
        ## increment vcoul 
        xorpd  %xmm4,%xmm4
        addsd  nb313_vctot(%rsp),%xmm5
        mulsd  nb313_rinvM(%rsp),%xmm3
        movlpd %xmm5,nb313_vctot(%rsp)
        mulsd  nb313_tsc(%rsp),%xmm3
        subsd  %xmm3,%xmm4

        movapd nb313_dxM(%rsp),%xmm0
        movapd nb313_dyM(%rsp),%xmm1
        movapd nb313_dzM(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb313_fixM(%rsp),%xmm3
        movapd nb313_fiyM(%rsp),%xmm4
        movapd nb313_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb313_fixM(%rsp)
        movlpd %xmm4,nb313_fiyM(%rsp)
        movlpd %xmm7,nb313_fizM(%rsp)

        movq nb313_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb313_fjx(%rsp),%xmm0
        addsd  nb313_fjy(%rsp),%xmm1
        addsd  nb313_fjz(%rsp),%xmm2

        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel313_x86_64_sse2.nb313_updateouterdata: 
        movl  nb313_ii3(%rsp),%ecx
        movq  nb313_faction(%rbp),%rdi
        movq  nb313_fshift(%rbp),%rsi
        movl  nb313_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb313_fixO(%rsp),%xmm0
        movapd nb313_fiyO(%rsp),%xmm1
        movapd nb313_fizO(%rsp),%xmm2

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
        movapd nb313_fixH1(%rsp),%xmm0
        movapd nb313_fiyH1(%rsp),%xmm1
        movapd nb313_fizH1(%rsp),%xmm2

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
        movapd nb313_fixH2(%rsp),%xmm0
        movapd nb313_fiyH2(%rsp),%xmm1
        movapd nb313_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb313_fixM(%rsp),%xmm0
        movapd nb313_fiyM(%rsp),%xmm1
        movapd nb313_fizM(%rsp),%xmm2

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
        movsd  72(%rdi,%rcx,8),%xmm3
        movsd  80(%rdi,%rcx,8),%xmm4
        movsd  88(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,72(%rdi,%rcx,8)
        movsd  %xmm4,80(%rdi,%rcx,8)
        movsd  %xmm5,88(%rdi,%rcx,8)

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
        movl nb313_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb313_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb313_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb313_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb313_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb313_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb313_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313_x86_64_sse2.nb313_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313_n(%rsp)
        jmp _nb_kernel313_x86_64_sse2.nb313_outer
_nb_kernel313_x86_64_sse2.nb313_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313_x86_64_sse2.nb313_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313_x86_64_sse2.nb313_threadloop
_nb_kernel313_x86_64_sse2.nb313_end: 
        movl nb313_nouter(%rsp),%eax
        movl nb313_ninner(%rsp),%ebx
        movq nb313_outeriter(%rbp),%rcx
        movq nb313_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1080,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel313nf_x86_64_sse2
.globl _nb_kernel313nf_x86_64_sse2
nb_kernel313nf_x86_64_sse2:     
_nb_kernel313nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb313nf_fshift, 16
.set nb313nf_gid, 24
.set nb313nf_pos, 32
.set nb313nf_faction, 40
.set nb313nf_charge, 48
.set nb313nf_p_facel, 56
.set nb313nf_argkrf, 64
.set nb313nf_argcrf, 72
.set nb313nf_Vc, 80
.set nb313nf_type, 88
.set nb313nf_p_ntype, 96
.set nb313nf_vdwparam, 104
.set nb313nf_Vvdw, 112
.set nb313nf_p_tabscale, 120
.set nb313nf_VFtab, 128
.set nb313nf_invsqrta, 136
.set nb313nf_dvda, 144
.set nb313nf_p_gbtabscale, 152
.set nb313nf_GBtab, 160
.set nb313nf_p_nthreads, 168
.set nb313nf_count, 176
.set nb313nf_mtx, 184
.set nb313nf_outeriter, 192
.set nb313nf_inneriter, 200
.set nb313nf_work, 208
        ## stack offsets for local variables 
        ## bottom of stack is cache-aligned for sse2 use 
.set nb313nf_ixO, 0
.set nb313nf_iyO, 16
.set nb313nf_izO, 32
.set nb313nf_ixH1, 48
.set nb313nf_iyH1, 64
.set nb313nf_izH1, 80
.set nb313nf_ixH2, 96
.set nb313nf_iyH2, 112
.set nb313nf_izH2, 128
.set nb313nf_ixM, 144
.set nb313nf_iyM, 160
.set nb313nf_izM, 176
.set nb313nf_iqM, 192
.set nb313nf_iqH, 208
.set nb313nf_qqM, 224
.set nb313nf_qqH, 240
.set nb313nf_rinvsqO, 256
.set nb313nf_rinvH1, 272
.set nb313nf_rinvH2, 288
.set nb313nf_rinvM, 304
.set nb313nf_rO, 320
.set nb313nf_rH1, 336
.set nb313nf_rH2, 352
.set nb313nf_rM, 368
.set nb313nf_tsc, 384
.set nb313nf_two, 400
.set nb313nf_c6, 416
.set nb313nf_c12, 432
.set nb313nf_vctot, 448
.set nb313nf_Vvdwtot, 464
.set nb313nf_half, 480
.set nb313nf_three, 496
.set nb313nf_is3, 512
.set nb313nf_ii3, 516
.set nb313nf_nri, 520
.set nb313nf_iinr, 528
.set nb313nf_jindex, 536
.set nb313nf_jjnr, 544
.set nb313nf_shift, 552
.set nb313nf_shiftvec, 560
.set nb313nf_facel, 568
.set nb313nf_innerjjnr, 576
.set nb313nf_ntia, 584
.set nb313nf_innerk, 588
.set nb313nf_n, 592
.set nb313nf_nn1, 596
.set nb313nf_nouter, 600
.set nb313nf_ninner, 604
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $616,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb313nf_nouter(%rsp)
        movl %eax,nb313nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb313nf_nri(%rsp)
        movq %rsi,nb313nf_iinr(%rsp)
        movq %rdx,nb313nf_jindex(%rsp)
        movq %rcx,nb313nf_jjnr(%rsp)
        movq %r8,nb313nf_shift(%rsp)
        movq %r9,nb313nf_shiftvec(%rsp)
        movq nb313nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb313nf_facel(%rsp)

        movq nb313nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb313nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb313nf_half(%rsp)
        movl %ebx,nb313nf_half+4(%rsp)
        movsd nb313nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb313nf_half(%rsp)
        movapd %xmm2,nb313nf_two(%rsp)
        movapd %xmm3,nb313nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb313nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb313nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb313nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb313nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb313nf_iqH(%rsp)
        movapd %xmm4,nb313nf_iqM(%rsp)

        movq  nb313nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb313nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313nf_ntia(%rsp)
_nb_kernel313nf_x86_64_sse2.nb313nf_threadloop: 
        movq  nb313nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel313nf_x86_64_sse2.nb313nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313nf_x86_64_sse2.nb313nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313nf_n(%rsp)
        movl %ebx,nb313nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313nf_x86_64_sse2.nb313nf_outerstart
        jmp _nb_kernel313nf_x86_64_sse2.nb313nf_end

_nb_kernel313nf_x86_64_sse2.nb313nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313nf_nouter(%rsp),%ebx
        movl %ebx,nb313nf_nouter(%rsp)

_nb_kernel313nf_x86_64_sse2.nb313nf_outer: 
        movq  nb313nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb313nf_is3(%rsp)            ## store is3 

        movq  nb313nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb313nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb313nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb313nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb313nf_ixO(%rsp)
        movapd %xmm4,nb313nf_iyO(%rsp)
        movapd %xmm5,nb313nf_izO(%rsp)
        movapd %xmm6,nb313nf_ixH1(%rsp)
        movapd %xmm7,nb313nf_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb313nf_izH1(%rsp)
        movapd %xmm0,nb313nf_ixH2(%rsp)
        movapd %xmm1,nb313nf_iyH2(%rsp)
        movapd %xmm2,nb313nf_izH2(%rsp)
        movapd %xmm3,nb313nf_ixM(%rsp)
        movapd %xmm4,nb313nf_iyM(%rsp)
        movapd %xmm5,nb313nf_izM(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb313nf_vctot(%rsp)
        movapd %xmm4,nb313nf_Vvdwtot(%rsp)

        movq  nb313nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb313nf_pos(%rbp),%rsi
        movq  nb313nf_faction(%rbp),%rdi
        movq  nb313nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb313nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb313nf_ninner(%rsp),%ecx
        movl  %ecx,nb313nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb313nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel313nf_x86_64_sse2.nb313nf_unroll_loop
        jmp   _nb_kernel313nf_x86_64_sse2.nb313nf_checksingle
_nb_kernel313nf_x86_64_sse2.nb313nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb313nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb313nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb313nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313nf_iqM(%rsp),%xmm3
        mulpd  nb313nf_iqH(%rsp),%xmm4
        movapd  %xmm3,nb313nf_qqM(%rsp)
        movapd  %xmm4,nb313nf_qqH(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movq nb313nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb313nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb313nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb313nf_c6(%rsp)
        movapd %xmm6,nb313nf_c12(%rsp)

        movq nb313nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb313nf_ixO(%rsp),%xmm4
        movapd nb313nf_iyO(%rsp),%xmm5
        movapd nb313nf_izO(%rsp),%xmm6

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
        movapd nb313nf_ixH1(%rsp),%xmm4
        movapd nb313nf_iyH1(%rsp),%xmm5
        movapd nb313nf_izH1(%rsp),%xmm6

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
        movapd nb313nf_ixH2(%rsp),%xmm3
        movapd nb313nf_iyH2(%rsp),%xmm4
        movapd nb313nf_izH2(%rsp),%xmm5

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
        ## move ixM-izM to xmm2-xmm4  
        movapd nb313nf_iyM(%rsp),%xmm3
        movapd nb313nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb313nf_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm3
        subpd   %xmm3,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,nb313nf_rinvsqO(%rsp)

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH1(%rsp)       ## rinvH1 
        mulpd  %xmm0,%xmm6
        movapd %xmm6,nb313nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH2(%rsp)   ## rinv 
        mulpd %xmm0,%xmm5
        movapd %xmm5,nb313nf_rH2(%rsp)   ## r 

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulpd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulpd %xmm0,%xmm0       ## lu*lu 
        mulpd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subpd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulpd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvM(%rsp)   ## rinv 
        mulpd %xmm0,%xmm4
        movapd %xmm4,nb313nf_rM(%rsp)   ## r 

        ## do O interactions 
        movapd nb313nf_rinvsqO(%rsp),%xmm0
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb313nf_c6(%rsp),%xmm1
        mulpd  nb313nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb313nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb313nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb313nf_rH1(%rsp),%xmm7
        mulpd nb313nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        movapd nb313nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb313nf_vctot(%rsp)

        ## H2 interactions 
        movapd nb313nf_rH2(%rsp),%xmm7
        mulpd   nb313nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        movapd nb313nf_qqH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb313nf_vctot(%rsp)

        ## M interactions 
        movapd nb313nf_rM(%rsp),%xmm7
        mulpd   nb313nf_tsc(%rsp),%xmm7
        cvttpd2pi %xmm7,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        movd %eax,%mm0
        movd %ebx,%mm1

        pslld $2,%mm6           ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

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
        movapd nb313nf_qqM(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addpd  nb313nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb313nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb313nf_innerk(%rsp)
        jl    _nb_kernel313nf_x86_64_sse2.nb313nf_checksingle
        jmp   _nb_kernel313nf_x86_64_sse2.nb313nf_unroll_loop
_nb_kernel313nf_x86_64_sse2.nb313nf_checksingle: 
        movl  nb313nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel313nf_x86_64_sse2.nb313nf_dosingle
        jmp   _nb_kernel313nf_x86_64_sse2.nb313nf_updateouterdata
_nb_kernel313nf_x86_64_sse2.nb313nf_dosingle: 
        movq  nb313nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb313nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb313nf_iqM(%rsp),%xmm3
        mulpd  nb313nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movapd  %xmm3,nb313nf_qqM(%rsp)
        movapd  %xmm4,nb313nf_qqH(%rsp)

        movq nb313nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb313nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb313nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb313nf_c6(%rsp)
        movapd %xmm6,nb313nf_c12(%rsp)

        movq nb313nf_pos(%rbp),%rsi        ## base of pos[] 
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coords to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb313nf_ixO(%rsp),%xmm4
        movapd nb313nf_iyO(%rsp),%xmm5
        movapd nb313nf_izO(%rsp),%xmm6

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
        movapd nb313nf_ixH1(%rsp),%xmm4
        movapd nb313nf_iyH1(%rsp),%xmm5
        movapd nb313nf_izH1(%rsp),%xmm6

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
        movapd nb313nf_ixH2(%rsp),%xmm3
        movapd nb313nf_iyH2(%rsp),%xmm4
        movapd nb313nf_izH2(%rsp),%xmm5

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
        ## move ixM-izM to xmm2-xmm4  
        movapd nb313nf_iyM(%rsp),%xmm3
        movapd nb313nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb313nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## 1/x for O - rsqO is in xmm7
        cvtsd2ss %xmm7,%xmm2
        movsd   %xmm7,%xmm3
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movsd   nb313nf_two(%rsp),%xmm1
        movsd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm3
        subsd   %xmm3,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movsd  %xmm0,nb313nf_rinvsqO(%rsp)

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm6,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH1(%rsp)       ## rinvH1 
        mulsd  %xmm0,%xmm6
        movapd %xmm6,nb313nf_rH1(%rsp)          ## rH1 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm5,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvH2(%rsp)   ## rinv 
        mulsd %xmm0,%xmm5
        movapd %xmm5,nb313nf_rH2(%rsp)   ## r 

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb313nf_three(%rsp),%xmm0
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulsd   nb313nf_half(%rsp),%xmm0   ## iter1 ( new lu) 

        movapd %xmm4,%xmm2
        movapd %xmm0,%xmm3
        mulsd %xmm0,%xmm0       ## lu*lu 
        mulsd %xmm0,%xmm2       ## rsq*lu*lu 
        movapd nb313nf_three(%rsp),%xmm0
        subsd %xmm2,%xmm0       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm0       ## lu*( 3-rsq*lu*lu) 
        mulsd nb313nf_half(%rsp),%xmm0   ## rinv 
        movapd %xmm0,nb313nf_rinvM(%rsp)   ## rinv 
        mulsd %xmm0,%xmm4
        movapd %xmm4,nb313nf_rM(%rsp)   ## r 

        ## do O interactions 
        movapd  nb313nf_rinvsqO(%rsp),%xmm0
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb313nf_c6(%rsp),%xmm1
        mulsd  nb313nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb313nf_Vvdwtot(%rsp),%xmm3
        movsd %xmm3,nb313nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movapd nb313nf_rH1(%rsp),%xmm7
        mulpd nb313nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi

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
        movapd nb313nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%rsp)

        ##  H2 interactions 
        movapd nb313nf_rH2(%rsp),%xmm7
        mulsd   nb313nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi

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
        movapd nb313nf_qqH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%rsp)

        ## M interactions 
        movapd nb313nf_rM(%rsp),%xmm7
        mulsd   nb313nf_tsc(%rsp),%xmm7
        cvttsd2si %xmm7,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm7
        movapd %xmm7,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb313nf_VFtab(%rbp),%rsi

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
        movapd nb313nf_qqM(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## increment vcoul 
        addsd  nb313nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb313nf_vctot(%rsp)

_nb_kernel313nf_x86_64_sse2.nb313nf_updateouterdata: 
        ## get n from stack
        movl nb313nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb313nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb313nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb313nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb313nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb313nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb313nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313nf_x86_64_sse2.nb313nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313nf_n(%rsp)
        jmp _nb_kernel313nf_x86_64_sse2.nb313nf_outer
_nb_kernel313nf_x86_64_sse2.nb313nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313nf_x86_64_sse2.nb313nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313nf_x86_64_sse2.nb313nf_threadloop
_nb_kernel313nf_x86_64_sse2.nb313nf_end: 
        movl nb313nf_nouter(%rsp),%eax
        movl nb313nf_ninner(%rsp),%ebx
        movq nb313nf_outeriter(%rbp),%rcx
        movq nb313nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $616,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



