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





.globl nb_kernel233_x86_64_sse2
.globl _nb_kernel233_x86_64_sse2
nb_kernel233_x86_64_sse2:       
_nb_kernel233_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb233_fshift, 16
.set nb233_gid, 24
.set nb233_pos, 32
.set nb233_faction, 40
.set nb233_charge, 48
.set nb233_p_facel, 56
.set nb233_argkrf, 64
.set nb233_argcrf, 72
.set nb233_Vc, 80
.set nb233_type, 88
.set nb233_p_ntype, 96
.set nb233_vdwparam, 104
.set nb233_Vvdw, 112
.set nb233_p_tabscale, 120
.set nb233_VFtab, 128
.set nb233_invsqrta, 136
.set nb233_dvda, 144
.set nb233_p_gbtabscale, 152
.set nb233_GBtab, 160
.set nb233_p_nthreads, 168
.set nb233_count, 176
.set nb233_mtx, 184
.set nb233_outeriter, 192
.set nb233_inneriter, 200
.set nb233_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb233_ixO, 0
.set nb233_iyO, 16
.set nb233_izO, 32
.set nb233_ixH1, 48
.set nb233_iyH1, 64
.set nb233_izH1, 80
.set nb233_ixH2, 96
.set nb233_iyH2, 112
.set nb233_izH2, 128
.set nb233_ixM, 144
.set nb233_iyM, 160
.set nb233_izM, 176
.set nb233_iqH, 192
.set nb233_iqM, 208
.set nb233_dxO, 224
.set nb233_dyO, 240
.set nb233_dzO, 256
.set nb233_dxH1, 272
.set nb233_dyH1, 288
.set nb233_dzH1, 304
.set nb233_dxH2, 320
.set nb233_dyH2, 336
.set nb233_dzH2, 352
.set nb233_dxM, 368
.set nb233_dyM, 384
.set nb233_dzM, 400
.set nb233_qqH, 416
.set nb233_qqM, 432
.set nb233_c6, 448
.set nb233_c12, 464
.set nb233_tsc, 480
.set nb233_fstmp, 496
.set nb233_vctot, 512
.set nb233_Vvdwtot, 528
.set nb233_fixO, 544
.set nb233_fiyO, 560
.set nb233_fizO, 576
.set nb233_fixH1, 592
.set nb233_fiyH1, 608
.set nb233_fizH1, 624
.set nb233_fixH2, 640
.set nb233_fiyH2, 656
.set nb233_fizH2, 672
.set nb233_fixM, 688
.set nb233_fiyM, 704
.set nb233_fizM, 720
.set nb233_fjx, 736
.set nb233_fjy, 752
.set nb233_fjz, 768
.set nb233_half, 784
.set nb233_three, 800
.set nb233_two, 816
.set nb233_rinvH1, 832
.set nb233_rinvH2, 848
.set nb233_rinvM, 864
.set nb233_krsqH1, 880
.set nb233_krsqH2, 896
.set nb233_krsqM, 912
.set nb233_krf, 928
.set nb233_crf, 944
.set nb233_rsqO, 960
.set nb233_facel, 976
.set nb233_innerjjnr, 984
.set nb233_iinr, 992
.set nb233_jindex, 1000
.set nb233_jjnr, 1008
.set nb233_shift, 1016
.set nb233_shiftvec, 1024
.set nb233_is3, 1032
.set nb233_ii3, 1036
.set nb233_nri, 1040
.set nb233_ntia, 1044
.set nb233_innerk, 1048
.set nb233_n, 1052
.set nb233_nn1, 1056
.set nb233_nouter, 1060
.set nb233_ninner, 1064
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
        movl %eax,nb233_nouter(%rsp)
        movl %eax,nb233_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb233_nri(%rsp)
        movq %rsi,nb233_iinr(%rsp)
        movq %rdx,nb233_jindex(%rsp)
        movq %rcx,nb233_jjnr(%rsp)
        movq %r8,nb233_shift(%rsp)
        movq %r9,nb233_shiftvec(%rsp)
        movq nb233_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb233_facel(%rsp)

        movq nb233_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb233_tsc(%rsp)

        movq nb233_argkrf(%rbp),%rsi
        movq nb233_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb233_krf(%rsp)
        movapd %xmm2,nb233_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb233_half(%rsp)
        movl %ebx,nb233_half+4(%rsp)
        movsd nb233_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb233_half(%rsp)
        movapd %xmm2,nb233_two(%rsp)
        movapd %xmm3,nb233_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb233_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb233_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb233_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb233_iqH(%rsp)
        movapd %xmm4,nb233_iqM(%rsp)

        movq  nb233_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb233_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233_ntia(%rsp)
_nb_kernel233_x86_64_sse2.nb233_threadloop: 
        movq  nb233_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel233_x86_64_sse2.nb233_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233_x86_64_sse2.nb233_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233_n(%rsp)
        movl %ebx,nb233_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233_x86_64_sse2.nb233_outerstart
        jmp _nb_kernel233_x86_64_sse2.nb233_end

_nb_kernel233_x86_64_sse2.nb233_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233_nouter(%rsp),%ebx
        movl %ebx,nb233_nouter(%rsp)

_nb_kernel233_x86_64_sse2.nb233_outer: 
        movq  nb233_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb233_is3(%rsp)      ## store is3 

        movq  nb233_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb233_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb233_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb233_ii3(%rsp)

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
        movapd %xmm3,nb233_ixO(%rsp)
        movapd %xmm4,nb233_iyO(%rsp)
        movapd %xmm5,nb233_izO(%rsp)
        movapd %xmm6,nb233_ixH1(%rsp)
        movapd %xmm7,nb233_iyH1(%rsp)

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
        movapd %xmm6,nb233_izH1(%rsp)
        movapd %xmm0,nb233_ixH2(%rsp)
        movapd %xmm1,nb233_iyH2(%rsp)
        movapd %xmm2,nb233_izH2(%rsp)
        movapd %xmm3,nb233_ixM(%rsp)
        movapd %xmm4,nb233_iyM(%rsp)
        movapd %xmm5,nb233_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb233_vctot(%rsp)
        movapd %xmm4,nb233_Vvdwtot(%rsp)
        movapd %xmm4,nb233_fixO(%rsp)
        movapd %xmm4,nb233_fiyO(%rsp)
        movapd %xmm4,nb233_fizO(%rsp)
        movapd %xmm4,nb233_fixH1(%rsp)
        movapd %xmm4,nb233_fiyH1(%rsp)
        movapd %xmm4,nb233_fizH1(%rsp)
        movapd %xmm4,nb233_fixH2(%rsp)
        movapd %xmm4,nb233_fiyH2(%rsp)
        movapd %xmm4,nb233_fizH2(%rsp)
        movapd %xmm4,nb233_fixM(%rsp)
        movapd %xmm4,nb233_fiyM(%rsp)
        movapd %xmm4,nb233_fizM(%rsp)

        movq  nb233_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb233_pos(%rbp),%rsi
        movq  nb233_faction(%rbp),%rdi
        movq  nb233_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb233_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb233_ninner(%rsp),%ecx
        movl  %ecx,nb233_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb233_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel233_x86_64_sse2.nb233_unroll_loop
        jmp   _nb_kernel233_x86_64_sse2.nb233_checksingle
_nb_kernel233_x86_64_sse2.nb233_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb233_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb233_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb233_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb233_iqM(%rsp),%xmm3
        mulpd  nb233_iqH(%rsp),%xmm4

        movapd  %xmm3,nb233_qqM(%rsp)
        movapd  %xmm4,nb233_qqH(%rsp)

        movq nb233_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb233_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb233_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movlpd (%rsi,%r9,8),%xmm7       ## c6b
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm7      ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb233_c6(%rsp)
        movapd %xmm6,nb233_c12(%rsp)

        movq nb233_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load j coordinates 
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

    subpd nb233_ixO(%rsp),%xmm3
    subpd nb233_iyO(%rsp),%xmm4
    subpd nb233_izO(%rsp),%xmm5

    movapd %xmm3,nb233_dxO(%rsp)
    movapd %xmm4,nb233_dyO(%rsp)
    movapd %xmm5,nb233_dzO(%rsp)

        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5

        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    ## xmm3=rsq

    cvtpd2ps %xmm3,%xmm5
    rsqrtps %xmm5,%xmm5
    cvtps2pd %xmm5,%xmm15    ## lu in low xmm2 

    ## lookup seed in xmm2 
    movapd %xmm15,%xmm5      ## copy of lu 
    mulpd %xmm15,%xmm15       ## lu*lu 
    movapd nb233_three(%rsp),%xmm7
    mulpd %xmm3,%xmm15       ## rsq*lu*lu                    
    movapd nb233_half(%rsp),%xmm6
    subpd %xmm15,%xmm7       ## 30-rsq*lu*lu 
    mulpd %xmm5,%xmm7
    mulpd %xmm6,%xmm7       ## xmm0=iter1 of rinv (new lu) 

    movapd %xmm7,%xmm5      ## copy of lu 
    mulpd %xmm7,%xmm7       ## lu*lu 
    movapd nb233_three(%rsp),%xmm15
    mulpd %xmm3,%xmm7       ## rsq*lu*lu                    
    movapd nb233_half(%rsp),%xmm6
    subpd %xmm7,%xmm15       ## 30-rsq*lu*lu 
    mulpd %xmm5,%xmm15
    mulpd %xmm6,%xmm15       ## xmm15=rinv

    mulpd %xmm15,%xmm3       ## xmm3=r 

    ## xmm15=rinv
    ## xmm3=r

    mulpd nb233_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttpd2pi %xmm3,%mm6

    ## convert back to float
    cvtpi2pd  %mm6,%xmm4

    ## multiply by 8
    pslld   $3,%mm6

    ## calculate eps
    subpd     %xmm4,%xmm3   ## xmm3=eps

    ## move to integer registers
    movd %mm6,%r10d
    psrlq $32,%mm6
    movd %mm6,%r11d

    ## xmm3=eps
    ## xmm15=rinv

        movq nb233_VFtab(%rbp),%rsi
    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movapd (%rsi,%r10,8),%xmm4          ## Y1d F1d  
    movapd (%rsi,%r11,8),%xmm12         ## Y2d F2d 
    movapd 32(%rsi,%r10,8),%xmm8        ## Y1r F1r  
    movapd 32(%rsi,%r11,8),%xmm13       ## Y2r F2r 
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm9
    unpcklpd %xmm12,%xmm4   ## Y1d Y2d 
    unpckhpd %xmm12,%xmm5   ## F1d F2d 
    unpcklpd %xmm13,%xmm8   ## Y1r Y2r 
    unpckhpd %xmm13,%xmm9   ## F1r F2r 

    movapd 16(%rsi,%r10,8),%xmm6        ## G1d H1d  
    movapd 16(%rsi,%r11,8),%xmm12           ## G2d H2d 
    movapd 48(%rsi,%r10,8),%xmm10           ## G1r H1r      
    movapd 48(%rsi,%r11,8),%xmm13           ## G2r H2r 
    movapd %xmm6,%xmm7
    movapd %xmm10,%xmm11
    unpcklpd %xmm12,%xmm6   ## G1d G2d 
    unpckhpd %xmm12,%xmm7   ## H1d H2d 
    unpcklpd %xmm13,%xmm10  ## G1r G2r 
    unpckhpd %xmm13,%xmm11  ## H1r H2r 
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  %xmm3,%xmm7   ## Heps
    mulpd  %xmm3,%xmm11
    mulpd  %xmm3,%xmm6  ## Geps
    mulpd  %xmm3,%xmm10
    mulpd  %xmm3,%xmm7  ## Heps2
    mulpd  %xmm3,%xmm11
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
    mulpd  %xmm3,%xmm5 ## eps*Fp
    mulpd  %xmm3,%xmm9
    movapd nb233_c6(%rsp),%xmm12
    movapd nb233_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb233_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb233_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    mulpd  nb233_tsc(%rsp),%xmm7
    mulpd  %xmm15,%xmm7  ## -fscal
    xorpd  %xmm9,%xmm9

    subpd  %xmm7,%xmm9    ## fscal
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    mulpd  nb233_dxO(%rsp),%xmm9    ## fx/fy/fz
    mulpd  nb233_dyO(%rsp),%xmm10
    mulpd  nb233_dzO(%rsp),%xmm11

    ## save j force temporarily
    movapd %xmm9,nb233_fjx(%rsp)
    movapd %xmm10,nb233_fjy(%rsp)
    movapd %xmm11,nb233_fjz(%rsp)

    ## increment i O force
    addpd nb233_fixO(%rsp),%xmm9
    addpd nb233_fiyO(%rsp),%xmm10
    addpd nb233_fizO(%rsp),%xmm11
    movapd %xmm9,nb233_fixO(%rsp)
    movapd %xmm10,nb233_fiyO(%rsp)
    movapd %xmm11,nb233_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb233_ixH1(%rsp),%xmm0
    subpd nb233_iyH1(%rsp),%xmm1
    subpd nb233_izH1(%rsp),%xmm2
    subpd nb233_ixH2(%rsp),%xmm3
    subpd nb233_iyH2(%rsp),%xmm4
    subpd nb233_izH2(%rsp),%xmm5
    subpd nb233_ixM(%rsp),%xmm6
    subpd nb233_iyM(%rsp),%xmm7
    subpd nb233_izM(%rsp),%xmm8

        movapd %xmm0,nb233_dxH1(%rsp)
        movapd %xmm1,nb233_dyH1(%rsp)
        movapd %xmm2,nb233_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb233_dxH2(%rsp)
        movapd %xmm4,nb233_dyH2(%rsp)
        movapd %xmm5,nb233_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb233_dxM(%rsp)
        movapd %xmm7,nb233_dyM(%rsp)
        movapd %xmm8,nb233_dzM(%rsp)
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

        movapd  nb233_three(%rsp),%xmm9
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

        movapd  nb233_half(%rsp),%xmm15
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

        movapd  nb233_three(%rsp),%xmm1
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

        movapd  nb233_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb233_krf(%rsp),%xmm2
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
    movapd nb233_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb233_qqH(%rsp),%xmm12
    movapd nb233_qqM(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb233_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb233_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j forces to xmm0-xmm2
    movq nb233_faction(%rbp),%rdi
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

    ## add forces from O interaction
    addpd nb233_fjx(%rsp),%xmm0
    addpd nb233_fjy(%rsp),%xmm1
    addpd nb233_fjz(%rsp),%xmm2

        mulpd nb233_dxH1(%rsp),%xmm7
        mulpd nb233_dyH1(%rsp),%xmm8
        mulpd nb233_dzH1(%rsp),%xmm9
        mulpd nb233_dxH2(%rsp),%xmm10
        mulpd nb233_dyH2(%rsp),%xmm11
        mulpd nb233_dzH2(%rsp),%xmm12
        mulpd nb233_dxM(%rsp),%xmm13
        mulpd nb233_dyM(%rsp),%xmm14
        mulpd nb233_dzM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb233_fixH1(%rsp),%xmm7
    addpd nb233_fiyH1(%rsp),%xmm8
    addpd nb233_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb233_fixH2(%rsp),%xmm10
    addpd nb233_fiyH2(%rsp),%xmm11
    addpd nb233_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb233_fixM(%rsp),%xmm13
    addpd nb233_fiyM(%rsp),%xmm14
    addpd nb233_fizM(%rsp),%xmm15

    movapd %xmm7,nb233_fixH1(%rsp)
    movapd %xmm8,nb233_fiyH1(%rsp)
    movapd %xmm9,nb233_fizH1(%rsp)
    movapd %xmm10,nb233_fixH2(%rsp)
    movapd %xmm11,nb233_fiyH2(%rsp)
    movapd %xmm12,nb233_fizH2(%rsp)
    movapd %xmm13,nb233_fixM(%rsp)
    movapd %xmm14,nb233_fiyM(%rsp)
    movapd %xmm15,nb233_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb233_innerk(%rsp)
        jl   _nb_kernel233_x86_64_sse2.nb233_checksingle
        jmp  _nb_kernel233_x86_64_sse2.nb233_unroll_loop
_nb_kernel233_x86_64_sse2.nb233_checksingle: 
        movl  nb233_innerk(%rsp),%edx
        andl  $1,%edx
        jnz  _nb_kernel233_x86_64_sse2.nb233_dosingle
        jmp  _nb_kernel233_x86_64_sse2.nb233_updateouterdata
_nb_kernel233_x86_64_sse2.nb233_dosingle: 
        movq  nb233_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb233_innerjjnr(%rsp)

        movq nb233_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb233_iqM(%rsp),%xmm3
        mulsd  nb233_iqH(%rsp),%xmm4

        movapd  %xmm3,nb233_qqM(%rsp)
        movapd  %xmm4,nb233_qqH(%rsp)

        movq nb233_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb233_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb233_ntia(%rsp),%edi
        addl %edi,%r8d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a
        movhpd 8(%rsi,%r8,8),%xmm6      ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb233_c6(%rsp)
        movapd %xmm6,nb233_c12(%rsp)

        movq nb233_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2  and xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb233_ixO(%rsp),%xmm4
        subsd nb233_iyO(%rsp),%xmm5
        subsd nb233_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxO(%rsp)
        movapd %xmm5,nb233_dyO(%rsp)
        movapd %xmm6,nb233_dzO(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 
        movapd %xmm7,nb233_rsqO(%rsp)

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb233_ixH1(%rsp),%xmm4
        subsd nb233_iyH1(%rsp),%xmm5
        subsd nb233_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb233_dxH1(%rsp)
        movapd %xmm5,nb233_dyH1(%rsp)
        movapd %xmm6,nb233_dzH1(%rsp)
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
        subsd nb233_ixH2(%rsp),%xmm3
        subsd nb233_iyH2(%rsp),%xmm4
        subsd nb233_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb233_dxH2(%rsp)
        movapd %xmm4,nb233_dyH2(%rsp)
        movapd %xmm5,nb233_dzH2(%rsp)
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
        subsd nb233_ixM(%rsp),%xmm4
        subsd nb233_iyM(%rsp),%xmm3
        subsd nb233_izM(%rsp),%xmm2

        ## store dr 
        movapd %xmm4,nb233_dxM(%rsp)
        movapd %xmm3,nb233_dyM(%rsp)
        movapd %xmm2,nb233_dzM(%rsp)

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb233_krf(%rsp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb233_krsqM(%rsp)
        movsd %xmm1,nb233_krsqH2(%rsp)
        movsd %xmm2,nb233_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233_rinvM(%rsp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb233_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb233_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233_half(%rsp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb233_rsqO(%rsp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb233_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movq nb233_VFtab(%rbp),%rsi

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
        mulsd  nb233_two(%rsp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addsd  nb233_Vvdwtot(%rsp),%xmm5
        xorpd  %xmm3,%xmm3
        mulsd  nb233_tsc(%rsp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm3,nb233_fstmp(%rsp)
        movsd %xmm5,nb233_Vvdwtot(%rsp)

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
        mulsd  nb233_two(%rsp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb233_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm7
        mulsd  %xmm4,%xmm5

        addsd  nb233_Vvdwtot(%rsp),%xmm5
        movsd nb233_fstmp(%rsp),%xmm3
        mulsd  nb233_tsc(%rsp),%xmm7
        subsd  %xmm7,%xmm3
        movsd %xmm5,nb233_Vvdwtot(%rsp)

        mulsd  %xmm0,%xmm3


        movsd nb233_dxO(%rsp),%xmm0
        movsd nb233_dyO(%rsp),%xmm1
        movsd nb233_dzO(%rsp),%xmm2

        movq   nb233_faction(%rbp),%rdi
        mulsd  %xmm3,%xmm0
        mulsd  %xmm3,%xmm1
        mulsd  %xmm3,%xmm2

        ## update O forces 
        movapd nb233_fixO(%rsp),%xmm3
        movapd nb233_fiyO(%rsp),%xmm4
        movapd nb233_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixO(%rsp)
        movsd %xmm4,nb233_fiyO(%rsp)
        movsd %xmm7,nb233_fizO(%rsp)
        ## update j forces with water O 
        movsd %xmm0,nb233_fjx(%rsp)
        movsd %xmm1,nb233_fjy(%rsp)
        movsd %xmm2,nb233_fjz(%rsp)

        ## H1 interactions
        movsd  nb233_rinvH1(%rsp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb233_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        mulsd   nb233_two(%rsp),%xmm0
        subsd   nb233_crf(%rsp),%xmm6
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqH(%rsp),%xmm6   ## vcoul 
        mulsd   nb233_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH1 in xmm4 

        addsd  nb233_vctot(%rsp),%xmm6

        movapd nb233_dxH1(%rsp),%xmm0
        movapd nb233_dyH1(%rsp),%xmm1
        movapd nb233_dzH1(%rsp),%xmm2
        movsd %xmm6,nb233_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb233_fixH1(%rsp),%xmm3
        movapd nb233_fiyH1(%rsp),%xmm4
        movapd nb233_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixH1(%rsp)
        movsd %xmm4,nb233_fiyH1(%rsp)
        movsd %xmm7,nb233_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb233_fjx(%rsp),%xmm0
        addsd  nb233_fjy(%rsp),%xmm1
        addsd  nb233_fjz(%rsp),%xmm2
        movsd %xmm0,nb233_fjx(%rsp)
        movsd %xmm1,nb233_fjy(%rsp)
        movsd %xmm2,nb233_fjz(%rsp)

        ## H2 interactions 
        movsd  nb233_rinvH2(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb233_two(%rsp),%xmm0
        subsd   nb233_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqH(%rsp),%xmm5   ## vcoul 
        mulsd   nb233_qqH(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb233_vctot(%rsp),%xmm5

        movapd nb233_dxH2(%rsp),%xmm0
        movapd nb233_dyH2(%rsp),%xmm1
        movapd nb233_dzH2(%rsp),%xmm2
        movsd %xmm5,nb233_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb233_fixH2(%rsp),%xmm3
        movapd nb233_fiyH2(%rsp),%xmm4
        movapd nb233_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixH2(%rsp)
        movsd %xmm4,nb233_fiyH2(%rsp)
        movsd %xmm7,nb233_fizH2(%rsp)
        ## update j forces with water H2 
        addsd  nb233_fjx(%rsp),%xmm0
        addsd  nb233_fjy(%rsp),%xmm1
        addsd  nb233_fjz(%rsp),%xmm2
        movsd %xmm0,nb233_fjx(%rsp)
        movsd %xmm1,nb233_fjy(%rsp)
        movsd %xmm2,nb233_fjz(%rsp)

        ## M interactions 
        movsd  nb233_rinvM(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233_krsqM(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        mulsd   nb233_two(%rsp),%xmm0
        subsd   nb233_crf(%rsp),%xmm5
        subsd   %xmm0,%xmm7     ## xmm7=rinv-2*krsq 
        mulsd   nb233_qqM(%rsp),%xmm5   ## vcoul 
        mulsd   nb233_qqM(%rsp),%xmm7
        mulsd  %xmm7,%xmm4              ## total fsH2 in xmm4 

        addsd  nb233_vctot(%rsp),%xmm5

        movapd nb233_dxM(%rsp),%xmm0
        movapd nb233_dyM(%rsp),%xmm1
        movapd nb233_dzM(%rsp),%xmm2
        movsd %xmm5,nb233_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb233_fixM(%rsp),%xmm3
        movapd nb233_fiyM(%rsp),%xmm4
        movapd nb233_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb233_fixM(%rsp)
        movsd %xmm4,nb233_fiyM(%rsp)
        movsd %xmm7,nb233_fizM(%rsp)

        movq nb233_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb233_fjx(%rsp),%xmm0
        addsd  nb233_fjy(%rsp),%xmm1
        addsd  nb233_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel233_x86_64_sse2.nb233_updateouterdata: 
        movl  nb233_ii3(%rsp),%ecx
        movq  nb233_faction(%rbp),%rdi
        movq  nb233_fshift(%rbp),%rsi
        movl  nb233_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb233_fixO(%rsp),%xmm0
        movapd nb233_fiyO(%rsp),%xmm1
        movapd nb233_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb233_fixH1(%rsp),%xmm0
        movapd nb233_fiyH1(%rsp),%xmm1
        movapd nb233_fizH1(%rsp),%xmm2

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
        movapd nb233_fixH2(%rsp),%xmm0
        movapd nb233_fiyH2(%rsp),%xmm1
        movapd nb233_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb233_fixM(%rsp),%xmm0
        movapd nb233_fiyM(%rsp),%xmm1
        movapd nb233_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movl nb233_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb233_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb233_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb233_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb233_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb233_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb233_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233_x86_64_sse2.nb233_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233_n(%rsp)
        jmp _nb_kernel233_x86_64_sse2.nb233_outer
_nb_kernel233_x86_64_sse2.nb233_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233_x86_64_sse2.nb233_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233_x86_64_sse2.nb233_threadloop
_nb_kernel233_x86_64_sse2.nb233_end: 
        movl nb233_nouter(%rsp),%eax
        movl nb233_ninner(%rsp),%ebx
        movq nb233_outeriter(%rbp),%rcx
        movq nb233_inneriter(%rbp),%rdx
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





.globl nb_kernel233nf_x86_64_sse2
.globl _nb_kernel233nf_x86_64_sse2
nb_kernel233nf_x86_64_sse2:     
_nb_kernel233nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb233nf_fshift, 16
.set nb233nf_gid, 24
.set nb233nf_pos, 32
.set nb233nf_faction, 40
.set nb233nf_charge, 48
.set nb233nf_p_facel, 56
.set nb233nf_argkrf, 64
.set nb233nf_argcrf, 72
.set nb233nf_Vc, 80
.set nb233nf_type, 88
.set nb233nf_p_ntype, 96
.set nb233nf_vdwparam, 104
.set nb233nf_Vvdw, 112
.set nb233nf_p_tabscale, 120
.set nb233nf_VFtab, 128
.set nb233nf_invsqrta, 136
.set nb233nf_dvda, 144
.set nb233nf_p_gbtabscale, 152
.set nb233nf_GBtab, 160
.set nb233nf_p_nthreads, 168
.set nb233nf_count, 176
.set nb233nf_mtx, 184
.set nb233nf_outeriter, 192
.set nb233nf_inneriter, 200
.set nb233nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb233nf_ixO, 0
.set nb233nf_iyO, 16
.set nb233nf_izO, 32
.set nb233nf_ixH1, 48
.set nb233nf_iyH1, 64
.set nb233nf_izH1, 80
.set nb233nf_ixH2, 96
.set nb233nf_iyH2, 112
.set nb233nf_izH2, 128
.set nb233nf_ixM, 144
.set nb233nf_iyM, 160
.set nb233nf_izM, 176
.set nb233nf_iqH, 192
.set nb233nf_iqM, 208
.set nb233nf_qqH, 224
.set nb233nf_qqM, 240
.set nb233nf_c6, 256
.set nb233nf_c12, 272
.set nb233nf_vctot, 288
.set nb233nf_Vvdwtot, 304
.set nb233nf_half, 320
.set nb233nf_three, 336
.set nb233nf_tsc, 352
.set nb233nf_rinvH1, 368
.set nb233nf_rinvH2, 384
.set nb233nf_rinvM, 400
.set nb233nf_krsqH1, 416
.set nb233nf_krsqH2, 432
.set nb233nf_krsqM, 448
.set nb233nf_krf, 464
.set nb233nf_crf, 480
.set nb233nf_rsqO, 496
.set nb233nf_facel, 512
.set nb233nf_iinr, 520
.set nb233nf_jindex, 528
.set nb233nf_jjnr, 536
.set nb233nf_shift, 544
.set nb233nf_shiftvec, 552
.set nb233nf_innerjjnr, 560
.set nb233nf_nri, 568
.set nb233nf_is3, 572
.set nb233nf_ii3, 576
.set nb233nf_ntia, 580
.set nb233nf_innerk, 584
.set nb233nf_n, 588
.set nb233nf_nn1, 592
.set nb233nf_nouter, 596
.set nb233nf_ninner, 600
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
        movl %eax,nb233nf_nouter(%rsp)
        movl %eax,nb233nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb233nf_nri(%rsp)
        movq %rsi,nb233nf_iinr(%rsp)
        movq %rdx,nb233nf_jindex(%rsp)
        movq %rcx,nb233nf_jjnr(%rsp)
        movq %r8,nb233nf_shift(%rsp)
        movq %r9,nb233nf_shiftvec(%rsp)
        movq nb233nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb233nf_facel(%rsp)

        movq nb233nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb233nf_tsc(%rsp)

        movq nb233nf_argkrf(%rbp),%rsi
        movq nb233nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb233nf_krf(%rsp)
        movapd %xmm2,nb233nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb233nf_half(%rsp)
        movl %ebx,nb233nf_half+4(%rsp)
        movsd nb233nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb233nf_half(%rsp)
        movapd %xmm3,nb233nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb233nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb233nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4

        movsd nb233nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb233nf_iqH(%rsp)
        movapd %xmm4,nb233nf_iqM(%rsp)

        movq  nb233nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb233nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233nf_ntia(%rsp)
_nb_kernel233nf_x86_64_sse2.nb233nf_threadloop: 
        movq  nb233nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel233nf_x86_64_sse2.nb233nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233nf_x86_64_sse2.nb233nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233nf_n(%rsp)
        movl %ebx,nb233nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233nf_x86_64_sse2.nb233nf_outerstart
        jmp _nb_kernel233nf_x86_64_sse2.nb233nf_end

_nb_kernel233nf_x86_64_sse2.nb233nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233nf_nouter(%rsp),%ebx
        movl %ebx,nb233nf_nouter(%rsp)

_nb_kernel233nf_x86_64_sse2.nb233nf_outer: 
        movq  nb233nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb233nf_is3(%rsp)            ## store is3 

        movq  nb233nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb233nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb233nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb233nf_ii3(%rsp)

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
        movapd %xmm3,nb233nf_ixO(%rsp)
        movapd %xmm4,nb233nf_iyO(%rsp)
        movapd %xmm5,nb233nf_izO(%rsp)
        movapd %xmm6,nb233nf_ixH1(%rsp)
        movapd %xmm7,nb233nf_iyH1(%rsp)

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
        movapd %xmm6,nb233nf_izH1(%rsp)
        movapd %xmm0,nb233nf_ixH2(%rsp)
        movapd %xmm1,nb233nf_iyH2(%rsp)
        movapd %xmm2,nb233nf_izH2(%rsp)
        movapd %xmm3,nb233nf_ixM(%rsp)
        movapd %xmm4,nb233nf_iyM(%rsp)
        movapd %xmm5,nb233nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb233nf_vctot(%rsp)
        movapd %xmm4,nb233nf_Vvdwtot(%rsp)

        movq  nb233nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb233nf_pos(%rbp),%rsi
        movq  nb233nf_faction(%rbp),%rdi
        movq  nb233nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb233nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb233nf_ninner(%rsp),%ecx
        movl  %ecx,nb233nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb233nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel233nf_x86_64_sse2.nb233nf_unroll_loop
        jmp   _nb_kernel233nf_x86_64_sse2.nb233nf_checksingle
_nb_kernel233nf_x86_64_sse2.nb233nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb233nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb233nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb233nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb233nf_iqM(%rsp),%xmm3
        mulpd  nb233nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb233nf_qqM(%rsp)
        movapd  %xmm4,nb233nf_qqH(%rsp)

        movq nb233nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb233nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb233nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb233nf_c6(%rsp)
        movapd %xmm6,nb233nf_c12(%rsp)

        movq nb233nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb233nf_ixO(%rsp),%xmm4
        movapd nb233nf_iyO(%rsp),%xmm5
        movapd nb233nf_izO(%rsp),%xmm6

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
        movapd nb233nf_ixH1(%rsp),%xmm4
        movapd nb233nf_iyH1(%rsp),%xmm5
        movapd nb233nf_izH1(%rsp),%xmm6

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
        movapd nb233nf_ixH2(%rsp),%xmm3
        movapd nb233nf_iyH2(%rsp),%xmm4
        movapd nb233nf_izH2(%rsp),%xmm5

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
        movapd nb233nf_iyM(%rsp),%xmm3
        movapd nb233nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2


        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movapd %xmm7,nb233nf_rsqO(%rsp)

        ## calculate krsq
        movapd nb233nf_krf(%rsp),%xmm0
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        mulpd %xmm4,%xmm0
        mulpd %xmm5,%xmm1
        mulpd %xmm6,%xmm2
        movapd %xmm0,nb233nf_krsqM(%rsp)
        movapd %xmm1,nb233nf_krsqH2(%rsp)
        movapd %xmm2,nb233nf_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,nb233nf_rinvM(%rsp)


        ## rsqO - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb233nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb233nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvO in xmm7 



        movapd nb233nf_rsqO(%rsp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulpd %xmm7,%xmm4       ## xmm4=r 
        mulpd nb233nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movq nb233nf_VFtab(%rbp),%rsi
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

        movapd nb233nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb233nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb233nf_Vvdwtot(%rsp)

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

        movapd nb233nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb233nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb233nf_Vvdwtot(%rsp)

        ## H1 interactions 
        movapd  nb233nf_rinvH1(%rsp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movapd  %xmm6,%xmm7
        movapd  nb233nf_krsqH1(%rsp),%xmm0
        addpd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subpd   nb233nf_crf(%rsp),%xmm6
        mulpd   nb233nf_qqH(%rsp),%xmm6   ## vcoul 
        addpd  nb233nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb233nf_vctot(%rsp)

        ## H2 interactions 
        movapd  nb233nf_rinvH2(%rsp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233nf_krsqH2(%rsp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb233nf_crf(%rsp),%xmm5
        mulpd   nb233nf_qqH(%rsp),%xmm5   ## vcoul 
        addpd  nb233nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb233nf_vctot(%rsp)

        ## M interactions 
        movapd  nb233nf_rinvM(%rsp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movapd  %xmm5,%xmm7
        movapd  nb233nf_krsqM(%rsp),%xmm0
        addpd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subpd   nb233nf_crf(%rsp),%xmm5
        mulpd   nb233nf_qqM(%rsp),%xmm5   ## vcoul 
        addpd  nb233nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb233nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb233nf_innerk(%rsp)
        jl   _nb_kernel233nf_x86_64_sse2.nb233nf_checksingle
        jmp  _nb_kernel233nf_x86_64_sse2.nb233nf_unroll_loop
_nb_kernel233nf_x86_64_sse2.nb233nf_checksingle: 
        movl  nb233nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz  _nb_kernel233nf_x86_64_sse2.nb233nf_dosingle
        jmp  _nb_kernel233nf_x86_64_sse2.nb233nf_updateouterdata
_nb_kernel233nf_x86_64_sse2.nb233nf_dosingle: 
        movq  nb233nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb233nf_innerjjnr(%rsp)

        movq nb233nf_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulsd  nb233nf_iqM(%rsp),%xmm3
        mulsd  nb233nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb233nf_qqM(%rsp)
        movapd  %xmm4,nb233nf_qqH(%rsp)

        movq nb233nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb233nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb233nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb233nf_c6(%rsp)
        movapd %xmm6,nb233nf_c12(%rsp)

        movq nb233nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb233nf_ixO(%rsp),%xmm4
        movapd nb233nf_iyO(%rsp),%xmm5
        movapd nb233nf_izO(%rsp),%xmm6

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
        movapd %xmm7,nb233nf_rsqO(%rsp)

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb233nf_ixH1(%rsp),%xmm4
        movapd nb233nf_iyH1(%rsp),%xmm5
        movapd nb233nf_izH1(%rsp),%xmm6

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
        movapd nb233nf_ixH2(%rsp),%xmm3
        movapd nb233nf_iyH2(%rsp),%xmm4
        movapd nb233nf_izH2(%rsp),%xmm5

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
        movapd nb233nf_iyM(%rsp),%xmm3
        movapd nb233nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb233nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## calculate krsq
        movsd nb233nf_krf(%rsp),%xmm0
        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2
        mulsd %xmm4,%xmm0
        mulsd %xmm5,%xmm1
        mulsd %xmm6,%xmm2
        movsd %xmm0,nb233nf_krsqM(%rsp)
        movsd %xmm1,nb233nf_krsqH2(%rsp)
        movsd %xmm2,nb233nf_krsqH1(%rsp)

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb233nf_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb233nf_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb233nf_rinvM(%rsp)

        ## rsqO - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movsd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movsd  nb233nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb233nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movsd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movsd nb233nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb233nf_half(%rsp),%xmm4   ## rinv 
        movsd  %xmm4,%xmm7      ## rinvO in xmm7 

        movsd nb233nf_rsqO(%rsp),%xmm4
        movapd %xmm7,%xmm0
        ## LJ table interaction.
        mulsd %xmm7,%xmm4       ## xmm4=r 
        mulsd nb233nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movq nb233nf_VFtab(%rbp),%rsi

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

        movsd nb233nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb233nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb233nf_Vvdwtot(%rsp)

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

        movsd nb233nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb233nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb233nf_Vvdwtot(%rsp)

        ## H1 interactions
        movsd  nb233nf_rinvH1(%rsp),%xmm6
        movsd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        movsd  %xmm6,%xmm7
        movsd  nb233nf_krsqH1(%rsp),%xmm0
        addsd   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subsd   nb233nf_crf(%rsp),%xmm6
        mulsd   nb233nf_qqH(%rsp),%xmm6   ## vcoul 
        addsd  nb233nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb233nf_vctot(%rsp)

        ## H2 interactions 
        movsd  nb233nf_rinvH2(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233nf_krsqH2(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb233nf_crf(%rsp),%xmm5
        mulsd   nb233nf_qqH(%rsp),%xmm5   ## vcoul 
        addsd  nb233nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb233nf_vctot(%rsp)

        ## M interactions 
        movsd  nb233nf_rinvM(%rsp),%xmm5
        movsd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        movsd  %xmm5,%xmm7
        movsd  nb233nf_krsqM(%rsp),%xmm0
        addsd   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subsd   nb233nf_crf(%rsp),%xmm5
        mulsd   nb233nf_qqM(%rsp),%xmm5   ## vcoul 
        addsd  nb233nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb233nf_vctot(%rsp)

_nb_kernel233nf_x86_64_sse2.nb233nf_updateouterdata: 
        ## get n from stack
        movl nb233nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb233nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb233nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb233nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb233nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb233nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb233nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233nf_x86_64_sse2.nb233nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233nf_n(%rsp)
        jmp _nb_kernel233nf_x86_64_sse2.nb233nf_outer
_nb_kernel233nf_x86_64_sse2.nb233nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233nf_x86_64_sse2.nb233nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233nf_x86_64_sse2.nb233nf_threadloop
_nb_kernel233nf_x86_64_sse2.nb233nf_end: 
        movl nb233nf_nouter(%rsp),%eax
        movl nb233nf_ninner(%rsp),%ebx
        movq nb233nf_outeriter(%rbp),%rcx
        movq nb233nf_inneriter(%rbp),%rdx
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


