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




.globl nb_kernel330_x86_64_sse2
.globl _nb_kernel330_x86_64_sse2
nb_kernel330_x86_64_sse2:       
_nb_kernel330_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb330_fshift, 16
.set nb330_gid, 24
.set nb330_pos, 32
.set nb330_faction, 40
.set nb330_charge, 48
.set nb330_p_facel, 56
.set nb330_argkrf, 64
.set nb330_argcrf, 72
.set nb330_Vc, 80
.set nb330_type, 88
.set nb330_p_ntype, 96
.set nb330_vdwparam, 104
.set nb330_Vvdw, 112
.set nb330_p_tabscale, 120
.set nb330_VFtab, 128
.set nb330_invsqrta, 136
.set nb330_dvda, 144
.set nb330_p_gbtabscale, 152
.set nb330_GBtab, 160
.set nb330_p_nthreads, 168
.set nb330_count, 176
.set nb330_mtx, 184
.set nb330_outeriter, 192
.set nb330_inneriter, 200
.set nb330_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb330_ix, 0
.set nb330_iy, 16
.set nb330_iz, 32
.set nb330_iq, 48
.set nb330_dx, 64
.set nb330_dy, 80
.set nb330_dz, 96
.set nb330_rinv, 112
.set nb330_tsc, 128
.set nb330_qq, 144
.set nb330_c6, 160
.set nb330_c12, 176
.set nb330_eps, 192
.set nb330_vctot, 208
.set nb330_Vvdwtot, 224
.set nb330_fix, 240
.set nb330_fiy, 256
.set nb330_fiz, 272
.set nb330_half, 288
.set nb330_three, 304
.set nb330_nri, 320
.set nb330_iinr, 328
.set nb330_jindex, 336
.set nb330_jjnr, 344
.set nb330_shift, 352
.set nb330_shiftvec, 360
.set nb330_facel, 368
.set nb330_innerjjnr, 376
.set nb330_is3, 384
.set nb330_ii3, 388
.set nb330_ntia, 392
.set nb330_innerk, 396
.set nb330_n, 400
.set nb330_nn1, 404
.set nb330_ntype, 408
.set nb330_nouter, 412
.set nb330_ninner, 416
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb330_nouter(%rsp)
        movl %eax,nb330_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb330_nri(%rsp)
        movq %rsi,nb330_iinr(%rsp)
        movq %rdx,nb330_jindex(%rsp)
        movq %rcx,nb330_jjnr(%rsp)
        movq %r8,nb330_shift(%rsp)
        movq %r9,nb330_shiftvec(%rsp)
        movq nb330_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb330_ntype(%rsp)
        movq nb330_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb330_facel(%rsp)

        movq nb330_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb330_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb330_half(%rsp)
        movl %ebx,nb330_half+4(%rsp)
        movsd nb330_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb330_half(%rsp)
        movapd %xmm3,nb330_three(%rsp)

_nb_kernel330_x86_64_sse2.nb330_threadloop: 
        movq  nb330_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel330_x86_64_sse2.nb330_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330_x86_64_sse2.nb330_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330_n(%rsp)
        movl %ebx,nb330_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330_x86_64_sse2.nb330_outerstart
        jmp _nb_kernel330_x86_64_sse2.nb330_end

_nb_kernel330_x86_64_sse2.nb330_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330_nouter(%rsp),%ebx
        movl %ebx,nb330_nouter(%rsp)

_nb_kernel330_x86_64_sse2.nb330_outer: 
        movq  nb330_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb330_is3(%rsp)      ## store is3 

        movq  nb330_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb330_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb330_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb330_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb330_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb330_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb330_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb330_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb330_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb330_ix(%rsp)
        movapd %xmm1,nb330_iy(%rsp)
        movapd %xmm2,nb330_iz(%rsp)

        movl  %ebx,nb330_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb330_vctot(%rsp)
        movapd %xmm4,nb330_Vvdwtot(%rsp)
        movapd %xmm4,nb330_fix(%rsp)
        movapd %xmm4,nb330_fiy(%rsp)
        movapd %xmm4,nb330_fiz(%rsp)

        movq  nb330_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb330_pos(%rbp),%rsi
        movq  nb330_faction(%rbp),%rdi
        movq  nb330_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb330_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330_ninner(%rsp),%ecx
        movl  %ecx,nb330_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb330_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel330_x86_64_sse2.nb330_unroll_loop
        jmp   _nb_kernel330_x86_64_sse2.nb330_checksingle
_nb_kernel330_x86_64_sse2.nb330_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb330_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%r12d
        movl  4(%rdx),%r13d
        addq $8,nb330_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 


        movq nb330_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%r12,%r12,2),%rax     ## replace jnr with j3 
        lea  (%r13,%r13,2),%rbx

        ## move two coordinates to xmm4-xmm6 
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb330_ix(%rsp),%xmm4
        subpd nb330_iy(%rsp),%xmm5
        subpd nb330_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb330_dx(%rsp)
        movapd %xmm5,nb330_dy(%rsp)
        movapd %xmm6,nb330_dz(%rsp)

        movq nb330_charge(%rbp),%rsi     ## base of charge[] 

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 
        movlpd (%rsi,%r12,8),%xmm3

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 
        movhpd (%rsi,%r13,8),%xmm3
        movq nb330_type(%rbp),%rsi

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb330_three(%rsp),%xmm1
        mulpd  nb330_iq(%rsp),%xmm3
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 
        movapd %xmm3,nb330_qq(%rsp)

        movl (%rsi,%r12,4),%r8d
        movl (%rsi,%r13,4),%r9d

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb330_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb330_tsc(%rsp),%xmm4
    movapd %xmm0,%xmm15 ## copy of rinv
        shll %r8d
        shll %r9d
        movl nb330_ntia(%rsp),%edi

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        addl %edi,%r8d
        addl %edi,%r9d
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,nb330_eps(%rsp)

        pslld $2,%mm6           ## idx *= 4 

        movq nb330_VFtab(%rbp),%rsi
        movd %mm6,%r10d
        psrlq $32,%mm6
        movd %mm6,%r11d         ## indices in r10,r11

        lea  (%r10,%r10,2),%r10        ## idx*=3 (12 total now) 
        lea  (%r11,%r11,2),%r11        ## idx*=3 (12 total now) 

    ## load Coulomb and LJ table data in parallel
    movapd (%rsi,%r10,8),%xmm0         ## Y1c F1c
    movapd (%rsi,%r11,8),%xmm12        ## Y2c F2c
    movapd 32(%rsi,%r10,8),%xmm4       ## Y1d F1d
    movapd 32(%rsi,%r11,8),%xmm13      ## Y2d F2d
    movapd 64(%rsi,%r10,8),%xmm8       ## Y1r F1r
    movapd 64(%rsi,%r11,8),%xmm14      ## Y2r F2r
        movapd %xmm0,%xmm1
        movapd %xmm4,%xmm5
        movapd %xmm8,%xmm9
        unpcklpd %xmm12,%xmm0   ## Y1c Y2c 
        unpckhpd %xmm12,%xmm1   ## F1c F2c 
        unpcklpd %xmm13,%xmm4   ## Y1d Y2d 
        unpckhpd %xmm13,%xmm5   ## F1d F2d 
        unpcklpd %xmm14,%xmm8   ## Y1r Y2r 
        unpckhpd %xmm14,%xmm9   ## F1r F2r 

    movapd 16(%rsi,%r10,8),%xmm2       ## G1c H1c
    movapd 16(%rsi,%r11,8),%xmm12      ## G2c H2c
    movapd 48(%rsi,%r10,8),%xmm6       ## G1d H1d
    movapd 48(%rsi,%r11,8),%xmm13      ## G2d H2d
    movapd 80(%rsi,%r10,8),%xmm10      ## G1r H1r
    movapd 80(%rsi,%r11,8),%xmm14      ## G2r H2r
        movapd %xmm2,%xmm3
        movapd %xmm6,%xmm7
        movapd %xmm10,%xmm11
        unpcklpd %xmm12,%xmm2   ## G1c G2c 
        unpckhpd %xmm12,%xmm3   ## H1c H2c 
        unpcklpd %xmm13,%xmm6   ## G1d G2d 
        unpckhpd %xmm13,%xmm7   ## H1d H2d 
        unpcklpd %xmm14,%xmm10  ## G1r G2r 
        unpckhpd %xmm14,%xmm11  ## H1r H2r 
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11
        movq nb330_vdwparam(%rbp),%rsi

    movapd nb330_eps(%rsp),%xmm12

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm12,%xmm7
    mulpd  %xmm12,%xmm11
    mulpd  %xmm12,%xmm2    ## Geps
    mulpd  %xmm12,%xmm6
    mulpd  %xmm12,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm12,%xmm7
    mulpd  %xmm12,%xmm11

    movlpd (%rsi,%r8,8),%xmm13   ## c6
    movlpd 8(%rsi,%r8,8),%xmm14     ## c12

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
    movhpd (%rsi,%r9,8),%xmm13   ## c6
    movhpd 8(%rsi,%r9,8),%xmm14     ## c12

    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm12,%xmm5
    mulpd  %xmm12,%xmm9
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  nb330_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulpd  %xmm13,%xmm5  ## vnb6
    mulpd  %xmm14,%xmm9  ## vnb12
    mulpd  nb330_qq(%rsp),%xmm3      ## FF*qq = fij
    mulpd  %xmm13,%xmm7  ## fijD
    mulpd  %xmm14,%xmm11  ##fijR

    ## accumulate vctot
    addpd  nb330_vctot(%rsp),%xmm1
    movapd %xmm1,nb330_vctot(%rsp)

    ## accumulate Vvdwtot
    addpd  nb330_Vvdwtot(%rsp),%xmm5
    addpd  %xmm9,%xmm5
    movapd %xmm5,nb330_Vvdwtot(%rsp)

    xorpd  %xmm9,%xmm9

    ## the fj's - start by accumulating forces from memory 
    movq nb330_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm4
        movlpd 8(%rdi,%rax,8),%xmm5
        movlpd 16(%rdi,%rax,8),%xmm6
        movhpd (%rdi,%rbx,8),%xmm4
        movhpd 8(%rdi,%rbx,8),%xmm5
        movhpd 16(%rdi,%rbx,8),%xmm6

    addpd  %xmm7,%xmm3
    addpd  %xmm11,%xmm3
    mulpd  %xmm15,%xmm3
    mulpd  nb330_tsc(%rsp),%xmm3    ## fscal

    subpd  %xmm3,%xmm9
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    movapd nb330_fix(%rsp),%xmm12
    movapd nb330_fiy(%rsp),%xmm13
    movapd nb330_fiz(%rsp),%xmm14

    mulpd  nb330_dx(%rsp),%xmm9
    mulpd  nb330_dy(%rsp),%xmm10
    mulpd  nb330_dz(%rsp),%xmm11

        addpd %xmm9,%xmm4
        addpd %xmm10,%xmm5
        addpd %xmm11,%xmm6

    ## accumulate i forces
    addpd %xmm9,%xmm12
    addpd %xmm10,%xmm13
    addpd %xmm11,%xmm14
        movlpd %xmm4,(%rdi,%rax,8)
        movlpd %xmm5,8(%rdi,%rax,8)
        movlpd %xmm6,16(%rdi,%rax,8)

    movapd %xmm12,nb330_fix(%rsp)
    movapd %xmm13,nb330_fiy(%rsp)
    movapd %xmm14,nb330_fiz(%rsp)


        movhpd %xmm4,(%rdi,%rbx,8)
        movhpd %xmm5,8(%rdi,%rbx,8)
        movhpd %xmm6,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb330_innerk(%rsp)
        jl    _nb_kernel330_x86_64_sse2.nb330_checksingle
        jmp   _nb_kernel330_x86_64_sse2.nb330_unroll_loop
_nb_kernel330_x86_64_sse2.nb330_checksingle: 
        movl  nb330_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330_x86_64_sse2.nb330_dosingle
        jmp    _nb_kernel330_x86_64_sse2.nb330_updateouterdata
_nb_kernel330_x86_64_sse2.nb330_dosingle: 
        movq nb330_charge(%rbp),%rsi
        movq nb330_pos(%rbp),%rdi
        movq  nb330_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movsd (%rsi,%rax,8),%xmm3
        mulsd  nb330_iq(%rsp),%xmm3
        movapd %xmm3,nb330_qq(%rsp)

        movq nb330_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb330_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb330_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb330_c6(%rsp)
        movapd %xmm6,nb330_c12(%rsp)

        movq nb330_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinate to xmm4-xmm6 
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        movq   nb330_faction(%rbp),%rdi

        ## calc dr 
        subsd nb330_ix(%rsp),%xmm4
        subsd nb330_iy(%rsp),%xmm5
        subsd nb330_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb330_dx(%rsp)
        movapd %xmm5,nb330_dy(%rsp)
        movapd %xmm6,nb330_dz(%rsp)

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb330_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb330_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb330_tsc(%rsp),%xmm4
    movapd %xmm0,%xmm15 ## copy of rinv

        cvttsd2si %xmm4,%r10d   ## mm6 = lu idx 
        cvtsi2sd %r10d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,nb330_eps(%rsp)

        shll  $2,%r10d          ## idx *= 4 

        movq nb330_VFtab(%rbp),%rsi

        lea  (%r10,%r10,2),%r10        ## idx*=3 (12 total now) 

    ## load Coulomb and LJ table data in parallel    
    movapd (%rsi,%r10,8),%xmm0         ## Y1c F1c
    movapd 32(%rsi,%r10,8),%xmm4       ## Y1d F1d
    movapd 64(%rsi,%r10,8),%xmm8       ## Y1r F1r
        movhlps %xmm0,%xmm1
        movhlps %xmm4,%xmm5
        movhlps %xmm8,%xmm9

    movapd 16(%rsi,%r10,8),%xmm2       ## G1c H1c
    movapd 48(%rsi,%r10,8),%xmm6       ## G1d H1d
    movapd 80(%rsi,%r10,8),%xmm10      ## G1r H1r
        movhlps %xmm2,%xmm3
        movhlps %xmm6,%xmm7
        movhlps %xmm10,%xmm11
    ## table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movapd nb330_eps(%rsp),%xmm12

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm12,%xmm7
    mulsd  %xmm12,%xmm11
    mulsd  %xmm12,%xmm2    ## Geps
    mulsd  %xmm12,%xmm6
    mulsd  %xmm12,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm12,%xmm7
    mulsd  %xmm12,%xmm11

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
    mulsd  nb330_eps(%rsp),%xmm1     ## eps*Fp
    mulsd  nb330_eps(%rsp),%xmm5
    mulsd  nb330_eps(%rsp),%xmm9
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  nb330_qq(%rsp),%xmm1     ## VV*qq = vcoul
    mulsd  nb330_c6(%rsp),%xmm5     ## vnb6
    mulsd  nb330_c12(%rsp),%xmm9     ## vnb12
    mulsd  nb330_qq(%rsp),%xmm3      ## FF*qq = fij
    mulsd  nb330_c6(%rsp),%xmm7     ## fijD
    mulsd  nb330_c12(%rsp),%xmm11     ##fijR

    ## accumulate vctot
    addsd  nb330_vctot(%rsp),%xmm1
    movsd %xmm1,nb330_vctot(%rsp)

    ## accumulate Vvdwtot
    addsd  nb330_Vvdwtot(%rsp),%xmm5
    addsd  %xmm9,%xmm5
    movsd %xmm5,nb330_Vvdwtot(%rsp)

    xorpd  %xmm9,%xmm9

    addsd  %xmm7,%xmm3
    addsd  %xmm11,%xmm3
    mulsd  %xmm15,%xmm3
    mulsd  nb330_tsc(%rsp),%xmm3    ## fscal

    subsd  %xmm3,%xmm9
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    movapd nb330_fix(%rsp),%xmm12
    movapd nb330_fiy(%rsp),%xmm13
    movapd nb330_fiz(%rsp),%xmm14

    mulsd  nb330_dx(%rsp),%xmm9
    mulsd  nb330_dy(%rsp),%xmm10
    mulsd  nb330_dz(%rsp),%xmm11

    ## accumulate i forces
    addsd %xmm9,%xmm12
    addsd %xmm10,%xmm13
    addsd %xmm11,%xmm14
    movsd %xmm12,nb330_fix(%rsp)
    movsd %xmm13,nb330_fiy(%rsp)
    movsd %xmm14,nb330_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb330_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel330_x86_64_sse2.nb330_updateouterdata: 
        movl  nb330_ii3(%rsp),%ecx
        movq  nb330_faction(%rbp),%rdi
        movq  nb330_fshift(%rbp),%rsi
        movl  nb330_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb330_fix(%rsp),%xmm0
        movapd nb330_fiy(%rsp),%xmm1
        movapd nb330_fiz(%rsp),%xmm2

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

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb330_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb330_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb330_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb330_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb330_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb330_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb330_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330_x86_64_sse2.nb330_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330_n(%rsp)
        jmp _nb_kernel330_x86_64_sse2.nb330_outer
_nb_kernel330_x86_64_sse2.nb330_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330_x86_64_sse2.nb330_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330_x86_64_sse2.nb330_threadloop
_nb_kernel330_x86_64_sse2.nb330_end: 
        movl nb330_nouter(%rsp),%eax
        movl nb330_ninner(%rsp),%ebx
        movq nb330_outeriter(%rbp),%rcx
        movq nb330_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel330nf_x86_64_sse2
.globl _nb_kernel330nf_x86_64_sse2
nb_kernel330nf_x86_64_sse2:     
_nb_kernel330nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb330nf_fshift, 16
.set nb330nf_gid, 24
.set nb330nf_pos, 32
.set nb330nf_faction, 40
.set nb330nf_charge, 48
.set nb330nf_p_facel, 56
.set nb330nf_argkrf, 64
.set nb330nf_argcrf, 72
.set nb330nf_Vc, 80
.set nb330nf_type, 88
.set nb330nf_p_ntype, 96
.set nb330nf_vdwparam, 104
.set nb330nf_Vvdw, 112
.set nb330nf_p_tabscale, 120
.set nb330nf_VFtab, 128
.set nb330nf_invsqrta, 136
.set nb330nf_dvda, 144
.set nb330nf_p_gbtabscale, 152
.set nb330nf_GBtab, 160
.set nb330nf_p_nthreads, 168
.set nb330nf_count, 176
.set nb330nf_mtx, 184
.set nb330nf_outeriter, 192
.set nb330nf_inneriter, 200
.set nb330nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330nf_ix, 0
.set nb330nf_iy, 16
.set nb330nf_iz, 32
.set nb330nf_iq, 48
.set nb330nf_tsc, 64
.set nb330nf_qq, 80
.set nb330nf_c6, 96
.set nb330nf_c12, 112
.set nb330nf_vctot, 128
.set nb330nf_Vvdwtot, 144
.set nb330nf_half, 160
.set nb330nf_three, 176
.set nb330nf_nri, 192
.set nb330nf_iinr, 200
.set nb330nf_jindex, 208
.set nb330nf_jjnr, 216
.set nb330nf_shift, 224
.set nb330nf_shiftvec, 232
.set nb330nf_facel, 240
.set nb330nf_innerjjnr, 248
.set nb330nf_is3, 256
.set nb330nf_ii3, 260
.set nb330nf_ntia, 264
.set nb330nf_innerk, 268
.set nb330nf_n, 272
.set nb330nf_nn1, 276
.set nb330nf_ntype, 280
.set nb330nf_nouter, 284
.set nb330nf_ninner, 288
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $312,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb330nf_nouter(%rsp)
        movl %eax,nb330nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb330nf_nri(%rsp)
        movq %rsi,nb330nf_iinr(%rsp)
        movq %rdx,nb330nf_jindex(%rsp)
        movq %rcx,nb330nf_jjnr(%rsp)
        movq %r8,nb330nf_shift(%rsp)
        movq %r9,nb330nf_shiftvec(%rsp)
        movq nb330nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb330nf_ntype(%rsp)
        movq nb330nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb330nf_facel(%rsp)

        movq nb330nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb330nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb330nf_half(%rsp)
        movl %ebx,nb330nf_half+4(%rsp)
        movsd nb330nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb330nf_half(%rsp)
        movapd %xmm3,nb330nf_three(%rsp)

_nb_kernel330nf_x86_64_sse2.nb330nf_threadloop: 
        movq  nb330nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel330nf_x86_64_sse2.nb330nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330nf_x86_64_sse2.nb330nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330nf_n(%rsp)
        movl %ebx,nb330nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330nf_x86_64_sse2.nb330nf_outerstart
        jmp _nb_kernel330nf_x86_64_sse2.nb330nf_end

_nb_kernel330nf_x86_64_sse2.nb330nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330nf_nouter(%rsp),%ebx
        movl %ebx,nb330nf_nouter(%rsp)

_nb_kernel330nf_x86_64_sse2.nb330nf_outer: 
        movq  nb330nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb330nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb330nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb330nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb330nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb330nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb330nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb330nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb330nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb330nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb330nf_ix(%rsp)
        movapd %xmm1,nb330nf_iy(%rsp)
        movapd %xmm2,nb330nf_iz(%rsp)

        movl  %ebx,nb330nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb330nf_vctot(%rsp)
        movapd %xmm4,nb330nf_Vvdwtot(%rsp)

        movq  nb330nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb330nf_pos(%rbp),%rsi
        movq  nb330nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb330nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb330nf_ninner(%rsp),%ecx
        movl  %ecx,nb330nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb330nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel330nf_x86_64_sse2.nb330nf_unroll_loop
        jmp   _nb_kernel330nf_x86_64_sse2.nb330nf_checksingle
_nb_kernel330nf_x86_64_sse2.nb330nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb330nf_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb330nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb330nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb330nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb330nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb330nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb330nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb330nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb330nf_c6(%rsp)
        movapd %xmm6,nb330nf_c12(%rsp)

        movq nb330nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move nb330nf_ix-iz to xmm4-xmm6 
        movapd nb330nf_ix(%rsp),%xmm4
        movapd nb330nf_iy(%rsp),%xmm5
        movapd nb330nf_iz(%rsp),%xmm6

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
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb330nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb330nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb330nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movq nb330nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        lea  (%rax,%rax,2),%rax        ## idx*=3 (12 total now) 
        lea  (%rbx,%rbx,2),%rbx        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        movapd nb330nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb330nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb330nf_vctot(%rsp)

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

        mulpd  nb330nf_c6(%rsp),%xmm5   ## Vvdw6 

        addpd  nb330nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb330nf_Vvdwtot(%rsp)

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

        mulpd  nb330nf_c12(%rsp),%xmm5   ## Vvdw12 

        addpd  nb330nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb330nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb330nf_innerk(%rsp)
        jl    _nb_kernel330nf_x86_64_sse2.nb330nf_checksingle
        jmp   _nb_kernel330nf_x86_64_sse2.nb330nf_unroll_loop
_nb_kernel330nf_x86_64_sse2.nb330nf_checksingle: 
        movl  nb330nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330nf_x86_64_sse2.nb330nf_dosingle
        jmp    _nb_kernel330nf_x86_64_sse2.nb330nf_updateouterdata
_nb_kernel330nf_x86_64_sse2.nb330nf_dosingle: 
        movq nb330nf_charge(%rbp),%rsi
        movq nb330nf_pos(%rbp),%rdi
        movq  nb330nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorpd  %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3      ## xmm6(0) has the charge       
        mulpd  nb330nf_iq(%rsp),%xmm3
        movapd %xmm3,nb330nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movq nb330nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb330nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb330nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb330nf_c6(%rsp)
        movapd %xmm6,nb330nf_c12(%rsp)

        movq nb330nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move nb330nf_ix-iz to xmm4-xmm6 
        movapd nb330nf_ix(%rsp),%xmm4
        movapd nb330nf_iy(%rsp),%xmm5
        movapd nb330nf_iz(%rsp),%xmm6

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
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb330nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb330nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb330nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb330nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb330nf_tsc(%rsp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb330nf_VFtab(%rbp),%rsi
        lea  (%rax,%rax,2),%rax        ## idx*=3 (12 total now) 

        ## Coulomb 
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
        movapd nb330nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb330nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb330nf_vctot(%rsp)

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
        movapd nb330nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb330nf_c6(%rsp),%xmm5    ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb330nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb330nf_Vvdwtot(%rsp)

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
        movapd nb330nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb330nf_c12(%rsp),%xmm5   ## Vvdw12 

        addsd  nb330nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb330nf_Vvdwtot(%rsp)

_nb_kernel330nf_x86_64_sse2.nb330nf_updateouterdata: 
        ## get n from stack
        movl nb330nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb330nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb330nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb330nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb330nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb330nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb330nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330nf_x86_64_sse2.nb330nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330nf_n(%rsp)
        jmp _nb_kernel330nf_x86_64_sse2.nb330nf_outer
_nb_kernel330nf_x86_64_sse2.nb330nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330nf_x86_64_sse2.nb330nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330nf_x86_64_sse2.nb330nf_threadloop
_nb_kernel330nf_x86_64_sse2.nb330nf_end: 
        movl nb330nf_nouter(%rsp),%eax
        movl nb330nf_ninner(%rsp),%ebx
        movq nb330nf_outeriter(%rbp),%rcx
        movq nb330nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $312,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

