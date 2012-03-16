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






.globl nb_kernel030_x86_64_sse2
.globl _nb_kernel030_x86_64_sse2
nb_kernel030_x86_64_sse2:       
_nb_kernel030_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb030_fshift, 16
.set nb030_gid, 24
.set nb030_pos, 32
.set nb030_faction, 40
.set nb030_charge, 48
.set nb030_p_facel, 56
.set nb030_argkrf, 64
.set nb030_argcrf, 72
.set nb030_Vc, 80
.set nb030_type, 88
.set nb030_p_ntype, 96
.set nb030_vdwparam, 104
.set nb030_Vvdw, 112
.set nb030_p_tabscale, 120
.set nb030_VFtab, 128
.set nb030_invsqrta, 136
.set nb030_dvda, 144
.set nb030_p_gbtabscale, 152
.set nb030_GBtab, 160
.set nb030_p_nthreads, 168
.set nb030_count, 176
.set nb030_mtx, 184
.set nb030_outeriter, 192
.set nb030_inneriter, 200
.set nb030_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb030_ix, 0
.set nb030_iy, 16
.set nb030_iz, 32
.set nb030_dx, 48
.set nb030_dy, 64
.set nb030_dz, 80
.set nb030_two, 96
.set nb030_tsc, 112
.set nb030_c6, 128
.set nb030_c12, 144
.set nb030_fscal, 160
.set nb030_Vvdwtot, 176
.set nb030_fix, 192
.set nb030_fiy, 208
.set nb030_fiz, 224
.set nb030_half, 240
.set nb030_three, 256
.set nb030_is3, 272
.set nb030_ii3, 276
.set nb030_nri, 280
.set nb030_iinr, 288
.set nb030_jindex, 296
.set nb030_jjnr, 304
.set nb030_shift, 312
.set nb030_shiftvec, 320
.set nb030_innerjjnr, 328
.set nb030_ntia, 336
.set nb030_innerk, 340
.set nb030_n, 344
.set nb030_nn1, 348
.set nb030_ntype, 352
.set nb030_nouter, 356
.set nb030_ninner, 360

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $376,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb030_nouter(%rsp)
        movl %eax,nb030_ninner(%rsp)



        movl (%rdi),%edi
        movl %edi,nb030_nri(%rsp)
        movq %rsi,nb030_iinr(%rsp)
        movq %rdx,nb030_jindex(%rsp)
        movq %rcx,nb030_jjnr(%rsp)
        movq %r8,nb030_shift(%rsp)
        movq %r9,nb030_shiftvec(%rsp)
        movq nb030_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb030_ntype(%rsp)

        movq nb030_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb030_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb030_half(%rsp)
        movl %ebx,nb030_half+4(%rsp)
        movsd nb030_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb030_half(%rsp)
        movapd %xmm2,nb030_two(%rsp)
        movapd %xmm3,nb030_three(%rsp)

_nb_kernel030_x86_64_sse2.nb030_threadloop: 
        movq  nb030_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel030_x86_64_sse2.nb030_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030_x86_64_sse2.nb030_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030_n(%rsp)
        movl %ebx,nb030_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030_x86_64_sse2.nb030_outerstart
        jmp _nb_kernel030_x86_64_sse2.nb030_end

_nb_kernel030_x86_64_sse2.nb030_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030_nouter(%rsp),%ebx
        movl %ebx,nb030_nouter(%rsp)

_nb_kernel030_x86_64_sse2.nb030_outer: 
        movq  nb030_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb030_is3(%rsp)      ## store is3 

        movq  nb030_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb030_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb030_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb030_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb030_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb030_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb030_ix(%rsp)
        movapd %xmm1,nb030_iy(%rsp)
        movapd %xmm2,nb030_iz(%rsp)

        movl  %ebx,nb030_ii3(%rsp)

        ## clear tot potential and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb030_Vvdwtot(%rsp)
        movapd %xmm4,nb030_fix(%rsp)
        movapd %xmm4,nb030_fiy(%rsp)
        movapd %xmm4,nb030_fiz(%rsp)

        movq  nb030_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb030_pos(%rbp),%rsi
        movq  nb030_faction(%rbp),%rdi
        movq  nb030_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb030_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030_ninner(%rsp),%ecx
        movl  %ecx,nb030_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb030_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel030_x86_64_sse2.nb030_unroll_loop
        jmp   _nb_kernel030_x86_64_sse2.nb030_checksingle
_nb_kernel030_x86_64_sse2.nb030_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb030_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        addq $8,nb030_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb030_pos(%rbp),%rsi               ## base of pos[] 
        lea  (%r8,%r8,2),%rax  ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        ## move two coordinates to xmm4-xmm6    
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb030_ix(%rsp),%xmm4
        subpd nb030_iy(%rsp),%xmm5
        subpd nb030_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb030_dx(%rsp)
        movapd %xmm5,nb030_dy(%rsp)
        movapd %xmm6,nb030_dz(%rsp)

        movq nb030_type(%rbp),%rsi

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        movl (%rsi,%r8,4),%r8d
        movl (%rsi,%r9,4),%r9d

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb030_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        shll %r8d
        shll %r9d

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb030_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm0,%xmm2       ## xmm2=iter2 of rinv (new lu) 

        movl nb030_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

    mulpd %xmm2,%xmm4   ## xmm4=r 
        mulpd nb030_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1
    ## xmm1=eps 
    ## xmm2=rinv
    movapd %xmm4,%xmm3  ## eps
        pslld $3,%mm6           ## idx *= 8 

        movq nb030_VFtab(%rbp),%rsi
        movd %mm6,%r10d
        psrlq $32,%mm6
        movd %mm6,%r11d

    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
        movapd (%rsi,%r10,8),%xmm4          ## Y1d F1d  
        movapd (%rsi,%r11,8),%xmm12         ## Y2d F2d 
        movapd 32(%rsi,%r10,8),%xmm8        ## Y1r F1r  
        movapd 32(%rsi,%r11,8),%xmm13           ## Y2r F2r 
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
        ## tables ready, in xmm4-xmm7 and xmm8-xmm11
        movq nb030_vdwparam(%rbp),%rsi

    mulpd  %xmm1,%xmm7   ## Heps
    mulpd  %xmm1,%xmm11
        movlpd (%rsi,%r8,8),%xmm12      ## c6a
        movlpd (%rsi,%r9,8),%xmm0       ## c6b
    mulpd  %xmm1,%xmm6  ## Geps
    mulpd  %xmm1,%xmm10
    mulpd  %xmm1,%xmm7  ## Heps2
    mulpd  %xmm1,%xmm11

        movhpd 8(%rsi,%r8,8),%xmm12     ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm0      ## c6b c12b 

    addpd  %xmm6,%xmm5 ## F+Geps
    addpd  %xmm10,%xmm9
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm11,%xmm9
    addpd  %xmm7,%xmm7   ## 2*Heps2
    addpd  %xmm11,%xmm11
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm10,%xmm11

        movapd %xmm12,%xmm13
        unpcklpd %xmm0,%xmm12
        unpckhpd %xmm0,%xmm13

    addpd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm9,%xmm11
    mulpd  %xmm1,%xmm5 ## eps*Fp
    mulpd  %xmm1,%xmm9
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb030_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb030_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7

    movq nb030_faction(%rbp),%rdi
        ## the fj's - start by combining forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5

    mulpd  nb030_tsc(%rsp),%xmm7
    mulpd  %xmm2,%xmm7
    xorpd  %xmm9,%xmm9

    subpd  %xmm7,%xmm9
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5

    movapd nb030_fix(%rsp),%xmm12
    movapd nb030_fiy(%rsp),%xmm13
    movapd nb030_fiz(%rsp),%xmm14

    mulpd  nb030_dx(%rsp),%xmm9
    mulpd  nb030_dy(%rsp),%xmm10
    mulpd  nb030_dz(%rsp),%xmm11

    ## accumulate i forces
    addpd %xmm9,%xmm12
    addpd %xmm10,%xmm13
    addpd %xmm11,%xmm14
    movapd %xmm12,nb030_fix(%rsp)
    movapd %xmm13,nb030_fiy(%rsp)
    movapd %xmm14,nb030_fiz(%rsp)

        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb030_innerk(%rsp)
        jl    _nb_kernel030_x86_64_sse2.nb030_checksingle
        jmp   _nb_kernel030_x86_64_sse2.nb030_unroll_loop

_nb_kernel030_x86_64_sse2.nb030_checksingle:    
        movl  nb030_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030_x86_64_sse2.nb030_dosingle
        jmp    _nb_kernel030_x86_64_sse2.nb030_updateouterdata
_nb_kernel030_x86_64_sse2.nb030_dosingle: 
        movq  nb030_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb030_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb030_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb030_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb030_c6(%rsp)
        movapd %xmm6,nb030_c12(%rsp)

        movq nb030_pos(%rbp),%rsi               ## base of pos[] 
        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6    
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb030_ix(%rsp),%xmm4
        subsd nb030_iy(%rsp),%xmm5
        subsd nb030_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb030_dx(%rsp)
        movapd %xmm5,nb030_dy(%rsp)
        movapd %xmm6,nb030_dz(%rsp)

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
        movapd nb030_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb030_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm0,%xmm2       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm2,%xmm4       ## xmm4=r 
        mulsd nb030_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r10d   ## mm6 = lu idx 
        cvtsi2sd %r10d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1
    ## xmm1=eps 
    ## xmm2=rinv
    movapd %xmm4,%xmm3  ## eps
        shll   $3,%r10d         ## idx *= 8 

        movq nb030_VFtab(%rbp),%rsi

    ## indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movsd  (%rsi,%r10,8),%xmm4
    movsd  8(%rsi,%r10,8),%xmm5
    movsd  16(%rsi,%r10,8),%xmm6
    movsd  24(%rsi,%r10,8),%xmm7
    movsd  32(%rsi,%r10,8),%xmm8
    movsd  40(%rsi,%r10,8),%xmm9
    movsd  48(%rsi,%r10,8),%xmm10
    movsd  56(%rsi,%r10,8),%xmm11
        ## tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulsd  %xmm1,%xmm7   ## Heps
    mulsd  %xmm1,%xmm11
    mulsd  %xmm1,%xmm6  ## Geps
    mulsd  %xmm1,%xmm10
    mulsd  %xmm1,%xmm7  ## Heps2
    mulsd  %xmm1,%xmm11
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
    mulsd  %xmm1,%xmm5 ## eps*Fp
    mulsd  %xmm1,%xmm9
    movapd nb030_c6(%rsp),%xmm12
    movapd nb030_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb030_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb030_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7

    mulsd  nb030_tsc(%rsp),%xmm7
    mulsd  %xmm2,%xmm7
    xorpd  %xmm9,%xmm9

    subsd  %xmm7,%xmm9
    movapd %xmm9,%xmm10
    movapd %xmm9,%xmm11

    movapd nb030_fix(%rsp),%xmm12
    movapd nb030_fiy(%rsp),%xmm13
    movapd nb030_fiz(%rsp),%xmm14

    mulsd  nb030_dx(%rsp),%xmm9
    mulsd  nb030_dy(%rsp),%xmm10
    mulsd  nb030_dz(%rsp),%xmm11

    ## accumulate i forces
    addsd %xmm9,%xmm12
    addsd %xmm10,%xmm13
    addsd %xmm11,%xmm14
    movsd %xmm12,nb030_fix(%rsp)
    movsd %xmm13,nb030_fiy(%rsp)
    movsd %xmm14,nb030_fiz(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb030_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel030_x86_64_sse2.nb030_updateouterdata: 
        movl  nb030_ii3(%rsp),%ecx
        movq  nb030_faction(%rbp),%rdi
        movq  nb030_fshift(%rbp),%rsi
        movl  nb030_is3(%rsp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb030_fix(%rsp),%xmm0
        movapd nb030_fiy(%rsp),%xmm1
        movapd nb030_fiz(%rsp),%xmm2

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
        subsd %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb030_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb030_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb030_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb030_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb030_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030_x86_64_sse2.nb030_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030_n(%rsp)
        jmp _nb_kernel030_x86_64_sse2.nb030_outer
_nb_kernel030_x86_64_sse2.nb030_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030_x86_64_sse2.nb030_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030_x86_64_sse2.nb030_threadloop
_nb_kernel030_x86_64_sse2.nb030_end: 

        movl nb030_nouter(%rsp),%eax
        movl nb030_ninner(%rsp),%ebx
        movq nb030_outeriter(%rbp),%rcx
        movq nb030_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $376,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret







.globl nb_kernel030nf_x86_64_sse2
.globl _nb_kernel030nf_x86_64_sse2
nb_kernel030nf_x86_64_sse2:     
_nb_kernel030nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb030nf_fshift, 16
.set nb030nf_gid, 24
.set nb030nf_pos, 32
.set nb030nf_faction, 40
.set nb030nf_charge, 48
.set nb030nf_p_facel, 56
.set nb030nf_argkrf, 64
.set nb030nf_argcrf, 72
.set nb030nf_Vc, 80
.set nb030nf_type, 88
.set nb030nf_p_ntype, 96
.set nb030nf_vdwparam, 104
.set nb030nf_Vvdw, 112
.set nb030nf_p_tabscale, 120
.set nb030nf_VFtab, 128
.set nb030nf_invsqrta, 136
.set nb030nf_dvda, 144
.set nb030nf_p_gbtabscale, 152
.set nb030nf_GBtab, 160
.set nb030nf_p_nthreads, 168
.set nb030nf_count, 176
.set nb030nf_mtx, 184
.set nb030nf_outeriter, 192
.set nb030nf_inneriter, 200
.set nb030nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030nf_ix, 0
.set nb030nf_iy, 16
.set nb030nf_iz, 32
.set nb030nf_tsc, 48
.set nb030nf_c6, 64
.set nb030nf_c12, 80
.set nb030nf_Vvdwtot, 96
.set nb030nf_half, 112
.set nb030nf_three, 128
.set nb030nf_is3, 144
.set nb030nf_ii3, 148
.set nb030nf_nri, 152
.set nb030nf_iinr, 160
.set nb030nf_jindex, 168
.set nb030nf_jjnr, 176
.set nb030nf_shift, 184
.set nb030nf_shiftvec, 192
.set nb030nf_innerjjnr, 200
.set nb030nf_ntia, 208
.set nb030nf_innerk, 212
.set nb030nf_n, 216
.set nb030nf_nn1, 220
.set nb030nf_ntype, 224
.set nb030nf_nouter, 228
.set nb030nf_ninner, 232

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $248,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb030nf_nouter(%rsp)
        movl %eax,nb030nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb030nf_nri(%rsp)
        movq %rsi,nb030nf_iinr(%rsp)
        movq %rdx,nb030nf_jindex(%rsp)
        movq %rcx,nb030nf_jjnr(%rsp)
        movq %r8,nb030nf_shift(%rsp)
        movq %r9,nb030nf_shiftvec(%rsp)
        movq nb030nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb030nf_ntype(%rsp)

        movq nb030nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb030nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb030nf_half(%rsp)
        movl %ebx,nb030nf_half+4(%rsp)
        movsd nb030nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb030nf_half(%rsp)
        movapd %xmm3,nb030nf_three(%rsp)

_nb_kernel030nf_x86_64_sse2.nb030nf_threadloop: 
        movq  nb030nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel030nf_x86_64_sse2.nb030nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030nf_x86_64_sse2.nb030nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030nf_n(%rsp)
        movl %ebx,nb030nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030nf_x86_64_sse2.nb030nf_outerstart
        jmp _nb_kernel030nf_x86_64_sse2.nb030nf_end

_nb_kernel030nf_x86_64_sse2.nb030nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030nf_nouter(%rsp),%ebx
        movl %ebx,nb030nf_nouter(%rsp)

_nb_kernel030nf_x86_64_sse2.nb030nf_outer: 
        movq  nb030nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb030nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb030nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb030nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb030nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb030nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb030nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb030nf_ix(%rsp)
        movapd %xmm1,nb030nf_iy(%rsp)
        movapd %xmm2,nb030nf_iz(%rsp)

        movl  %ebx,nb030nf_ii3(%rsp)

        ## clear tot potential 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb030nf_Vvdwtot(%rsp)

        movq  nb030nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb030nf_pos(%rbp),%rsi
        movq  nb030nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb030nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb030nf_ninner(%rsp),%ecx
        movl  %ecx,nb030nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb030nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel030nf_x86_64_sse2.nb030nf_unroll_loop
        jmp   _nb_kernel030nf_x86_64_sse2.nb030nf_checksingle
_nb_kernel030nf_x86_64_sse2.nb030nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb030nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb030nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb030nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb030nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb030nf_ntia(%rsp),%edi
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

        movapd %xmm4,nb030nf_c6(%rsp)
        movapd %xmm6,nb030nf_c12(%rsp)

        movq nb030nf_pos(%rbp),%rsi             ## base of pos[] 
        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movapd nb030nf_ix(%rsp),%xmm4
        movapd nb030nf_iy(%rsp),%xmm5
        movapd nb030nf_iz(%rsp),%xmm6

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
        movapd nb030nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb030nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb030nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movd %eax,%mm0
        movd %ebx,%mm1

        movq nb030nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
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
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb030nf_c6(%rsp),%xmm5   ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb030nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb030nf_Vvdwtot(%rsp)

        ## repulsion 
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

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb030nf_c12(%rsp),%xmm5

        addpd  nb030nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb030nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb030nf_innerk(%rsp)
        jl    _nb_kernel030nf_x86_64_sse2.nb030nf_checksingle
        jmp   _nb_kernel030nf_x86_64_sse2.nb030nf_unroll_loop

_nb_kernel030nf_x86_64_sse2.nb030nf_checksingle: 
        movl  nb030nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030nf_x86_64_sse2.nb030nf_dosingle
        jmp    _nb_kernel030nf_x86_64_sse2.nb030nf_updateouterdata
_nb_kernel030nf_x86_64_sse2.nb030nf_dosingle: 
        movq  nb030nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movq nb030nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb030nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb030nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb030nf_c6(%rsp)
        movapd %xmm6,nb030nf_c12(%rsp)

        movq nb030nf_pos(%rbp),%rsi             ## base of pos[] 
        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movapd nb030nf_ix(%rsp),%xmm4
        movapd nb030nf_iy(%rsp),%xmm5
        movapd nb030nf_iz(%rsp),%xmm6

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
        movapd nb030nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb030nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb030nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb030nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb030nf_tsc(%rsp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 
        shll $3,%eax

        movq nb030nf_VFtab(%rbp),%rsi

        ## dispersion 
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
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb030nf_c6(%rsp),%xmm5  ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb030nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb030nf_Vvdwtot(%rsp)

        ## repulsion 
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

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb030nf_c12(%rsp),%xmm5

        addsd  nb030nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb030nf_Vvdwtot(%rsp)

_nb_kernel030nf_x86_64_sse2.nb030nf_updateouterdata: 
        ## get n from stack
        movl nb030nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb030nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb030nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb030nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb030nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030nf_x86_64_sse2.nb030nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030nf_n(%rsp)
        jmp _nb_kernel030nf_x86_64_sse2.nb030nf_outer
_nb_kernel030nf_x86_64_sse2.nb030nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030nf_x86_64_sse2.nb030nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030nf_x86_64_sse2.nb030nf_threadloop
_nb_kernel030nf_x86_64_sse2.nb030nf_end: 
        movl nb030nf_nouter(%rsp),%eax
        movl nb030nf_ninner(%rsp),%ebx
        movq nb030nf_outeriter(%rbp),%rcx
        movq nb030nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $248,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

