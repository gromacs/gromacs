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





.globl nb_kernel010_x86_64_sse2
.globl _nb_kernel010_x86_64_sse2
nb_kernel010_x86_64_sse2:       
_nb_kernel010_x86_64_sse2:      
.set nb010_fshift, 16
.set nb010_gid, 24
.set nb010_pos, 32
.set nb010_faction, 40
.set nb010_charge, 48
.set nb010_p_facel, 56
.set nb010_argkrf, 64
.set nb010_argcrf, 72
.set nb010_Vc, 80
.set nb010_type, 88
.set nb010_p_ntype, 96
.set nb010_vdwparam, 104
.set nb010_Vvdw, 112
.set nb010_p_tabscale, 120
.set nb010_VFtab, 128
.set nb010_invsqrta, 136
.set nb010_dvda, 144
.set nb010_p_gbtabscale, 152
.set nb010_GBtab, 160
.set nb010_p_nthreads, 168
.set nb010_count, 176
.set nb010_mtx, 184
.set nb010_outeriter, 192
.set nb010_inneriter, 200
.set nb010_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010_ix, 0
.set nb010_iy, 16
.set nb010_iz, 32
.set nb010_dx, 48
.set nb010_dy, 64
.set nb010_dz, 80
.set nb010_two, 96
.set nb010_c6, 112
.set nb010_c12, 128
.set nb010_six, 144
.set nb010_twelve, 160
.set nb010_Vvdwtot, 176
.set nb010_fix, 192
.set nb010_fiy, 208
.set nb010_fiz, 224
.set nb010_half, 240
.set nb010_three, 256
.set nb010_nri, 272
.set nb010_iinr, 280
.set nb010_jindex, 288
.set nb010_jjnr, 296
.set nb010_shift, 304
.set nb010_shiftvec, 312
.set nb010_facel, 320
.set nb010_innerjjnr, 328
.set nb010_is3, 336
.set nb010_ii3, 340
.set nb010_ntia, 344
.set nb010_innerk, 348
.set nb010_n, 352
.set nb010_nn1, 356
.set nb010_ntype, 360
.set nb010_nouter, 364
.set nb010_ninner, 368
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $392,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb010_nouter(%rsp)
        movl %eax,nb010_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb010_nri(%rsp)
        movq %rsi,nb010_iinr(%rsp)
        movq %rdx,nb010_jindex(%rsp)
        movq %rcx,nb010_jjnr(%rsp)
        movq %r8,nb010_shift(%rsp)
        movq %r9,nb010_shiftvec(%rsp)
        movq nb010_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb010_ntype(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double two IEEE (hex)
        movl $0x40000000,%ebx
        movl %eax,nb010_two(%rsp)
        movl %ebx,nb010_two+4(%rsp)
        movsd nb010_two(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm2
        addpd  %xmm1,%xmm2      ## 4.0
        addpd  %xmm1,%xmm2      ## six
        movapd %xmm2,%xmm3
        addpd  %xmm3,%xmm3      ## twelve
        movapd %xmm1,nb010_two(%rsp)
        movapd %xmm2,nb010_six(%rsp)
        movapd %xmm3,nb010_twelve(%rsp)

_nb_kernel010_x86_64_sse2.nb010_threadloop: 
        movq  nb010_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel010_x86_64_sse2.nb010_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010_x86_64_sse2.nb010_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010_n(%rsp)
        movl %ebx,nb010_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010_x86_64_sse2.nb010_outerstart
        jmp _nb_kernel010_x86_64_sse2.nb010_end

_nb_kernel010_x86_64_sse2.nb010_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010_nouter(%rsp),%ebx
        movl %ebx,nb010_nouter(%rsp)

_nb_kernel010_x86_64_sse2.nb010_outer: 
        movq  nb010_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb010_is3(%rsp)      ## store is3 

        movq  nb010_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb010_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb010_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb010_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb010_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb010_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb010_ix(%rsp)
        movapd %xmm1,nb010_iy(%rsp)
        movapd %xmm2,nb010_iz(%rsp)

        movl  %ebx,nb010_ii3(%rsp)

        ## clear vvdwtot (xmm12) and i forces (xmm13-xmm15)
        xorpd %xmm12,%xmm12
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15

        movq  nb010_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb010_pos(%rbp),%rsi
        movq  nb010_faction(%rbp),%rdi
        movq  nb010_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb010_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010_ninner(%rsp),%ecx
        movl  %ecx,nb010_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb010_innerk(%rsp)      ## number of innerloop atoms 

        jge   _nb_kernel010_x86_64_sse2.nb010_unroll_loop
        jmp   _nb_kernel010_x86_64_sse2.nb010_checksingle
_nb_kernel010_x86_64_sse2.nb010_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb010_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        addq  $8,nb010_innerjjnr(%rsp)              ## advance pointer (unrolled 2) 

        movq nb010_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        ## move two coordinates to xmm0-xmm2    

        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb010_ix(%rsp),%xmm4
        subpd nb010_iy(%rsp),%xmm5
        subpd nb010_iz(%rsp),%xmm6

        movq nb010_type(%rbp),%rsi

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movl (%rsi,%r8,4),%r8d
        movl (%rsi,%r9,4),%r9d

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4       ## rsq in xmm4 

        movq nb010_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d

        cvtpd2ps %xmm4,%xmm6
        rcpps %xmm6,%xmm6
        cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

        movl nb010_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        ## 1/x lookup seed in xmm6 
        movapd nb010_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulpd %xmm6,%xmm4       ## lu*rsq 
        subpd %xmm4,%xmm0       ## 2-lu*rsq 
        mulpd %xmm0,%xmm6       ## (new lu) 

        movlpd (%rsi,%r8,8),%xmm8       ## c6a
        movlpd 8(%rsi,%r8,8),%xmm4      ## c612b

        movapd nb010_two(%rsp),%xmm0
        mulpd %xmm6,%xmm5       ## lu*rsq 
        subpd %xmm5,%xmm0       ## 2-lu*rsq 
        mulpd %xmm6,%xmm0       ## xmm0=rinvsq 

        movhpd (%rsi,%r9,8),%xmm8       ## c6a c12a 
        movhpd 8(%rsi,%r9,8),%xmm4      ## c6b c12b 

        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

    ## c6 in xmm7, c12 in xmm8

        mulpd  %xmm8,%xmm1  ## mult by c6
        mulpd  %xmm4,%xmm2  ## mult by c12
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulpd  nb010_six(%rsp),%xmm1
        mulpd  nb010_twelve(%rsp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm2,%xmm0      ## xmm4=total fscal 
        movapd %xmm0,%xmm4

        ## the fj's - start by accumulating forces from memory 
        movq   nb010_faction(%rbp),%rdi
    movlpd (%rdi,%rax,8),%xmm6
        movlpd 8(%rdi,%rax,8),%xmm7
        movlpd 16(%rdi,%rax,8),%xmm8
        movhpd (%rdi,%rbx,8),%xmm6
        movhpd 8(%rdi,%rbx,8),%xmm7
        movhpd 16(%rdi,%rbx,8),%xmm8

    ## increment potential
    addpd  %xmm5,%xmm12

        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        addpd %xmm9,%xmm6
        addpd %xmm10,%xmm7
        addpd %xmm11,%xmm8
        movlpd %xmm6,(%rdi,%rax,8)
        movlpd %xmm7,8(%rdi,%rax,8)
        movlpd %xmm8,16(%rdi,%rax,8)
        movhpd %xmm6,(%rdi,%rbx,8)
        movhpd %xmm7,8(%rdi,%rbx,8)
        movhpd %xmm8,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl  $2,nb010_innerk(%rsp)
        jl    _nb_kernel010_x86_64_sse2.nb010_checksingle
        jmp   _nb_kernel010_x86_64_sse2.nb010_unroll_loop
_nb_kernel010_x86_64_sse2.nb010_checksingle:    
        movl  nb010_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010_x86_64_sse2.nb010_dosingle
        jmp    _nb_kernel010_x86_64_sse2.nb010_updateouterdata
_nb_kernel010_x86_64_sse2.nb010_dosingle: 
        movq nb010_pos(%rbp),%rdi
        movq  nb010_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb010_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb010_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb010_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm7            ## c6
        movsd 8(%rsi,%r8,8),%xmm8       ## c12  
    ## c6 in xmm7, c12 in xmm8

        movq nb010_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb010_ix(%rsp),%xmm4
        subsd nb010_iy(%rsp),%xmm5
        subsd nb010_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4       ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm6
        rcpss %xmm6,%xmm6
        cvtss2sd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulsd %xmm6,%xmm4       ## lu*rsq 
        subsd %xmm4,%xmm0       ## 2-lu*rsq 
        mulsd %xmm0,%xmm6       ## (new lu) 

        movapd nb010_two(%rsp),%xmm0
        mulsd %xmm6,%xmm5       ## lu*rsq 
        subsd %xmm5,%xmm0       ## 2-lu*rsq 
        mulsd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulsd  %xmm7,%xmm1  ## mult by c6
        mulsd  %xmm8,%xmm2  ## mult by c12
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulsd  nb010_six(%rsp),%xmm1
        mulsd  nb010_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm2,%xmm0      ## xmm4=total fscal 
        movapd %xmm0,%xmm4

    ## increment potential
    addsd  %xmm5,%xmm12

        movq   nb010_faction(%rbp),%rdi
        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15

        ## the fj's - start by accumulating forces from memory 
    movsd (%rdi,%rax,8),%xmm3
        movsd 8(%rdi,%rax,8),%xmm4
        movsd 16(%rdi,%rax,8),%xmm5
        addsd %xmm9,%xmm3
        addsd %xmm10,%xmm4
        addsd %xmm11,%xmm5
        movsd %xmm3,(%rdi,%rax,8)
        movsd %xmm4,8(%rdi,%rax,8)
        movsd %xmm5,16(%rdi,%rax,8)

_nb_kernel010_x86_64_sse2.nb010_updateouterdata: 
        movl  nb010_ii3(%rsp),%ecx
        movq  nb010_faction(%rbp),%rdi
        movq  nb010_fshift(%rbp),%rsi
        movl  nb010_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm3
        movhlps %xmm14,%xmm4
        movhlps %xmm15,%xmm5
        addsd  %xmm3,%xmm13
        addsd  %xmm4,%xmm14
        addsd  %xmm5,%xmm15 ## sum is in low xmm13-xmm15

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb010_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb010_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 have the sum now 

        ## add earlier value from mem 
        movq  nb010_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## finish if last 
        movl nb010_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010_x86_64_sse2.nb010_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010_n(%rsp)
        jmp _nb_kernel010_x86_64_sse2.nb010_outer
_nb_kernel010_x86_64_sse2.nb010_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010_x86_64_sse2.nb010_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010_x86_64_sse2.nb010_threadloop
_nb_kernel010_x86_64_sse2.nb010_end: 

        movl nb010_nouter(%rsp),%eax
        movl nb010_ninner(%rsp),%ebx
        movq nb010_outeriter(%rbp),%rcx
        movq nb010_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $392,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret







.globl nb_kernel010nf_x86_64_sse2
.globl _nb_kernel010nf_x86_64_sse2
nb_kernel010nf_x86_64_sse2:     
_nb_kernel010nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb010nf_fshift, 16
.set nb010nf_gid, 24
.set nb010nf_pos, 32
.set nb010nf_faction, 40
.set nb010nf_charge, 48
.set nb010nf_p_facel, 56
.set nb010nf_argkrf, 64
.set nb010nf_argcrf, 72
.set nb010nf_Vc, 80
.set nb010nf_type, 88
.set nb010nf_p_ntype, 96
.set nb010nf_vdwparam, 104
.set nb010nf_Vvdw, 112
.set nb010nf_p_tabscale, 120
.set nb010nf_VFtab, 128
.set nb010nf_invsqrta, 136
.set nb010nf_dvda, 144
.set nb010nf_p_gbtabscale, 152
.set nb010nf_GBtab, 160
.set nb010nf_p_nthreads, 168
.set nb010nf_count, 176
.set nb010nf_mtx, 184
.set nb010nf_outeriter, 192
.set nb010nf_inneriter, 200
.set nb010nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010nf_ix, 0
.set nb010nf_iy, 16
.set nb010nf_iz, 32
.set nb010nf_two, 48
.set nb010nf_c6, 64
.set nb010nf_c12, 80
.set nb010nf_Vvdwtot, 96
.set nb010nf_half, 112
.set nb010nf_three, 128
.set nb010nf_nri, 144
.set nb010nf_iinr, 152
.set nb010nf_jindex, 160
.set nb010nf_jjnr, 168
.set nb010nf_shift, 176
.set nb010nf_shiftvec, 184
.set nb010nf_innerjjnr, 192
.set nb010nf_is3, 200
.set nb010nf_ii3, 204
.set nb010nf_ntia, 208
.set nb010nf_innerk, 212
.set nb010nf_n, 216
.set nb010nf_nn1, 220
.set nb010nf_ntype, 224
.set nb010nf_nouter, 228
.set nb010nf_ninner, 232
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
        movl %eax,nb010nf_nouter(%rsp)
        movl %eax,nb010nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb010nf_nri(%rsp)
        movq %rsi,nb010nf_iinr(%rsp)
        movq %rdx,nb010nf_jindex(%rsp)
        movq %rcx,nb010nf_jjnr(%rsp)
        movq %r8,nb010nf_shift(%rsp)
        movq %r9,nb010nf_shiftvec(%rsp)
        movq nb010nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb010nf_ntype(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double two IEEE (hex)
        movl $0x40000000,%ebx
        movl %eax,nb010nf_two(%rsp)
        movl %ebx,nb010nf_two+4(%rsp)
        movsd nb010nf_two(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,nb010nf_two(%rsp)

_nb_kernel010nf_x86_64_sse2.nb010nf_threadloop: 
        movq  nb010nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel010nf_x86_64_sse2.nb010nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010nf_x86_64_sse2.nb010nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010nf_n(%rsp)
        movl %ebx,nb010nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010nf_x86_64_sse2.nb010nf_outerstart
        jmp _nb_kernel010nf_x86_64_sse2.nb010nf_end

_nb_kernel010nf_x86_64_sse2.nb010nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010nf_nouter(%rsp),%ebx
        movl %ebx,nb010nf_nouter(%rsp)

_nb_kernel010nf_x86_64_sse2.nb010nf_outer: 
        movq  nb010nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb010nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb010nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb010nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb010nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb010nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb010nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb010nf_ix(%rsp)
        movapd %xmm1,nb010nf_iy(%rsp)
        movapd %xmm2,nb010nf_iz(%rsp)

        movl  %ebx,nb010nf_ii3(%rsp)

        ## clear Vvdwtot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb010nf_Vvdwtot(%rsp)

        movq  nb010nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb010nf_pos(%rbp),%rsi
        movq  nb010nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb010nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010nf_ninner(%rsp),%ecx
        movl  %ecx,nb010nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb010nf_innerk(%rsp)      ## number of innerloop atoms 

        jge   _nb_kernel010nf_x86_64_sse2.nb010nf_unroll_loop
        jmp   _nb_kernel010nf_x86_64_sse2.nb010nf_checksingle
_nb_kernel010nf_x86_64_sse2.nb010nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb010nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq  $8,nb010nf_innerjjnr(%rsp)              ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb010nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb010nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb010nf_ntia(%rsp),%edi
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

        movapd %xmm4,nb010nf_c6(%rsp)
        movapd %xmm6,nb010nf_c12(%rsp)

        movq nb010nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    

        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010nf_ix(%rsp),%xmm4
        movapd nb010nf_iy(%rsp),%xmm5
        movapd nb010nf_iz(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4       ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm6
        rcpps %xmm6,%xmm6
        cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010nf_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulpd %xmm6,%xmm4       ## lu*rsq 
        subpd %xmm4,%xmm0       ## 2-lu*rsq 
        mulpd %xmm0,%xmm6       ## (new lu) 

        movapd nb010nf_two(%rsp),%xmm0
        mulpd %xmm6,%xmm5       ## lu*rsq 
        subpd %xmm5,%xmm0       ## 2-lu*rsq 
        mulpd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulpd  nb010nf_c6(%rsp),%xmm1
        mulpd  nb010nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb010nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb010nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl  $2,nb010nf_innerk(%rsp)
        jl    _nb_kernel010nf_x86_64_sse2.nb010nf_checksingle
        jmp   _nb_kernel010nf_x86_64_sse2.nb010nf_unroll_loop
_nb_kernel010nf_x86_64_sse2.nb010nf_checksingle: 
        movl  nb010nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010nf_x86_64_sse2.nb010nf_dosingle
        jmp    _nb_kernel010nf_x86_64_sse2.nb010nf_updateouterdata
_nb_kernel010nf_x86_64_sse2.nb010nf_dosingle: 
        movq nb010nf_pos(%rbp),%rdi
        movq  nb010nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage    
        movq nb010nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb010nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb010nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb010nf_c6(%rsp)
        movapd %xmm6,nb010nf_c12(%rsp)

        movq nb010nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        

        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010nf_ix(%rsp),%xmm4
        movapd nb010nf_iy(%rsp),%xmm5
        movapd nb010nf_iz(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4       ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm6
        rcpss %xmm6,%xmm6
        cvtss2sd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010nf_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulsd %xmm6,%xmm4       ## lu*rsq 
        subsd %xmm4,%xmm0       ## 2-lu*rsq 
        mulsd %xmm0,%xmm6       ## (new lu) 

        movapd nb010nf_two(%rsp),%xmm0
        mulsd %xmm6,%xmm5       ## lu*rsq 
        subsd %xmm5,%xmm0       ## 2-lu*rsq 
        mulsd %xmm6,%xmm0       ## xmm0=rinvsq 
        movapd %xmm0,%xmm4

        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulsd  nb010nf_c6(%rsp),%xmm1
        mulsd  nb010nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb010nf_Vvdwtot(%rsp),%xmm5
        movlpd %xmm5,nb010nf_Vvdwtot(%rsp)

_nb_kernel010nf_x86_64_sse2.nb010nf_updateouterdata: 
        ## get n from stack
        movl nb010nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb010nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb010nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb010nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb010nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010nf_x86_64_sse2.nb010nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010nf_n(%rsp)
        jmp _nb_kernel010nf_x86_64_sse2.nb010nf_outer
_nb_kernel010nf_x86_64_sse2.nb010nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010nf_x86_64_sse2.nb010nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010nf_x86_64_sse2.nb010nf_threadloop
_nb_kernel010nf_x86_64_sse2.nb010nf_end: 
        movl nb010nf_nouter(%rsp),%eax
        movl nb010nf_ninner(%rsp),%ebx
        movq nb010nf_outeriter(%rbp),%rcx
        movq nb010nf_inneriter(%rbp),%rdx
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

