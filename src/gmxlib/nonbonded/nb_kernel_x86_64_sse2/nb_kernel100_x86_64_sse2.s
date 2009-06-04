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







.globl nb_kernel100_x86_64_sse2
.globl _nb_kernel100_x86_64_sse2
nb_kernel100_x86_64_sse2:       
_nb_kernel100_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb100_fshift, 16
.set nb100_gid, 24
.set nb100_pos, 32
.set nb100_faction, 40
.set nb100_charge, 48
.set nb100_p_facel, 56
.set nb100_argkrf, 64
.set nb100_argcrf, 72
.set nb100_Vc, 80
.set nb100_type, 88
.set nb100_p_ntype, 96
.set nb100_vdwparam, 104
.set nb100_Vvdw, 112
.set nb100_p_tabscale, 120
.set nb100_VFtab, 128
.set nb100_invsqrta, 136
.set nb100_dvda, 144
.set nb100_p_gbtabscale, 152
.set nb100_GBtab, 160
.set nb100_p_nthreads, 168
.set nb100_count, 176
.set nb100_mtx, 184
.set nb100_outeriter, 192
.set nb100_inneriter, 200
.set nb100_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb100_ix, 0
.set nb100_iy, 16
.set nb100_iz, 32
.set nb100_iq, 48
.set nb100_dx, 64
.set nb100_dy, 80
.set nb100_dz, 96
.set nb100_vctot, 112
.set nb100_fix, 128
.set nb100_fiy, 144
.set nb100_fiz, 160
.set nb100_half, 176
.set nb100_three, 192
.set nb100_is3, 208
.set nb100_ii3, 212
.set nb100_nri, 216
.set nb100_iinr, 224
.set nb100_jindex, 232
.set nb100_jjnr, 240
.set nb100_shift, 248
.set nb100_shiftvec, 256
.set nb100_facel, 264
.set nb100_innerjjnr, 272
.set nb100_innerk, 280
.set nb100_n, 284
.set nb100_nn1, 288
.set nb100_nouter, 292
.set nb100_ninner, 296
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
        movl %eax,nb100_nouter(%rsp)
        movl %eax,nb100_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb100_nri(%rsp)
        movq %rsi,nb100_iinr(%rsp)
        movq %rdx,nb100_jindex(%rsp)
        movq %rcx,nb100_jjnr(%rsp)
        movq %r8,nb100_shift(%rsp)
        movq %r9,nb100_shiftvec(%rsp)
        movq nb100_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb100_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb100_half(%rsp)
        movsd nb100_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movapd %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movapd %xmm1,nb100_half(%rsp)
        movapd %xmm3,nb100_three(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb100_half(%rsp)
        movl %ebx,nb100_half+4(%rsp)
        movsd nb100_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb100_half(%rsp)
        movapd %xmm3,nb100_three(%rsp)

_nb_kernel100_x86_64_sse2.nb100_threadloop: 
        movq  nb100_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel100_x86_64_sse2.nb100_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100_x86_64_sse2.nb100_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100_n(%rsp)
        movl %ebx,nb100_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100_x86_64_sse2.nb100_outerstart
        jmp _nb_kernel100_x86_64_sse2.nb100_end

_nb_kernel100_x86_64_sse2.nb100_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100_nouter(%rsp),%ebx
        movl %ebx,nb100_nouter(%rsp)

_nb_kernel100_x86_64_sse2.nb100_outer: 
        movq  nb100_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb100_is3(%rsp)      ## store is3 

        movq  nb100_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb100_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb100_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb100_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3
        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb100_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb100_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb100_ix(%rsp)
        movapd %xmm1,nb100_iy(%rsp)
        movapd %xmm2,nb100_iz(%rsp)

        movl  %ebx,nb100_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorpd %xmm12,%xmm12
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15

        movq  nb100_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb100_pos(%rbp),%rsi
        movq  nb100_faction(%rbp),%rdi
        movq  nb100_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb100_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100_ninner(%rsp),%ecx
        movl  %ecx,nb100_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb100_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel100_x86_64_sse2.nb100_unroll_loop
        jmp   _nb_kernel100_x86_64_sse2.nb100_checksingle
_nb_kernel100_x86_64_sse2.nb100_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb100_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        addq $8,nb100_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb100_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        movq   nb100_faction(%rbp),%rdi

        ## calc dr 
        subpd nb100_ix(%rsp),%xmm4
        subpd nb100_iy(%rsp),%xmm5
        subpd nb100_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb100_charge(%rbp),%rsi     ## base of charge[] 

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

        movlpd (%rsi,%r8,8),%xmm3       ## jq A 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb100_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movhpd (%rsi,%r9,8),%xmm3       ## jq B 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb100_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100_half(%rsp),%xmm0
        mulpd nb100_iq(%rsp),%xmm3              ## qq 
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 

        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  %xmm3,%xmm4      ## xmm4=fscal 

        ## the fj's - start by combining forces from memory 
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2

    ## increment vctot
        addpd  %xmm3,%xmm12

        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

        addpd %xmm9,%xmm0
        addpd %xmm10,%xmm1
        addpd %xmm11,%xmm2

        ## xmm9-xmm11 contains tx-tz (partial force) 
        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb100_innerk(%rsp)
        jl    _nb_kernel100_x86_64_sse2.nb100_checksingle
        jmp   _nb_kernel100_x86_64_sse2.nb100_unroll_loop

_nb_kernel100_x86_64_sse2.nb100_checksingle:    
        movl  nb100_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100_x86_64_sse2.nb100_dosingle
        jmp    _nb_kernel100_x86_64_sse2.nb100_updateouterdata
_nb_kernel100_x86_64_sse2.nb100_dosingle: 
        movq nb100_charge(%rbp),%rsi
        movq nb100_pos(%rbp),%rdi

        movq nb100_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl (%rdx),%eax

        movq nb100_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3       ## jq A 
        mulsd nb100_iq(%rsp),%xmm3              ## qq 

        movq nb100_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        movq   nb100_faction(%rbp),%rdi

        ## calc dr 
        subsd nb100_ix(%rsp),%xmm4
        subsd nb100_iy(%rsp),%xmm5
        subsd nb100_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

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
        movapd nb100_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb100_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 

        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  %xmm3,%xmm4      ## xmm4=fscal 

    ## increment vctot
        addsd  %xmm3,%xmm12

        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11
        ## xmm9-xmm11 contains tx-tz (partial force) 
        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15

        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel100_x86_64_sse2.nb100_updateouterdata: 
        movl  nb100_ii3(%rsp),%ecx
        movq  nb100_faction(%rbp),%rdi
        movq  nb100_fshift(%rbp),%rsi
        movl  nb100_is3(%rsp),%edx

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
        movl nb100_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb100_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total coulomb energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 have the sum now 

        ## add earlier value from mem 
        movq  nb100_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## finish if last 
        movl nb100_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100_x86_64_sse2.nb100_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100_n(%rsp)
        jmp _nb_kernel100_x86_64_sse2.nb100_outer
_nb_kernel100_x86_64_sse2.nb100_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100_x86_64_sse2.nb100_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100_x86_64_sse2.nb100_threadloop
_nb_kernel100_x86_64_sse2.nb100_end: 
        movl nb100_nouter(%rsp),%eax
        movl nb100_ninner(%rsp),%ebx
        movq nb100_outeriter(%rbp),%rcx
        movq nb100_inneriter(%rbp),%rdx
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




.globl nb_kernel100nf_x86_64_sse2
.globl _nb_kernel100nf_x86_64_sse2
nb_kernel100nf_x86_64_sse2:     
_nb_kernel100nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb100nf_fshift, 16
.set nb100nf_gid, 24
.set nb100nf_pos, 32
.set nb100nf_faction, 40
.set nb100nf_charge, 48
.set nb100nf_p_facel, 56
.set nb100nf_argkrf, 64
.set nb100nf_argcrf, 72
.set nb100nf_Vc, 80
.set nb100nf_type, 88
.set nb100nf_p_ntype, 96
.set nb100nf_vdwparam, 104
.set nb100nf_Vvdw, 112
.set nb100nf_p_tabscale, 120
.set nb100nf_VFtab, 128
.set nb100nf_invsqrta, 136
.set nb100nf_dvda, 144
.set nb100nf_p_gbtabscale, 152
.set nb100nf_GBtab, 160
.set nb100nf_p_nthreads, 168
.set nb100nf_count, 176
.set nb100nf_mtx, 184
.set nb100nf_outeriter, 192
.set nb100nf_inneriter, 200
.set nb100nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100nf_ix, 0
.set nb100nf_iy, 16
.set nb100nf_iz, 32
.set nb100nf_iq, 48
.set nb100nf_vctot, 64
.set nb100nf_half, 80
.set nb100nf_three, 96
.set nb100nf_is3, 112
.set nb100nf_ii3, 116
.set nb100nf_nri, 120
.set nb100nf_iinr, 128
.set nb100nf_jindex, 136
.set nb100nf_jjnr, 144
.set nb100nf_shift, 156
.set nb100nf_shiftvec, 164
.set nb100nf_facel, 172
.set nb100nf_innerjjnr, 180
.set nb100nf_innerk, 188
.set nb100nf_n, 192
.set nb100nf_nn1, 196
.set nb100nf_nouter, 200
.set nb100nf_ninner, 204
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $216,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb100nf_nouter(%rsp)
        movl %eax,nb100nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb100nf_nri(%rsp)
        movq %rsi,nb100nf_iinr(%rsp)
        movq %rdx,nb100nf_jindex(%rsp)
        movq %rcx,nb100nf_jjnr(%rsp)
        movq %r8,nb100nf_shift(%rsp)
        movq %r9,nb100nf_shiftvec(%rsp)
        movq nb100nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb100nf_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb100nf_half(%rsp)
        movl %ebx,nb100nf_half+4(%rsp)
        movsd nb100nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb100nf_half(%rsp)
        movapd %xmm3,nb100nf_three(%rsp)

_nb_kernel100nf_x86_64_sse2.nb100nf_threadloop: 
        movq  nb100nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel100nf_x86_64_sse2.nb100nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100nf_x86_64_sse2.nb100nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100nf_n(%rsp)
        movl %ebx,nb100nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100nf_x86_64_sse2.nb100nf_outerstart
        jmp _nb_kernel100nf_x86_64_sse2.nb100nf_end

_nb_kernel100nf_x86_64_sse2.nb100nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100nf_nouter(%rsp),%ebx
        movl %ebx,nb100nf_nouter(%rsp)

_nb_kernel100nf_x86_64_sse2.nb100nf_outer: 
        movq  nb100nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb100nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb100nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb100nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb100nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3
        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb100nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb100nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb100nf_ix(%rsp)
        movapd %xmm1,nb100nf_iy(%rsp)
        movapd %xmm2,nb100nf_iz(%rsp)

        movl  %ebx,nb100nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb100nf_vctot(%rsp)

        movq  nb100nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb100nf_pos(%rbp),%rsi
        movq  nb100nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb100nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb100nf_ninner(%rsp),%ecx
        movl  %ecx,nb100nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb100nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel100nf_x86_64_sse2.nb100nf_unroll_loop
        jmp   _nb_kernel100nf_x86_64_sse2.nb100nf_checksingle
_nb_kernel100nf_x86_64_sse2.nb100nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb100nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb100nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb100nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3      ## jq A 
        movhpd (%rsi,%rbx,8),%xmm3      ## jq B 

        movapd nb100nf_iq(%rsp),%xmm5

        mulpd %xmm5,%xmm3               ## qq 

        movq nb100nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movapd nb100nf_ix(%rsp),%xmm4
        movapd nb100nf_iy(%rsp),%xmm5
        movapd nb100nf_iz(%rsp),%xmm6

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
        movapd nb100nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb100nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        movapd nb100nf_vctot(%rsp),%xmm5
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        addpd  %xmm3,%xmm5
        movapd %xmm5,nb100nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb100nf_innerk(%rsp)
        jl    _nb_kernel100nf_x86_64_sse2.nb100nf_checksingle
        jmp   _nb_kernel100nf_x86_64_sse2.nb100nf_unroll_loop

_nb_kernel100nf_x86_64_sse2.nb100nf_checksingle: 
        movl  nb100nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100nf_x86_64_sse2.nb100nf_dosingle
        jmp    _nb_kernel100nf_x86_64_sse2.nb100nf_updateouterdata
_nb_kernel100nf_x86_64_sse2.nb100nf_dosingle: 
        movq nb100nf_charge(%rbp),%rsi
        movq nb100nf_pos(%rbp),%rdi

        movq nb100nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl (%rdx),%eax

        xorpd %xmm3,%xmm3
        movsd (%rsi,%rax,8),%xmm3       ## jq A 
        movapd nb100nf_iq(%rsp),%xmm5
        unpcklpd %xmm6,%xmm3
        mulpd %xmm5,%xmm3               ## qq 

        movq nb100nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movapd nb100nf_ix(%rsp),%xmm4
        movapd nb100nf_iy(%rsp),%xmm5
        movapd nb100nf_iz(%rsp),%xmm6

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
        movapd nb100nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb100nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb100nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb100nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        movlpd nb100nf_vctot(%rsp),%xmm5
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        addsd  %xmm3,%xmm5
        movlpd %xmm5,nb100nf_vctot(%rsp)

_nb_kernel100nf_x86_64_sse2.nb100nf_updateouterdata: 
        ## get n from stack
        movl nb100nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb100nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb100nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb100nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb100nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100nf_x86_64_sse2.nb100nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100nf_n(%rsp)
        jmp _nb_kernel100nf_x86_64_sse2.nb100nf_outer
_nb_kernel100nf_x86_64_sse2.nb100nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100nf_x86_64_sse2.nb100nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100nf_x86_64_sse2.nb100nf_threadloop
_nb_kernel100nf_x86_64_sse2.nb100nf_end: 
        movl nb100nf_nouter(%rsp),%eax
        movl nb100nf_ninner(%rsp),%ebx
        movq nb100nf_outeriter(%rbp),%rcx
        movq nb100nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $216,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





