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



.globl nb_kernel200_ia32_sse
.globl _nb_kernel200_ia32_sse
nb_kernel200_ia32_sse:  
_nb_kernel200_ia32_sse: 
.set nb200_p_nri, 8
.set nb200_iinr, 12
.set nb200_jindex, 16
.set nb200_jjnr, 20
.set nb200_shift, 24
.set nb200_shiftvec, 28
.set nb200_fshift, 32
.set nb200_gid, 36
.set nb200_pos, 40
.set nb200_faction, 44
.set nb200_charge, 48
.set nb200_p_facel, 52
.set nb200_argkrf, 56
.set nb200_argcrf, 60
.set nb200_Vc, 64
.set nb200_type, 68
.set nb200_p_ntype, 72
.set nb200_vdwparam, 76
.set nb200_Vvdw, 80
.set nb200_p_tabscale, 84
.set nb200_VFtab, 88
.set nb200_invsqrta, 92
.set nb200_dvda, 96
.set nb200_p_gbtabscale, 100
.set nb200_GBtab, 104
.set nb200_p_nthreads, 108
.set nb200_count, 112
.set nb200_mtx, 116
.set nb200_outeriter, 120
.set nb200_inneriter, 124
.set nb200_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200_ix, 0
.set nb200_iy, 16
.set nb200_iz, 32
.set nb200_iq, 48
.set nb200_dx, 64
.set nb200_dy, 80
.set nb200_dz, 96
.set nb200_vctot, 112
.set nb200_fix, 128
.set nb200_fiy, 144
.set nb200_fiz, 160
.set nb200_half, 176
.set nb200_three, 192
.set nb200_two, 208
.set nb200_krf, 224
.set nb200_crf, 240
.set nb200_is3, 256
.set nb200_ii3, 260
.set nb200_innerjjnr, 264
.set nb200_innerk, 268
.set nb200_n, 272
.set nb200_nn1, 276
.set nb200_nri, 280
.set nb200_facel, 284
.set nb200_nouter, 288
.set nb200_ninner, 292
.set nb200_salign, 296
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $300,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb200_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb200_p_nri(%ebp),%ecx
        movl nb200_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb200_nri(%esp)
        movl %esi,nb200_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb200_nouter(%esp)
        movl %eax,nb200_ninner(%esp)


        movl nb200_argkrf(%ebp),%esi
        movl nb200_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb200_krf(%esp)
        shufps $0,%xmm6,%xmm6
        movaps %xmm6,nb200_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb200_half(%esp)
        movss nb200_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb200_half(%esp)
        movaps %xmm2,nb200_two(%esp)
        movaps %xmm3,nb200_three(%esp)

_nb_kernel200_ia32_sse.nb200_threadloop: 
        movl  nb200_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel200_ia32_sse.nb200_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200_ia32_sse.nb200_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200_n(%esp)
        movl %ebx,nb200_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200_ia32_sse.nb200_outerstart
        jmp _nb_kernel200_ia32_sse.nb200_end

        ## assume we have at least one i particle - start directly      
_nb_kernel200_ia32_sse.nb200_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200_nouter(%esp),%ebx
        movl %ebx,nb200_nouter(%esp)

_nb_kernel200_ia32_sse.nb200_outer: 
        movl  nb200_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb200_is3(%esp)      ## store is3 

        movl  nb200_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb200_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb200_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb200_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb200_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb200_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb200_ix(%esp)
        movaps %xmm1,nb200_iy(%esp)
        movaps %xmm2,nb200_iz(%esp)

        movl  %ebx,nb200_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb200_vctot(%esp)
        movaps %xmm4,nb200_fix(%esp)
        movaps %xmm4,nb200_fiy(%esp)
        movaps %xmm4,nb200_fiz(%esp)

        movl  nb200_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb200_pos(%ebp),%esi
        movl  nb200_faction(%ebp),%edi
        movl  nb200_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb200_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb200_ninner(%esp),%ecx
        movl  %ecx,nb200_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb200_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel200_ia32_sse.nb200_unroll_loop
        jmp   _nb_kernel200_ia32_sse.nb200_finish_inner
_nb_kernel200_ia32_sse.nb200_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb200_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb200_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb200_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb200_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  

        movl nb200_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb200_ix(%esp),%xmm4
        movaps nb200_iy(%esp),%xmm5
        movaps nb200_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb200_dx(%esp)
        movaps %xmm5,nb200_dy(%esp)
        movaps %xmm6,nb200_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subps  nb200_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulps  nb200_two(%esp),%xmm7

        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3
        mulps  %xmm3,%xmm4      ## xmm4=total fscal 
        addps  nb200_vctot(%esp),%xmm6

        movaps nb200_dx(%esp),%xmm0
        movaps nb200_dy(%esp),%xmm1
        movaps nb200_dz(%esp),%xmm2

        movaps %xmm6,nb200_vctot(%esp)

        movl   nb200_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb200_fix(%esp),%xmm3
        movaps nb200_fiy(%esp),%xmm4
        movaps nb200_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb200_fix(%esp)
        movaps %xmm4,nb200_fiy(%esp)
        movaps %xmm5,nb200_fiz(%esp)
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm6
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm6

        movaps %xmm4,%xmm3
        shufps $136,%xmm6,%xmm3 ## constant 10001000
        shufps $221,%xmm6,%xmm4 ## constant 11011101                          

        ## now xmm3-xmm5 contains fjx, fjy, fjz 
        subps  %xmm0,%xmm3
        subps  %xmm1,%xmm4

        ## unpack them back so we can store them - first x & y in xmm3/xmm4 

        movaps %xmm3,%xmm6
        unpcklps %xmm4,%xmm6
        unpckhps %xmm4,%xmm3
        ## xmm6(l)=x & y for j1, (h) for j2 
        ## xmm3(l)=x & y for j3, (h) for j4 
        movlps %xmm6,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)

        movhps %xmm6,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)

        ## and the z forces 
        movss  8(%edi,%eax,4),%xmm4
        movss  8(%edi,%ebx,4),%xmm5
        movss  8(%edi,%ecx,4),%xmm6
        movss  8(%edi,%edx,4),%xmm7
        subss  %xmm2,%xmm4
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm5
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm6
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm7
        movss  %xmm4,8(%edi,%eax,4)
        movss  %xmm5,8(%edi,%ebx,4)
        movss  %xmm6,8(%edi,%ecx,4)
        movss  %xmm7,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb200_innerk(%esp)
        jl    _nb_kernel200_ia32_sse.nb200_finish_inner
        jmp   _nb_kernel200_ia32_sse.nb200_unroll_loop
_nb_kernel200_ia32_sse.nb200_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb200_innerk(%esp)
        movl  nb200_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel200_ia32_sse.nb200_dopair
        jmp   _nb_kernel200_ia32_sse.nb200_checksingle
_nb_kernel200_ia32_sse.nb200_dopair: 
        movl nb200_charge(%ebp),%esi

        movl  nb200_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb200_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges         

        movl nb200_pos(%ebp),%edi

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb200_iq(%esp),%xmm3

        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb200_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 

        movaps nb200_ix(%esp),%xmm4
        movaps nb200_iy(%esp),%xmm5
        movaps nb200_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb200_dx(%esp)
        movaps %xmm5,nb200_dy(%esp)
        movaps %xmm6,nb200_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subps  nb200_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 
        mulps  nb200_two(%esp),%xmm7

        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3

        mulps  %xmm3,%xmm4      ## xmm4=total fscal 
        addps  nb200_vctot(%esp),%xmm6

        movaps nb200_dx(%esp),%xmm0
        movaps nb200_dy(%esp),%xmm1
        movaps nb200_dz(%esp),%xmm2

        movaps %xmm6,nb200_vctot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb200_fix(%esp),%xmm3
        movaps nb200_fiy(%esp),%xmm4
        movaps nb200_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb200_fix(%esp)
        movaps %xmm4,nb200_fiy(%esp)
        movaps %xmm5,nb200_fiz(%esp)
        ## update the fj's 
        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)

        shufps $225,%xmm0,%xmm0 ## constant 11100001
        shufps $225,%xmm1,%xmm1 ## constant 11100001
        shufps $225,%xmm2,%xmm2 ## constant 11100001

        movss   (%edi,%ebx,4),%xmm3
        movss   4(%edi,%ebx,4),%xmm4
        movss   8(%edi,%ebx,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%ebx,4)
        movss   %xmm4,4(%edi,%ebx,4)
        movss   %xmm5,8(%edi,%ebx,4)

_nb_kernel200_ia32_sse.nb200_checksingle:       
        movl  nb200_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200_ia32_sse.nb200_dosingle
        jmp    _nb_kernel200_ia32_sse.nb200_updateouterdata
_nb_kernel200_ia32_sse.nb200_dosingle:  
        movl nb200_charge(%ebp),%esi
        movl nb200_pos(%ebp),%edi
        movl  nb200_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge               

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb200_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb200_ix(%esp),%xmm4
        movaps nb200_iy(%esp),%xmm5
        movaps nb200_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb200_dx(%esp)
        movaps %xmm5,nb200_dy(%esp)
        movaps %xmm6,nb200_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subps  nb200_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        mulps  nb200_two(%esp),%xmm7

        subps  %xmm7,%xmm0
        mulps  %xmm0,%xmm3
        mulps  %xmm3,%xmm4      ## xmm4=total fscal 
        addss  nb200_vctot(%esp),%xmm6

        movl   nb200_faction(%ebp),%edi

        movaps nb200_dx(%esp),%xmm0
        movaps nb200_dy(%esp),%xmm1
        movaps nb200_dz(%esp),%xmm2

        movss %xmm6,nb200_vctot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb200_fix(%esp),%xmm3
        movaps nb200_fiy(%esp),%xmm4
        movaps nb200_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb200_fix(%esp)
        movaps %xmm4,nb200_fiy(%esp)
        movaps %xmm5,nb200_fiz(%esp)
        ## update fj 

        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)
_nb_kernel200_ia32_sse.nb200_updateouterdata: 
        movl  nb200_ii3(%esp),%ecx
        movl  nb200_faction(%ebp),%edi
        movl  nb200_fshift(%ebp),%esi
        movl  nb200_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb200_fix(%esp),%xmm0
        movaps nb200_fiy(%esp),%xmm1
        movaps nb200_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## increment fshift force  
        movss  (%esi,%edx,4),%xmm3
        movss  4(%esi,%edx,4),%xmm4
        movss  8(%esi,%edx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%esi,%edx,4)
        movss  %xmm4,4(%esi,%edx,4)
        movss  %xmm5,8(%esi,%edx,4)

        ## get n from stack
        movl nb200_n(%esp),%esi
        ## get group index for i particle 
        movl  nb200_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb200_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb200_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb200_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200_ia32_sse.nb200_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200_n(%esp)
        jmp _nb_kernel200_ia32_sse.nb200_outer
_nb_kernel200_ia32_sse.nb200_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200_ia32_sse.nb200_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200_ia32_sse.nb200_threadloop
_nb_kernel200_ia32_sse.nb200_end: 
        emms

        movl nb200_nouter(%esp),%eax
        movl nb200_ninner(%esp),%ebx
        movl nb200_outeriter(%ebp),%ecx
        movl nb200_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb200_salign(%esp),%eax
        addl %eax,%esp
        addl $300,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel200nf_ia32_sse
.globl _nb_kernel200nf_ia32_sse
nb_kernel200nf_ia32_sse:        
_nb_kernel200nf_ia32_sse:       
.set nb200nf_p_nri, 8
.set nb200nf_iinr, 12
.set nb200nf_jindex, 16
.set nb200nf_jjnr, 20
.set nb200nf_shift, 24
.set nb200nf_shiftvec, 28
.set nb200nf_fshift, 32
.set nb200nf_gid, 36
.set nb200nf_pos, 40
.set nb200nf_faction, 44
.set nb200nf_charge, 48
.set nb200nf_p_facel, 52
.set nb200nf_argkrf, 56
.set nb200nf_argcrf, 60
.set nb200nf_Vc, 64
.set nb200nf_type, 68
.set nb200nf_p_ntype, 72
.set nb200nf_vdwparam, 76
.set nb200nf_Vvdw, 80
.set nb200nf_p_tabscale, 84
.set nb200nf_VFtab, 88
.set nb200nf_invsqrta, 92
.set nb200nf_dvda, 96
.set nb200nf_p_gbtabscale, 100
.set nb200nf_GBtab, 104
.set nb200nf_p_nthreads, 108
.set nb200nf_count, 112
.set nb200nf_mtx, 116
.set nb200nf_outeriter, 120
.set nb200nf_inneriter, 124
.set nb200nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200nf_ix, 0
.set nb200nf_iy, 16
.set nb200nf_iz, 32
.set nb200nf_iq, 48
.set nb200nf_vctot, 64
.set nb200nf_half, 80
.set nb200nf_three, 96
.set nb200nf_krf, 112
.set nb200nf_crf, 128
.set nb200nf_is3, 144
.set nb200nf_ii3, 148
.set nb200nf_innerjjnr, 152
.set nb200nf_innerk, 156
.set nb200nf_n, 160
.set nb200nf_nn1, 164
.set nb200nf_nri, 168
.set nb200nf_facel, 172
.set nb200nf_nouter, 176
.set nb200nf_ninner, 180
.set nb200nf_salign, 184
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $188,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb200nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb200nf_p_nri(%ebp),%ecx
        movl nb200nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb200nf_nri(%esp)
        movl %esi,nb200nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb200nf_nouter(%esp)
        movl %eax,nb200nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb200nf_half(%esp)
        movss nb200nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb200nf_half(%esp)
        movaps %xmm3,nb200nf_three(%esp)

        movl nb200nf_argkrf(%ebp),%esi
        movl nb200nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb200nf_krf(%esp)
        shufps $0,%xmm6,%xmm6
        movaps %xmm6,nb200nf_crf(%esp)

_nb_kernel200nf_ia32_sse.nb200nf_threadloop: 
        movl  nb200nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel200nf_ia32_sse.nb200nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200nf_ia32_sse.nb200nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200nf_n(%esp)
        movl %ebx,nb200nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200nf_ia32_sse.nb200nf_outerstart
        jmp _nb_kernel200nf_ia32_sse.nb200nf_end

_nb_kernel200nf_ia32_sse.nb200nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200nf_nouter(%esp),%ebx
        movl %ebx,nb200nf_nouter(%esp)

_nb_kernel200nf_ia32_sse.nb200nf_outer: 
        movl  nb200nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb200nf_is3(%esp)            ## store is3 

        movl  nb200nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb200nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb200nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb200nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb200nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb200nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb200nf_ix(%esp)
        movaps %xmm1,nb200nf_iy(%esp)
        movaps %xmm2,nb200nf_iz(%esp)

        movl  %ebx,nb200nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb200nf_vctot(%esp)

        movl  nb200nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb200nf_pos(%ebp),%esi
        movl  nb200nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb200nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb200nf_ninner(%esp),%ecx
        movl  %ecx,nb200nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb200nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel200nf_ia32_sse.nb200nf_unroll_loop
        jmp   _nb_kernel200nf_ia32_sse.nb200nf_finish_inner
_nb_kernel200nf_ia32_sse.nb200nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb200nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb200nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb200nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb200nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  

        movl nb200nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb200nf_ix(%esp),%xmm4
        movaps nb200nf_iy(%esp),%xmm5
        movaps nb200nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addps  nb200nf_vctot(%esp),%xmm6
        movaps %xmm6,nb200nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb200nf_innerk(%esp)
        jl    _nb_kernel200nf_ia32_sse.nb200nf_finish_inner
        jmp   _nb_kernel200nf_ia32_sse.nb200nf_unroll_loop
_nb_kernel200nf_ia32_sse.nb200nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb200nf_innerk(%esp)
        movl  nb200nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel200nf_ia32_sse.nb200nf_dopair
        jmp   _nb_kernel200nf_ia32_sse.nb200nf_checksingle
_nb_kernel200nf_ia32_sse.nb200nf_dopair: 
        movl nb200nf_charge(%ebp),%esi

        movl  nb200nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb200nf_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges         

        movl nb200nf_pos(%ebp),%edi

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb200nf_iq(%esp),%xmm3

        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move ix-iz to xmm4-xmm6 

        movaps nb200nf_ix(%esp),%xmm4
        movaps nb200nf_iy(%esp),%xmm5
        movaps nb200nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 

        addps  nb200nf_vctot(%esp),%xmm6
        movaps %xmm6,nb200nf_vctot(%esp)

_nb_kernel200nf_ia32_sse.nb200nf_checksingle:   
        movl  nb200nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200nf_ia32_sse.nb200nf_dosingle
        jmp    _nb_kernel200nf_ia32_sse.nb200nf_updateouterdata
_nb_kernel200nf_ia32_sse.nb200nf_dosingle: 
        movl nb200nf_charge(%ebp),%esi
        movl nb200nf_pos(%ebp),%edi
        movl  nb200nf_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge 
        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb200nf_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb200nf_ix(%esp),%xmm4
        movaps nb200nf_iy(%esp),%xmm5
        movaps nb200nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb200nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        addss  nb200nf_vctot(%esp),%xmm6
        movss %xmm6,nb200nf_vctot(%esp)

_nb_kernel200nf_ia32_sse.nb200nf_updateouterdata: 
        ## get n from stack
        movl nb200nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb200nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb200nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb200nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

       ## finish if last 
        movl nb200nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200nf_ia32_sse.nb200nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200nf_n(%esp)
        jmp _nb_kernel200nf_ia32_sse.nb200nf_outer
_nb_kernel200nf_ia32_sse.nb200nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200nf_ia32_sse.nb200nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200nf_ia32_sse.nb200nf_threadloop
_nb_kernel200nf_ia32_sse.nb200nf_end: 
        emms

        movl nb200nf_nouter(%esp),%eax
        movl nb200nf_ninner(%esp),%ebx
        movl nb200nf_outeriter(%ebp),%ecx
        movl nb200nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb200nf_salign(%esp),%eax
        addl %eax,%esp
        addl $188,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


