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


## nb010 - forces are calculated
.globl nb_kernel010_ia32_sse
.globl _nb_kernel010_ia32_sse
nb_kernel010_ia32_sse:  
_nb_kernel010_ia32_sse: 
.set nb010_p_nri, 8
.set nb010_iinr, 12
.set nb010_jindex, 16
.set nb010_jjnr, 20
.set nb010_shift, 24
.set nb010_shiftvec, 28
.set nb010_fshift, 32
.set nb010_gid, 36
.set nb010_pos, 40
.set nb010_faction, 44
.set nb010_charge, 48
.set nb010_p_facel, 52
.set nb010_p_krf, 56
.set nb010_p_crf, 60
.set nb010_Vc, 64
.set nb010_type, 68
.set nb010_p_ntype, 72
.set nb010_vdwparam, 76
.set nb010_Vvdw, 80
.set nb010_p_tabscale, 84
.set nb010_VFtab, 88
.set nb010_invsqrta, 92
.set nb010_dvda, 96
.set nb010_p_gbtabscale, 100
.set nb010_GBtab, 104
.set nb010_p_nthreads, 108
.set nb010_count, 112
.set nb010_mtx, 116
.set nb010_outeriter, 120
.set nb010_inneriter, 124
.set nb010_work, 128
        ## The mutex (last arg) is not used in assembly.
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
.set nb010_is3, 272
.set nb010_ii3, 276
.set nb010_ntia, 280
.set nb010_innerjjnr, 284
.set nb010_innerk, 288
.set nb010_n, 292
.set nb010_nn1, 296
.set nb010_nri, 300
.set nb010_ntype, 304
.set nb010_nouter, 308
.set nb010_ninner, 312
.set nb010_salign, 316
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $320,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax          ## constant 16-byte align bottom of stack 
        subl %eax,%esp
        movl %eax,nb010_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb010_p_nri(%ebp),%ecx
        movl nb010_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb010_nri(%esp)
        movl %edi,nb010_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010_nouter(%esp)
        movl %eax,nb010_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x40000000,%eax   ## constant 2.0 in IEEE (hex)
        movl %eax,nb010_two(%esp)
        movss nb010_two(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm1,%xmm2      ## 4.0
        addps  %xmm1,%xmm2      ## 6.0
        movaps %xmm2,%xmm3
        addps  %xmm3,%xmm3      ## constant 12.0
        movaps %xmm1,nb010_two(%esp)
        movaps %xmm2,nb010_six(%esp)
        movaps %xmm3,nb010_twelve(%esp)

_nb_kernel010_ia32_sse.nb010_threadloop: 
        movl  nb010_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010_ia32_sse.nb010_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010_ia32_sse.nb010_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010_n(%esp)
        movl %ebx,nb010_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010_ia32_sse.nb010_outerstart
        jmp _nb_kernel010_ia32_sse.nb010_end

_nb_kernel010_ia32_sse.nb010_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010_nouter(%esp),%ebx
        movl %ebx,nb010_nouter(%esp)

_nb_kernel010_ia32_sse.nb010_outer: 
        movl  nb010_shift(%ebp),%eax            ## eax = base of shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx                ## ebx=3*is 
        movl  %ebx,nb010_is3(%esp)              ## store is3 

        movl  nb010_shiftvec(%ebp),%eax         ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb010_iinr(%ebp),%ecx             ## ecx = base of iinr[] 
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movl nb010_type(%ebp),%edx
        movl (%edx,%ebx,4),%edx
        imull nb010_ntype(%esp),%edx
        shll %edx
        movl %edx,nb010_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx                ## ebx = 3*ii=ii3 
        movl  nb010_pos(%ebp),%eax              ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb010_ix(%esp)
        movaps %xmm1,nb010_iy(%esp)
        movaps %xmm2,nb010_iz(%esp)

        movl  %ebx,nb010_ii3(%esp)

        ## clear Vvdwtot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb010_Vvdwtot(%esp)
        movaps %xmm4,nb010_fix(%esp)
        movaps %xmm4,nb010_fiy(%esp)
        movaps %xmm4,nb010_fiz(%esp)

        movl  nb010_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movl  nb010_pos(%ebp),%esi
        movl  nb010_faction(%ebp),%edi
        movl  nb010_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb010_ninner(%esp),%ecx
        movl  %ecx,nb010_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb010_innerk(%esp)           ## number of innerloop atoms 

        jge   _nb_kernel010_ia32_sse.nb010_unroll_loop
        jmp   _nb_kernel010_ia32_sse.nb010_finish_inner
_nb_kernel010_ia32_sse.nb010_unroll_loop: 
        ## quad-unrolled innerloop starts here 
        movl  nb010_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx                     ## eax-edx=jnr1-4 
        ## advance pointer (unrolled 4) 
        addl  $16,nb010_innerjjnr(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb010_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb010_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb010_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb010_c6(%esp)
        movaps %xmm6,nb010_c12(%esp)

        movl nb010_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb010_ix(%esp),%xmm4
        movaps nb010_iy(%esp),%xmm5
        movaps nb010_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb010_dx(%esp)
        movaps %xmm5,nb010_dy(%esp)
        movaps %xmm6,nb010_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4

        ## rsq in xmm4 
        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq

        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010_c6(%esp),%xmm1
        mulps  nb010_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010_Vvdwtot(%esp),%xmm5
        mulps  nb010_six(%esp),%xmm1
        mulps  nb010_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 

        movaps nb010_dx(%esp),%xmm0
        movaps nb010_dy(%esp),%xmm1
        movaps nb010_dz(%esp),%xmm2

        movaps %xmm5,nb010_Vvdwtot(%esp)

        movl   nb010_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb010_fix(%esp),%xmm3
        movaps nb010_fiy(%esp),%xmm4
        movaps nb010_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb010_fix(%esp)
        movaps %xmm4,nb010_fiy(%esp)
        movaps %xmm5,nb010_fiz(%esp)
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
        subl  $4,nb010_innerk(%esp)
        jl    _nb_kernel010_ia32_sse.nb010_finish_inner
        jmp   _nb_kernel010_ia32_sse.nb010_unroll_loop
_nb_kernel010_ia32_sse.nb010_finish_inner: 
        ## check if at least two particles remain 
        addl  $4,nb010_innerk(%esp)
        movl  nb010_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel010_ia32_sse.nb010_dopair
        jmp   _nb_kernel010_ia32_sse.nb010_checksingle
_nb_kernel010_ia32_sse.nb010_dopair: 
        movl  nb010_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl  $8,nb010_innerjjnr(%esp)

        movl nb010_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb010_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb010_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb010_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb010_c6(%esp)
        movaps %xmm6,nb010_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0


        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb010_faction(%ebp),%edi
        ## move nb010_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb010_ix(%esp),%xmm4
        movaps nb010_iy(%esp),%xmm5
        movaps nb010_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb010_dx(%esp)
        movaps %xmm5,nb010_dy(%esp)
        movaps %xmm6,nb010_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 


        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010_c6(%esp),%xmm1
        mulps  nb010_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010_Vvdwtot(%esp),%xmm5
        mulps  nb010_six(%esp),%xmm1
        mulps  nb010_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 

        movaps nb010_dx(%esp),%xmm0
        movaps nb010_dy(%esp),%xmm1
        movaps nb010_dz(%esp),%xmm2

        movaps %xmm5,nb010_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb010_fix(%esp),%xmm3
        movaps nb010_fiy(%esp),%xmm4
        movaps nb010_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb010_fix(%esp)
        movaps %xmm4,nb010_fiy(%esp)
        movaps %xmm5,nb010_fiz(%esp)
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

_nb_kernel010_ia32_sse.nb010_checksingle:       
        movl  nb010_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010_ia32_sse.nb010_dosingle
        jmp    _nb_kernel010_ia32_sse.nb010_updateouterdata
_nb_kernel010_ia32_sse.nb010_dosingle: 
        movl nb010_pos(%ebp),%edi
        movl  nb010_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        movl nb010_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb010_vdwparam(%ebp),%esi
        shll %ecx
        addl nb010_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb010_c6(%esp)
        movaps %xmm6,nb010_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        xorps   %xmm7,%xmm7

        movaps nb010_ix(%esp),%xmm4
        movaps nb010_iy(%esp),%xmm5
        movaps nb010_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb010_dx(%esp)
        movaps %xmm5,nb010_dy(%esp)
        movaps %xmm6,nb010_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010_c6(%esp),%xmm1
        mulps  nb010_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb010_Vvdwtot(%esp),%xmm5
        mulps  nb010_six(%esp),%xmm1
        mulps  nb010_twelve(%esp),%xmm2
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 

        movl   nb010_faction(%ebp),%edi

        movaps nb010_dx(%esp),%xmm0
        movaps nb010_dy(%esp),%xmm1
        movaps nb010_dz(%esp),%xmm2

        movss %xmm5,nb010_Vvdwtot(%esp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb010_fix(%esp),%xmm3
        movaps nb010_fiy(%esp),%xmm4
        movaps nb010_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb010_fix(%esp)
        movaps %xmm4,nb010_fiy(%esp)
        movaps %xmm5,nb010_fiz(%esp)
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

_nb_kernel010_ia32_sse.nb010_updateouterdata: 
        movl  nb010_ii3(%esp),%ecx
        movl  nb010_faction(%ebp),%edi
        movl  nb010_fshift(%ebp),%esi
        movl  nb010_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb010_fix(%esp),%xmm0
        movaps nb010_fiy(%esp),%xmm1
        movaps nb010_fiz(%esp),%xmm2

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
        movl nb010_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

         ## accumulate total lj energy and update it 
        movaps nb010_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb010_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb010_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010_ia32_sse.nb010_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010_n(%esp)
        jmp _nb_kernel010_ia32_sse.nb010_outer
_nb_kernel010_ia32_sse.nb010_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010_ia32_sse.nb010_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010_ia32_sse.nb010_threadloop
_nb_kernel010_ia32_sse.nb010_end: 
        emms

        movl nb010_nouter(%esp),%eax
        movl nb010_ninner(%esp),%ebx
        movl nb010_outeriter(%ebp),%ecx
        movl nb010_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb010_salign(%esp),%eax
        addl %eax,%esp          ## account for stack alignment 
        addl $320,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel010nf_ia32_sse
.globl _nb_kernel010nf_ia32_sse
nb_kernel010nf_ia32_sse:        
_nb_kernel010nf_ia32_sse:       
.set nb010nf_p_nri, 8
.set nb010nf_iinr, 12
.set nb010nf_jindex, 16
.set nb010nf_jjnr, 20
.set nb010nf_shift, 24
.set nb010nf_shiftvec, 28
.set nb010nf_fshift, 32
.set nb010nf_gid, 36
.set nb010nf_pos, 40
.set nb010nf_faction, 44
.set nb010nf_charge, 48
.set nb010nf_p_facel, 52
.set nb010nf_p_krf, 56
.set nb010nf_p_crf, 60
.set nb010nf_Vc, 64
.set nb010nf_type, 68
.set nb010nf_p_ntype, 72
.set nb010nf_vdwparam, 76
.set nb010nf_Vvdw, 80
.set nb010nf_p_tabscale, 84
.set nb010nf_VFtab, 88
.set nb010nf_invsqrta, 92
.set nb010nf_dvda, 96
.set nb010nf_p_gbtabscale, 100
.set nb010nf_GBtab, 104
.set nb010nf_p_nthreads, 108
.set nb010nf_count, 112
.set nb010nf_mtx, 116
.set nb010nf_outeriter, 120
.set nb010nf_inneriter, 124
.set nb010nf_work, 128
        ## The mutex (last arg) is not used in assembly.
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
.set nb010nf_is3, 144
.set nb010nf_ii3, 148
.set nb010nf_ntia, 152
.set nb010nf_innerjjnr, 156
.set nb010nf_innerk, 160
.set nb010nf_n, 164
.set nb010nf_nn1, 168
.set nb010nf_nri, 172
.set nb010nf_ntype, 176
.set nb010nf_nouter, 180
.set nb010nf_ninner, 184
.set nb010nf_salign, 188
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $192,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb010nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb010nf_p_nri(%ebp),%ecx
        movl nb010nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb010nf_nri(%esp)
        movl %edi,nb010nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010nf_nouter(%esp)
        movl %eax,nb010nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x40000000,%eax   ## constant 2.0 in IEEE (hex)
        movl %eax,nb010nf_two(%esp)
        movss nb010nf_two(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,nb010nf_two(%esp)

_nb_kernel010nf_ia32_sse.nb010nf_threadloop: 
        movl  nb010nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010nf_ia32_sse.nb010nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010nf_ia32_sse.nb010nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010nf_n(%esp)
        movl %ebx,nb010nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010nf_ia32_sse.nb010nf_outerstart
        jmp _nb_kernel010nf_ia32_sse.nb010nf_end

_nb_kernel010nf_ia32_sse.nb010nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010nf_nouter(%esp),%ebx
        movl %ebx,nb010nf_nouter(%esp)

_nb_kernel010nf_ia32_sse.nb010nf_outer: 
        movl  nb010nf_shift(%ebp),%eax          ## eax = base of shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx                ## ebx=3*is 
        movl  %ebx,nb010nf_is3(%esp)            ## store is3 

        movl  nb010nf_shiftvec(%ebp),%eax       ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb010nf_iinr(%ebp),%ecx           ## ecx = base of iinr[] 
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movl  nb010nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb010nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb010nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb010nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb010nf_ix(%esp)
        movaps %xmm1,nb010nf_iy(%esp)
        movaps %xmm2,nb010nf_iz(%esp)

        movl  %ebx,nb010nf_ii3(%esp)

        ## clear Vvdwtot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb010nf_Vvdwtot(%esp)

        movl  nb010nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb010nf_pos(%ebp),%esi
        movl  nb010nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb010nf_ninner(%esp),%ecx
        movl  %ecx,nb010nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb010nf_innerk(%esp)         ## number of innerloop atoms 

        jge   _nb_kernel010nf_ia32_sse.nb010nf_unroll_loop
        jmp   _nb_kernel010nf_ia32_sse.nb010nf_finish_inner
_nb_kernel010nf_ia32_sse.nb010nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb010nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        ## advance pointer (unrolled 4) 
        addl  $16,nb010nf_innerjjnr(%esp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb010nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb010nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb010nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb010nf_c6(%esp)
        movaps %xmm6,nb010nf_c12(%esp)

        movl nb010nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb010nf_ix(%esp),%xmm4
        movaps nb010nf_iy(%esp),%xmm5
        movaps nb010nf_iz(%esp),%xmm6

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
        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010nf_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%esp),%xmm1
        mulps  nb010nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb010nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl  $4,nb010nf_innerk(%esp)
        jl    _nb_kernel010nf_ia32_sse.nb010nf_finish_inner
        jmp   _nb_kernel010nf_ia32_sse.nb010nf_unroll_loop
_nb_kernel010nf_ia32_sse.nb010nf_finish_inner: 
        ## check if at least two particles remain 
        addl  $4,nb010nf_innerk(%esp)
        movl  nb010nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel010nf_ia32_sse.nb010nf_dopair
        jmp   _nb_kernel010nf_ia32_sse.nb010nf_checksingle
_nb_kernel010nf_ia32_sse.nb010nf_dopair: 
        movl  nb010nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl  $8,nb010nf_innerjjnr(%esp)

        movl nb010nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb010nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb010nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb010nf_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb010nf_c6(%esp)
        movaps %xmm6,nb010nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move nb010nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb010nf_ix(%esp),%xmm4
        movaps nb010nf_iy(%esp),%xmm5
        movaps nb010nf_iz(%esp),%xmm6

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


        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010nf_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%esp),%xmm1
        mulps  nb010nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb010nf_Vvdwtot(%esp)

_nb_kernel010nf_ia32_sse.nb010nf_checksingle:   
        movl  nb010nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010nf_ia32_sse.nb010nf_dosingle
        jmp    _nb_kernel010nf_ia32_sse.nb010nf_updateouterdata
_nb_kernel010nf_ia32_sse.nb010nf_dosingle: 
        movl nb010nf_pos(%ebp),%edi
        movl  nb010nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        movl nb010nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb010nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb010nf_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb010nf_c6(%esp)
        movaps %xmm6,nb010nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        xorps   %xmm7,%xmm7

        movaps nb010nf_ix(%esp),%xmm4
        movaps nb010nf_iy(%esp),%xmm5
        movaps nb010nf_iz(%esp),%xmm6

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

        rcpps %xmm4,%xmm5
        ## constant 1/x lookup seed in xmm5 
        movaps nb010nf_two(%esp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%esp),%xmm1
        mulps  nb010nf_c12(%esp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb010nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb010nf_Vvdwtot(%esp)

_nb_kernel010nf_ia32_sse.nb010nf_updateouterdata: 
        ## get n from stack
        movl nb010nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb010nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb010nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb010nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010nf_ia32_sse.nb010nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010nf_n(%esp)
        jmp _nb_kernel010nf_ia32_sse.nb010nf_outer
_nb_kernel010nf_ia32_sse.nb010nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010nf_ia32_sse.nb010nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010nf_ia32_sse.nb010nf_threadloop
_nb_kernel010nf_ia32_sse.nb010nf_end: 
        emms

        movl nb010nf_nouter(%esp),%eax
        movl nb010nf_ninner(%esp),%ebx
        movl nb010nf_outeriter(%ebp),%ecx
        movl nb010nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb010nf_salign(%esp),%eax
        addl %eax,%esp
        addl $192,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

